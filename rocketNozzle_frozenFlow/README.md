# Rocket Nozzle CFD — LOX/CH4 Thrust Chamber (10 kN Design Point)

CFD study of a LOX/CH4 rocket engine convergent-divergent nozzle, sized for a
**10 kN** design thrust and simulated in OpenFOAM using `rhoCentralFoam`.

The project is split into three parts: a Python pipeline that sizes the
nozzle and generates its mesh, the OpenFOAM case itself, and a set of
Python post-processing scripts that compute thrust and check convergence
directly from the CFD results.

## Repository structure

```
rocketNozzle_coldFlow/
├── case_study/       OpenFOAM case: mesh, boundary conditions, solver setup
├── modules/           Python scripts: sizing, geometry, post-processing
└── Output/            Saved logs / results from each run (log.thrust, log.convergence, ...)
```

---

## 1. Design & geometry pipeline (`modules/`)

Run **before** the OpenFOAM case, in this order, to go from design
parameters (target thrust, chamber pressure, mixture ratio) to a mesh:

| Script | Purpose |
|---|---|
| `thermochemistry.py` | Computes equilibrium combustion properties (`Tc`, `gamma`, `MW`) for a given chamber pressure and O/F mixture ratio, via Cantera's GRI-30 mechanism. Falls back to preliminary placeholder values if Cantera isn't installed — **always confirm Cantera is actually installed** (`pip install cantera --break-system-packages`) before sizing a real case, otherwise the geometry will be built from generic fallback numbers instead of your actual propellant chemistry. |
| `nozzle_sizing.py` | 1-D isentropic gas dynamics. Given target thrust, chamber pressure, ambient pressure, and O/F ratio, solves for exit Mach number, area ratio, throat/exit diameters, mass flow, and Isp. |
| `contour_generator.py` | Generates a Rao 80% bell nozzle contour from the throat/exit geometry, exports it as a CSV of (x, r) wall coordinates, and plots it. |
| `main.py` | Driver script — ties the three modules above together and runs the full sizing pipeline for a given design point. |
| `<csv_to_blockmesh>.py` | *(rename to match your actual script)* Converts the exported contour CSV into an OpenFOAM `blockMeshDict` (vertices, blocks, splines, wedge/wall/inlet/outlet patches). |

> **Note:** The `blockMeshDict` currently committed in `case_study/system/`
> was generated for the **10 kN** design point. If you change the target
> thrust, chamber pressure, or mixture ratio, you must re-run the sizing
> pipeline and regenerate `blockMeshDict` — don't hand-edit the mesh file
> directly, since its vertices/splines are derived from the sizing
> calculation.

---

## 2. OpenFOAM case (`case_study/`)

- **Solver**: `rhoCentralFoam` (density-based, compressible — appropriate for the transonic/supersonic flow through the nozzle)
- **Mesh**: axisymmetric 5° wedge representing the full engine (scaled by 72× for full 360° results)
- **Thermophysical model**: `hePsiThermo` / `pureMixture` / `perfectGas`, using a single equivalent gas with the equilibrium combustion properties from `thermochemistry.py` (frozen-flow assumption — the gas composition is fixed at chamber equilibrium values; no further chemical reaction is modeled as it expands through the nozzle)
- **Turbulence**: RAS `kOmegaSST`
- **Boundary conditions**:
  - `inlet`: `totalPressure` / `totalTemperature`, set to chamber conditions
  - `outlet`: `zeroGradient`
  - `nozzleWall`: no-slip, adiabatic wall
  - `wedgeMinus` / `wedgePlus`: `wedge` (axisymmetric representation)
- **Parallel execution**: domain decomposed 6-way (`scotch` method) and run via `mpirun`

### Live function objects (`system/controlDict`)

Two function objects run automatically during the solve (not as a separate post-processing step):

```
functions
{
    exitProperties
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             true;
        writeFields     false;
        regionType      patch;
        name            outlet;
        operation       areaAverage;
        fields          (p U rho);
    }

    forces
    {
        type            forces;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   10;
        patches         (nozzleWall);
        rho             rho;
        rhoInf          1;
        CofR            (0 0 0);
        pRef            101325;
        log             true;
    }
}
```

- **`exitProperties`** area-averages `p`, `U`, `rho` on the outlet patch each write interval — this is what the actual thrust calculation is based on.
- **`forces`** integrates pressure + viscous load on `nozzleWall`. This is used to independently verify the simulation has reached **steady state**, not as a direct thrust measurement — since the domain's `inlet` is an open boundary (not a modeled closed chamber dome), the wall-force-only integral is *not* directly comparable to total engine thrust (see [Design notes](#design-notes--known-issues) below).

### Running the case

```bash
cd case_study
chmod +x Allrun Allclean   # first time only
./Allrun
```

`Allrun` runs, in order: `blockMesh` → `checkMesh` → `decomposePar` →
`runParallel rhoCentralFoam` (6 processes) → `reconstructPar` → the two
post-processing scripts below.

To reset the case to a clean state before a fresh run:
```bash
./Allclean
```

---

## 3. Post-processing (`modules/`)

Both scripts run automatically at the end of `Allrun` and write their
output into `../Output/`. Their default case directory resolves relative
to their own file location (a sibling `case_study/` folder), so they can
also be run manually from anywhere without extra flags.

### `calculate_thrust.py`

Reads the full `exitProperties` time history and computes total thrust via
the standard exit-plane momentum balance, scaled from the modeled 5°
wedge slice to the full 360° engine:

```
F = ṁ · Ve + (Pe − Pa) · Ae
```

Reports the final (steady-state) value along with a quick check of how
much the result has moved over the trailing samples. For this design
point, the converged result is **≈9.76 kN** against the 10 kN target
(≈2.4% low — consistent with real viscous and multi-dimensional losses
the 1-D design theory doesn't capture).

```bash
python3 modules/calculate_thrust.py --case-dir case_study --json-out Output/log.thrust.json
```

### `check_convergence.py`

Two independent convergence checks:

- **Physical steady-state check** — from `postProcessing/forces/*/force.dat`, finds the longest trailing time window where the wall force stays within tolerance of its own mean.
- **Numerical residual check** — parses the raw solver log (`log.rhoCentralFoam`) for each equation's `Initial`/`Final residual` per timestep, confirming the final timestep is within tolerance and showing how far residuals dropped over the run.

```bash
python3 modules/check_convergence.py --case-dir case_study
```

Both scripts exit with status `0` if converged and `1` otherwise, so they
can be used in scripted checks if desired.

---

## Design notes / known issues

A few non-obvious things worth knowing if you're picking this project back up:

- **`forces` ≠ total engine thrust.** Because the CFD domain starts at an open `inlet` boundary rather than a closed chamber dome, integrating force over `nozzleWall` alone misses a large compensating pressure term. Use `exitProperties`-based `calculate_thrust.py` for the actual thrust number; use `forces`-based `check_convergence.py` only to confirm steady state.
- **Nozzle contour convergent-section fix.** An earlier version of `contour_generator.py` sampled the convergent-section arc beyond its valid geometric domain, creating a near-90° corner just upstream of the throat. This inflated the (non-thrust-relevant) radial wall force significantly but had negligible effect on the axial thrust number. Fixed by constraining the arc sampling range to match its radius of curvature.
- **`decomposeParDict` must be named exactly that** (not `decomposerParDict` or similar) — OpenFOAM fails silently-ish (`Error getting 'numberOfSubdomains'`) rather than reporting a missing file if it's misnamed.

---

## Requirements

- OpenFOAM (ESI/openfoam.com, v2412 or compatible)
- Python 3
- `pip install cantera scipy numpy pandas matplotlib --break-system-packages`
