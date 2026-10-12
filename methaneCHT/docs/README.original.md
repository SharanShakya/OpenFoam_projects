# Cryogenic methane CHT in a rectangular cooling channel (OpenFOAM v2412)

Replication set-up for Wang, Chen & Meng, *Appl. Therm. Eng.* 54 (2013) 237-246.
Solver: `chtMultiRegionSimpleFoam` (ESI OpenFOAM v2412). Half-width model (symmetry at z = 0).

> **Status:** written without access to an OpenFOAM installation, so the case has **not been run**.
> `tools/makeCase.py` was tested (mesh topology, file generation); the OpenFOAM dictionaries were
> checked against my knowledge of v2412 only. See "Things to verify first" below.

## Quick start

```bash
pip install CoolProp numpy            # property tables
source /path/to/OpenFOAM-v2412/etc/bashrc
./Allrun --preset coarse              # ~0.3 M cells, quick test
./Allrun.parallel --np 64             # paper mesh (1000x110x55 = 6.05 M cells)
```
`Allrun` first runs `tools/makeCase.py "$@"`, so every option below can be passed straight to it.

| Preset | Cells (x y z) | first-cell height | use |
|---|---|---|---|
| `coarse` | 240 x 52 x 25 | 5 um | debugging (wall-function regime) |
| `medium` | 480 x 80 x 40 | 1 um | grid-independence step |
| `fine` (default) | 1000 x 110 x 55 | 0.25 um | paper mesh, y+ < 1 expected |

## Case matrix from the paper

```bash
# Section 3.1 - solid conductivity (square channel, U = 15 m/s)
for l in 20 50 100 200 400; do ./Allrun.parallel --lambda-s $l; done      # use a copy of the case per run

# Section 3.2 - aspect ratio (constant mass flow, lambda_s = 100)
./Allrun.parallel --H 2e-3 --U 20      # R = 1
./Allrun.parallel --H 4e-3 --U 10      # R = 2
./Allrun.parallel --H 6e-3 --U 6.67    # R = 3
./Allrun.parallel --H 8e-3 --U 5       # R = 4
```
Fixed for all runs: T_in = 120 K, p = 8 MPa, q = 3 MW/m2 on the top exterior surface for x = 0.15-0.45 m,
all other exterior surfaces adiabatic, W = 2 mm, wall = 0.5 mm, L = 0.6 m.

## Layout

```
Allrun / Allrun.parallel / Allclean
tools/makeCase.py                 generator (edit options, not the generated files)
constant/regionProperties         fluid, solid
constant/fluid/                   thermophysicalProperties (generated), turbulenceProperties, g, pRef
constant/solid/                   thermophysicalProperties (generated; only kappa matters)
system/controlDict                function objects (mass flow, heat flux, y+, min/max, bulk planes)
system/{fluid,solid}/             fvSchemes, fvSolution (+ decomposeParDict, topoSetDict generated)
0.orig/{fluid,solid}/             fields (generated)
```
Patches: fluid `inlet outlet symmetryFluid fluid_to_solid`; solid `topHeated solidInsulated symmetrySolid solid_to_fluid`.

## Modelling choices and deviations from the paper

| Item | Paper | This case | Comment |
|---|---|---|---|
| Equation of state / transport | SRK + extended corresponding states | CoolProp (Helmholtz EOS) tables, **function of T at p = 8 MPa** | Channel pressure drop is only ~0.01-0.18 MPa; shifts of the pseudo-critical point are small. Expect a few % difference in rho, mu, lambda vs the paper. Upgrade path: 2-D (p,T) custom thermo. |
| Turbulence | k-eps, enhanced wall treatment, y+ < 1 | k-omega SST, resolved wall (`kLowReWallFunction`, `omegaWallFunction`, `nutUSpaldingWallFunction`) | Not identical; results near the pseudo-critical wall region can differ. |
| Inlet turbulence | not given | I = 5 %, l = 0.07 Dh, k = 1.5 (U I)2, omega = sqrt(k)/(C_mu^0.25 l) | assumption |
| Outlet | extrapolation | fixed p = 8 MPa | assumption (paper: "operation pressure 8 MPa") |
| Solid rho, Cp | not needed | 8000 kg/m3, 500 J/kg/K | irrelevant in steady state |
| Gravity | - | g = 0 (p_rgh = p) | negligible |
| Bulk temperature | Eq. 7 (rho u Cp weighted) | mass-flux weighted T and h at x-planes | invert h -> T with the Hs table for a true mixing-cup temperature; Eq. 7 differs slightly |

## Validation targets (paper)

Energy check: the imposed power in the half model is q x 0.3 m x (W/2 + t_w) = **1350 W**
(`wallHeatFlux` integral on `fluid_to_solid` must converge to this).

*Section 3.1 (U = 15 m/s, H = 2 mm)*
- T_b rises from 120 K to about 150 K at x = 0.45 m, independent of lambda_s (Fig. 2a).
- Top wall exterior max T about 360 K (lambda_s = 20) vs 195 K (400); exterior about 75 K above interior at 20 (Fig. 2b).
- Interior heat flux, top / side / bottom: 2.8 / 0.8 / <0.1 MW/m2 at lambda_s = 20; 1.6 / 1.1 / 0.85 MW/m2 at 400 (Fig. 5).
- Nu at top wall within about 10-20 % of the modified Jackson & Hall correlation, Eq. 8 (exponent n is not given in the paper: see Ruan & Meng 2012, Ref. [18]).

*Section 3.2 (lambda_s = 100, constant mass flow)*

| H (mm) | U (m/s) | Re | dp over heated section (MPa) | eta |
|---|---|---|---|---|
| 2 | 20 | 1.551e5 | 0.180 | 1 |
| 4 | 10 | 1.034e5 | 0.036 | 0.777 |
| 6 | 6.67 | 0.776e5 | 0.015 | 0.656 |
| 8 | 5 | 0.620e5 | 0.009 | 0.557 |

T_b differs by < 3 K between aspect ratios; max exterior top wall T about 200 K (R = 1) to 360 K (R = 4).
`makeCase.py` prints the inlet Re from the CoolProp tables: compare with the first row as a property check.

## Post-processing

- `postProcessing/bulk_planeN/.../surfaceFieldValue.dat`: mass-weighted T and h at plane N (positions in `system/planes.txt`).
- `postProcessing/pAvg_planeN/...`: area-averaged p; dp over the heated section = pAvg_plane0 - pAvg_plane6.
  Friction factor f = (dp/L) (2 Dh / (rho u^2)); eta from Eqs. 10-11.
- `postProcessing/topHeatedT/...`: mean exterior top-wall temperature (T_w,twe in Eq. 11). h_eq = q / (T_w,twe - T_b).
- `postProcessing/wallHeatFlux`, `yPlus` (log + written fields), `inletMassFlow`, `outletMassFlow` (mass balance).
- Profiles along x for top/side/bottom interior walls (Figs. 2-5, 7): sample in ParaView (`paraFoam -region fluid`).

## Things to verify first (most likely places for small syntax fixes)

1. **Tabulated thermo keywords** in `constant/fluid/thermophysicalProperties`
   (`icoTabulated`, `hTabulated`, `tabulated` transport with `rho`, `Cp`/`Hs`, `mu`/`kappa` tables).
   Check the class headers in `$WM_PROJECT_DIR/src/thermophysicalModels/specie/` for your exact release.
2. `externalWallHeatFluxTemperature` entries on `solid/T` (`mode flux; q uniform ...`).
3. `turbulenceProperties`: `model kOmegaSST;` (older releases used `RASModel`).
4. The interface patches must be called `fluid_to_solid` / `solid_to_fluid` after `splitMeshRegions`
   (check `constant/fluid/polyMesh/boundary`); if they differ, rename them in `0.orig` and `controlDict`.
5. `boxToFace`/`setToFaceZone` in `system/fluid/topoSetDict`; if zones come out empty, check `log.topoSet.fluid`.
6. After the first iterations, check `yPlus` (target max < 1). Adjust `--y1` and regenerate if needed.

## Convergence tips

- Start with `--preset coarse`, then `medium`, then `fine`; use the coarse result as a sanity check of signs and magnitudes.
- Near the pseudo-critical temperature (printed by `makeCase.py`, roughly 200 K at 8 MPa) lower relaxation factors (U, h, k, omega 0.5; p_rgh 0.2)
  and switch `div(phi,U)`, `div(phi,h)` to `bounded Gauss upwind` for the first few thousand iterations.
- Ramp the heat flux: generate with `--q 3e5`, run ~2000 iterations, then edit `q` in `<latestTime>/solid/T`
  (`sed -i 's/uniform 3e+05/uniform 1e+06/' <t>/solid/T`), restart, and repeat up to 3e6.
- If steady state will not converge, switch to `chtMultiRegionFoam` with a large `deltaT` (and set solid rho, Cp).
