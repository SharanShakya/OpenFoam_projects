# methaneCHT: preliminary methane cooling-channel CHT study

An OpenFOAM v2412 learning and reproduction study inspired by:

Leilei Wang, Zunjing Chen and Hua Meng (2013), **Numerical study of conjugate heat
transfer of cryogenic methane in rectangular engine cooling channels at supercritical
pressures**, Applied Thermal Engineering 54(1), 237–246.
https://doi.org/10.1016/j.applthermaleng.2013.02.007

This is an independent implementation with modelling differences, not the authors'
original code, an exact reproduction, or an endorsed project.

## Current status

A coarse-mesh baseline reportedly completed 30,000 steady iterations on six MPI
processes. Medium and fine results are not available because of computational
limitations. Mesh independence and full validation have **not** been demonstrated.
No structural analysis has been completed in this repository.

| Log observation at iteration 30000 | Value |
|---|---:|
| Maximum fluid temperature | 210.33775 K |
| Maximum solid temperature | 224.43278 K |
| Reported fluid-interface heat-transfer integral | 1211.7357 W |
| Wall y+ minimum / average / maximum | 1.518 / 8.912 / 22.145 |

The documented default heating is 1350 W for the half-width model. The roughly
10.2% difference from the reported interface integral remains unresolved; it must
not be presented as a validated energy balance. The y+ values do not meet the
fine-grid target of y+ < 1. Stable output histories do not establish accuracy.
Confirm these observations against the archived logs and actual run settings.

## Model and limitations

- Steady `chtMultiRegionSimpleFoam`, fluid and solid regions, half-width symmetry.
- Default geometry: 0.6 m long, full channel width 2 mm, inner height 2 mm,
  wall thickness 0.5 mm; heating from x = 0.15 to 0.45 m.
- Default inlet: methane at 120 K and 15 m/s; outlet pressure 8 MPa.
- Default imposed flux: 3 MW/m2; solid conductivity: 100 W/(m K).
- CoolProp properties tabulated versus temperature at fixed pressure, rather than
  the paper's SRK/ECS treatment; k-omega SST rather than the paper's k-epsilon
  treatment. These differences must accompany any comparison with the paper.
- No combustion calculation, complete thrust chamber, or structural solver here.

## Reproducing the baseline

Use a fresh case directory. Load OpenFOAM v2412, install CoolProp (and numpy if
required by your environment), then use an explicit preset:

```bash
./Allrun.parallel --preset coarse --np 6
```

The generator defaults to fine without `--preset coarse`. The coarse preset
contains 312,000 total cells (192,000 fluid and 120,000 solid). Generated
configurations in `output/coarse/case_settings/` preserve the collected run inputs.
Original usage and model notes are retained in `docs/README.original.md`.

For an existing decomposed run, keep the processor folders and resume the solver
directly with `startFrom latestTime`; do not rerun the preparation script.
Reconstruct completed results using `reconstructPar -allRegions -latestTime`.

## Results and reproducibility

See `output/coarse/` for case settings, logs, monitor histories, exported tables
and figures as available. Figures must identify field, units, section location,
mesh, and iteration. The original processor checkpoints and reconstructed fields
are retained separately from ordinary Git history. If shared, add the release URL
and checksums here. No download of full fields is claimed until that is done.

The raw postProcessing data may contain overlapping restart histories; retain
provenance when combining them. Record OpenFOAM, Python and CoolProp versions,
actual input settings and generation command before publishing.

## Planned work

1. Check inlet/outlet mass balance and independently evaluate the thermal balance.
2. Investigate property treatment, interface coupling and heat-flux diagnostics.
3. Perform feasible mesh and near-wall sensitivity studies.
4. Compare with published results, reporting modelling differences and uncertainty.
5. Explore a separate thermal-to-structural workflow as a proposed master's topic.


