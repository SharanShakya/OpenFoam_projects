"""
main.py
Ties together thermochemistry, sizing, contour generation, and OpenFOAM
blockMeshDict generation for a LOX/CH4 rocket nozzle design.
"""

from pathlib import Path

from nozzle_sizing import compute_nozzle_sizing
from contour_generator import generate_rao_contour, export_contour_csv, plot_contour
from csv_to_blockmesh import generate_blockmesh

# ==========================================
# DESIGN INPUTS
# ==========================================
F_target = 10_000.0     # Target thrust [N]
Pc = 30.0e5              # Chamber pressure [Pa]
Pa = 101_325.0           # Ambient pressure [Pa]
OF = 3.4                 # Oxidizer/Fuel mass ratio
OUTPUT_FOLDER = str(Path.home() / "Desktop" / "projectFiles" / "rocketNozzle" / "rocketNozzle_frozenFlow" / "output") # change this 
# OUTPUT_FOLDER = str(Path.home() / "Desktop" / "projectfile" / "case1" / "output")
# ==========================================
# BLOCKMESH INPUTS
# ==========================================
GENERATE_BLOCKMESH = True                       # Set False to skip CFD mesh step
BLOCKMESH_OUTPUT = Path("/home/mr-shakya/Desktop/projectFiles/rocketNozzle/rocketNozzle_frozenFlow/Case_Study/system/blockMeshDict") # change this
WEDGE_ANGLE = 5.0        # Total wedge angle [deg]
RADIAL_CELLS = 30        # Radial cells per axial block
AXIAL_CELLS = 180        # Approximate total axial cells
RADIAL_GRADING = 0.15    # <1 concentrates cells toward the wall

# ==========================================
# SIZING
# ==========================================
results = compute_nozzle_sizing(F_target, Pc, Pa, OF)

# ==========================================
# CONTOUR
# ==========================================
x_wall, r_wall, Ln = generate_rao_contour(results["Rt"], results["Re"], results["epsilon"])

csv_path = export_contour_csv(x_wall, r_wall, OUTPUT_FOLDER)
plot_path = plot_contour(
    x_wall, r_wall,
    results["Rt"], results["Re"], results["Dt"], results["De"], Ln,
    OUTPUT_FOLDER,
    title="10 kN LOX/CH4 Rocket Nozzle Geometry (80% Bell)"
)

# ==========================================
# OPENFOAM BLOCKMESH
# ==========================================
blockmesh_path = None
if GENERATE_BLOCKMESH:
    blockmesh_path = generate_blockmesh(
        Path(csv_path),
        Path(BLOCKMESH_OUTPUT),
        wedge_angle=WEDGE_ANGLE,
        radial_cells=RADIAL_CELLS,
        axial_cells=AXIAL_CELLS,
        radial_grading=RADIAL_GRADING,
    )

# ==========================================
# SUMMARY
# ==========================================
print("=== 1-D ROCKET NOZZLE DESIGN RESULTS ===")
print(f"Thrust                   : {F_target/1000:.3f} kN")
print(f"Throat Diameter (Dt)     : {results['Dt']*1000:.2f} mm")
print(f"Exit Diameter (De)       : {results['De']*1000:.2f} mm")
print(f"Nozzle Length (Ln)       : {Ln*1000:.2f} mm")
print(f"Specific Impulse (Isp)   : {results['Isp']:.2f} s")
print(f"Chemistry source         : {results['chem_source']}")
print("-" * 45)
print(f"CSV file saved to       : {csv_path}")
print(f"Plot image saved to     : {plot_path}")
if blockmesh_path is not None:
    print(f"blockMeshDict saved to  : {blockmesh_path}")
