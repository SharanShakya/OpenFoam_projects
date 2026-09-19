#!/usr/bin/env python3
"""
Convert an axisymmetric rocket-nozzle (x,r) CSV contour to an OpenFOAM
blockMeshDict using a small-angle wedge sector.

Input CSV columns:
    x_m, r_m

The generated mesh is an axisymmetric wedge. The nozzle wall is represented
by blockMesh spline edges so the CFD mesh follows the supplied contour.

Usage (CLI):
    python csv_to_blockmesh.py nozzle_contour.csv system/blockMeshDict

Usage (as a module):
    from csv_to_blockmesh import generate_blockmesh
    generate_blockmesh(Path("nozzle_contour.csv"), Path("system/blockMeshDict"))

The script assumes:
    - x is the axial direction
    - r is the radius from the axis
    - the smallest-r contour point is the throat
    - the domain extends from the axis (r=0) to the nozzle wall

The output is intended for a standard OpenFOAM blockMesh workflow.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np


def read_contour(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read and clean x,r data from a CSV file."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.dtype.names is None:
        raise ValueError("CSV must have a header row containing x_m and r_m.")

    names = {n.lower(): n for n in data.dtype.names}
    if "x_m" not in names or "r_m" not in names:
        raise ValueError("CSV must contain columns named x_m and r_m.")

    x = np.asarray(data[names["x_m"]], dtype=float)
    r = np.asarray(data[names["r_m"]], dtype=float)

    mask = np.isfinite(x) & np.isfinite(r) & (r >= 0.0)
    x, r = x[mask], r[mask]
    if len(x) < 4:
        raise ValueError("Not enough valid contour points.")

    order = np.argsort(x)
    x, r = x[order], r[order]

    unique_x = []
    unique_r = []
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and math.isclose(x[j], x[i], rel_tol=0.0, abs_tol=1e-12):
            j += 1
        unique_x.append(x[i])
        unique_r.append(np.min(r[i:j]))
        i = j

    x = np.asarray(unique_x)
    r = np.asarray(unique_r)

    if np.any(np.diff(x) <= 0):
        raise ValueError("x coordinates must be strictly increasing after cleanup.")

    return x, r


def interp_radius(xp: np.ndarray, x: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Linear interpolation of contour radius at requested x locations."""
    return np.interp(xp, x, r)


def choose_breakpoints(x: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Choose sensible multi-block break locations around the throat."""
    throat_i = int(np.argmin(r))
    xt = float(x[throat_i])
    xmin = float(x[0])
    xmax = float(x[-1])

    left_fracs = [0.0, 0.30, 0.60, 0.82, 1.0]
    right_fracs = [0.0, 0.18, 0.40, 0.70, 1.0]

    left = [xmin + f * (xt - xmin) for f in left_fracs]
    right = [xt + f * (xmax - xt) for f in right_fracs]

    bp = np.array(left[:-1] + right, dtype=float)
    bp[0] = xmin
    bp[-1] = xmax

    bp = np.unique(np.round(bp, 14))
    return bp


def distribute_cells(lengths: np.ndarray, total: int, minimum: int = 4) -> np.ndarray:
    """Distribute axial cells approximately in proportion to block length."""
    n = len(lengths)
    if total < n * minimum:
        total = n * minimum

    remaining = total - n * minimum
    if remaining <= 0:
        return np.full(n, minimum, dtype=int)

    weights = lengths / np.sum(lengths)
    raw = remaining * weights
    extra = np.floor(raw).astype(int)
    leftover = remaining - int(np.sum(extra))

    frac = raw - extra
    for idx in np.argsort(-frac)[:leftover]:
        extra[idx] += 1

    return extra + minimum


def point3(x: float, r: float, theta: float) -> tuple[float, float, float]:
    """Map x-r coordinates into a small 3-D wedge."""
    return (x, r * math.cos(theta), r * math.sin(theta))


def fmt_point(p: tuple[float, float, float]) -> str:
    return f"({p[0]:.12g} {p[1]:.12g} {p[2]:.12g})"


def make_spline_points(
    x0: float,
    x1: float,
    x: np.ndarray,
    r: np.ndarray,
    theta: float,
) -> list[tuple[float, float, float]]:
    """Return intermediate wall points strictly between block endpoints."""
    mask = (x > x0 + 1e-12) & (x < x1 - 1e-12)
    return [point3(float(xi), float(ri), theta) for xi, ri in zip(x[mask], r[mask])]


def build_dict(
    x: np.ndarray,
    r: np.ndarray,
    breakpoints: np.ndarray,
    wedge_angle_deg: float,
    radial_cells: int,
    total_axial_cells: int,
    radial_grading: float,
) -> str:
    half = math.radians(wedge_angle_deg / 2.0)

    rb = interp_radius(breakpoints, x, r)
    nblocks = len(breakpoints) - 1
    lengths = np.diff(breakpoints)
    axial_cells = distribute_cells(lengths, total_axial_cells, minimum=4)

    lines: list[str] = []
    add = lines.append

    add("FoamFile")
    add("{")
    add("    version     2.0;")
    add("    format      ascii;")
    add("    class       dictionary;")
    add('    object      blockMeshDict;')
    add("}")
    add("")
    add(f"// Generated from x-r contour CSV")
    add(f"// Wedge angle = {wedge_angle_deg:g} degrees")
    add(f"// Blocks = {nblocks}, radial cells/block = {radial_cells}")
    add("")
    add("scale 1.0;")
    add("")

    add("vertices")
    add("(")
    for xi, ri in zip(breakpoints, rb):
        add(f"    // x = {xi:.12g} m, r = {ri:.12g} m")
        add(f"    {fmt_point(point3(float(xi), 0.0, -half))}")
        add(f"    {fmt_point(point3(float(xi), float(ri), -half))}")
        add(f"    {fmt_point(point3(float(xi), 0.0, +half))}")
        add(f"    {fmt_point(point3(float(xi), float(ri), +half))}")
    add(");")
    add("")

    def v(i: int, kind: str) -> int:
        base = 4 * i
        return {
            "am": base + 0,
            "wm": base + 1,
            "ap": base + 2,
            "wp": base + 3,
        }[kind]

    add("blocks")
    add("(")
    for i in range(nblocks):
        ids = [
            v(i, "am"),
            v(i + 1, "am"),
            v(i + 1, "wm"),
            v(i, "wm"),
            v(i, "ap"),
            v(i + 1, "ap"),
            v(i + 1, "wp"),
            v(i, "wp"),
        ]
        add(
            f"    hex ({' '.join(map(str, ids))}) "
            f"({int(axial_cells[i])} {int(radial_cells)} 1) "
            f"simpleGrading (1 {radial_grading:g} 1)"
        )
    add(");")
    add("")

    add("edges")
    add("(")
    for i in range(nblocks):
        x0, x1 = float(breakpoints[i]), float(breakpoints[i + 1])
        for theta, vm0, vm1 in [
            (-half, v(i, "wm"), v(i + 1, "wm")),
            (+half, v(i, "wp"), v(i + 1, "wp")),
        ]:
            pts = make_spline_points(x0, x1, x, r, theta)
            if pts:
                add(f"    spline {vm0} {vm1}")
                add("    (")
                for p in pts:
                    add(f"        {fmt_point(p)}")
                add("    )")
    add(");")
    add("")

    add("boundary")
    add("(")
    add("    inlet")
    add("    {")
    add("        type patch;")
    add("        faces")
    add("        (")
    add(f"            ({v(0,'am')} {v(0,'ap')} {v(0,'wp')} {v(0,'wm')})")
    add("        );")
    add("    }")
    add("")

    add("    outlet")
    add("    {")
    add("        type patch;")
    add("        faces")
    add("        (")
    add(f"            ({v(nblocks,'am')} {v(nblocks,'wm')} {v(nblocks,'wp')} {v(nblocks,'ap')})")
    add("        );")
    add("    }")
    add("")

    add("    nozzleWall")
    add("    {")
    add("        type wall;")
    add("        faces")
    add("        (")
    for i in range(nblocks):
        add(f"            ({v(i,'wm')} {v(i+1,'wm')} {v(i+1,'wp')} {v(i,'wp')})")
    add("        );")
    add("    }")
    add("")

    add("    wedgeMinus")
    add("    {")
    add("        type wedge;")
    add("        faces")
    add("        (")
    for i in range(nblocks):
        add(f"            ({v(i,'am')} {v(i,'wm')} {v(i+1,'wm')} {v(i+1,'am')})")
    add("        );")
    add("    }")
    add("")

    add("    wedgePlus")
    add("    {")
    add("        type wedge;")
    add("        faces")
    add("        (")
    for i in range(nblocks):
        add(f"            ({v(i,'ap')} {v(i+1,'ap')} {v(i+1,'wp')} {v(i,'wp')})")
    add("        );")
    add("    }")
    add(");")
    add("")

    add("mergePatchPairs")
    add("(")
    add(");")
    add("")

    return "\n".join(lines)


def generate_blockmesh(
    csv_path: Path,
    output_path: Path,
    wedge_angle: float = 5.0,
    radial_cells: int = 30,
    axial_cells: int = 180,
    radial_grading: float = 0.15,
    verbose: bool = True,
) -> Path:
    """
    Read a nozzle contour CSV and write an OpenFOAM blockMeshDict.

    This is the importable equivalent of running the script's CLI.

    Parameters
    ----------
    csv_path : Path
        Input contour CSV with x_m, r_m columns.
    output_path : Path
        Output path for blockMeshDict (typically system/blockMeshDict).
    wedge_angle : float
        Total wedge angle in degrees.
    radial_cells : int
        Radial cells per axial block.
    axial_cells : int
        Approximate total axial cells.
    radial_grading : float
        Simple radial grading; <1 concentrates cells toward the wall.
    verbose : bool
        Print a summary to stdout.

    Returns
    -------
    Path
        The path the blockMeshDict was written to.
    """
    if wedge_angle <= 0 or wedge_angle >= 30:
        raise ValueError("Use a small wedge angle between 0 and 30 degrees.")
    if radial_cells < 4:
        raise ValueError("radial_cells must be >= 4")
    if axial_cells < 20:
        raise ValueError("axial_cells must be >= 20")
    if radial_grading <= 0:
        raise ValueError("radial_grading must be > 0")

    csv_path = Path(csv_path)
    output_path = Path(output_path)

    x, r = read_contour(csv_path)
    bp = choose_breakpoints(x, r)
    text = build_dict(x, r, bp, wedge_angle, radial_cells, axial_cells, radial_grading)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text)

    if verbose:
        throat_i = int(np.argmin(r))
        throat_x = x[throat_i]
        throat_r = r[throat_i]

        print("Generated OpenFOAM blockMeshDict")
        print(f"  Input contour : {csv_path}")
        print(f"  Output        : {output_path}")
        print(f"  Points        : {len(x)}")
        print(f"  Throat        : x={throat_x:.8g} m, r={throat_r:.8g} m")
        print(f"  Wedge angle   : {wedge_angle:g} deg")
        print(f"  Blocks        : {len(bp)-1}")
        print(f"  Radial cells  : {radial_cells}/block")
        print(f"  Axial cells   : ~{axial_cells} total")
        print("  Breakpoints   :")
        for value in bp:
            print(f"      {value:.8g} m")

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="Input contour CSV")
    parser.add_argument("output", type=Path, help="Output system/blockMeshDict")
    parser.add_argument("--wedge-angle", type=float, default=5.0,
                        help="Total wedge angle in degrees (default: 5)")
    parser.add_argument("--radial-cells", type=int, default=30,
                        help="Radial cells per axial block (default: 30)")
    parser.add_argument("--axial-cells", type=int, default=180,
                        help="Approximate total axial cells (default: 180)")
    parser.add_argument("--radial-grading", type=float, default=0.15,
                        help="Simple radial grading; <1 concentrates cells toward wall (default: 0.15)")
    args = parser.parse_args()

    generate_blockmesh(
        args.csv,
        args.output,
        wedge_angle=args.wedge_angle,
        radial_cells=args.radial_cells,
        axial_cells=args.axial_cells,
        radial_grading=args.radial_grading,
    )


if __name__ == "__main__":
    main()
