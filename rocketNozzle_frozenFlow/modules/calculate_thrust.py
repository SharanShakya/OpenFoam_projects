#!/usr/bin/env python3
"""
calculate_thrust.py

Reads the exit-plane control-volume data written by the OpenFOAM
'exitProperties' (surfaceFieldValue, areaAverage) function object and
computes total engine thrust via:

    F = m_dot * Ve + (Pe - Pa) * Ae

scaled from the modeled wedge slice up to the full 360-degree engine.

Works whether exitProperties has one row (a single postProcess snapshot)
or many rows (the normal case when exitProperties runs live during the
solve, since it's listed in system/controlDict) - it always uses the row
with the LATEST time as the representative (steady-state) value, and
reports how much that value has moved over the trailing samples so you
can see at a glance whether it's actually settled.

Usage:
    python3 scripts/calculate_thrust.py
    python3 scripts/calculate_thrust.py --pa 101325 --wedge-angle 5
    python3 scripts/calculate_thrust.py --json-out log.thrust.json
"""

import argparse
import glob
import json
import os
import pathlib
import re
import sys

# This script now lives in a fixed location (e.g. project_root/modules/),
# with the OpenFOAM case as a sibling folder (e.g. project_root/Case_Study/).
# Resolving the default case directory relative to *this file's* location -
# rather than relative to the current working directory - means the script
# finds the right case whether it's run from inside the case folder (as
# Allrun does), from modules/, or from anywhere else.
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
DEFAULT_CASE_DIR = str(SCRIPT_DIR.parent / "Case_Study")

ROW_RE = re.compile(
    r"""^\s*
    ([\d.eE+\-]+)\s+          # time
    ([\d.eE+\-]+)\s+          # p (scalar)
    \(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)\s+  # U (vector)
    ([\d.eE+\-]+)\s*$         # rho (scalar)
    """,
    re.VERBOSE,
)
AREA_RE = re.compile(r"#\s*Area\s*:\s*([\d.eE+\-]+)")


def find_data_files(case_dir: str, function_name: str):
    pattern = os.path.join(case_dir, "postProcessing", function_name, "*", "surfaceFieldValue.dat")
    files = sorted(glob.glob(pattern))
    return files


def parse_file(path):
    """Return (area, [ (time, p, Ux, Uy, Uz, rho), ... ])."""
    area = None
    rows = []
    with open(path, "r") as f:
        for line in f:
            if area is None:
                m = AREA_RE.search(line)
                if m:
                    area = float(m.group(1))
                    continue
            if line.startswith("#") or not line.strip():
                continue
            m = ROW_RE.match(line)
            if m:
                t, p, ux, uy, uz, rho = (float(x) for x in m.groups())
                rows.append((t, p, ux, uy, uz, rho))
    return area, rows


def gather_all_rows(case_dir: str, function_name: str):
    files = find_data_files(case_dir, function_name)
    if not files:
        print(
            f"Error: no '{function_name}' surfaceFieldValue.dat found under "
            f"'{os.path.join(case_dir, 'postProcessing', function_name)}'.\n"
            f"Did the solver run with '{function_name}' defined in system/controlDict, "
            f"or did you mean to run postProcess first?"
        )
        sys.exit(1)

    area = None
    all_rows = []
    for path in files:
        file_area, rows = parse_file(path)
        if file_area is not None:
            area = file_area  # geometry doesn't change mid-run; last seen value is fine
        all_rows.extend(rows)

    if area is None or not all_rows:
        print(f"Error: failed to parse area or any data rows from {files}")
        sys.exit(1)

    # De-duplicate by time (keep the LAST occurrence in file order, in case
    # Allrun was re-run without cleaning postProcessing/) and sort by time.
    by_time = {}
    for row in all_rows:
        by_time[row[0]] = row
    sorted_rows = sorted(by_time.values(), key=lambda r: r[0])
    return area, sorted_rows


def compute_thrust(area, p, ux, rho, pa, wedge_angle):
    scale = 360.0 / wedge_angle
    mdot_slice = rho * ux * area
    mdot_total = mdot_slice * scale
    f_mom = mdot_slice * ux * scale
    f_press = (p - pa) * area * scale
    f_total = f_mom + f_press
    return mdot_total, f_mom, f_press, f_total


def steadiness_hint(rows, pa, wedge_angle, area, window=10):
    """Return a short human-readable line describing how much the computed
    thrust has moved over the trailing `window` samples."""
    if len(rows) < 2:
        return "Only one data point available - cannot assess steadiness from this file alone."
    tail = rows[-window:] if len(rows) >= window else rows
    totals = []
    for (t, p, ux, uy, uz, rho) in tail:
        _, _, _, f_total = compute_thrust(area, p, ux, rho, pa, wedge_angle)
        totals.append(f_total)
    spread = (max(totals) - min(totals))
    mean = sum(totals) / len(totals)
    pct = 100.0 * spread / abs(mean) if mean else float("nan")
    return (
        f"Last {len(tail)} samples (t={tail[0][0]:.4f}s to t={tail[-1][0]:.4f}s): "
        f"thrust varied by {pct:.3f}% (spread {spread/1000:.4f} kN around a "
        f"{mean/1000:.3f} kN mean)."
    )


def main():
    ap = argparse.ArgumentParser(description="Compute rocket nozzle thrust from OpenFOAM exit-plane data.")
    ap.add_argument("--case-dir", default=DEFAULT_CASE_DIR,
                     help=f"Case root directory (default: {DEFAULT_CASE_DIR}, "
                          f"i.e. the 'Case_Study' folder next to this script)")
    ap.add_argument("--function-name", default="exitProperties",
                     help="Name of the surfaceFieldValue function object (default: exitProperties)")
    ap.add_argument("--pa", type=float, default=101325.0, help="Ambient pressure [Pa] (default: 101325)")
    ap.add_argument("--wedge-angle", type=float, default=5.0, help="Mesh wedge angle in degrees (default: 5)")
    ap.add_argument("--json-out", default=None, help="Optional path to also write results as JSON")
    args = ap.parse_args()

    area, rows = gather_all_rows(args.case_dir, args.function_name)
    t_final, p_final, ux_final, uy_final, uz_final, rho_final = rows[-1]

    mdot_total, f_mom, f_press, f_total = compute_thrust(
        area, p_final, ux_final, rho_final, args.pa, args.wedge_angle
    )
    hint = steadiness_hint(rows, args.pa, args.wedge_angle, area)

    report = f"""
=================================================
  AUTOMATED THRUST CALCULATION REPORT (t = {t_final:.6g}s)
=================================================
Data points found    : {len(rows)} (t = {rows[0][0]:.6g}s to {rows[-1][0]:.6g}s)
Wedge slice area      : {area:.6e} m^2
Exit Pressure (p_e)   : {p_final:.2f} Pa
Axial Velocity (Ux)   : {ux_final:.2f} m/s
Radial Velocity (Uy)  : {uy_final:.2f} m/s
Exit Density (rho)    : {rho_final:.6f} kg/m^3
Ambient pressure (Pa) : {args.pa:.2f} Pa
Wedge angle           : {args.wedge_angle:.3f} deg  (scale factor {360.0/args.wedge_angle:.3f}x)
-------------------------------------------------
Total Mass Flow       : {mdot_total:.3f} kg/s
Momentum Thrust       : {f_mom / 1000.0:.3f} kN
Pressure Thrust       : {f_press / 1000.0:.3f} kN
TOTAL NET THRUST      : {f_total / 1000.0:.3f} kN
-------------------------------------------------
Steadiness check      : {hint}
=================================================
"""
    print(report)

    if args.json_out:
        payload = {
            "time": t_final,
            "n_samples": len(rows),
            "area_m2": area,
            "p_pa": p_final,
            "Ux_mps": ux_final,
            "Uy_mps": uy_final,
            "rho_kgm3": rho_final,
            "pa_pa": args.pa,
            "wedge_angle_deg": args.wedge_angle,
            "mdot_total_kgps": mdot_total,
            "F_momentum_N": f_mom,
            "F_pressure_N": f_press,
            "F_total_N": f_total,
            "F_total_kN": f_total / 1000.0,
        }
        with open(args.json_out, "w") as f:
            json.dump(payload, f, indent=2)


if __name__ == "__main__":
    main()
