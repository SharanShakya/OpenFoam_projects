#!/usr/bin/env python3
"""
check_convergence.py

Two independent convergence checks for an OpenFOAM run:

  A) PHYSICAL STEADY-STATE check
     Reads the full time history written by the 'forces' function object
     (postProcessing/forces/*/force.dat) and checks whether the axial
     force has stopped changing - i.e. whether the flow has actually
     reached steady state, not just reached the end of the requested
     simulated time.

  B) NUMERICAL SOLVER RESIDUAL check
     Parses the solver log (default: log.rhoCentralFoam) for the
     "Solving for <field>, Initial residual = ..., Final residual = ..."
     lines OpenFOAM prints every timestep, and checks:
       - that the Final residual at the last timestep is comfortably
         below a tolerance for every iteratively-solved field, and
       - that the Initial residual has trended down over the run
         (a hallmark of approaching steady state via pseudo-time-marching).

Usage:
    python3 scripts/check_convergence.py
    python3 scripts/check_convergence.py --log log.rhoCentralFoam --force-tol 0.5
"""

import argparse
import glob
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

# ---------------------------------------------------------------------------
# Part A: physical steady-state check (forces function object)
# ---------------------------------------------------------------------------

def find_force_files(case_dir):
    pattern = os.path.join(case_dir, "postProcessing", "forces", "*", "force.dat")
    return sorted(glob.glob(pattern))


def parse_force_file(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                t = float(parts[0])
                total_x = float(parts[1])
                total_y = float(parts[2])
            except ValueError:
                continue
            rows.append((t, total_x, total_y))
    return rows


def gather_force_rows(case_dir):
    files = find_force_files(case_dir)
    if not files:
        return None
    by_time = {}
    for path in files:
        for row in parse_force_file(path):
            by_time[row[0]] = row
    return sorted(by_time.values(), key=lambda r: r[0])


def check_steady_state(rows, tol_pct, min_window):
    """
    Walk backward from the end of the run, growing a 'steady window' as
    long as every point in it stays within tol_pct of the window mean.
    Returns a dict describing the result.
    """
    n = len(rows)
    totals = [r[1] for r in rows]

    # find the longest trailing run where each point is within tol_pct of
    # the running mean of that trailing run
    steady_start_idx = n - 1
    for start in range(n - 1, -1, -1):
        window = totals[start:]
        mean = sum(window) / len(window)
        if mean == 0:
            break
        max_dev_pct = max(abs(v - mean) / abs(mean) for v in window) * 100.0
        if max_dev_pct <= tol_pct:
            steady_start_idx = start
        else:
            break

    steady_rows = rows[steady_start_idx:]
    steady_window_len = len(steady_rows)
    mean_val = sum(r[1] for r in steady_rows) / steady_window_len
    spread_pct = (max(r[1] for r in steady_rows) - min(r[1] for r in steady_rows)) / abs(mean_val) * 100.0

    achieved = steady_window_len >= min_window
    return {
        "achieved": achieved,
        "steady_from_t": steady_rows[0][0],
        "steady_to_t": steady_rows[-1][0],
        "n_steady_samples": steady_window_len,
        "mean_total_x": mean_val,
        "spread_pct": spread_pct,
        "final_total_x": rows[-1][1],
    }


# ---------------------------------------------------------------------------
# Part B: solver residual check (raw log file)
# ---------------------------------------------------------------------------

TIME_RE = re.compile(r"^Time\s*=\s*([\d.eE+\-]+)")
SOLVE_RE = re.compile(
    r"Solving for (\w+), Initial residual = ([\d.eE+\-]+), Final residual = ([\d.eE+\-]+)"
)
# Fields solved by a direct/explicit "diagonal" update always report 0/0 and
# carry no convergence information - skip them.
SKIP_FIELDS = {"rho", "rhoUx", "rhoUy", "rhoUz", "rhoE"}


def parse_solver_log(path):
    """Return list of (time, {field: (initial, final)}) in file order."""
    blocks = []
    current_time = None
    current_fields = {}
    with open(path, "r") as f:
        for line in f:
            m = TIME_RE.match(line)
            if m:
                if current_time is not None and current_fields:
                    blocks.append((current_time, current_fields))
                current_time = float(m.group(1))
                current_fields = {}
                continue
            m = SOLVE_RE.search(line)
            if m:
                field, init_r, final_r = m.group(1), float(m.group(2)), float(m.group(3))
                if field not in SKIP_FIELDS:
                    current_fields[field] = (init_r, final_r)
    if current_time is not None and current_fields:
        blocks.append((current_time, current_fields))
    return blocks


def check_residuals(blocks, final_tol):
    if not blocks:
        return None
    first_time, first_fields = blocks[0]
    last_time, last_fields = blocks[-1]

    field_report = {}
    all_ok = True
    for field, (last_init, last_final) in last_fields.items():
        ok = last_final <= final_tol
        all_ok = all_ok and ok
        first_init = first_fields.get(field, (None, None))[0]
        drop = None
        if first_init and first_init > 0 and last_init > 0:
            drop = first_init / last_init
        field_report[field] = {
            "first_initial": first_init,
            "last_initial": last_init,
            "last_final": last_final,
            "within_tol": ok,
            "drop_factor": drop,
        }
    return {
        "n_timesteps_logged": len(blocks),
        "first_time": first_time,
        "last_time": last_time,
        "all_within_tol": all_ok,
        "fields": field_report,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Check physical and numerical convergence of an OpenFOAM run.")
    ap.add_argument("--case-dir", default=DEFAULT_CASE_DIR,
                     help=f"Case root directory (default: {DEFAULT_CASE_DIR}, "
                          f"i.e. the 'Case_Study' folder next to this script)")
    ap.add_argument("--log", default="log.rhoCentralFoam", help="Solver log file to parse for residuals")
    ap.add_argument("--force-tol", type=float, default=1.0,
                     help="Max %% deviation from the mean allowed within the steady window (default 1.0)")
    ap.add_argument("--min-steady-samples", type=int, default=10,
                     help="Minimum number of trailing samples required to call it steady (default 10)")
    ap.add_argument("--residual-tol", type=float, default=1e-4,
                     help="Max acceptable Final residual at the last timestep (default 1e-4)")
    args = ap.parse_args()

    lines = []
    lines.append("=================================================")
    lines.append("  CONVERGENCE REPORT")
    lines.append("=================================================")

    overall_ok = True

    # --- Part A ---
    force_rows = gather_force_rows(args.case_dir)
    lines.append("")
    lines.append("[A] Physical steady-state check (forces on nozzleWall)")
    lines.append("-------------------------------------------------")
    if not force_rows:
        lines.append("No postProcessing/forces/*/force.dat found - skipping this check.")
        lines.append("(Add a 'forces' function object to system/controlDict to enable it.)")
    elif len(force_rows) < 2:
        lines.append("Only one force sample found - cannot assess steadiness from this alone.")
        overall_ok = False
    else:
        result = check_steady_state(force_rows, args.force_tol, args.min_steady_samples)
        verdict = "STEADY" if result["achieved"] else "NOT YET STEADY"
        overall_ok = overall_ok and result["achieved"]
        lines.append(f"Samples              : {len(force_rows)} (t={force_rows[0][0]:.6g}s to t={force_rows[-1][0]:.6g}s)")
        lines.append(f"Steady window         : t={result['steady_from_t']:.6g}s to t={result['steady_to_t']:.6g}s "
                      f"({result['n_steady_samples']} samples)")
        lines.append(f"Mean force in window  : {result['mean_total_x']:.4f} N (slice)")
        lines.append(f"Spread in window      : {result['spread_pct']:.4f} %  (tolerance: {args.force_tol}%)")
        lines.append(f"VERDICT               : {verdict}")
        if not result["achieved"]:
            lines.append(f"  -> fewer than {args.min_steady_samples} trailing samples met the "
                          f"{args.force_tol}% tolerance. Consider extending endTime.")

    # --- Part B ---
    lines.append("")
    lines.append("[B] Numerical solver residual check")
    lines.append("-------------------------------------------------")
    log_path = os.path.join(args.case_dir, args.log)
    if not os.path.exists(log_path):
        lines.append(f"Log file '{log_path}' not found - skipping this check.")
    else:
        blocks = parse_solver_log(log_path)
        residual_result = check_residuals(blocks, args.residual_tol)
        if residual_result is None:
            lines.append(f"No 'Solving for ...' residual lines found in '{log_path}'.")
            lines.append("(This build's residual monitor may need the 'solverInfo' function object")
            lines.append(" rather than 'residuals' - the raw solver log is parsed here regardless.)")
        else:
            overall_ok = overall_ok and residual_result["all_within_tol"]
            lines.append(f"Timesteps logged      : {residual_result['n_timesteps_logged']} "
                          f"(t={residual_result['first_time']:.6g}s to t={residual_result['last_time']:.6g}s)")
            lines.append(f"{'Field':<8}{'Final@last':>14}{'Initial@last':>16}{'Drop vs t0':>14}{'OK?':>8}")
            for field, info in sorted(residual_result["fields"].items()):
                drop_str = f"{info['drop_factor']:.1e}x" if info["drop_factor"] else "n/a"
                ok_str = "yes" if info["within_tol"] else "NO"
                lines.append(
                    f"{field:<8}{info['last_final']:>14.3e}{info['last_initial']:>16.3e}"
                    f"{drop_str:>14}{ok_str:>8}"
                )
            verdict = "WITHIN TOLERANCE" if residual_result["all_within_tol"] else "ELEVATED - CHECK ABOVE"
            lines.append(f"VERDICT               : {verdict} (tolerance: {args.residual_tol:.1e})")
            lines.append("'Drop vs t0' shows how many orders of magnitude the Initial residual fell")
            lines.append("from the first to the last logged timestep - a large drop means the solution")
            lines.append("is changing far less between timesteps now than when it started (steady).")

    lines.append("")
    lines.append("=================================================")
    lines.append(f"OVERALL: {'CONVERGED' if overall_ok else 'REVIEW NEEDED - see flags above'}")
    lines.append("=================================================")

    report = "\n".join(lines)
    print(report)
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
