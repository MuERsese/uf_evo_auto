#!/usr/bin/env python3
import argparse
import os
import subprocess
import numpy as np


def load_tum_file(input_path, dedup_output_path=None):
    """
    Load a standard TUM trajectory file.
    Expected format: timestamp tx ty tz qx qy qz qw
    If duplicate timestamps exist, keep the first row and drop later duplicates.
    """
    try:
        data = np.loadtxt(input_path, comments="#")
    except OSError as exc:
        raise RuntimeError(f"Cannot read TUM file: {input_path}") from exc
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid TUM file: {input_path}. "
            "Expected numeric whitespace-separated columns."
        ) from exc

    if data.size == 0:
        raise RuntimeError(f"No TUM rows found in {input_path}")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] != 8:
        raise RuntimeError(
            f"Invalid TUM file: {input_path}. "
            f"Expected exactly 8 columns, got {data.shape[1]}."
        )

    if not np.all(np.isfinite(data)):
        raise RuntimeError(f"Invalid TUM file: {input_path}. Contains NaN or Inf.")

    out_of_order_count = 0
    duplicate_count = 0
    conflict_count = 0
    evo_input_path = input_path

    if data.shape[0] > 1:
        original_dt = np.diff(data[:, 0])
        out_of_order_count = int(np.sum(original_dt < 0))

        if out_of_order_count > 0:
            order = np.argsort(data[:, 0], kind="mergesort")
            data = data[order]

        dt = np.diff(data[:, 0])
        duplicate_idx = np.where(dt == 0)[0] + 1
        duplicate_count = int(len(duplicate_idx))

        if duplicate_count > 0:
            keep = np.ones(data.shape[0], dtype=bool)
            keep[duplicate_idx] = False

            for idx in duplicate_idx:
                if not np.allclose(data[idx, 1:], data[idx - 1, 1:], atol=1e-9, rtol=0.0):
                    conflict_count += 1

            data = data[keep]

        if (out_of_order_count > 0 or duplicate_count > 0) and dedup_output_path is not None:
            np.savetxt(
                dedup_output_path,
                data,
                fmt="%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f"
            )
            evo_input_path = dedup_output_path

    return data, evo_input_path, out_of_order_count, duplicate_count, conflict_count


def median_dt(t):
    dt = np.diff(t)
    dt = dt[dt > 0]
    if len(dt) == 0:
        return 0.05
    return float(np.median(dt))


def count_matches(gt_t, res_t_shifted, threshold):
    """
    Count how many shifted result timestamps have a nearest GT timestamp within threshold.
    """
    idx = np.searchsorted(gt_t, res_t_shifted)

    best = np.full(len(res_t_shifted), np.inf)

    valid = idx < len(gt_t)
    best[valid] = np.minimum(
        best[valid],
        np.abs(gt_t[idx[valid]] - res_t_shifted[valid])
    )

    valid = idx > 0
    best[valid] = np.minimum(
        best[valid],
        np.abs(gt_t[idx[valid] - 1] - res_t_shifted[valid])
    )

    finite = best[np.isfinite(best)]
    median_best = float(np.median(finite)) if len(finite) > 0 else float("inf")

    return int(np.sum(best <= threshold)), median_best


def run_cmd(cmd, log_path):
    print("\nRunning:")
    print(" ".join(cmd))

    with open(log_path, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)

    if p.returncode != 0:
        print(f"[ERROR] command failed, see log: {log_path}")
    else:
        print(f"[OK] log saved: {log_path}")

    return p.returncode


def remove_existing_file(path):
    if os.path.exists(path):
        os.remove(path)
        print(f"[INFO] removed existing file: {path}")


def parse_evo_stats(log_path):
    """
    Parse max/mean/rmse/std and related metrics from an evo terminal log.
    """
    stats = {}
    if not os.path.exists(log_path):
        return stats

    keys = {"max", "mean", "median", "min", "rmse", "sse", "std"}

    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in keys:
                try:
                    stats[parts[0]] = float(parts[1])
                except ValueError:
                    pass

    return stats


def check_track_lost(
    gt,
    res,
    best_offset,
    best_count,
    gap_threshold,
    jump_threshold,
    coverage_threshold,
    match_ratio_threshold,
    ate_stats=None,
    ate_lost_threshold=None,
):
    """
    Estimate possible tracking loss from trajectory continuity and time coverage.
    This is an automatic heuristic, not Ultra-Fusion's internal tracking state.
    """
    gt_t = gt[:, 0]
    res_t = res[:, 0]
    res_xyz = res[:, 1:4]
    res_t_shifted = res_t + best_offset

    gt_duration = gt_t[-1] - gt_t[0]
    res_duration = res_t[-1] - res_t[0]

    overlap_start = max(gt_t[0], res_t_shifted[0])
    overlap_end = min(gt_t[-1], res_t_shifted[-1])
    overlap_duration = max(0.0, overlap_end - overlap_start)

    coverage_ratio = overlap_duration / gt_duration if gt_duration > 0 else 0.0
    match_ratio = best_count / len(res_t) if len(res_t) > 0 else 0.0

    dt = np.diff(res_t)
    dxyz = np.diff(res_xyz, axis=0)
    step_dist = np.linalg.norm(dxyz, axis=1)

    large_gap_idx = np.where(dt > gap_threshold)[0]
    large_jump_idx = np.where(step_dist > jump_threshold)[0]

    reasons = []

    if coverage_ratio < coverage_threshold:
        reasons.append(
            f"coverage_ratio {coverage_ratio:.3f} < threshold {coverage_threshold:.3f}"
        )

    if len(large_gap_idx) > 0:
        reasons.append(
            f"large time gaps: {len(large_gap_idx)} gaps > {gap_threshold:.3f}s"
        )

    if len(large_jump_idx) > 0:
        reasons.append(
            f"large position jumps: {len(large_jump_idx)} jumps > {jump_threshold:.3f}m"
        )

    if match_ratio < match_ratio_threshold:
        reasons.append(
            f"match_ratio {match_ratio:.3f} < threshold {match_ratio_threshold:.3f}"
        )

    if ate_stats is not None and ate_lost_threshold is not None:
        if "max" in ate_stats and ate_stats["max"] > ate_lost_threshold:
            reasons.append(
                f"ATE max {ate_stats['max']:.3f}m > threshold {ate_lost_threshold:.3f}m"
            )

    lost_candidate = len(reasons) > 0

    report = {
        "lost_candidate": lost_candidate,
        "reasons": reasons,
        "gt_duration": gt_duration,
        "res_duration": res_duration,
        "overlap_duration": overlap_duration,
        "coverage_ratio": coverage_ratio,
        "match_ratio": match_ratio,
        "gap_threshold": gap_threshold,
        "jump_threshold": jump_threshold,
        "num_large_gaps": len(large_gap_idx),
        "num_large_jumps": len(large_jump_idx),
        "max_time_gap": float(np.max(dt)) if len(dt) > 0 else 0.0,
        "mean_time_gap": float(np.mean(dt)) if len(dt) > 0 else 0.0,
        "max_step_distance": float(np.max(step_dist)) if len(step_dist) > 0 else 0.0,
        "mean_step_distance": float(np.mean(step_dist)) if len(step_dist) > 0 else 0.0,
        "large_gap_idx": large_gap_idx,
        "large_jump_idx": large_jump_idx,
        "res_t": res_t,
        "dt": dt,
        "step_dist": step_dist,
    }

    return report


def write_track_report(report, out_path):
    with open(out_path, "w") as f:
        f.write("===== Track Lost Check =====\n")
        f.write(f"TRACK_LOST_CANDIDATE: {'YES' if report['lost_candidate'] else 'NO'}\n\n")

        f.write("===== Reasons =====\n")
        if report["reasons"]:
            for r in report["reasons"]:
                f.write(f"- {r}\n")
        else:
            f.write("- None\n")

        f.write("\n===== Metrics =====\n")
        f.write(f"gt_duration       : {report['gt_duration']:.6f} s\n")
        f.write(f"result_duration   : {report['res_duration']:.6f} s\n")
        f.write(f"overlap_duration  : {report['overlap_duration']:.6f} s\n")
        f.write(f"coverage_ratio    : {report['coverage_ratio']:.6f}\n")
        f.write(f"match_ratio       : {report['match_ratio']:.6f}\n")
        f.write(f"max_time_gap      : {report['max_time_gap']:.6f} s\n")
        f.write(f"mean_time_gap     : {report['mean_time_gap']:.6f} s\n")
        f.write(f"max_step_distance : {report['max_step_distance']:.6f} m\n")
        f.write(f"mean_step_distance: {report['mean_step_distance']:.6f} m\n")
        f.write(f"num_large_gaps    : {report['num_large_gaps']}\n")
        f.write(f"num_large_jumps   : {report['num_large_jumps']}\n")
        f.write(f"gap_threshold     : {report['gap_threshold']:.6f} s\n")
        f.write(f"jump_threshold    : {report['jump_threshold']:.6f} m\n")

        f.write("\n===== Large Gaps, first 10 =====\n")
        for i in report["large_gap_idx"][:10]:
            t0 = report["res_t"][i]
            t1 = report["res_t"][i + 1]
            gap = report["dt"][i]
            f.write(f"idx {i}: {t0:.9f} -> {t1:.9f}, gap = {gap:.6f} s\n")

        f.write("\n===== Large Jumps, first 10 =====\n")
        for i in report["large_jump_idx"][:10]:
            t0 = report["res_t"][i]
            t1 = report["res_t"][i + 1]
            dist = report["step_dist"][i]
            f.write(f"idx {i}: {t0:.9f} -> {t1:.9f}, jump = {dist:.6f} m\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run evo APE/RPE for TUM trajectories and generate a simple "
            "tracking-loss candidate report."
        )
    )
    parser.add_argument("gt", help="GT TUM txt")
    parser.add_argument("result", help="Ultra-Fusion result TUM txt")
    parser.add_argument("--out_dir", default="evo_results")
    parser.add_argument("--plot", action="store_true", help="show evo plots")
    parser.add_argument("--rpe_delta", type=float, default=1.0)
    parser.add_argument("--rpe_delta_unit", default="m", choices=["m", "s", "f", "deg", "rad"])

    # Tracking-loss heuristic thresholds.
    parser.add_argument(
        "--gap_threshold",
        type=float,
        default=None,
        help="If result pose time gap is larger than this, mark as possible lost. Default: max(1.0, 5 * result_median_dt)"
    )
    parser.add_argument(
        "--jump_threshold",
        type=float,
        default=5.0,
        help="If consecutive pose jump is larger than this in meters, mark as possible lost. Default: 5.0m"
    )
    parser.add_argument(
        "--coverage_threshold",
        type=float,
        default=0.90,
        help="If result/GT time overlap ratio is lower than this, mark as possible lost. Default: 0.90"
    )
    parser.add_argument(
        "--match_ratio_threshold",
        type=float,
        default=0.30,
        help="If matched result pose ratio is lower than this, mark as possible lost. Default: 0.30"
    )
    parser.add_argument(
        "--ate_lost_threshold",
        type=float,
        default=None,
        help="Optional. If ATE max is larger than this in meters, mark as possible lost. Example: --ate_lost_threshold 20"
    )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    gt_dedup_path = os.path.join(args.out_dir, "gt_evo_input.txt")
    res_dedup_path = os.path.join(args.out_dir, "result_evo_input.txt")

    gt, gt_evo_input, gt_out_of_order, gt_duplicates, gt_conflicts = load_tum_file(
        args.gt,
        dedup_output_path=gt_dedup_path,
    )
    res, res_evo_input, res_out_of_order, res_duplicates, res_conflicts = load_tum_file(
        args.result,
        dedup_output_path=res_dedup_path,
    )

    if (
        gt_out_of_order > 0 or res_out_of_order > 0
        or gt_duplicates > 0 or res_duplicates > 0
    ):
        print("\n===== TUM Input Preparation =====")
        print(f"GT out-of-order timestamp pairs    : {gt_out_of_order}")
        print(f"Result out-of-order timestamp pairs: {res_out_of_order}")
        print(f"GT duplicate timestamps removed    : {gt_duplicates}")
        print(f"Result duplicate timestamps removed: {res_duplicates}")
        if gt_out_of_order > 0 or gt_duplicates > 0:
            print(f"GT evo input                       : {gt_evo_input}")
        if res_out_of_order > 0 or res_duplicates > 0:
            print(f"Result evo input                   : {res_evo_input}")
        if gt_conflicts > 0 or res_conflicts > 0:
            print("[WARNING] Duplicate timestamps with different poses were found.")
            print(f"GT conflicting duplicates          : {gt_conflicts}")
            print(f"Result conflicting duplicates      : {res_conflicts}")
            print("The first pose for each duplicated timestamp was kept.")

    gt_t = gt[:, 0]
    res_t = res[:, 0]

    gt_dt = median_dt(gt_t)
    res_dt = median_dt(res_t)

    # Set the timestamp matching threshold from the result trajectory frequency.
    # For example, for a 10 Hz result trajectory, median dt is about 0.1 s and
    # t_max_diff becomes about 0.05 s.
    t_max_diff = max(0.02, min(0.10, 0.5 * res_dt + 1e-6))

    # Time-gap threshold used by the tracking-loss heuristic.
    # Default: at least 1 second, or 5 times the result trajectory period.
    gap_threshold = args.gap_threshold
    if gap_threshold is None:
        gap_threshold = max(1.0, 5.0 * res_dt)

    # Timestamp offset candidates:
    # 1. 0: use this when both trajectories already share the same clock.
    # 2. gt_start - res_start: handle a global start-time shift.
    # 3. gt_mid - res_mid: more stable when initialization makes the start noisy.
    candidates = [
        0.0,
        gt_t[0] - res_t[0],
        np.median(gt_t) - np.median(res_t),
    ]

    best_offset = 0.0
    best_count = -1
    best_median_diff = float("inf")

    print("\n===== Time Info =====")
    print(f"GT     start/end/duration: {gt_t[0]:.9f} / {gt_t[-1]:.9f} / {gt_t[-1] - gt_t[0]:.3f}s")
    print(f"Result start/end/duration: {res_t[0]:.9f} / {res_t[-1]:.9f} / {res_t[-1] - res_t[0]:.3f}s")
    print(f"GT median dt    : {gt_dt:.6f}s")
    print(f"Result median dt: {res_dt:.6f}s")
    print(f"Auto t_max_diff : {t_max_diff:.6f}s")
    print(f"Gap threshold   : {gap_threshold:.6f}s")
    print(f"Jump threshold  : {args.jump_threshold:.6f}m")

    print("\n===== Offset Candidates =====")
    for off in candidates:
        count, med = count_matches(gt_t, res_t + off, t_max_diff)
        print(f"offset {off:+.9f}s -> matches {count}/{len(res_t)}, median nearest diff {med:.6f}s")

        # Selection rule:
        # 1. Prefer the offset with more timestamp matches.
        # 2. If match counts are within 5%, prefer the smaller median time error.
        if count > best_count:
            best_count = count
            best_offset = float(off)
            best_median_diff = med
        elif best_count > 0 and count >= 0.95 * best_count and med < best_median_diff:
            best_count = count
            best_offset = float(off)
            best_median_diff = med

    print("\n===== Selected =====")
    print(f"selected t_offset : {best_offset:+.9f}s")
    print(f"selected matches  : {best_count}/{len(res_t)}")
    print(f"median time diff  : {best_median_diff:.6f}s")

    if best_count == 0:
        print("\n[WARNING] No timestamp matches found even after auto offset.")
        print("Please check whether GT and result cover the same time range.")
        return

    ape_zip = os.path.join(args.out_dir, "ape.zip")
    rpe_zip = os.path.join(args.out_dir, "rpe.zip")
    ape_log = os.path.join(args.out_dir, "ape.txt")
    rpe_log = os.path.join(args.out_dir, "rpe.txt")
    track_log = os.path.join(args.out_dir, "track_check.txt")

    remove_existing_file(ape_zip)
    remove_existing_file(rpe_zip)

    ape_cmd = [
        "evo_ape", "tum", gt_evo_input, res_evo_input,
        "-a",
        "--t_offset", f"{best_offset:.9f}",
        "--t_max_diff", f"{t_max_diff:.6f}",
        "--save_results", ape_zip,
    ]

    rpe_cmd = [
        "evo_rpe", "tum", gt_evo_input, res_evo_input,
        "-a",
        "--t_offset", f"{best_offset:.9f}",
        "--t_max_diff", f"{t_max_diff:.6f}",
        "--delta", str(args.rpe_delta),
        "--delta_unit", args.rpe_delta_unit,
        "--save_results", rpe_zip,
    ]

    if args.plot:
        ape_cmd.append("-p")
        rpe_cmd.append("-p")

    run_cmd(ape_cmd, ape_log)
    run_cmd(rpe_cmd, rpe_log)

    ape_stats = parse_evo_stats(ape_log)

    track_report = check_track_lost(
        gt=gt,
        res=res,
        best_offset=best_offset,
        best_count=best_count,
        gap_threshold=gap_threshold,
        jump_threshold=args.jump_threshold,
        coverage_threshold=args.coverage_threshold,
        match_ratio_threshold=args.match_ratio_threshold,
        ate_stats=ape_stats,
        ate_lost_threshold=args.ate_lost_threshold,
    )

    write_track_report(track_report, track_log)

    print("\n===== Track Lost Check =====")
    print(f"TRACK_LOST_CANDIDATE: {'YES' if track_report['lost_candidate'] else 'NO'}")
    if track_report["reasons"]:
        for r in track_report["reasons"]:
            print(f"- {r}")
    else:
        print("- None")

    print("\n===== Results Files =====")
    print(f"APE zip     : {ape_zip}")
    print(f"RPE zip     : {rpe_zip}")
    print(f"APE log     : {ape_log}")
    print(f"RPE log     : {rpe_log}")
    print(f"Track check : {track_log}")

    print("\nCheck key numbers:")
    print(f"grep -E 'rmse|max|mean|std' {ape_log}")
    print(f"grep -E 'rmse|max|mean|std' {rpe_log}")
    print(f"cat {track_log}")


if __name__ == "__main__":
    main()
