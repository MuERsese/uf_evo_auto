# UF Evo Auto

`uf_evo_auto.py` is a small helper script for evaluating Ultra-Fusion trajectories with
[`evo`](https://github.com/MichaelGrupp/evo). It runs APE and RPE, automatically chooses a
reasonable timestamp offset, and writes a simple track-lost candidate report.

## Features

- Runs `evo_ape tum` and `evo_rpe tum`.
- Tests several timestamp offset candidates and selects the one with the most matches.
- Automatically sets `t_max_diff` from the result trajectory frequency.
- Writes APE/RPE logs and `.zip` result files.
- Reports possible tracking loss based on coverage, timestamp gaps, position jumps, match ratio,
  and optional ATE max threshold.

## Requirements

```bash
pip install numpy evo
```

Make sure `evo_ape` and `evo_rpe` are available:

```bash
evo_ape --help
evo_rpe --help
```

## Input Format

Both input files must be standard TUM trajectory files:

```text
timestamp tx ty tz qx qy qz qw
```

Example:

```text
1735888008.319916725 0.000084527 -0.000031948 -0.000094470 -0.014919607 0.233888848 -0.000040759 0.972148863
1735888008.419853210 -0.038656056 -0.020440916 0.001053196 -0.014993476 0.234208084 0.000345852 0.972070801
```

The files should already be valid TUM:

- exactly 8 numeric columns
- whitespace-separated values
- quaternion order: `qx qy qz qw`

Rows may be out of timestamp order, and duplicate timestamps are allowed. When either case is
found, the script writes a sorted and deduplicated file in the output directory and passes that
file to evo. For duplicate timestamps, the first pose for that timestamp is kept. The original
input files are not modified.

## Input Preparation

Before running evo, the script loads both trajectories and checks that all rows are valid numeric
TUM rows. If a trajectory has out-of-order or duplicate timestamps, the script creates an evo-ready
copy:

```text
gt_evo_input.txt
result_evo_input.txt
```

These files are written only when needed. The console output reports how many out-of-order
timestamp pairs were found and how many duplicate timestamps were removed.

## Usage

```bash
./uf_evo_auto.py groundtruth.txt result.txt --out_dir evo_results
```

Example:

```bash
./uf_evo_auto.py \
  dataset/Outdoor01/groundtruth.txt \
  output/Outdoor01/LIO/result_LIO.txt \
  --out_dir output/Outdoor01/LIO
```

Show evo plots:

```bash
./uf_evo_auto.py gt.txt result.txt --out_dir evo_results --plot
```

Use RPE over 1 meter, which is the default:

```bash
./uf_evo_auto.py gt.txt result.txt --rpe_delta 1.0 --rpe_delta_unit m
```

Use RPE over 1 second:

```bash
./uf_evo_auto.py gt.txt result.txt --rpe_delta 1.0 --rpe_delta_unit s
```

## Options

```text
--out_dir DIR                 Output directory. Default: evo_results
--plot                        Pass -p to evo_ape and evo_rpe
--rpe_delta VALUE             RPE delta value. Default: 1.0
--rpe_delta_unit UNIT         RPE delta unit: m, s, f, deg, or rad. Default: m
--gap_threshold SECONDS       Tracking-loss time-gap threshold
--jump_threshold METERS       Tracking-loss position-jump threshold. Default: 5.0
--coverage_threshold RATIO    Minimum time-overlap ratio. Default: 0.90
--match_ratio_threshold RATIO Minimum timestamp match ratio. Default: 0.30
--ate_lost_threshold METERS   Optional APE max threshold for possible tracking loss
```

## Output Files

The output directory contains:

```text
ape.zip           evo APE result archive
rpe.zip           evo RPE result archive
ape.txt           evo APE terminal log
rpe.txt           evo RPE terminal log
track_check.txt   automatic track-lost candidate report
gt_evo_input.txt  sorted/deduplicated GT file, only written if needed
result_evo_input.txt
                  sorted/deduplicated result file, only written if needed
```

When `ape.zip` or `rpe.zip` already exists, the script removes the old archive before calling evo.
This avoids evo's interactive overwrite prompt during repeated runs.

Quickly inspect key metrics:

```bash
grep -E 'rmse|max|mean|std' evo_results/ape.txt
grep -E 'rmse|max|mean|std' evo_results/rpe.txt
cat evo_results/track_check.txt
```

## Track-Lost Heuristic

`track_check.txt` is an automatic warning report, not an internal Ultra-Fusion tracking state.
It marks `TRACK_LOST_CANDIDATE: YES` if one or more checks fail:

- trajectory time overlap is lower than `--coverage_threshold`
- result timestamp gap is larger than `--gap_threshold`
- consecutive position jump is larger than `--jump_threshold`
- timestamp match ratio is lower than `--match_ratio_threshold`
- optional APE max is larger than `--ate_lost_threshold`

Read `track_check.txt` first when checking for possible tracking loss:

```text
TRACK_LOST_CANDIDATE: YES
```

means the trajectory should be inspected. The `Reasons` section explains which heuristic triggered
the warning. `NO` means none of the configured heuristic thresholds were crossed.

Useful options:

```bash
--gap_threshold 1.0
--jump_threshold 5.0
--coverage_threshold 0.90
--match_ratio_threshold 0.30
--ate_lost_threshold 20
```

## Converting ROS Path to TUM

If your result is stored as a `nav_msgs/Path` topic in a rosbag, convert it first:

```bash
./path_bag_to_tum.py \
  -i output/Outdoor01/result.bag \
  -t /result_path \
  -o output/Outdoor01/result_path.txt
```

Then run:

```bash
./uf_evo_auto.py gt.txt output/Outdoor01/result_path.txt --out_dir output/Outdoor01/evo
```
