#!/usr/bin/env python3
"""Sensitivity analysis of exported RNAscope detections, without changing profiles.

Reuses accepted nuclei and all detected peaks from a complete analysis run.
Sweeps counting-region radius and independent Pnoc/Cre positivity cutoffs;
does not optimize detection accuracy or re-estimate segmentation/background.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path

import numpy as np
import tifffile

import rnascope as rn


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def read_scene(scene, radii, baseline_radius):
    """Reassign original peaks; verify the saved baseline before sweeping."""
    metadata = json.loads((scene / "metadata.json").read_text())
    spacing = np.asarray(metadata["pixel_size_yx_um"], dtype=float)
    if spacing.shape != (2,) or not np.all(np.isfinite(spacing) & (spacing > 0)):
        raise ValueError(f"Invalid physical calibration: {scene}")
    labels = tifffile.imread(scene / "nuclei_labels.tif")
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError(f"Expected a two-dimensional integer label mask: {scene}")
    rows = read_csv(scene / "nuclei.csv")
    nuclei = [{"nucleus_id": int(r["nucleus_id"]), "accepted": r["accepted"] == "True"} for r in rows]
    ids = np.array([n["nucleus_id"] for n in nuclei], dtype=int)
    if len(set(ids)) != len(ids) or set(ids) != set(np.unique(labels)) - {0}:
        raise ValueError(f"Nucleus IDs do not match the label mask: {scene}")
    accepted = np.array([n["accepted"] for n in nuclei], dtype=bool)
    spots = read_csv(scene / "spots.csv")
    if any(s["gene"] not in ("pnoc", "cre") for s in spots):
        raise ValueError(f"Unexpected gene in spot table: {scene}")
    coordinates = {}
    for gene in ("pnoc", "cre"):
        coordinates[gene] = np.array([[int(s["y_px"]), int(s["x_px"])]
                                      for s in spots if s["gene"] == gene], dtype=int).reshape(-1, 2)
        xy = coordinates[gene]
        if np.any(xy < 0) or np.any(xy >= np.array(labels.shape)):
            raise ValueError(f"Spot coordinate outside the label mask: {scene}")

    counts = {}
    areas = {}
    for radius in sorted(set([0.0, baseline_radius, *radii])):
        assignments = rn.build_assignment_labels(labels, nuclei, spacing, radius)
        counts[radius] = {}
        for gene, xy in coordinates.items():
            values = np.bincount(assignments[xy[:, 0], xy[:, 1]], minlength=int(labels.max()) + 1)
            counts[radius][gene] = values[ids]
        area = np.bincount(assignments.ravel(), minlength=int(labels.max()) + 1)[ids]
        areas[radius] = float(area[accepted].sum() * np.prod(spacing))
    for gene in ("pnoc", "cre"):
        nuclear = np.array([int(r[f"{gene}_nuclear_spots"]) for r in rows])
        saved = np.array([int(r[f"{gene}_spots"]) for r in rows])
        if not np.array_equal(counts[0.0][gene], nuclear):
            raise ValueError(f"Nuclear {gene} counts do not reproduce: {scene}")
        # Expanded runs intentionally discard rejected objects from assignments.
        if not np.array_equal(counts[baseline_radius][gene][accepted], saved[accepted]):
            raise ValueError(f"Baseline assigned {gene} counts do not reproduce: {scene}")
    return ids[accepted], {r: {g: v[accepted] for g, v in d.items()} for r, d in counts.items()}, areas


def count_setting(counts, baseline_cre, pnoc_cutoff, cre_cutoff):
    p = counts["pnoc"] >= pnoc_cutoff
    c = counts["cre"] >= cre_cutoff
    return {"total_nuclei": len(p), "pnoc_positive": int(p.sum()), "cre_positive": int(c.sum()),
            "double_positive": int((p & c).sum()), "pnoc_only": int((p & ~c).sum()),
            "cre_only": int((~p & c).sum()), "neither": int((~p & ~c).sum()),
            "baseline_cre_positive": int(baseline_cre.sum()),
            "baseline_cre_pnoc_positive": int((baseline_cre & p).sum()),
            "cre_retained_vs_baseline": int((baseline_cre & c).sum()),
            "cre_added_vs_baseline": int((~baseline_cre & c).sum()),
            "cre_lost_vs_baseline": int((baseline_cre & ~c).sum())}


def add_percentages(row):
    for name, numerator, denominator in [
            ("pct_cre_expressing_pnoc", "double_positive", "cre_positive"),
            ("pct_pnoc_expressing_cre", "double_positive", "pnoc_positive"),
            ("pct_baseline_cre_expressing_pnoc", "baseline_cre_pnoc_positive", "baseline_cre_positive")]:
        row[name] = 100 * row[numerator] / row[denominator] if row[denominator] else None
    return row


def plot_sweep(rows, radii, cutoffs, baseline, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    lookup = {(r["expand_um"], r["pnoc_min_spots"], r["cre_min_spots"]): r for r in rows}
    columns = min(3, len(radii)); nrows = math.ceil(len(radii) / columns)
    fig, axes = plt.subplots(nrows, columns, figsize=(5.4 * columns, 4.8 * nrows), squeeze=False,
                             layout="constrained")
    for ax, radius in zip(axes.flat, radii):
        values = [[lookup[radius, p, c]["pct_cre_expressing_pnoc"] for c in cutoffs] for p in cutoffs]
        values = np.array([[np.nan if v is None else v for v in row] for row in values])
        heat = ax.imshow(values, vmin=0, vmax=100, cmap="viridis", aspect="auto")
        for i, p in enumerate(cutoffs):
            for j, c in enumerate(cutoffs):
                row = lookup[radius, p, c]
                pct = row["pct_cre_expressing_pnoc"]
                label = f"{pct:.1f}%\n{row['double_positive']}/{row['cre_positive']}" if pct is not None else "NA\n0/0"
                ax.text(j, i, label, ha="center", va="center", fontsize=7,
                        color="black" if pct is not None and pct > 65 else "white")
        if radius == baseline[0] and baseline[1] in cutoffs and baseline[2] in cutoffs:
            ax.add_patch(Rectangle((cutoffs.index(baseline[2]) - .5, cutoffs.index(baseline[1]) - .5),
                                   1, 1, fill=False, edgecolor="#ff3333", linewidth=3))
        ax.set(xticks=range(len(cutoffs)), xticklabels=cutoffs, yticks=range(len(cutoffs)), yticklabels=cutoffs,
               xlabel="Minimum Cre peaks", ylabel="Minimum Pnoc peaks", title=f"Nucleus + {radius:g} µm")
    for ax in list(axes.flat)[len(radii):]:
        ax.set_visible(False)
    fig.colorbar(heat, ax=list(axes.flat), label="Detected Cre-positive objects also Pnoc-positive (%)", shrink=.7)
    fig.suptitle("Sensitivity analysis: fixed nuclei and peak detections\nEach box: percentage and double-positive / Cre-positive count; red = saved baseline")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"threshold_heatmaps.{suffix}", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), layout="constrained")
    for cutoff in cutoffs:
        selected = [lookup[r, cutoff, cutoff] for r in radii]
        axes[0].plot(radii, [r["pct_cre_expressing_pnoc"] for r in selected], marker="o", label=f"Both ≥{cutoff}")
        axes[1].plot(radii, [r["cre_positive"] for r in selected], marker="o", label=f"Cre ≥{cutoff}")
    p, c = baseline[1:]
    if p in cutoffs and c in cutoffs:
        selected = [lookup[r, p, c] for r in radii]
        axes[2].plot(radii, [r["pct_cre_expressing_pnoc"] for r in selected], marker="o", label="Reclassify Cre at each radius")
        axes[2].plot(radii, [r["pct_baseline_cre_expressing_pnoc"] for r in selected], marker="s", label="Keep saved Cre-positive group")
    for ax in axes:
        ax.set_xlabel("Expansion around nucleus (µm)"); ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set(ylabel="Double-positive / Cre-positive (%)", title="Both cutoffs changed together", ylim=(0, 105))
    axes[1].set(ylabel="Detected Cre-positive objects", title="Changing the denominator")
    axes[2].set(ylabel="Cre-positive with detectable Pnoc (%)", title=f"Pnoc ≥{p}; Cre ≥{c}; fixed-group comparison", ylim=(0, 105))
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"radius_sensitivity.{suffix}", dpi=180)
    plt.close(fig)


def run_sweep(results, output, radii, cutoffs):
    results, output = Path(results).resolve(), Path(output).resolve()
    radii, cutoffs = sorted(set(radii)), sorted(set(cutoffs))
    if not radii or any(not math.isfinite(r) or r < 0 for r in radii):
        raise ValueError("Radii must be finite and nonnegative")
    if not cutoffs or any(isinstance(c, bool) or not isinstance(c, int) or c < 1 for c in cutoffs):
        raise ValueError("Peak cutoffs must be positive integers; zero would count signal-free objects")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output folder is not empty; choose a new folder")
    batch = json.loads((results / "batch_summary.json").read_text())
    if not batch["complete"]:
        raise ValueError("Sweep requires a complete analysis run; resolve its failures first")
    config = json.loads((results / "config_used.json").read_text())
    source_manifest = json.loads((results / "run_manifest.json").read_text())
    if rn.sha256(results / "config_used.json") != source_manifest["config_sha256"]:
        raise ValueError("Saved configuration does not match the run manifest")
    baseline = (float(config["assignment"]["expand_um"]), config["spots"]["pnoc"]["positive_min_spots"],
                config["spots"]["cre"]["positive_min_spots"])
    summaries = read_csv(results / "summary.csv")
    if len(summaries) != batch["successful_scenes"] or any(r["status"] != "ok" for r in summaries):
        raise ValueError("Scene summary and batch completeness disagree")
    all_rows, cells, hashes = [], [], {}
    baseline_totals = {"total_nuclei": 0, "pnoc_positive": 0, "cre_positive": 0, "double_positive": 0}
    for summary in summaries:
        scene = (results / "images" / summary["file"] / f"scene_{summary['scene']}").resolve()
        if not scene.is_relative_to(results / "images"):
            raise ValueError("Scene path escapes the analysis results")
        ids, counts, areas = read_scene(scene, radii, baseline[0])
        baseline_cre = counts[baseline[0]]["cre"] >= baseline[2]
        base_counts = count_setting(counts[baseline[0]], baseline_cre, *baseline[1:])
        for key in baseline_totals:
            if base_counts[key] != int(summary[key]):
                raise ValueError(f"Baseline summary does not reproduce: {scene} / {key}")
            baseline_totals[key] += base_counts[key]
        for name in ("metadata.json", "nuclei.csv", "spots.csv", "nuclei_labels.tif"):
            hashes[(scene / name).relative_to(results).as_posix()] = rn.sha256(scene / name)
        for radius in radii:
            for idx, nucleus_id in enumerate(ids):
                cells.append({"file": summary["file"], "scene": summary["scene"], "nucleus_id": int(nucleus_id),
                              "expand_um": radius, "pnoc_spots": int(counts[radius]["pnoc"][idx]),
                              "cre_spots": int(counts[radius]["cre"][idx]), "baseline_cre_positive": bool(baseline_cre[idx])})
            for p, c in itertools.product(cutoffs, repeat=2):
                row = {"file": summary["file"], "scene": summary["scene"], "expand_um": radius,
                       "pnoc_min_spots": p, "cre_min_spots": c,
                       **count_setting(counts[radius], baseline_cre, p, c), "assigned_area_um2": areas[radius]}
                all_rows.append(add_percentages(row))
    if not all_rows:
        raise ValueError("No successful scenes to sweep")
    for key, value in baseline_totals.items():
        if value != batch[key]:
            raise ValueError(f"Baseline batch count does not reproduce: {key}")
    totals = []
    keys = [k for k in all_rows[0] if k not in ("file", "scene", "expand_um", "pnoc_min_spots", "cre_min_spots") and not k.startswith("pct_")]
    for radius, p, c in itertools.product(radii, cutoffs, cutoffs):
        selected = [r for r in all_rows if (r["expand_um"], r["pnoc_min_spots"], r["cre_min_spots"]) == (radius, p, c)]
        totals.append(add_percentages({"expand_um": radius, "pnoc_min_spots": p, "cre_min_spots": c,
                                      **{k: sum(r[k] for r in selected) for k in keys}}))
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in [("sweep.csv", totals), ("per_scene_sweep.csv", all_rows)]:
        rn.write_csv(output / name, rows, list(rows[0]))
    rn.write_csv(output / "counts_by_nucleus.csv", cells,
                 ["file", "scene", "nucleus_id", "expand_um", "pnoc_spots", "cre_spots", "baseline_cre_positive"])
    rn.write_json(output / "sweep_manifest.json", {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "source_results": str(results),
        "source_analysis_manifest": source_manifest, "sweep_script_sha256": rn.sha256(__file__),
        "assignment_helper_source_sha256": rn.sha256(rn.__file__), "source_export_sha256": hashes,
        "source_summary_sha256": rn.sha256(results / "summary.csv"),
        "packages": {p: rn.package_version(p) for p in ("numpy", "scikit-image", "tifffile", "matplotlib")},
        "radii_um": radii, "pnoc_and_cre_cutoffs": cutoffs, "baseline_radius_pnoc_cre": baseline,
        "baseline_counts": baseline_totals, "settings_tested": len(totals), "scenes": len(summaries),
        "interpretation": "Exploratory sensitivity analysis. Nuclei/peaks fixed; no accuracy optimum, negative-control error rate or independent biological replication established. Fixed-group endpoint retains saved baseline Cre-positive objects regardless of current radius/Cre cutoff."})
    plot_sweep(totals, radii, cutoffs, baseline, output)
    return totals


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Complete output folder from rnascope.py")
    parser.add_argument("--output", type=Path, required=True, help="New/empty sweep output folder")
    parser.add_argument("--radii", nargs="+", type=float, default=[0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
    parser.add_argument("--cutoffs", nargs="+", type=int, default=list(range(1, 7)))
    args = parser.parse_args(argv)
    try:
        rows = run_sweep(args.results, args.output, args.radii, args.cutoffs)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"Sweep failed: {exc}\n")
    print(f"Wrote {len(rows)} settings to {args.output}; analysis profiles unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
