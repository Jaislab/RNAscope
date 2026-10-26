#!/usr/bin/env python3
"""Compare frozen RNAscope detections under StarDist and Cellpose masks.

This is an exploratory comparison and blinded review generator. StarDist is a
reference method, not ground truth; no biological accuracy is inferred from IoU
or a larger coexpression percentage. Pnoc-guided cell outlines are estimates.
"""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from skimage import segmentation
import tifffile

import rnascope as rn


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def overlap_table(reference, query):
    if reference.shape != query.shape or reference.ndim != 2:
        raise ValueError("Instance masks must be two-dimensional with identical shapes")
    if any(not np.issubdtype(x.dtype, np.integer) or np.any(x < 0) for x in (reference, query)):
        raise ValueError("Instance labels must be nonnegative integers")
    a, b = (np.union1d([0], np.unique(x)) for x in (reference, query))
    ai, bi = np.searchsorted(a, reference.ravel()), np.searchsorted(b, query.ravel())
    table = np.bincount(ai * len(b) + bi, minlength=len(a) * len(b)).reshape(len(a), len(b))
    return a[1:], b[1:], table[1:, 1:], table[1:, :].sum(axis=1), table[:, 1:].sum(axis=0)


def match_instances(reference, query, threshold=.5):
    """Maximum-cardinality one-to-one matching, then maximize IoU among ties."""
    if not 0 < threshold <= 1:
        raise ValueError("IoU threshold must lie in (0, 1]")
    ref_ids, query_ids, inter, ref_area, query_area = overlap_table(reference, query)
    if not len(ref_ids) or not len(query_ids):
        return [], {"reference_objects": len(ref_ids), "query_objects": len(query_ids),
                    "matched_at_iou_threshold": 0, "reference_unmatched": len(ref_ids),
                    "query_unmatched": len(query_ids), "median_matched_iou": None}
    iou = inter / (ref_area[:, None] + query_area[None, :] - inter)
    eligible = iou >= threshold
    rows, cols = linear_sum_assignment(eligible * (min(iou.shape) + 1) + iou, maximize=True)
    pairs = [{"reference_id": int(ref_ids[i]), "query_id": int(query_ids[j]), "iou": float(iou[i, j]),
              "reference_coverage": float(inter[i, j] / ref_area[i]),
              "query_coverage": float(inter[i, j] / query_area[j])}
             for i, j in zip(rows, cols) if eligible[i, j]]
    return pairs, {"reference_objects": len(ref_ids), "query_objects": len(query_ids),
                   "matched_at_iou_threshold": len(pairs), "reference_unmatched": len(ref_ids) - len(pairs),
                   "query_unmatched": len(query_ids) - len(pairs),
                   "median_matched_iou": float(np.median([p["iou"] for p in pairs])) if pairs else None,
                   "reference_with_multiple_20pct_overlaps": int(((inter / ref_area[:, None] >= .2).sum(axis=1) > 1).sum()),
                   "query_with_multiple_20pct_overlaps": int(((inter / query_area[None, :] >= .2).sum(axis=0) > 1).sum())}


def accepted_labels(labels, nuclei):
    accepted = np.zeros(max(int(labels.max()), max((n["nucleus_id"] for n in nuclei), default=0)) + 1, bool)
    for n in nuclei:
        accepted[n["nucleus_id"]] = n["accepted"]
    return np.where(accepted[labels], labels, 0)


def link_cell_masks(nuclear_labels, nuclei, cell_labels, min_coverage=.5):
    """Require one cell outline per accepted nucleus; report every rejected link.

    Any candidate nucleus with >=min_coverage of its area in a cell counts as a
    nuclear occupant, including rejected reference objects. Border cells and
    non-unique links are excluded; other nuclear pixels are blocked afterward.
    """
    if not 0 < min_coverage <= 1:
        raise ValueError("Nuclear coverage must lie in (0, 1]")
    ref_ids, cell_ids, inter, ref_area, cell_area = overlap_table(nuclear_labels, cell_labels)
    coverage = inter / ref_area[:, None] if len(ref_ids) else np.zeros((0, len(cell_ids)))
    lookup = {n["nucleus_id"]: n for n in nuclei}
    border = set(np.unique(np.concatenate([cell_labels[0], cell_labels[-1], cell_labels[:, 0], cell_labels[:, -1]])))
    links = []
    for j, cid in enumerate(cell_ids):
        occupants = np.flatnonzero(coverage[:, j] >= min_coverage)
        anchor = int(ref_ids[occupants[0]]) if len(occupants) == 1 else 0
        if not len(occupants):
            status = "no_majority_nucleus"
        elif len(occupants) > 1:
            status = "multiple_nuclear_occupants"
        elif not lookup[anchor]["accepted"]:
            status = "rejected_reference_nucleus"
        elif int(cid) in border:
            status = "image_border_cell"
        else:
            status = "eligible"
        links.append({"cell_id": int(cid), "reference_nucleus_id": anchor, "status": status,
                      "majority_nuclear_occupants": len(occupants),
                      "nuclear_coverage": float(coverage[occupants[0], j]) if anchor else None,
                      "cell_area_px": int(cell_area[j]),
                      "other_nuclear_overlap_px": int(inter[:, j].sum() - (inter[occupants[0], j] if anchor else 0))})
    anchors = [r["reference_nucleus_id"] for r in links if r["status"] == "eligible"]
    for r in links:
        if r["status"] == "eligible" and anchors.count(r["reference_nucleus_id"]) > 1:
            r["status"] = "nucleus_split_across_cells"
    mapping = np.zeros(int(cell_labels.max()) + 1, np.uint32)
    for r in links:
        if r["status"] == "eligible":
            mapping[r["cell_id"]] = r["reference_nucleus_id"]
    assignments = mapping[cell_labels]
    assignments[(nuclear_labels > 0) & (nuclear_labels != assignments)] = 0
    return assignments, links


def load_reference(scene, image):
    metadata = json.loads((scene / "metadata.json").read_text())
    nuclei = []
    for row in read_csv(scene / "nuclei.csv"):
        n = dict(row)
        for field in ("nucleus_id", "pnoc_spots", "cre_spots", "pnoc_nuclear_spots", "cre_nuclear_spots"):
            n[field] = int(n[field])
        for field in ("accepted", "pnoc_positive", "cre_positive"):
            n[field] = n[field] == "True"
        for field in ("x_px", "y_px", "x_um", "y_um", "area_um2", "solidity", "eccentricity"):
            n[field] = float(n[field])
        nuclei.append(n)
    spots = {"pnoc": [], "cre": []}
    for row in read_csv(scene / "spots.csv"):
        for field in ("spot_id", "x_px", "y_px", "nucleus_id", "assignment_id"):
            row[field] = int(row[field])
        spots[row["gene"]].append(row)
    labels = tifffile.imread(scene / "nuclei_labels.tif")
    if image.shape[1:] != labels.shape:
        raise ValueError("Input image and reference mask shapes disagree")
    return metadata, nuclei, spots, labels


def analyze_masks(labels, nuclei, assignments, spots, image, metadata, config, output, name):
    output.mkdir(parents=True, exist_ok=False)
    nuclei, spots = copy.deepcopy(nuclei), copy.deepcopy(spots)
    for n in nuclei:
        # Recompute assignment-dependent flags for this method, preserving only
        # the nucleus-shape reason inherited from the selected nuclear mask.
        n["qc_flags"] = ";".join(f for f in n.get("qc_flags", "").split(";") if f == "irregular_nuclear_shape")
    masks = {}
    for gene, c in metadata["resolved_channels"].items():
        bits = metadata["channels"][c].get("component_bit_count")
        saturation = rn.saturation_level(image[c], config["qc"]["saturation_value"], bits)
        if saturation is not None:
            masks[f"{gene}_saturation"] = image[c] >= saturation
    rn.assign_and_classify(labels, nuclei, spots, config["spots"], masks, assignments)
    summary = {"method": name, **rn.summarize(nuclei)}
    for n in nuclei:
        n.setdefault("file", "reference.czi"); n.setdefault("scene", "0")
    rn.write_csv(output / "nuclei.csv", nuclei, rn.NUCLEUS_FIELDS)
    rn.write_json(output / "summary.json", summary)
    spacing = metadata["pixel_size_yx_um"]
    for filename, values in [("nuclei_labels", labels), ("assignment_labels", assignments)]:
        tifffile.imwrite(output / f"{filename}.tif", values.astype(np.uint32), photometric="minisblack",
                         resolution=(10000 / spacing[1], 10000 / spacing[0]), resolutionunit="CENTIMETER",
                         metadata={"axes": "YX", "unit": "um"})
        rn.write_imagej_labels(output / f"{filename}_imagej.tif", values, spacing)
    return summary, nuclei


def run_comparison(reference, input_path, inference, output):
    reference, inference, output = map(Path, (reference, inference, output))
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output folder must be new or empty")
    root = rn.analysis_root(reference)
    config_path = root / "config_used.json"
    config = json.loads(config_path.read_text())
    reference_manifest = json.loads((root / "run_manifest.json").read_text())
    if rn.sha256(config_path) != reference_manifest["config_sha256"] or config["assignment"]["expand_um"] != 2:
        raise ValueError("Use a verified reference run with the current 2 um counting region")
    if any(config["spots"][gene]["positive_min_spots"] != 3 for gene in ("pnoc", "cre")):
        raise ValueError("This comparison preserves the requested >=3 cutoffs for both genes")
    with np.load(input_path, allow_pickle=False) as archive:
        image = archive["image"]
    metadata, nuclei, spots, labels = load_reference(reference, image)
    spacing = metadata["pixel_size_yx_um"]
    baseline = json.loads((reference / "summary.json").read_text())
    inference_manifest = json.loads((inference / "inference_manifest.json").read_text())
    if inference_manifest["input_sha256"] != rn.sha256(input_path):
        raise ValueError("Cellpose input does not match its inference manifest")
    source_hashes = [r["sha256"] for r in reference_manifest["inputs"] if r["file"] == baseline["file"]]
    if source_hashes != [inference_manifest["input_metadata"]["czi_sha256"]]:
        raise ValueError("Reference and Cellpose did not use the same source CZI")
    variants = ["nuclei_dapi", "cyto3_dapi", "cyto3_pnoc_dapi"]
    variants.extend(v for v in ("nuclei_dapi_clahe", "cyto3_dapi_clahe") if v in inference_manifest["profiles"])
    cp = {}
    for variant in variants:
        path = inference / f"{variant}_labels.tif"
        if rn.sha256(path) != inference_manifest["profiles"][variant]["mask_sha256"]:
            raise ValueError(f"Cellpose mask hash mismatch: {variant}")
        cp[variant] = tifffile.imread(path)
        if cp[variant].shape != labels.shape:
            raise ValueError(f"Cellpose shape mismatch: {variant}")
    output.mkdir(parents=True, exist_ok=True)
    results, agreements, display_masks = [], {}, {}
    ref_accepted = accepted_labels(labels, nuclei)
    reference_by_id = {n["nucleus_id"]: n for n in nuclei}
    for variant in ["stardist", *[v for v in variants if v != "cyto3_pnoc_dapi"]]:
        variant_labels = labels if variant == "stardist" else cp[variant]
        records = nuclei if variant == "stardist" else rn.nucleus_records(variant_labels, spacing, config["nuclei"], image[metadata["resolved_channels"]["dapi"]])
        display_masks[variant] = accepted_labels(variant_labels, records)
        pairs, agreement = match_instances(ref_accepted, display_masks[variant])
        agreements[variant] = agreement
        if variant != "stardist":
            rn.write_csv(output / f"{variant}_nuclear_matches.csv", pairs,
                         ["reference_id", "query_id", "iou", "reference_coverage", "query_coverage"])
        for radius in [0., 2.]:
            name = f"{variant}_{'nuclear' if radius == 0 else 'plus2um'}"
            assignments = rn.build_assignment_labels(variant_labels, records, spacing, radius)
            summary, classified = analyze_masks(variant_labels, records, assignments, spots, image, metadata, config, output / name, name)
            results.append(summary)
            if variant == "stardist" and radius == 2:
                for key in ("total_nuclei", "pnoc_positive", "cre_positive", "double_positive", "pct_cre_expressing_pnoc"):
                    if summary[key] != baseline[key]:
                        raise ValueError(f"Frozen baseline does not reproduce: {key}")
                display_masks["stardist_plus2um"] = accepted_labels(assignments, records)
            if variant != "stardist" and radius == 2:
                query_by_id = {n["nucleus_id"]: n for n in classified}
                changed = []
                for pair in pairs:
                    a, b = reference_by_id[pair["reference_id"]], query_by_id[pair["query_id"]]
                    changed.append({**pair, "reference_pnoc_spots": a["pnoc_spots"], "query_pnoc_spots": b["pnoc_spots"],
                                    "reference_cre_spots": a["cre_spots"], "query_cre_spots": b["cre_spots"],
                                    "reference_class": a["classification"], "query_class": b["classification"],
                                    "class_changed": a["classification"] != b["classification"]})
                rn.write_csv(output / f"{variant}_matched_classifications.csv", changed,
                             ["reference_id", "query_id", "iou", "reference_coverage", "query_coverage",
                              "reference_pnoc_spots", "query_pnoc_spots", "reference_cre_spots", "query_cre_spots",
                              "reference_class", "query_class", "class_changed"])
                agreements[variant]["matched_class_changes_plus2um"] = sum(r["class_changed"] for r in changed)
    assignments, links = link_cell_masks(labels, nuclei, cp["cyto3_pnoc_dapi"])
    rn.write_csv(output / "cyto3_cell_links.csv", links, ["cell_id", "reference_nucleus_id", "status", "majority_nuclear_occupants",
                                                        "nuclear_coverage", "cell_area_px", "other_nuclear_overlap_px"])
    eligible_ids = {r["reference_nucleus_id"] for r in links if r["status"] == "eligible"}
    records = copy.deepcopy(nuclei)
    for n in records:
        if n["accepted"] and n["nucleus_id"] not in eligible_ids:
            n["accepted"] = False
            n["exclusion_reason"] = "no_eligible_cell_outline"
    summary, classified = analyze_masks(labels, records, assignments, spots, image, metadata, config,
                                         output / "cyto3_pnoc_cell_regions", "cyto3_pnoc_cell_regions")
    results.append(summary)
    display_masks["cyto3_pnoc_cell_regions"] = assignments
    display_masks["cyto3_raw_cells"] = cp["cyto3_pnoc_dapi"]
    baseline_cre = {n["nucleus_id"] for n in nuclei if n["accepted"] and n["cre_positive"]}
    mapped_baseline_cre = [n for n in classified if n["accepted"] and n["nucleus_id"] in baseline_cre]
    positive = sum(n["pnoc_positive"] for n in mapped_baseline_cre)
    unknown = len(baseline_cre) - len(mapped_baseline_cre)
    fixed = {"baseline_cre_positive": len(baseline_cre), "mapped_baseline_cre": len(mapped_baseline_cre),
             "unmapped_baseline_cre_unknown": unknown, "pnoc_positive_mapped_baseline_cre": positive,
             "pct_pnoc_among_mapped_baseline_cre": 100 * positive / len(mapped_baseline_cre) if mapped_baseline_cre else None,
             "whole_baseline_cre_fraction_bounds_pct": [100 * positive / len(baseline_cre), 100 * (positive + unknown) / len(baseline_cre)] if baseline_cre else [None, None]}
    rn.write_json(output / "cell_region_coverage.json", {**fixed, "link_status_counts": {
        status: sum(r["status"] == status for r in links) for status in sorted({r["status"] for r in links})},
        "interpretation": "Unknown/unmapped objects are not negatives. The mapped subset may be selected by Pnoc-guided segmentation; these percentages are not whole-cohort accuracy estimates."})
    rn.write_csv(output / "comparison.csv", results, list(results[0]))
    rn.write_json(output / "nuclear_agreement.json", agreements)
    rn.write_json(output / "comparison_manifest.json", {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": rn.sha256(__file__),
        "analysis_helper_sha256": rn.sha256(rn.__file__), "source_reference": str(reference.resolve()),
        "source_export_sha256": {n: rn.sha256(reference / n) for n in ("metadata.json", "nuclei.csv", "spots.csv", "nuclei_labels.tif", "summary.json")},
        "inference_manifest": inference_manifest, "cutoffs": {"pnoc": 3, "cre": 3},
        "cell_link_nuclear_coverage_min": .5, "iou_match_threshold": .5,
        "nuclear_acceptance_filters": config["nuclei"], "packages": {p: rn.package_version(p) for p in ("numpy", "scipy", "scikit-image", "tifffile")},
        "cell_link_policy": "One majority-covered candidate nucleus, accepted by reference filters; unique mapping; no image-border cell; block other nuclei. No cell-area cutoff is inferred from nuclear area limits.",
        "annotations_available": False, "interpretation": "Method agreement only. No expert segmentation ground truth supplied; Pnoc guidance can bias Pnoc association."})
    save_figures(output, image, metadata, display_masks)
    save_blinded_review(output, image, metadata, display_masks)
    return results


def save_figures(output, image, metadata, masks):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    channels = metadata["resolved_channels"]
    rgb = np.stack([rn.display_scale(image[channels[g]]) for g in ("cre", "pnoc", "dapi")], axis=-1)
    panels = [("stardist", "StarDist: accepted nuclei"), ("nuclei_dapi", "Cellpose nuclei: accepted DAPI objects"),
              ("cyto3_dapi", "Cellpose cyto3: accepted DAPI objects"),
              *[(k, f"Cellpose {k.split('_')[0]}: reference-preprocessed DAPI")
                for k in ("nuclei_dapi_clahe", "cyto3_dapi_clahe") if k in masks],
              ("stardist_plus2um", "StarDist: nucleus +2 µm regions"),
              ("cyto3_raw_cells", "Cellpose cyto3: all Pnoc/DAPI cell estimates"),
              ("cyto3_pnoc_cell_regions", "Cellpose cell regions with eligible nuclear links")]
    for crop_name, sl in [("overview", np.s_[:, :]), ("detail", np.s_[256:512, 256:512])]:
        columns = 4 if len(panels) > 6 else 3
        fig, axes = plt.subplots(2, columns, figsize=(6 * columns, 12), layout="constrained")
        for ax, (key, title) in zip(axes.flat, panels):
            ax.imshow(rgb[sl], interpolation="nearest")
            boundary = segmentation.find_boundaries(masks[key], mode="inner")[sl]
            overlay = np.zeros((*boundary.shape, 4)); overlay[boundary] = [1, 1, 0, .9]
            ax.imshow(overlay, interpolation="nearest")
            ax.set_title(title, fontsize=11); ax.axis("off")
        for ax in list(axes.flat)[len(panels):]:
            ax.axis("off")
        fig.suptitle("Same image and frozen peak detections; ≥3 peaks per gene\nCre red / Pnoc green / DAPI blue; yellow = estimated boundaries; display contrast adjusted")
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"segmentation_{crop_name}.{suffix}", dpi=170)
        plt.close(fig)


def save_blinded_review(output, image, metadata, masks, seed=20261023):
    """Eight spatial tiles picked before considering any expression results.

    The model key stays outside the reviewer directory. Draw raw-image expert
    labels first; prediction overlays are a separate second-stage review.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder = output / "blinded_review"
    raw_folder, prediction_folder = folder / "raw_images", folder / "prediction_review"
    raw_folder.mkdir(parents=True); prediction_folder.mkdir()
    spacing, channels = metadata["pixel_size_yx_um"], metadata["resolved_channels"]
    rgb = np.stack([rn.display_scale(image[channels[g]]) for g in ("cre", "pnoc", "dapi")], axis=-1)
    candidates = [(x, y, min(image.shape[2], x + 256), min(image.shape[1], y + 256))
                  for y in range(0, image.shape[1], 256) for x in range(0, image.shape[2], 256)]
    rng = np.random.default_rng(seed)
    selection = rng.permutation(len(candidates))[:min(8, len(candidates))]
    tiles, review_rows, key_rows = [], [], []
    sequence = 0
    tasks = [("nuclei", [k for k in ("stardist", "nuclei_dapi", "cyto3_dapi", "nuclei_dapi_clahe", "cyto3_dapi_clahe") if k in masks]),
             ("region_boundary_estimate", [k for k in ("stardist_plus2um", "cyto3_raw_cells", "cyto3_pnoc_cell_regions") if k in masks])]
    for number, index in enumerate(selection, 1):
        x0, y0, x1, y1 = candidates[index]
        hx0, hy0, hx1, hy1 = max(0, x0 - 32), max(0, y0 - 32), min(image.shape[2], x1 + 32), min(image.shape[1], y1 + 32)
        tile_id = f"T{number:02d}"
        split = "calibration" if number <= min(4, len(selection) // 2) else "held_out"
        tiles.append({"tile_id": tile_id, "split": split, "core_x0": x0, "core_y0": y0, "core_x1": x1, "core_y1": y1,
                      "image_x0": hx0, "image_y0": hy0, "image_x1": hx1, "image_y1": hy1})
        sl = np.s_[hy0:hy1, hx0:hx1]
        tifffile.imwrite(raw_folder / f"{tile_id}_raw_channels.tif", image[:, sl[0], sl[1]], imagej=True,
                         resolution=(1 / spacing[1], 1 / spacing[0]),
                         metadata={"axes": "CYX", "unit": "um", "Labels": [c["name"] for c in metadata["channels"]]})
        fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
        axes[0].imshow(rn.display_scale(image[channels["dapi"]])[sl], cmap="gray")
        axes[1].imshow(rn.display_scale(image[channels["pnoc"]])[sl], cmap="gray")
        axes[2].imshow(rgb[sl])
        for ax, title in zip(axes, ["DAPI", "Pnoc", "Cre red / Pnoc green / DAPI blue"]):
            ax.set_title(title); ax.axis("off")
        fig.suptitle(f"{tile_id}: raw image; no model predictions")
        fig.savefig(raw_folder / f"{tile_id}_raw.png", dpi=120); plt.close(fig)
        for task, variants in tasks:
            for variant in rng.permutation(variants):
                sequence += 1
                review_id = f"R{sequence:03d}"
                fig, ax = plt.subplots(figsize=(5, 5), layout="constrained")
                ax.imshow(rgb[sl], interpolation="nearest")
                boundary = segmentation.find_boundaries(masks[variant], mode="inner")[sl]
                overlay = np.zeros((*boundary.shape, 4)); overlay[boundary] = [1, 1, 0, 1]
                ax.imshow(overlay, interpolation="nearest")
                ax.set_title(f"{review_id} | {tile_id} | {task.replace('_', ' ')}", fontsize=9); ax.axis("off")
                fig.savefig(prediction_folder / f"{review_id}.png", dpi=130); plt.close(fig)
                rn.write_imagej_labels(prediction_folder / f"{review_id}_labels.tif", masks[variant][sl], spacing)
                review_rows.append({"review_id": review_id, "tile_id": tile_id, "split": split, "task": task,
                                    "reviewer": "", "missed_objects": "", "false_objects": "", "split_errors": "",
                                    "merge_errors": "", "boundary_quality_1_to_5_or_unknown": "",
                                    "cell_boundary_not_discernible": "", "comments": ""})
                key_rows.append({"review_id": review_id, "tile_id": tile_id, "task": task, "method": variant})
    rn.write_csv(folder / "tiles.csv", tiles, list(tiles[0]))
    rn.write_csv(folder / "review_scores.csv", review_rows, list(review_rows[0]))
    rn.write_csv(output / "blinding_key_DO_NOT_SHARE.csv", key_rows, list(key_rows[0]))
    tifffile.imwrite(raw_folder / "full_image_raw_channels.tif", image, imagej=True,
                     resolution=(1 / spacing[1], 1 / spacing[0]),
                     metadata={"axes": "CYX", "unit": "um", "Labels": [c["name"] for c in metadata["channels"]]})
    for name in ("manual_nuclei", "manual_cells", "annotated_roi"):
        rn.write_imagej_labels(folder / f"TEMPLATE_{name}.tif", np.zeros(image.shape[1:], np.uint16), spacing)
    (folder / "README.md").write_text(
        "# Blinded segmentation review\n\n"
        "No expert annotations are included. TEMPLATE files contain only zeros and are not ground truth.\n\n"
        "1. Work from raw_images first. Channels are Cy5/Cre, Cy3 (probe identity unconfirmed), ATTO520/Pnoc, DAPI.\n"
        "2. Eight spatial tiles were selected without using RNA classifications, with a 32-pixel context margin. "
        "tiles.csv gives core/image coordinates in the full image. Four tiles are calibration and four held out. "
        "Use the held-out tiles only after review criteria or any training choices are fixed.\n"
        "3. In full-image coordinates, draw instance IDs 1,2,... for manually annotated nuclei; "
        "save a completed copy of TEMPLATE_manual_nuclei.tif. Mark fully reviewed pixels in TEMPLATE_annotated_roi.tif. "
        "Distinguish unreviewed regions from background. Do not reuse model masks as expert labels.\n"
        "4. If genuine cell boundaries are visible, draw a separate manual-cell mask. "
        "Use uncertain/not discernible when cytoplasmic or membrane information is insufficient; "
        "do not invent whole-cell borders from isolated Pnoc dots. Keep nuclear and cell annotations separate.\n"
        "5. Only after independent annotation, inspect prediction_review and complete review_scores.csv. "
        "Nuclear and assignment-region tasks are distinct. Each task's model order is randomized independently. "
        "A score of 1 means very poor boundary agreement and 5 means very good; use unknown if unassessable.\n"
        "6. Keep the model key, numerical method results and coexpression percentages away from reviewers. "
        "The administrator's blinding key is outside this folder. Algorithm identity is hidden, but distinctive "
        "shapes can still make methods recognizable; this does not guarantee successful blinding.\n"
        "7. Return completed label masks, reviewed ROI, scores and reviewer IDs. "
        "No accuracy comparison should be calculated from untouched zero templates or unreviewed pixels.\n",
        encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True, help="Reference scene output folder")
    parser.add_argument("--input", type=Path, required=True, help="Raw CYX input_channels.npz")
    parser.add_argument("--inference", type=Path, required=True, help="Completed cellpose_infer.py output")
    parser.add_argument("--output", type=Path, required=True, help="New or empty comparison folder")
    args = parser.parse_args(argv)
    results = run_comparison(args.reference, args.input, args.inference, args.output)
    for row in results:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
