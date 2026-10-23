#!/usr/bin/env python3
"""Batch DAPI-associated Pnoc/Cre analysis of single-section CZI images.

Counts detected spot centers within nuclear masks or a configured surrounding
region. The primary endpoint is detectable Pnoc association among Cre-positive
objects; this does not estimate the total cytoplasmic Pnoc-positive population.
Parameters are provisional; see README.md before interpreting.
"""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import logging
import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import czifile
import numpy as np
from scipy import ndimage as ndi
from skimage import exposure, filters, measure, morphology, segmentation
from skimage.feature import peak_local_max
import tifffile

VERSION = "0.2.3"
DEFAULT_CONFIG = {
    "channels": {
        "dapi": {"names": ["DAPI-T4", "DAPI"], "index": None, "fallback_index": 3},
        "pnoc": {"names": ["At520-T3", "ATTO 520", "At520"], "index": None, "fallback_index": 2},
        "cre": {"names": ["Cy5-T1", "Cy5"], "index": None, "fallback_index": 0},
    },
    "pixel_size_um": None,
    "assignment": {"expand_um": 0.0},
    "nuclei": {
        "method": "watershed",
        "threshold_method": "li",
        "smooth_sigma_um": 0.7,
        "background_sigma_um": 12.0,
        "threshold_factor": 1.0,
        "min_area_um2": 15.0,
        "max_area_um2": 250.0,
        "fill_holes_um2": 8.0,
        "seed_min_distance_um": 3.0,
        "exclude_border": True,
        "review_min_solidity": 0.8,
        "review_max_eccentricity": 0.95,
        "local_block_um": 25.0,
        "local_otsu_floor": 0.6,
        "local_min_solidity": 0.8,
        "stardist_model": "2D_versatile_fluo",
        "stardist_smooth_um": 1.4,
        "stardist_target_pixel_um": 0.92,
        "stardist_tiles": [2, 2],
        "stardist_clahe": False,
        "stardist_clahe_kernel_um": 25.0,
        "stardist_clahe_clip": 0.01,
        "stardist_probability_threshold": None,
        "min_dapi_contrast": None,
        "contrast_ring_inner_um": 1.5,
        "contrast_ring_outer_um": 3.5,
    },
    "spots": {
        "pnoc": {"method": "log", "sigma_um": [0.45, 0.7, 1.0], "background_sigma_um": 4.0,
                 "min_distance_um": 1.4, "snr_threshold": 5.0,
                 "min_response": 0.0, "min_contrast": 0.0,
                 "intensity_background_radius_um": 8.0, "intensity_smooth_um": 0.55,
                 "cluster_area_um2": 12.0, "positive_min_spots": 3},
        "cre": {"method": "log", "sigma_um": [0.45, 0.7, 1.0], "background_sigma_um": 4.0,
                "min_distance_um": 1.4, "snr_threshold": 5.0,
                "min_response": 0.0, "min_contrast": 0.0,
                "intensity_background_radius_um": 8.0, "intensity_smooth_um": 0.55,
                "cluster_area_um2": 12.0, "positive_min_spots": 3},
    },
    "qc": {"saturation_value": None, "review_crops": 12, "crop_width_um": 45.0},
}

QC_FLAG_NAMES = ["irregular_nuclear_shape", "dapi_saturation", "pnoc_saturation", "cre_saturation",
                 "pnoc_saturated_or_clustered_region", "cre_saturated_or_clustered_region",
                 "saturated_signal", "large_bright_cluster", "near_expression_threshold",
                 "near_pnoc_cutoff", "near_cre_cutoff", "saturated_pnoc_below_cutoff", "saturated_cre_below_cutoff"]
NUCLEUS_FIELDS = ["file", "scene", "nucleus_id", "x_px", "y_px", "x_um", "y_um",
                  "area_um2", "solidity", "eccentricity", "dapi_contrast_ratio", "accepted", "exclusion_reason",
                  "pnoc_nuclear_spots", "cre_nuclear_spots", "pnoc_spots", "cre_spots",
                  "pnoc_positive", "cre_positive", "classification", "needs_review", "review_priority", "review_reasons", "qc_flags",
                  *[f"flag_{name}" for name in QC_FLAG_NAMES]]
SPOT_FIELDS = ["file", "scene", "gene", "spot_id", "x_px", "y_px", "x_um", "y_um",
               "sigma_um", "response", "detection_snr", "contrast", "nucleus_id", "assignment_id", "counted", "qc_flags"]
COUNT_FIELDS = ["total_nuclei", "excluded_nuclei", "pnoc_positive", "cre_positive",
                "double_positive", "pnoc_only", "cre_only", "neither", "nuclei_needing_review", "nuclei_cre_endpoint_review"]
SUMMARY_FIELDS = ["file", "scene", "status", *COUNT_FIELDS,
                  "pct_pnoc_expressing_cre", "pct_cre_expressing_pnoc", "assignment_expand_um", "analysis_profile", "warnings", "error"]
COLORS = {"neither": "#a0a0a0", "pnoc_only": "#00dd66", "cre_only": "#ff5577",
          "double_positive": "#ffe14a", "excluded": "#7080a0"}


def merge_config(overrides: dict, base: dict | None = None, prefix: str = "") -> dict:
    """Merge partial JSON configuration, rejecting misspelled/unknown options."""
    result = copy.deepcopy(DEFAULT_CONFIG if base is None else base)
    if not isinstance(overrides, dict):
        raise ValueError(f"{prefix or 'config'} must be an object")
    for key, value in overrides.items():
        if key not in result:
            raise ValueError(f"Unknown setting: {prefix}{key}")
        if isinstance(result[key], dict):
            result[key] = merge_config(value, result[key], f"{prefix}{key}.")
        else:
            result[key] = value
    return result


def validate_config(config: dict) -> None:
    def number(value, name, *, zero=False, integer=False):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{name} must be a finite number")
        if value < 0 or (not zero and value == 0) or (integer and int(value) != value):
            raise ValueError(f"Invalid value for {name}: {value}")

    if config["pixel_size_um"] is not None:
        p = config["pixel_size_um"]
        if not isinstance(p, list) or len(p) != 2:
            raise ValueError("pixel_size_um must be null or [y_um, x_um]")
        for v in p:
            number(v, "pixel_size_um")
    number(config["assignment"]["expand_um"], "assignment.expand_um", zero=True)
    for role, channel in config["channels"].items():
        if not isinstance(channel["names"], list) or not channel["names"] or not all(
                isinstance(n, str) and n for n in channel["names"]):
            raise ValueError(f"channels.{role}.names must be a nonempty list of strings")
        for name in ["index", "fallback_index"]:
            if channel[name] is not None:
                number(channel[name], f"channels.{role}.{name}", zero=True, integer=True)
    n = config["nuclei"]
    for key, value in n.items():
        if key == "method":
            if value not in ["watershed", "local_watershed", "stardist"]:
                raise ValueError("nuclei.method must be watershed, local_watershed, or stardist")
        elif key in ["min_dapi_contrast", "stardist_probability_threshold"]:
            if value is not None:
                number(value, f"nuclei.{key}")
                if key == "stardist_probability_threshold" and value >= 1:
                    raise ValueError("stardist_probability_threshold must be between 0 and 1")
        elif key == "stardist_model":
            if not isinstance(value, str) or not value:
                raise ValueError("stardist_model must name a registered StarDist model")
        elif key == "stardist_tiles":
            if not isinstance(value, list) or len(value) != 2:
                raise ValueError("stardist_tiles must be [y_tiles, x_tiles]")
            for v in value:
                number(v, "stardist_tiles", integer=True)
        elif key == "threshold_method":
            if value not in ["li", "otsu"]:
                raise ValueError("threshold_method must be 'li' or 'otsu'")
        elif key in ["exclude_border", "stardist_clahe"]:
            if not isinstance(value, bool):
                raise ValueError(f"{key} must be true or false")
        else:
            number(value, f"nuclei.{key}", zero=(key == "fill_holes_um2"))
            if key in ["review_min_solidity", "review_max_eccentricity", "local_min_solidity", "stardist_clahe_clip"] and value > 1:
                raise ValueError(f"nuclei.{key} must be <= 1")
    if n["max_area_um2"] < n["min_area_um2"]:
        raise ValueError("max_area_um2 must be >= min_area_um2")
    if n["background_sigma_um"] <= n["smooth_sigma_um"]:
        raise ValueError("Nuclear background scale must exceed smoothing scale")
    if n["contrast_ring_outer_um"] <= n["contrast_ring_inner_um"]:
        raise ValueError("DAPI contrast outer ring must be larger than inner ring")
    for gene, s in config["spots"].items():
        if not isinstance(s["sigma_um"], list) or not s["sigma_um"]:
            raise ValueError(f"spots.{gene}.sigma_um must be a nonempty list")
        for v in s["sigma_um"]:
            number(v, f"spots.{gene}.sigma_um")
        if s["sigma_um"] != sorted(set(s["sigma_um"])):
            raise ValueError("sigma_um values must be increasing and unique")
        for key, value in s.items():
            if key == "method":
                if value not in ["log", "intensity"]:
                    raise ValueError(f"spots.{gene}.method must be log or intensity")
            elif key != "sigma_um":
                number(value, f"spots.{gene}.{key}", zero=key in ["min_response", "min_contrast"],
                       integer=key == "positive_min_spots")
        if s["background_sigma_um"] <= max(s["sigma_um"]):
            raise ValueError("Spot background scale must exceed all spot scales")
    q = config["qc"]
    number(q["review_crops"], "review_crops", zero=True, integer=True)
    number(q["crop_width_um"], "crop_width_um")
    if q["saturation_value"] is not None:
        number(q["saturation_value"], "saturation_value")


def normalize_name(name: str) -> str:
    return "".join(c.lower() for c in name if c.isalnum())


def resolve_channels(channel_metadata: list[dict], count: int, config: dict):
    if channel_metadata and len(channel_metadata) != count:
        raise ValueError("Channel metadata count differs from scene channel count; scene-specific mapping is required")
    selected, warnings = {}, []
    for role, settings in config.items():
        index = settings["index"]
        if index is None:
            aliases = {normalize_name(n) for n in settings["names"]}
            matches = [i for i, m in enumerate(channel_metadata)
                       if any(normalize_name(m.get(k, "")) in aliases for k in ["name", "fluor"])]
            if len(matches) > 1:
                raise ValueError(f"Ambiguous metadata for {role}; set its channel index explicitly")
            if matches:
                index = matches[0]
            elif not channel_metadata or not any(m.get("name") or m.get("fluor") for m in channel_metadata):
                index = settings["fallback_index"]
                warnings.append(f"{role}: missing channel metadata; using fallback index {index}")
            else:
                raise ValueError(f"No metadata match for {role}; verify channel identity and set index explicitly")
        if index is None or not 0 <= index < count:
            raise ValueError(f"Channel index {index} for {role} is out of range")
        selected[role] = int(index)
    if len(set(selected.values())) != len(selected):
        raise ValueError("DAPI, Pnoc and Cre must use different channels")
    return selected, warnings


def metadata_from_xml(xml: str):
    root = ET.fromstring(xml)
    channels = []
    for n in root.findall(".//Information/Image/Dimensions/Channels/Channel"):
        bit_text = n.findtext("ComponentBitCount") or root.findtext(".//Information/Image/ComponentBitCount")
        try:
            bit_count = int(bit_text)
        except (ValueError, TypeError):
            bit_count = None
        channels.append({"name": n.get("Name", ""), "fluor": n.findtext("Fluor", ""),
                         "id": n.get("Id", ""),
                         "digital_gain": n.findtext("DetectorSettings/DigitalGain"),
                         "detector_voltage": n.findtext("DetectorSettings/Voltage"),
                         "component_bit_count": bit_count})
    sizes = {}
    for n in root.findall(".//Scaling/Items/Distance"):
        if n.get("Id") in ["X", "Y"]:
            sizes[n.get("Id")] = float(n.findtext("Value")) * 1e6
    return channels, sizes


def as_cyx(array: np.ndarray, dims) -> np.ndarray:
    dims = list(dims)
    if len(dims) != array.ndim or len(set(dims)) != len(dims):
        raise ValueError(f"Invalid CZI dimensions: {dims} / {array.shape}")
    for axis in range(len(dims) - 1, -1, -1):
        if dims[axis] not in ["C", "Y", "X"]:
            if array.shape[axis] != 1:
                raise ValueError(f"Only single 2D sections are supported; {dims[axis]} has size {array.shape[axis]}")
            array = np.take(array, 0, axis=axis)
            dims.pop(axis)
    if set(dims) != {"C", "Y", "X"}:
        raise ValueError(f"Expected C, Y, X dimensions, found {dims}")
    return np.transpose(array, [dims.index(k) for k in "CYX"])


def physical_footprint(radius_um: float, spacing) -> np.ndarray:
    ry, rx = [max(1, math.ceil(radius_um / p)) for p in spacing]
    y, x = np.mgrid[-ry:ry + 1, -rx:rx + 1]
    return (y * spacing[0]) ** 2 + (x * spacing[1]) ** 2 <= radius_um ** 2


_STARDIST_MODELS = {}


def package_version(distribution, module=None):
    """Best-effort provenance, including tensorflow-cpu / tensorflow-macos installs."""
    try:
        return str(importlib.metadata.version(distribution))
    except Exception:
        # Missing or damaged distribution metadata does not imply a broken import.
        pass
    module = module or {"scikit-image": "skimage"}.get(distribution, distribution.replace("-", "_"))
    try:
        version = getattr(importlib.import_module(module), "__version__", None)
        return str(version) if version is not None else "unknown"
    except Exception:
        return "unknown"


def describe_profile(config):
    n = config["nuclei"]
    method = n["method"]
    if method == "watershed":
        method += "/" + n["threshold_method"]
    elif method == "stardist" and n["stardist_clahe"]:
        method += "/CLAHE"
    return (f"{method}; Pnoc={config['spots']['pnoc']['method']}; Cre={config['spots']['cre']['method']}; "
            f"expansion={config['assignment']['expand_um']:g} um")


def preflight(config):
    """Fail once for an unusable runtime before hashing or processing any images."""
    if not hasattr(czifile.CziFile, "scenes"):
        raise RuntimeError("The installed czifile lacks the required scenes API. "
                           "Use Python 3.12 and install this repository's requirements.")
    if config["nuclei"]["method"] == "stardist":
        try:
            load_stardist_model(config["nuclei"]["stardist_model"])
        except Exception as exc:
            raise RuntimeError(f"StarDist initialization failed before analysis: {exc}") from exc


def load_stardist_model(name):
    """Use the publisher's registered archive hash and normal TLS verification."""
    if name not in _STARDIST_MODELS:
        os.environ.setdefault("KERAS_HOME", str(Path(__file__).resolve().parent / ".cache" / "keras"))
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
        os.environ.setdefault("TF_NUM_INTRAOP_THREADS", "2")
        os.environ.setdefault("TF_NUM_INTEROP_THREADS", "2")
        try:
            from stardist.models import StarDist2D
            from csbdeep.models.pretrained import get_model_folder
        except ImportError as exc:
            raise RuntimeError("StarDist selected but unavailable. Install requirements-stardist.txt; "
                               "no automatic fallback is used.") from exc
        model = StarDist2D.from_pretrained(name)
        folder = get_model_folder(StarDist2D, name)
        info = {"model": name,
                "packages": {p: package_version(p) for p in ["stardist", "csbdeep", "tensorflow"]},
                "model_file_sha256": {p.name: sha256(p) for p in Path(folder).iterdir() if p.is_file()}}
        _STARDIST_MODELS[name] = model, info
    return _STARDIST_MODELS[name]


def stardist_labels(image, spacing, config):
    model, info = load_stardist_model(config["stardist_model"])
    from csbdeep.utils import normalize
    smooth = ndi.gaussian_filter(image.astype(np.float32), [config["stardist_smooth_um"] / p for p in spacing])
    normalized = np.clip(normalize(smooth, 1, 99.8), 0, 1)
    if config["stardist_clahe"]:
        normalized = exposure.equalize_adapthist(
            normalized, kernel_size=tuple(max(1, int(round(config["stardist_clahe_kernel_um"] / p))) for p in spacing),
            clip_limit=config["stardist_clahe_clip"])
    kwargs = {} if config["stardist_probability_threshold"] is None else {
        "prob_thresh": config["stardist_probability_threshold"]}
    labels, _ = model.predict_instances(
        normalized,
        scale=tuple(p / config["stardist_target_pixel_um"] for p in spacing),
        n_tiles=tuple(config["stardist_tiles"]), show_tile_progress=False, **kwargs)
    return labels.astype(np.int32), {"method": "stardist", **info}


def local_watershed_labels(image, spacing, config):
    """Local DAPI thresholding inspired by the supplied rnascope_pnoc_cre.py."""
    smooth = ndi.gaussian_filter(image.astype(np.float32), [config["smooth_sigma_um"] / p for p in spacing])
    block = tuple(max(3, int(round(config["local_block_um"] / p)) | 1) for p in spacing)
    threshold = filters.threshold_local(smooth, block_size=block, method="gaussian")
    floor = float(config["local_otsu_floor"] * filters.threshold_otsu(smooth))
    mask = (smooth > threshold) & (smooth > floor)
    mask = ndi.binary_fill_holes(mask)
    mask = ndi.binary_opening(mask, structure=morphology.disk(1))
    distance = ndi.distance_transform_edt(mask, sampling=spacing)
    distance = ndi.gaussian_filter(distance, 1.0)
    # Rectangular footprint matches the reference's Chebyshev peak suppression.
    footprint = np.ones(tuple(2 * max(2, int(round(config["seed_min_distance_um"] / p))) + 1 for p in spacing), bool)
    coords = peak_local_max(distance, footprint=footprint, labels=measure.label(mask), exclude_border=False)
    markers = np.zeros(mask.shape, np.int32)
    if len(coords):
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    labels = segmentation.watershed(-distance, markers, mask=mask).astype(np.int32)
    return labels, {"method": "local_watershed", "block_size_yx_px": block, "global_floor": floor}


def segment_nuclei(image: np.ndarray, spacing, config: dict):
    if config["method"] != "watershed":
        if not np.ptp(image):
            return np.zeros(image.shape, np.int32), [], {"method": config["method"], "flat_dapi": True}
        labels, info = (stardist_labels(image, spacing, config) if config["method"] == "stardist" else
                        local_watershed_labels(image, spacing, config))
        return labels, nucleus_records(labels, spacing, config, image), info
    f = image.astype(np.float32)
    smooth = ndi.gaussian_filter(f, [config["smooth_sigma_um"] / p for p in spacing])
    background = ndi.gaussian_filter(f, [config["background_sigma_um"] / p for p in spacing])
    corrected = np.maximum(smooth - background, 0)
    if not np.any(corrected):
        return np.zeros(image.shape, np.int32), [], {"threshold": None}
    threshold_function = filters.threshold_li if config["threshold_method"] == "li" else filters.threshold_otsu
    threshold = float(threshold_function(corrected) * config["threshold_factor"])
    mask = corrected > threshold
    if config["fill_holes_um2"]:
        # Fill only enclosed holes, never background regions touching image edges.
        holes, num = ndi.label(~mask)
        sizes = np.bincount(holes.ravel())
        border = np.unique(np.concatenate([holes[0], holes[-1], holes[:, 0], holes[:, -1]]))
        fill = sizes * np.prod(spacing) <= config["fill_holes_um2"]
        fill[border] = False
        fill[0] = False
        mask |= fill[holes]
    # Remove tiny debris before seeding; final area checks happen after splitting.
    components, _ = ndi.label(mask)
    areas = np.bincount(components.ravel()) * np.prod(spacing)
    keep = areas >= config["min_area_um2"] * 0.35
    keep[0] = False
    mask = keep[components]
    distance = ndi.distance_transform_edt(mask, sampling=spacing)
    coords = peak_local_max(distance, labels=mask.astype(np.uint8), exclude_border=False,
                            footprint=physical_footprint(config["seed_min_distance_um"], spacing))
    markers = np.zeros(mask.shape, np.int32)
    if len(coords):
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    labels = segmentation.watershed(-distance, markers, mask=mask).astype(np.int32)
    return labels, nucleus_records(labels, spacing, config, image), {
        "threshold": threshold, "method": "watershed", "threshold_method": config["threshold_method"]}


def dapi_contrast_ratios(labels, image, spacing, config):
    ids = np.arange(1, int(labels.max()) + 1)
    if not len(ids):
        return {}
    smooth = ndi.gaussian_filter(image.astype(np.float32), [config["stardist_smooth_um"] / p for p in spacing])
    near = segmentation.expand_labels(labels, config["contrast_ring_inner_um"], spacing=spacing)
    ring = segmentation.expand_labels(labels, config["contrast_ring_outer_um"], spacing=spacing)
    ring[near > 0] = 0
    # Empty rings are unknown, rather than evidence for rejecting a crowded nucleus.
    sums = np.bincount(ring.ravel(), weights=smooth.ravel(), minlength=len(ids) + 1)
    counts = np.bincount(ring.ravel(), minlength=len(ids) + 1)
    inner_sum = np.bincount(labels.ravel(), weights=smooth.ravel(), minlength=len(ids) + 1)
    inner_count = np.bincount(labels.ravel(), minlength=len(ids) + 1)
    result = {}
    for label in ids:
        if counts[label] and sums[label] > 0 and inner_count[label]:
            result[int(label)] = float((inner_sum[label] / inner_count[label]) / (sums[label] / counts[label]))
        else:
            result[int(label)] = None
    return result


def nucleus_records(labels, spacing, config, image=None):
    contrast = dapi_contrast_ratios(labels, image, spacing, config) if config["min_dapi_contrast"] is not None else {}
    records = []
    for r in measure.regionprops(labels):
        y, x = r.centroid
        area = float(r.area * np.prod(spacing))
        reasons = []
        if area < config["min_area_um2"]:
            reasons.append("too_small")
        if area > config["max_area_um2"]:
            reasons.append("too_large")
        if config["exclude_border"] and (r.bbox[0] == 0 or r.bbox[1] == 0 or
                r.bbox[2] == labels.shape[0] or r.bbox[3] == labels.shape[1]):
            reasons.append("image_border")
        if config["method"] == "local_watershed" and r.solidity < config["local_min_solidity"]:
            reasons.append("low_solidity")
        ratio = contrast.get(int(r.label))
        if ratio is not None and ratio < config["min_dapi_contrast"]:
            reasons.append("low_dapi_contrast")
        records.append({"nucleus_id": int(r.label), "x_px": float(x), "y_px": float(y),
                        "x_um": float(x * spacing[1]), "y_um": float(y * spacing[0]),
                        "area_um2": area, "solidity": float(r.solidity), "eccentricity": float(r.eccentricity),
                        "dapi_contrast_ratio": ratio,
                        "qc_flags": "irregular_nuclear_shape" if r.solidity < config["review_min_solidity"] or
                        r.eccentricity > config["review_max_eccentricity"] else "",
                        "accepted": not reasons,
                        "exclusion_reason": ";".join(reasons)})
    return records


def robust_noise(a: np.ndarray) -> tuple[float, float]:
    median = float(np.median(a))
    noise = float(1.4826 * np.median(np.abs(a - median)))
    return median, max(noise, np.finfo(np.float32).eps)


def saturation_level(image, configured, component_bit_count=None):
    if configured is not None:
        return configured
    if np.issubdtype(image.dtype, np.integer):
        storage = np.iinfo(image.dtype)
        if isinstance(component_bit_count, int) and not isinstance(component_bit_count, bool) and 1 <= component_bit_count <= storage.bits:
            return min((1 << component_bit_count) - 1, int(storage.max))
        return int(storage.max)
    return None


def detect_spots(image: np.ndarray, spacing, config: dict, saturation=None):
    """Scale-normalized LoG peaks; no conversion of cluster intensity to molecules."""
    if config["method"] == "intensity":
        return detect_intensity_spots(image, spacing, config, saturation)
    f = image.astype(np.float32)
    residual = f - ndi.gaussian_filter(f, [config["background_sigma_um"] / p for p in spacing])
    responses = []
    for scale in config["sigma_um"]:
        sigma = [scale / p for p in spacing]
        response = sum(-ndi.gaussian_filter(residual, sigma, order=order) * sigma[axis] ** 2
                       for axis, order in enumerate([(2, 0), (0, 2)]))
        responses.append(response)
    stack = np.stack(responses)
    # Estimate noise on each SIGNED filter response. Estimating it after taking
    # the maximum across scales compresses the MAD and overcalls noise peaks.
    scale_stats = [robust_noise(r) for r in responses]
    scores = np.stack([(r - baseline) / noise for r, (baseline, noise) in zip(responses, scale_stats)])
    scores[stack < config["min_response"]] = -np.inf
    score = np.max(scores, axis=0)
    best_scale = np.argmax(stack, axis=0)
    response = np.max(stack, axis=0)
    # A truncated filter footprint cannot establish a spot's shape at an edge.
    excluded_edge = tuple(math.ceil(3 * max(config["sigma_um"]) / p) for p in spacing)
    coords = peak_local_max(score, threshold_abs=config["snr_threshold"], exclude_border=excluded_edge,
                            footprint=physical_footprint(config["min_distance_um"], spacing))
    raw_baseline, raw_noise = robust_noise(residual)
    bright, _ = ndi.label(residual > raw_baseline + 3 * raw_noise)
    areas = np.bincount(bright.ravel()) * np.prod(spacing)
    clusters = areas > config["cluster_area_um2"]
    clusters[0] = False
    cluster_mask = clusters[bright]
    saturated = np.zeros(image.shape, bool) if saturation is None else image >= saturation
    saturated_nearby = ndi.maximum_filter(saturated, size=3)
    spots = []
    for y, x in coords:
        if residual[y, x] < config["min_contrast"]:
            continue
        flags = []
        if cluster_mask[y, x]:
            flags.append("large_bright_cluster")
        if saturated_nearby[y, x]:
            flags.append("saturated_signal")
        spots.append({"spot_id": len(spots) + 1, "x_px": int(x), "y_px": int(y),
                      "x_um": float(x * spacing[1]), "y_um": float(y * spacing[0]),
                      "sigma_um": config["sigma_um"][int(best_scale[y, x])],
                      "response": float(response[y, x]), "detection_snr": float(score[y, x]),
                      "contrast": float(residual[y, x]),
                      "qc_flags": ";".join(flags)})
    diagnostics = {"method": "log", "snr_threshold": config["snr_threshold"],
                   "scale_response_thresholds": [max(config["min_response"], b + config["snr_threshold"] * n)
                                                 for b, n in scale_stats],
                   "scale_response_noise_mad": [n for _, n in scale_stats],
                   "excluded_image_edge_yx_px": excluded_edge,
                   "residual_noise_mad": raw_noise, "detected_spots": len(spots),
                   "saturation_value": saturation, "saturated_pixels": int(saturated.sum()),
                   "cluster_pixels": int(cluster_mask.sum())}
    return spots, diagnostics, cluster_mask | saturated


def detect_intensity_spots(image, spacing, config, saturation=None):
    """Morphological background removal and intensity peaks from the supplied script."""
    f = image.astype(np.float32)
    radius = tuple(max(3, int(round(config["intensity_background_radius_um"] / p))) for p in spacing)
    background = ndi.uniform_filter(ndi.grey_opening(f, size=tuple(2 * r + 1 for r in radius)), size=radius)
    residual = np.maximum(f - background, 0)
    smooth = ndi.gaussian_filter(residual, [config["intensity_smooth_um"] / p for p in spacing])
    baseline, noise = robust_noise(smooth)
    threshold = max(config["min_response"], baseline + config["snr_threshold"] * noise)
    coords = peak_local_max(smooth, threshold_abs=threshold, exclude_border=False,
                            footprint=physical_footprint(config["min_distance_um"], spacing))
    raw_baseline, raw_noise = robust_noise(residual)
    bright, _ = ndi.label(residual > raw_baseline + 3 * raw_noise)
    areas = np.bincount(bright.ravel()) * np.prod(spacing)
    clusters = areas > config["cluster_area_um2"]; clusters[0] = False
    cluster_mask = clusters[bright]
    saturated = np.zeros(image.shape, bool) if saturation is None else image >= saturation
    nearby = ndi.maximum_filter(saturated, size=3)
    spots = []
    for y, x in coords:
        if residual[y, x] < config["min_contrast"]:
            continue
        flags = []
        if cluster_mask[y, x]:
            flags.append("large_bright_cluster")
        if nearby[y, x]:
            flags.append("saturated_signal")
        spots.append({"spot_id": len(spots) + 1, "x_px": int(x), "y_px": int(y),
                      "x_um": float(x * spacing[1]), "y_um": float(y * spacing[0]),
                      "sigma_um": None, "response": float(smooth[y, x]),
                      "detection_snr": float((smooth[y, x] - baseline) / noise),
                      "contrast": float(residual[y, x]), "qc_flags": ";".join(flags)})
    return spots, {"method": "intensity", "response_threshold": threshold, "noise_mad": noise,
                   "detected_spots": len(spots), "saturation_value": saturation,
                   "saturated_pixels": int(saturated.sum()), "cluster_pixels": int(cluster_mask.sum())}, cluster_mask | saturated


def build_assignment_labels(labels, nuclei, spacing, expand_um):
    if not expand_um:
        return labels
    keep = np.zeros(int(labels.max()) + 1, bool)
    for n in nuclei:
        keep[n["nucleus_id"]] = n["accepted"]
    accepted = np.where(keep[labels], labels, 0)
    expanded = segmentation.expand_labels(accepted, distance=expand_um, spacing=spacing)
    # A rejected DAPI object is not cytoplasm belonging to a nearby accepted one.
    expanded[(labels > 0) & ~keep[labels]] = 0
    return expanded


def assign_and_classify(labels, nuclei, spots_by_gene, config, review_masks=None, assignment_labels=None):
    assignment_labels = labels if assignment_labels is None else assignment_labels
    lookup = {n["nucleus_id"]: n for n in nuclei}
    flags = {k: set(filter(None, n.get("qc_flags", "").split(";"))) for k, n in lookup.items()}
    for n in nuclei:
        n.update(pnoc_spots=0, cre_spots=0, pnoc_nuclear_spots=0, cre_nuclear_spots=0)
    for gene, spots in spots_by_gene.items():
        for spot in spots:
            label = int(labels[spot["y_px"], spot["x_px"]])
            owner = int(assignment_labels[spot["y_px"], spot["x_px"]])
            spot["nucleus_id"] = label
            spot["assignment_id"] = owner
            spot["counted"] = bool(owner and lookup[owner]["accepted"])
            if label:
                lookup[label][f"{gene}_nuclear_spots"] += 1
            if owner:
                lookup[owner][f"{gene}_spots"] += 1
                flags[owner].update(filter(None, spot["qc_flags"].split(";")))
                if "saturated_signal" in spot["qc_flags"].split(";"):
                    flags[owner].add(f"{gene}_saturation")
    for key, mask in (review_masks or {}).items():
        for label in np.unique(assignment_labels[mask]):
            if label:
                flags[int(label)].add(key)
    for n in nuclei:
        p = n["pnoc_spots"] >= config["pnoc"]["positive_min_spots"]
        c = n["cre_spots"] >= config["cre"]["positive_min_spots"]
        n["pnoc_positive"], n["cre_positive"] = p, c
        n["classification"] = ("double_positive" if p and c else "pnoc_only" if p else
                               "cre_only" if c else "neither") if n["accepted"] else "excluded"
        cell_flags = flags[n["nucleus_id"]]
        reasons = []
        for gene in ["pnoc", "cre"]:
            count, cutoff = n[f"{gene}_spots"], config[gene]["positive_min_spots"]
            if abs(count - cutoff) <= 1:
                reasons.append(f"near_{gene}_cutoff")
                cell_flags.add("near_expression_threshold")
            if count < cutoff and f"{gene}_saturation" in cell_flags:
                reasons.append(f"saturated_{gene}_below_cutoff")
        cell_flags.update(reasons)
        n["qc_flags"] = ";".join(sorted(cell_flags))
        n["review_reasons"] = ";".join(reasons)
        n["needs_review"] = bool(n["accepted"] and reasons)
        n["review_priority"] = ("cre_endpoint" if c else "general") if n["needs_review"] else "none"
        for flag in QC_FLAG_NAMES:
            n[f"flag_{flag}"] = flag in cell_flags


def review_sort_key(nucleus):
    """Triage observed Cre-positive objects first without changing their calls."""
    tier = {"cre_endpoint": 0, "general": 1, "none": 2}[nucleus["review_priority"]]
    return tier, -nucleus["pnoc_spots"] - nucleus["cre_spots"], nucleus["nucleus_id"]


def summarize(nuclei):
    accepted = [n for n in nuclei if n["accepted"]]
    result = {"total_nuclei": len(accepted), "excluded_nuclei": len(nuclei) - len(accepted),
              "pnoc_positive": sum(n["pnoc_positive"] for n in accepted),
              "cre_positive": sum(n["cre_positive"] for n in accepted),
              "nuclei_needing_review": sum(n["needs_review"] for n in accepted),
              "nuclei_cre_endpoint_review": sum(n["review_priority"] == "cre_endpoint" for n in accepted)}
    for group in ["double_positive", "pnoc_only", "cre_only", "neither"]:
        result[group] = sum(n["classification"] == group for n in accepted)
    return add_percentages(result)


def add_percentages(result):
    for key, denominator in [("pct_pnoc_expressing_cre", "pnoc_positive"),
                             ("pct_cre_expressing_pnoc", "cre_positive")]:
        result[key] = 100.0 * result["double_positive"] / result[denominator] if result[denominator] else None
    return result


def aggregate(summaries):
    return add_percentages({k: sum(s[k] for s in summaries if s["status"] == "ok") for k in COUNT_FIELDS})


def write_csv(path, records, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        for record in records:
            writer.writerow({k: "NA" if record.get(k) is None else record[k] for k in fields})


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def display_scale(image):
    lo, hi = np.percentile(image, [1, 99.7])
    return np.clip((image.astype(np.float32) - lo) / max(float(hi - lo), 1), 0, 1)


def write_imagej_labels(path, labels, spacing):
    """Lossless Fiji companion: ImageJ accepts uint16 or float32, not uint32."""
    maximum = int(labels.max(initial=0))
    if maximum > 2 ** 24:
        logging.getLogger("rnascope").warning("ImageJ companion omitted: label IDs exceed exact float32 range; native uint32 labels preserved")
        return {"written": False, "reason": "IDs exceed exact ImageJ float32 range"}
    data = labels.astype(np.uint16 if maximum <= 65535 else np.float32)
    tifffile.imwrite(path, data, imagej=True, photometric="minisblack", resolutionunit="NONE",
                     resolution=(1 / spacing[1], 1 / spacing[0]), metadata={"unit": "um", "axes": "YX"})
    return {"written": True, "dtype": str(data.dtype), "unit": "um"}


def save_annotated_overlay(out, rgb, boundaries, nuclei, spots):
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    height, width = rgb.shape[:2]
    scale = 2 if max(height, width) < 2048 else 1
    fig = plt.figure(figsize=(width / 100, height / 100), dpi=100 * scale)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, interpolation="nearest")
    ax.imshow(boundaries, interpolation="nearest")
    for gene, color in [("pnoc", "lime"), ("cre", "#ff5588")]:
        selected = [s for s in spots[gene] if s["counted"]]
        if selected:
            ax.scatter([s["x_px"] for s in selected], [s["y_px"] for s in selected],
                       s=4, facecolors="none", edgecolors=color, linewidths=0.3)
    for n in nuclei:
        if n["accepted"]:
            ax.text(n["x_px"], n["y_px"], str(n["nucleus_id"]), color=COLORS[n["classification"]],
                    fontsize=6 / scale, ha="center", va="center",
                    path_effects=[patheffects.withStroke(linewidth=0.65, foreground="black")])
    ax.set_xlim(-0.5, width - 0.5); ax.set_ylim(height - 0.5, -0.5); ax.set_axis_off()
    fig.savefig(out / "annotated_overlay.png", dpi=100 * scale)
    fig.savefig(out / "annotated_overlay.pdf", dpi=100 * scale)
    plt.close(fig)
    write_json(out / "annotated_overlay_legend.json", {"source_size_yx": [height, width], "render_scale": scale,
               "IDs": "accepted nucleus_id values in nuclei.csv", "label_colors": COLORS,
               "channels": {"red": "Cre", "green": "Pnoc", "blue": "DAPI"},
               "boundaries": {"white": "nuclei", "cyan": "expanded assignment, when enabled"}})


def save_qc(out, image, channels, labels, nuclei, spots, spacing, cy3_index=None, crop_count=12, crop_um=45,
            assignment_labels=None, expand_um=0):
    # Import only when rendering; matplotlib's Agg backend supports headless runs.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.colors import to_rgb

    rgb = np.zeros((*labels.shape, 3), np.float32)
    rgb[..., 0] = display_scale(image[channels["cre"]])
    rgb[..., 1] = display_scale(image[channels["pnoc"]])
    rgb[..., 2] = display_scale(image[channels["dapi"]])
    palette = np.zeros((int(labels.max()) + 1, 3), np.float32)
    for n in nuclei:
        palette[n["nucleus_id"]] = to_rgb(COLORS[n["classification"]])
    class_rgb = palette[labels]
    boundary = segmentation.find_boundaries(labels, mode="inner")
    overlay = np.zeros((*labels.shape, 4), np.float32)
    overlay[boundary] = [1, 1, 1, 0.65]
    if expand_um:
        assignment_boundary = segmentation.find_boundaries(assignment_labels, mode="inner")
        overlay[assignment_boundary] = [0, 1, 1, 0.7]
    save_annotated_overlay(out, rgb, overlay, nuclei, spots)

    def draw(ax, extent=None):
        ax.imshow(rgb)
        ax.imshow(overlay)
        for gene, color in [("pnoc", "lime"), ("cre", "#ff5588")]:
            selected = [s for s in spots[gene] if s["counted"] and
                        (extent is None or extent[0] <= s["x_px"] <= extent[1] and extent[2] <= s["y_px"] <= extent[3])]
            if selected:
                ax.scatter([s["x_px"] for s in selected], [s["y_px"] for s in selected],
                           s=18 if extent else 5, facecolors="none", edgecolors=color, linewidths=0.6)
        if extent:
            ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[3], extent[2])
        ax.set_axis_off()

    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    axs[0, 0].imshow(image[channels["dapi"]], cmap="gray",
                     vmin=np.percentile(image[channels["dapi"]], 1),
                     vmax=np.percentile(image[channels["dapi"]], 99.7))
    boundary_rgb = np.zeros((*labels.shape, 4)); boundary_rgb[boundary] = [1, 0.7, 0, 0.8]
    axs[0, 0].imshow(boundary_rgb); axs[0, 0].set_title("DAPI and candidate nuclear boundaries")
    for ax, gene in [(axs[0, 1], "pnoc"), (axs[0, 2], "cre")]:
        ax.imshow(display_scale(image[channels[gene]]), cmap="gray", vmin=0, vmax=1)
        for counted, color in [(False, "#888888"), (True, "lime" if gene == "pnoc" else "#ff5588")]:
            selected = [s for s in spots[gene] if s["counted"] == counted]
            if selected:
                ax.scatter([s["x_px"] for s in selected], [s["y_px"] for s in selected],
                           s=5, facecolors="none", edgecolors=color, linewidths=0.4)
        ax.set_title(f"{gene}: detected spots (gray = not counted)")
    draw(axs[1, 0]); axs[1, 0].set_title("Cre red / Pnoc green / DAPI blue + counted spots")
    axs[1, 1].imshow(class_rgb); axs[1, 1].set_title("Provisional nuclear classifications")
    axs[1, 1].legend(handles=[Patch(color=c, label=k.replace("_", " ")) for k, c in COLORS.items()],
                     fontsize=7, loc="upper right")
    if cy3_index is not None:
        axs[1, 2].imshow(display_scale(image[cy3_index]), cmap="gray", vmin=0, vmax=1)
        axs[1, 2].set_title("Cy3 reference: inspect shared bright artifacts")
    else:
        axs[1, 2].text(0.1, 0.5, "No Cy3 reference channel identified")
    for ax in axs.flat:
        ax.set_axis_off()
    region = f"Nuclei + {expand_um:g} µm (cyan assignment / white nuclear boundary)" if expand_um else "Nuclear-associated RNA only"
    fig.suptitle(region + " | display contrast adjusted; analysis uses original intensities")
    fig.tight_layout(); fig.savefig(out / "qc_overview.png", dpi=160); plt.close(fig)

    candidates = [n for n in nuclei if n["accepted"]]
    # Include each class, then prioritized ambiguous objects; deterministic selection.
    chosen = []
    ranked = sorted(candidates, key=review_sort_key)
    for group in ["double_positive", "pnoc_only", "cre_only", "neither"]:
        representative = next((n for n in ranked if n["classification"] == group), None)
        if representative is not None and len(chosen) < crop_count:
            chosen.append(representative)
    chosen.extend(n for n in ranked if n not in chosen)
    chosen = sorted(chosen[:crop_count], key=review_sort_key)
    if chosen:
        columns = min(4, len(chosen)); rows = math.ceil(len(chosen) / columns)
        fig, axs = plt.subplots(rows, columns, figsize=(4 * columns, 4.4 * rows), squeeze=False)
        for ax, n in zip(axs.flat, chosen):
            x, y = n["x_px"], n["y_px"]
            extent = [max(0, x - crop_um / spacing[1] / 2), min(labels.shape[1] - 1, x + crop_um / spacing[1] / 2),
                      max(0, y - crop_um / spacing[0] / 2), min(labels.shape[0] - 1, y + crop_um / spacing[0] / 2)]
            draw(ax, extent)
            ax.text(x, y, str(n["nucleus_id"]), color="yellow", fontsize=9, ha="center")
            ax.set_title(f"ID {n['nucleus_id']} | Pnoc {n['pnoc_spots']} / Cre {n['cre_spots']}\n"
                         f"{n['classification']} | review: {n['review_priority']}", fontsize=9)
        for ax in list(axs.flat)[len(chosen):]:
            ax.set_axis_off()
        fig.tight_layout(); fig.savefig(out / "review_crops.png", dpi=160); plt.close(fig)
        write_json(out / "review_crop_ids.json", [n["nucleus_id"] for n in chosen])
    save_cre_only_audit(out, image, channels, labels, nuclei, spots, spacing,
                        assignment_labels if assignment_labels is not None else labels, expand_um, crop_um)


def save_cre_only_audit(out, image, channels, labels, nuclei, spots, spacing, assignments, expand_um, crop_um):
    """Every accepted Cre+/Pnoc- object, not a selected sample of disagreements."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    cells = [n for n in nuclei if n["classification"] == "cre_only"]
    write_csv(out / "cre_only_review.csv", cells, NUCLEUS_FIELDS)
    if not cells:
        return
    scaled = {g: display_scale(image[c]) for g, c in channels.items()}
    rgb = np.stack([scaled["cre"], scaled["pnoc"], scaled["dapi"]], axis=-1)
    with PdfPages(out / "cre_only_review.pdf") as pdf:
        for page, start in enumerate(range(0, len(cells), 4), 1):
            group = cells[start:start + 4]
            fig, axs = plt.subplots(len(group), 4, figsize=(12, 3.2 * len(group)), squeeze=False)
            for row, n in enumerate(group):
                x, y = n["x_px"], n["y_px"]
                y0, y1 = max(0, int(y - crop_um / spacing[0] / 2)), min(labels.shape[0], int(y + crop_um / spacing[0] / 2) + 1)
                x0, x1 = max(0, int(x - crop_um / spacing[1] / 2)), min(labels.shape[1], int(x + crop_um / spacing[1] / 2) + 1)
                sl = np.s_[y0:y1, x0:x1]
                for col, gene in enumerate(["dapi", "pnoc", "cre", "merge"]):
                    ax = axs[row, col]
                    ax.imshow(rgb[sl] if gene == "merge" else scaled[gene][sl], cmap="gray", vmin=0, vmax=1,
                              interpolation="nearest")
                    for lab, color in [(labels, "yellow"), (assignments, "cyan")][:2 if expand_um else 1]:
                        mask = lab[sl] == n["nucleus_id"]
                        if mask.any() and not mask.all():
                            ax.contour(mask, levels=[0.5], colors=[color], linewidths=0.7)
                    if gene in ["pnoc", "cre"]:
                        selected = [s for s in spots[gene] if x0 <= s["x_px"] < x1 and y0 <= s["y_px"] < y1]
                        for assigned, color in [(False, "#aaaaaa"), (True, "lime" if gene == "pnoc" else "#ff6688")]:
                            candidates = [s for s in selected if (s["assignment_id"] == n["nucleus_id"]) == assigned]
                            if candidates:
                                ax.scatter([s["x_px"] - x0 for s in candidates], [s["y_px"] - y0 for s in candidates],
                                           s=16, facecolors="none", edgecolors=color, linewidths=0.6)
                    title = gene.upper() if gene == "dapi" else gene
                    if gene in ["pnoc", "cre"]:
                        title += f": {n[gene + '_spots']} assigned / {n[gene + '_nuclear_spots']} nuclear"
                    ax.set_title(f"ID {n['nucleus_id']} | {title}", fontsize=8)
                    ax.axis("off")
            fig.suptitle(f"All Cre+/Pnoc- objects, page {page} | yellow nucleus / cyan +{expand_um:g} µm region\n"
                         "Review segmentation, background peaks, saturation and near-threshold counts", fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig); fig.savefig(out / f"cre_only_review_{page:02d}.png", dpi=140); plt.close(fig)


def analyze_scene(image, spacing, config, channel_metadata, out, source, scene):
    channels, warnings = resolve_channels(channel_metadata, image.shape[0], config["channels"])
    if config["pixel_size_um"] is not None:
        warnings.append("Pixel calibration overridden by configuration")
    spacing = tuple(spacing)
    if len(spacing) != 2 or not all(np.isfinite(s) and s > 0 for s in spacing):
        raise ValueError("Missing or invalid pixel calibration; set pixel_size_um to [y_um, x_um]")
    dapi = image[channels["dapi"]]
    labels, nuclei, segmentation_info = segment_nuclei(dapi, spacing, config["nuclei"])
    expand_um = config["assignment"]["expand_um"]
    assignments = build_assignment_labels(labels, nuclei, spacing, expand_um)
    spots, detection_info, masks = {}, {}, {}
    for gene in ["pnoc", "cre"]:
        raw = image[channels[gene]]
        bits = channel_metadata[channels[gene]].get("component_bit_count") if channel_metadata else None
        level = saturation_level(raw, config["qc"]["saturation_value"], bits)
        spots[gene], detection_info[gene], mask = detect_spots(raw, spacing, config["spots"][gene], level)
        masks[f"{gene}_saturated_or_clustered_region"] = mask
        if level is not None:
            masks[f"{gene}_saturation"] = raw >= level
        if detection_info[gene]["saturated_pixels"]:
            warnings.append(f"{gene}: saturated pixels present; inspect clusters")
    bits = channel_metadata[channels["dapi"]].get("component_bit_count") if channel_metadata else None
    level = saturation_level(dapi, config["qc"]["saturation_value"], bits)
    if level is not None:
        masks["dapi_saturation"] = dapi >= level
    assign_and_classify(labels, nuclei, spots, config["spots"], masks, assignments)
    summary = summarize(nuclei)
    if not summary["total_nuclei"]:
        warnings.append("No accepted nuclei; inspect DAPI segmentation and size limits")
    summary.update(file=source, scene=str(scene), status="ok", assignment_expand_um=expand_um,
                   analysis_profile=describe_profile(config),
                   warnings="; ".join(warnings), error="")
    for n in nuclei:
        n.update(file=source, scene=str(scene))
    for gene in spots:
        for s in spots[gene]:
            s.update(file=source, scene=str(scene), gene=gene)
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / "nuclei.csv", nuclei, NUCLEUS_FIELDS)
    write_csv(out / "review_queue.csv", sorted((n for n in nuclei if n["needs_review"]), key=review_sort_key), NUCLEUS_FIELDS)
    write_csv(out / "spots.csv", [s for gene in spots.values() for s in gene], SPOT_FIELDS)
    write_json(out / "summary.json", summary)
    scene_metadata = {"pixel_size_yx_um": spacing, "channels": channel_metadata,
                                       "resolved_channels": channels, "segmentation": segmentation_info,
                                       "spot_detection": detection_info,
                                       "assignment": {"expand_um": expand_um, "blocked_excluded_nuclei": True},
                      "review_policy": "accepted objects within +/-1 peak of either cutoff, or saturated RNA below that gene cutoff; other flags informational",
                      "review_priority_policy": "cre_endpoint: queued Cre-positive objects; general: other queued objects; none: not queued; triage does not alter expression calls",
                      "coordinate_origin": "top-left of reconstructed scene; x=column, y=row"}
    tifffile.imwrite(out / "nuclei_labels.tif", labels.astype(np.uint32), photometric="minisblack",
                     resolution=(10000 / spacing[1], 10000 / spacing[0]), resolutionunit="CENTIMETER",
                     metadata={"axes": "YX", "unit": "um", "pixel_size_yx_um": spacing})
    scene_metadata["imagej_labels"] = {"nuclei": write_imagej_labels(out / "nuclei_labels_imagej.tif", labels, spacing)}
    if expand_um:
        tifffile.imwrite(out / "assignment_labels.tif", assignments.astype(np.uint32), photometric="minisblack",
                         resolution=(10000 / spacing[1], 10000 / spacing[0]), resolutionunit="CENTIMETER",
                         metadata={"axes": "YX", "unit": "um", "expand_um": expand_um})
        scene_metadata["imagej_labels"]["assignment"] = write_imagej_labels(out / "assignment_labels_imagej.tif", assignments, spacing)
    write_json(out / "metadata.json", scene_metadata)
    cy3 = next((i for i, m in enumerate(channel_metadata)
                if normalize_name(m["fluor"]) == "cy3" or normalize_name(m["name"]) == "cy3t2"), None)
    save_qc(out, image, channels, labels, nuclei, spots, spacing, cy3,
            config["qc"]["review_crops"], config["qc"]["crop_width_um"], assignments, expand_um)
    return summary


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_batch(input_dir, output_dir, config, recursive=False):
    input_dir, output_dir = Path(input_dir).resolve(), Path(output_dir).resolve()
    if not input_dir.is_dir():
        raise ValueError(f"Input folder does not exist: {input_dir}")
    files = sorted(p for p in (input_dir.rglob("*") if recursive else input_dir.iterdir())
                   if p.is_file() and p.suffix.lower() == ".czi")
    if not files:
        raise ValueError(f"No CZI files found in {input_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Output folder is not empty; choose a new folder to preserve previous results")
    profile = describe_profile(config)
    print(f"Active profile: {profile}", flush=True)
    preflight(config)
    output_dir.mkdir(parents=True, exist_ok=True)
    log = logging.getLogger("rnascope")
    log.setLevel(logging.INFO)
    handler = logging.FileHandler(output_dir / "run.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    summaries, inputs = [], []
    write_json(output_dir / "config_used.json", config)
    started_at = datetime.now(timezone.utc).isoformat()
    packages = {name: package_version(name) for name in
                ["numpy", "scipy", "scikit-image", "czifile", "matplotlib", "tifffile"]}
    try:
        log.info("Active profile: %s", profile)
        for path in files:
            relative = path.relative_to(input_dir).as_posix()
            log.info("Reading %s", relative)
            print(f"Analyzing {relative}", flush=True)
            try:
                inputs.append({"file": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
                with czifile.CziFile(path) as czi:
                    metadata, calibration = metadata_from_xml(czi.metadata())
                    if not len(czi.scenes):
                        raise ValueError("No image scenes found")
                    for scene_id, scene in czi.scenes.items():
                        try:
                            # Check BEFORE decoding potentially large stacks.
                            if any(size != 1 for dim, size in scene.sizes.items() if dim not in "CYX"):
                                raise ValueError(f"Only single sections supported; scene dimensions: {scene.sizes}")
                            spacing = config["pixel_size_um"] or [calibration.get("Y"), calibration.get("X")]
                            if any(s is None for s in spacing):
                                raise ValueError("Missing pixel size; supply pixel_size_um in configuration")
                            image = as_cyx(scene.asarray(maxworkers=2), scene.dims)
                            dest = output_dir / "images" / relative / f"scene_{scene_id}"
                            result = analyze_scene(image, spacing, config, metadata, dest, relative, scene_id)
                            summaries.append(result)
                            log.info("%s scene %s: %s", relative, scene_id, json.dumps(result))
                            print(f"  Scene {scene_id}: {result['total_nuclei']} nuclei; "
                                  f"Pnoc+ {result['pnoc_positive']}, Cre+ {result['cre_positive']}, "
                                  f"double+ {result['double_positive']} (provisional)", flush=True)
                        except Exception as exc:
                            log.exception("Scene failed: %s / %s", relative, scene_id)
                            summaries.append({"file": relative, "scene": str(scene_id), "status": "failed", "error": str(exc),
                                              "analysis_profile": profile, "assignment_expand_um": config["assignment"]["expand_um"]})
                            print(f"  FAILED scene {scene_id}: {exc}", file=sys.stderr)
            except Exception as exc:
                log.exception("File failed: %s", relative)
                summaries.append({"file": relative, "scene": "", "status": "failed", "error": str(exc),
                                  "analysis_profile": profile, "assignment_expand_um": config["assignment"]["expand_um"]})
                print(f"  FAILED {relative}: {exc}", file=sys.stderr)
            write_csv(output_dir / "summary.csv", summaries, SUMMARY_FIELDS)
        total = aggregate(summaries)
        failed = sum(s["status"] == "failed" for s in summaries)
        total.update(successful_scenes=sum(s["status"] == "ok" for s in summaries), failed_items=failed,
                     assignment_expand_um=config["assignment"]["expand_um"],
                     analysis_profile=profile,
                     primary_endpoint="pct_cre_expressing_pnoc",
                     complete=not failed, interpretation=f"Provisional RNA peak counts; nucleus expansion {config['assignment']['expand_um']:g} um; review QC")
        write_json(output_dir / "batch_summary.json", total)
        write_csv(output_dir / "batch_summary.csv", [total], [*COUNT_FIELDS, "pct_pnoc_expressing_cre",
                  "pct_cre_expressing_pnoc", "assignment_expand_um", "analysis_profile", "successful_scenes", "failed_items", "complete"])
        write_json(output_dir / "run_manifest.json", {"script_version": VERSION, "script_sha256": sha256(__file__),
                   "analysis_profile": profile,
                   "started_at_utc": started_at, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                   "python_version": sys.version, "packages": packages, "inputs": inputs,
                   "config_sha256": sha256(output_dir / "config_used.json")})
        print(f"Results: {output_dir} | {total['successful_scenes']} successful scenes, {failed} failures", flush=True)
        return 1 if failed else 0
    finally:
        handler.close()
        log.removeHandler(handler)


def load_profile(path):
    config = merge_config(json.loads(Path(path).read_text(encoding="utf-8")))
    validate_config(config)
    return config


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, help="Folder containing CZI files")
    parser.add_argument("--output", type=Path, help="New or empty results folder")
    parser.add_argument("--config", type=Path, help="Required for analysis: choose an explicit JSON profile")
    parser.add_argument("--recursive", action="store_true", help="Include CZI files in subfolders")
    parser.add_argument("--write-default-config", type=Path, metavar="PATH", help="Copy the shipped StarDist + 2 um profile and exit")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    if args.write_default_config:
        try:
            config = load_profile(Path(__file__).with_name("settings.json"))
            with args.write_default_config.open("x", encoding="utf-8") as f:
                json.dump(config, f, indent=2); f.write("\n")
        except (ValueError, OSError) as exc:
            parser.error(str(exc))
        return 0
    if args.input is None or args.output is None:
        parser.error("--input and --output are required")
    if args.config is None:
        parser.error("--config is required: use settings_nucleus_only.json for nuclear overlap "
                     "or settings.json for nuclei plus 2 um; no method is selected silently")
    try:
        config = load_profile(args.config)
        return run_batch(args.input, args.output, config, args.recursive)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    sys.exit(main())
