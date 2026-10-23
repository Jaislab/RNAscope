#!/usr/bin/env python3
"""Run legacy Cellpose nuclei/cyto3 as an isolated, exploratory second opinion.

Input is a CYX array in input_channels.npz plus input_metadata.json exported
from the same CZI used for the reference analysis. No RNA detections or
production settings are changed. Use the separate Cellpose environment.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import time

import numpy as np
import tifffile


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def model_input(image, channels, variant, processed_dapi=None):
    """Cellpose channel indices are 1-based; zero requests grayscale."""
    if variant in ("nuclei_dapi", "cyto3_dapi"):
        return image[channels["dapi"]], [0, 0], None
    if variant in ("nuclei_dapi_clahe", "cyto3_dapi_clahe"):
        if processed_dapi is None or processed_dapi.shape != image.shape[1:]:
            raise ValueError("Preprocessing control requires the reference-preprocessed DAPI array")
        return processed_dapi, [0, 0], None
    if variant == "cyto3_pnoc_dapi":
        return np.stack([image[channels["pnoc"]], image[channels["dapi"]]], axis=-1), [1, 2], -1
    raise ValueError(f"Unknown Cellpose variant: {variant}")


def run_inference(input_path, metadata_path, output, nuclear_diameter=17.0, cell_diameter=30.0, threads=2, preprocessing_control=False):
    if any(not np.isfinite(d) or d <= 0 for d in (nuclear_diameter, cell_diameter)) or threads < 1:
        raise ValueError("Diameters and thread count must be positive")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output folder must be new or empty")
    if importlib.metadata.version("cellpose") != "3.1.1.3":
        raise RuntimeError("Use the separate pinned Cellpose 3.1.1.3 environment; newer SAM models are a different comparison")
    os.environ.setdefault("CELLPOSE_LOCAL_MODELS_PATH", str(Path(__file__).resolve().parent / ".cache" / "cellpose"))
    os.environ.setdefault("NUMBA_NUM_THREADS", str(threads))
    os.environ.setdefault("OMP_NUM_THREADS", str(threads))
    import torch
    from cellpose import models
    # The legacy /models route returns HTTP 500; the publisher's static route
    # serves the same named legacy files over normally verified HTTPS.
    models._MODEL_URL = "https://www.cellpose.org/static/models"
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    metadata = json.loads(Path(metadata_path).read_text())
    with np.load(input_path, allow_pickle=False) as archive:
        image = archive["image"]
        processed_dapi = archive["dapi_stardist_preprocessed"] if preprocessing_control else None
    if image.ndim != 3 or any(not 0 <= c < image.shape[0] for c in metadata["channel_indices"].values()):
        raise ValueError("Expected CYX input with valid biological channel indices")
    output.mkdir(parents=True, exist_ok=True)
    manifest = {"started_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_sha256": sha256(input_path), "input_metadata": metadata,
                "script_sha256": sha256(__file__), "threads": threads, "gpu": False,
                "packages": {p: importlib.metadata.version(p) for p in ("cellpose", "torch", "numpy", "scipy", "numba", "tifffile")},
                "normalization": "Cellpose per-channel 1st/99th percentile normalization; raw profiles have no smoothing/CLAHE; *_clahe profiles use reference DAPI preprocessing first",
                "profiles": {}}
    cache = {}
    profiles = [("nuclei_dapi", "nuclei", nuclear_diameter), ("cyto3_dapi", "cyto3", nuclear_diameter),
                ("cyto3_pnoc_dapi", "cyto3", cell_diameter)]
    if preprocessing_control:
        profiles.extend([("nuclei_dapi_clahe", "nuclei", nuclear_diameter), ("cyto3_dapi_clahe", "cyto3", nuclear_diameter)])
    for variant, model_name, diameter in profiles:
        if model_name not in cache:
            weight_path = Path(models.model_path(model_name))
            expected = "nucleitorch_0" if model_name == "nuclei" else "cyto3"
            if weight_path.name != expected:
                raise RuntimeError("Unexpected model path; refusing an implicit fallback")
            cache[model_name] = (models.CellposeModel(gpu=False, model_type=model_name), weight_path)
        model, weights = cache[model_name]
        x, channel_spec, channel_axis = model_input(image, metadata["channel_indices"], variant, processed_dapi)
        started = time.monotonic()
        print(f"Running {variant}: {model_name}, diameter {diameter:g} px, channels {channel_spec}", flush=True)
        masks, _, _ = model.eval(x, channels=channel_spec, channel_axis=channel_axis,
                                 diameter=diameter, normalize={"percentile": [1., 99.]},
                                 flow_threshold=.4, cellprob_threshold=0., min_size=15,
                                 resample=True, augment=False, do_3D=False)
        spacing = metadata["pixel_size_yx_um"]
        path = output / f"{variant}_labels.tif"
        tifffile.imwrite(path, masks.astype(np.uint32), photometric="minisblack",
                         resolution=(10000 / spacing[1], 10000 / spacing[0]), resolutionunit="CENTIMETER",
                         metadata={"axes": "YX", "unit": "um", "pixel_size_yx_um": spacing})
        manifest["profiles"][variant] = {"model": model_name, "diameter_px": diameter,
                                        "channels": channel_spec, "channel_axis": channel_axis,
                                        "reference_dapi_preprocessing": variant.endswith("_clahe"),
                                        "flow_threshold": .4, "cellprob_threshold": 0., "min_size_px": 15,
                                        "candidates": int(np.unique(masks[masks > 0]).size),
                                        "elapsed_seconds": time.monotonic() - started,
                                        "weights_file": weights.name, "weights_sha256": sha256(weights),
                                        "weights_url": f"{models._MODEL_URL}/{weights.name}",
                                        "mask_sha256": sha256(path)}
        print(f"{variant}: {manifest['profiles'][variant]['candidates']} masks", flush=True)
        manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        (output / "inference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nuclear-diameter-px", type=float, default=17.)
    parser.add_argument("--cell-diameter-px", type=float, default=30.)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--preprocessing-control", action="store_true", help="Also run both DAPI models with the reference smoothing/CLAHE")
    args = parser.parse_args(argv)
    run_inference(args.input, args.metadata, args.output, args.nuclear_diameter_px, args.cell_diameter_px, args.threads, args.preprocessing_control)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
