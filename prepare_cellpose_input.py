#!/usr/bin/env python3
"""Export one verified reference CZI scene for isolated Cellpose inference."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage import exposure

import rnascope as rn


def prepare(czi_path, reference, output):
    czi_path, reference, output = map(Path, (czi_path, reference, output))
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output folder must be new or empty")
    config_path = reference.parents[2] / "config_used.json"
    manifest = json.loads((reference.parents[2] / "run_manifest.json").read_text())
    summary = json.loads((reference / "summary.json").read_text())
    metadata = json.loads((reference / "metadata.json").read_text())
    matching = [item for item in manifest["inputs"] if item["file"] == summary["file"]]
    source_hash = rn.sha256(czi_path)
    if len(matching) != 1 or matching[0]["sha256"] != source_hash:
        raise ValueError("CZI hash does not match the reference analysis input")
    if rn.sha256(config_path) != manifest["config_sha256"]:
        raise ValueError("Reference configuration hash mismatch")
    config = rn.load_profile(config_path)
    spacing = metadata["pixel_size_yx_um"]
    if not np.isclose(spacing[0], spacing[1], rtol=.01):
        raise ValueError("This comparison currently requires square XY pixels; anisotropic resampling needs separate validation")
    with rn.czifile.CziFile(czi_path) as czi:
        scene = czi.scenes[int(summary["scene"])]
        if any(size != 1 for dim, size in scene.sizes.items() if dim not in "CYX"):
            raise ValueError("Only single optical sections are supported")
        image = rn.as_cyx(scene.asarray(maxworkers=2), scene.dims)
    # Match stardist_labels preprocessing exactly, without running StarDist.
    from csbdeep.utils import normalize
    n = config["nuclei"]
    dapi = image[metadata["resolved_channels"]["dapi"]]
    smooth = ndi.gaussian_filter(dapi.astype(np.float32), [n["stardist_smooth_um"] / p for p in spacing])
    processed = np.clip(normalize(smooth, 1, 99.8), 0, 1)
    if n["stardist_clahe"]:
        processed = exposure.equalize_adapthist(processed,
                    kernel_size=tuple(max(1, int(round(n["stardist_clahe_kernel_um"] / p))) for p in spacing),
                    clip_limit=n["stardist_clahe_clip"])
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "input_channels.npz", image=image,
                        dapi_stardist_preprocessed=processed.astype(np.float32), pixel_size_yx_um=spacing)
    rn.write_json(output / "input_metadata.json", {
        "source_czi": str(czi_path.resolve()), "czi_sha256": source_hash,
        "channel_indices": metadata["resolved_channels"], "pixel_size_yx_um": spacing,
        "reference_scene": str(reference.resolve()), "export_script_sha256": rn.sha256(__file__),
        "reference_analysis_script_sha256": manifest["script_sha256"],
        "dapi_preprocessing_control": {"smooth_um": n["stardist_smooth_um"], "percentile_normalization": [1, 99.8],
                                       "clahe": n["stardist_clahe"], "clahe_kernel_um": n["stardist_clahe_kernel_um"],
                                       "clahe_clip": n["stardist_clahe_clip"]}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--czi", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.czi, args.reference, args.output)


if __name__ == "__main__":
    main()
