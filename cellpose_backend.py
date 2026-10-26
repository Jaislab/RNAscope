"""Isolated pinned Cellpose CPU worker. Called by rnascope.py, not a batch entry point.

The worker applies exactly the ``model.eval`` keyword arguments it receives in a
JSON parameter file and echoes them into its metadata output, so the recorded
provenance in the main analysis is the parameter set that actually ran.
"""
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


_CELLPOSE_MODELS = {}
CELLPOSE_WEIGHTS_URL = "https://www.cellpose.org/static/models/nucleitorch_0"
# Reference fingerprint recorded in docs/CELLPOSE_COMPARISON.md, not a publisher checksum.
CELLPOSE_WEIGHTS_SHA256 = "89ca45e4a45048d5010d29621466b1274b91b3dcb9714dce9f9e90a9a8671303"
# Single definition of the fixed model.eval keyword arguments, shared by the
# batch runner, this worker and the exploratory cellpose_infer.py. Callers add
# "diameter" (and may override "channels"/"channel_axis" for multi-channel input).
CELLPOSE_INFERENCE = {"channels": [0, 0], "channel_axis": None, "normalize": {"percentile": [1.0, 99.0]},
                      "flow_threshold": 0.4, "cellprob_threshold": 0.0, "min_size": 15,
                      "resample": True, "augment": False, "do_3D": False}
INFERENCE_KEYS = frozenset({"diameter", *CELLPOSE_INFERENCE})


def load_inference_parameters(path):
    """Require the complete keyword set; never fill gaps from Cellpose defaults."""
    parameters = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(parameters, dict) or set(parameters) != INFERENCE_KEYS:
        raise ValueError(f"Inference parameters must define exactly: {sorted(INFERENCE_KEYS)}")
    diameter = parameters["diameter"]
    if isinstance(diameter, bool) or not isinstance(diameter, (int, float)) or not math.isfinite(diameter) or diameter <= 0:
        raise ValueError("diameter must be a positive finite number")
    return parameters


def load_cellpose_model(threads):
    """Load only the pinned legacy nuclei model; never substitute another model."""
    folder = Path(os.environ.get("CELLPOSE_LOCAL_MODELS_PATH",
                  str(Path(__file__).resolve().parent / ".cache" / "cellpose")))
    key = (str(folder.resolve()), threads)
    if key not in _CELLPOSE_MODELS:
        if importlib.metadata.version("cellpose") != "3.1.1.3":
            raise RuntimeError("Install requirements-cellpose.txt in an isolated environment; "
                               "Cellpose 3.1.1.3 is required, with no fallback")
        import torch
        from cellpose import models
        from urllib.request import urlopen
        import tempfile
        folder.mkdir(parents=True, exist_ok=True)
        weights = folder / "nucleitorch_0"
        if not weights.exists():
            # The batch runner kills this process on timeout, which skips the cleanup
            # below, so remove partial downloads from earlier attempts first.
            for stale in folder.glob(f"{weights.name}.partial-*"):
                stale.unlink(missing_ok=True)
            # Download atomically; partial or mismatched artifacts are never loaded.
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=folder, prefix=f"{weights.name}.partial-", delete=False) as stream:
                    temporary = Path(stream.name)
                    with urlopen(CELLPOSE_WEIGHTS_URL, timeout=120) as response:
                        while block := response.read(1024 * 1024):
                            stream.write(block)
                if sha256(temporary) != CELLPOSE_WEIGHTS_SHA256:
                    raise RuntimeError("Cellpose nuclei weight fingerprint mismatch")
                temporary.replace(weights)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        digest = sha256(weights)
        if digest != CELLPOSE_WEIGHTS_SHA256:
            raise RuntimeError("Cellpose nuclei weight fingerprint mismatch; refusing to load")
        torch.set_num_threads(threads)
        model = models.CellposeModel(gpu=False, pretrained_model=str(weights), model_type=None)
        # Report the digest of the file actually loaded; the caller compares it with the constant.
        info = {"model": "nuclei", "weights_url": CELLPOSE_WEIGHTS_URL,
                "weights_sha256": digest, "gpu": False, "threads": threads,
                "packages": {p: importlib.metadata.version(p) for p in ("cellpose", "torch", "numpy")}}
        _CELLPOSE_MODELS[key] = model, info
    return _CELLPOSE_MODELS[key]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--parameters', type=Path, help='JSON object holding every model.eval keyword argument')
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error('threads must be positive')
    # Thread limits (OMP/OpenBLAS/MKL/Numba) arrive in the environment from the
    # batch runner; setting them here would be too late for NumPy's BLAS pool.
    import numpy as np
    model, info = load_cellpose_model(args.threads)
    info = {**info, "backend_sha256": sha256(__file__)}
    if args.input is not None:
        if args.output is None or args.parameters is None:
            parser.error('inference requires --output and --parameters')
        try:
            parameters = load_inference_parameters(args.parameters)
        except ValueError as exc:
            parser.error(str(exc))
        image = np.load(args.input, allow_pickle=False)
        if image.ndim != 2:
            parser.error('input must be a 2D DAPI image')
        labels, _, _ = model.eval(image, **parameters)
        np.save(args.output, labels.astype(np.int32), allow_pickle=False)
        # Echo the exact keyword arguments the model received, only after inference succeeded.
        info["inference"] = parameters
        info["input_array"] = {"shape": list(image.shape), "dtype": str(image.dtype)}
    args.metadata.write_text(json.dumps(info, indent=2) + '\n', encoding="utf-8")


if __name__ == '__main__':
    main()
