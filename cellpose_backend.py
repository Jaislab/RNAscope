"""Isolated pinned Cellpose CPU worker. Called by rnascope.py, not a batch entry point."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import numpy as np


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
            # Download atomically; partial or mismatched artifacts are never loaded.
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=folder, delete=False) as stream:
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
        if sha256(weights) != CELLPOSE_WEIGHTS_SHA256:
            raise RuntimeError("Cellpose nuclei weight fingerprint mismatch; refusing to load")
        torch.set_num_threads(threads)
        model = models.CellposeModel(gpu=False, pretrained_model=str(weights), model_type=None)
        info = {"model": "nuclei", "weights_url": CELLPOSE_WEIGHTS_URL,
                "weights_sha256": CELLPOSE_WEIGHTS_SHA256, "gpu": False, "threads": threads,
                "packages": {p: importlib.metadata.version(p) for p in ("cellpose", "torch", "numpy")}}
        _CELLPOSE_MODELS[key] = model, info
    return _CELLPOSE_MODELS[key]



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--diameter', type=float)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('threads must be positive')
    os.environ.setdefault('NUMBA_NUM_THREADS', str(args.threads))
    os.environ.setdefault('OMP_NUM_THREADS', str(args.threads))
    model, info = load_cellpose_model(args.threads)
    if args.input is not None:
        if args.output is None or args.diameter is None or not np.isfinite(args.diameter) or args.diameter <= 0:
            parser.error('inference requires output and a positive diameter')
        image = np.load(args.input, allow_pickle=False)
        if image.ndim != 2:
            parser.error('input must be a 2D DAPI image')
        labels, _, _ = model.eval(image, channels=[0, 0], channel_axis=None,
                                 diameter=args.diameter, normalize={'percentile': [1., 99.]},
                                 flow_threshold=.4, cellprob_threshold=0., min_size=15,
                                 resample=True, augment=False, do_3D=False)
        np.save(args.output, labels.astype(np.int32), allow_pickle=False)
    args.metadata.write_text(json.dumps(info, indent=2) + '\n')


if __name__ == '__main__':
    main()
