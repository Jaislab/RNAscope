"""Run actual pinned Cellpose and verify main-pipeline exports on synthetic data.

From repository root: python tests/smoke_cellpose.py
Requires CELLPOSE_PYTHON when the isolated environment is not .venv-cellpose.
This is a software smoke test, not a biological validation benchmark.
"""
import json
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import rnascope as rn


def main():
    config = rn.load_profile(Path(rn.__file__).with_name('settings_cellpose.json'))
    rn.preflight(config)
    y, x = np.mgrid[:128, :128]
    centers = [(32, 32), (32, 90), (90, 32), (90, 90)]
    dapi = sum(np.exp(-((y-cy)**2+(x-cx)**2)/(2*6**2)) for cy, cx in centers)
    image = np.zeros((4, 128, 128), np.uint16)
    image[3] = (100 + 8000*dapi).astype(np.uint16)
    puncta = sum(np.exp(-((y-cy-dy)**2+(x-cx-dx)**2)/(2*.8**2))
                 for cy, cx in centers for dy, dx in [(-4,-4),(-4,4),(4,0)])
    image[0] = image[2] = (100 + 4000*puncta).astype(np.uint16)
    metadata = [{'name': n, 'fluor': ''} for n in ['Cy5-T1', 'Cy3-T2', 'At520-T3', 'DAPI-T4']]
    with tempfile.TemporaryDirectory() as directory:
        out = Path(directory) / 'scene'
        summary = rn.analyze_scene(image, (.69, .69), config, metadata, out, 'synthetic.czi', 0)
        assert summary['total_nuclei'] == 4, summary
        assert summary['double_positive'] == 4, summary
        assert summary['assignment_expand_um'] == 2
        assert 'cellpose/nuclei/raw-DAPI' in summary['analysis_profile']
        info = json.loads((out / 'metadata.json').read_text())
        assert info['segmentation']['weights_sha256'] == rn.CELLPOSE_WEIGHTS_SHA256
        assert info['segmentation']['packages']['cellpose'] == '3.1.1.3'
        assert info['segmentation']['diameter_px'] == config['nuclei']['cellpose_diameter_um'] / .69
        nuclei = rn.tifffile.imread(out / 'nuclei_labels.tif')
        expanded = rn.tifffile.imread(out / 'assignment_labels.tif')
        assert np.count_nonzero(expanded) > np.count_nonzero(nuclei)
        for name in ['nuclei.csv', 'spots.csv', 'review_queue.csv', 'annotated_overlay.pdf',
                     'assignment_labels_imagej.tif']:
            assert (out / name).stat().st_size > 0
        print('PASS: actual Cellpose CPU inference, 4/4 nuclei and double-positive objects, '
              '2 um expansion, weight provenance, CSV/TIFF/PDF exports')


if __name__ == '__main__':
    main()
