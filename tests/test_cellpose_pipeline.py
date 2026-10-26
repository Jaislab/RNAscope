"""Main-pipeline routing, calibration, provenance and export regression checks."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
from types import SimpleNamespace
import cellpose_backend as backend

import numpy as np
import rnascope as rn
from test_rnascope import gaussian_spots


class CellposePipelineTests(unittest.TestCase):
    def setUp(self):
        self.config = rn.load_profile(Path(rn.__file__).with_name('settings_cellpose.json'))

    def test_diameter_physical_units_and_raw_dapi_only(self):
        model = Mock()
        image = np.arange(1600, dtype=np.uint16).reshape(40, 40)
        model.eval.return_value = (np.zeros((40, 40), np.int32), None, None)
        with patch.object(rn, 'load_cellpose_model', return_value=(model, {'model': 'nuclei'})):
            _, info = rn.cellpose_labels(image, (.5, .5), self.config['nuclei'])
        args, kwargs = model.eval.call_args
        self.assertIs(args[0], image)
        self.assertEqual(kwargs['channels'], [0, 0])
        self.assertAlmostEqual(kwargs['diameter'], 23.48)
        self.assertEqual(info['diameter_um'], 11.74)
        self.assertEqual(info['input'], 'raw DAPI only')

    def test_anisotropic_sampling_rejected_before_loading(self):
        with patch.object(rn, 'load_cellpose_model') as loader:
            with self.assertRaisesRegex(ValueError, 'square XY'):
                rn.cellpose_labels(np.ones((10, 10)), (.5, 1), self.config['nuclei'])
            loader.assert_not_called()

    def test_unavailable_model_never_falls_back(self):
        image = np.arange(1600).reshape(40, 40)
        with patch.object(rn, 'load_cellpose_model', side_effect=RuntimeError('missing model')), \
             patch.object(rn, 'local_watershed_labels') as fallback:
            with self.assertRaisesRegex(RuntimeError, 'missing model'):
                rn.segment_nuclei(image, (.5, .5), self.config['nuclei'])
            fallback.assert_not_called()
        with patch.object(rn, 'load_cellpose_model', side_effect=RuntimeError('missing model')):
            with self.assertRaisesRegex(RuntimeError, 'before analysis'):
                rn.preflight(self.config)

    def test_worker_failure_is_reported(self):
        worker = rn.CellposeWorker('/does/not/exist/python', 2)
        with self.assertRaisesRegex(RuntimeError, 'Cannot run isolated'):
            worker.run()

    def test_backend_rejects_wrong_version(self):
        with patch.object(backend.importlib.metadata, 'version', return_value='4.0'):
            with self.assertRaisesRegex(RuntimeError, '3.1.1.3'):
                backend.load_cellpose_model(2)

    def test_backend_refuses_corrupted_cached_weights_before_model_load(self):
        import sys
        model_module = Mock()
        with tempfile.TemporaryDirectory() as directory, \
             patch.dict(backend.os.environ, {'CELLPOSE_LOCAL_MODELS_PATH': directory}), \
             patch.object(backend.importlib.metadata, 'version', return_value='3.1.1.3'), \
             patch.dict(sys.modules, {'torch': Mock(), 'cellpose': SimpleNamespace(models=model_module)}):
            (Path(directory) / 'nucleitorch_0').write_bytes(b'corrupted')
            with self.assertRaisesRegex(RuntimeError, 'fingerprint mismatch'):
                backend.load_cellpose_model(2)
            model_module.CellposeModel.assert_not_called()

    def test_invalid_parameters(self):
        for key, value in [('cellpose_diameter_um', 0), ('cellpose_threads', 1.5),
                           ('cellpose_threads', True)]:
            config = copy.deepcopy(self.config)
            config['nuclei'][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                rn.validate_config(config)

    def test_cellpose_exports_expanded_counts_and_model_identity(self):
        y, x = np.indices((100, 100))
        labels = np.zeros((100, 100), np.int32)
        labels[(y - 50)**2 + (x - 50)**2 < 10**2] = 1
        image = np.zeros((4, 100, 100), np.uint16)
        image[3] = 100 + 5000 * (labels > 0)
        # All three RNA spots are outside the nucleus but inside its 2 um region.
        positions = [(50, 62), (44, 60), (56, 60)]
        image[0] = image[2] = gaussian_spots((100, 100), positions, sigma=.6)
        metadata = [{'name': n, 'fluor': ''} for n in ['Cy5-T1', 'Cy3-T2', 'At520-T3', 'DAPI-T4']]
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(rn, 'cellpose_labels', return_value=(labels, {'method': 'cellpose', 'model': 'nuclei',
                              'weights_sha256': rn.CELLPOSE_WEIGHTS_SHA256})):
            out = Path(directory) / 'scene'
            summary = rn.analyze_scene(image, (.5, .5), self.config, metadata, out, 'known.czi', 0)
            self.assertEqual(summary['double_positive'], 1)
            self.assertIn('cellpose/nuclei/raw-DAPI', summary['analysis_profile'])
            import csv
            with (out / 'nuclei.csv').open() as stream:
                nucleus = next(csv.DictReader(stream))
            self.assertEqual(int(nucleus['pnoc_nuclear_spots']), 0)
            self.assertEqual(int(nucleus['pnoc_spots']), 3)
            info = json.loads((out / 'metadata.json').read_text())
            self.assertEqual(info['segmentation']['weights_sha256'], rn.CELLPOSE_WEIGHTS_SHA256)
            self.assertEqual(info['assignment']['expand_um'], 2)
            self.assertTrue((out / 'assignment_labels_imagej.tif').is_file())


if __name__ == '__main__':
    unittest.main()
