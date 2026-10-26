"""Main-pipeline routing, calibration, provenance and export regression checks."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
from types import SimpleNamespace
import cellpose_backend as backend

import numpy as np
import rnascope as rn
from test_rnascope import gaussian_spots


def echoing_worker(labels):
    """A worker double that echoes the request through a JSON round trip, like the real backend."""
    worker = Mock()
    worker.run.side_effect = lambda image, parameters: (
        labels, {'model': 'nuclei', 'weights_sha256': rn.CELLPOSE_WEIGHTS_SHA256,
                 'inference': json.loads(json.dumps(parameters))})
    return worker


class CellposePipelineTests(unittest.TestCase):
    def setUp(self):
        self.config = rn.load_profile(Path(rn.__file__).with_name('settings_cellpose.json'))

    def test_exact_pixel_diameter_and_raw_dapi_only(self):
        worker = echoing_worker(np.zeros((40, 40), np.int32))
        image = np.arange(1600, dtype=np.uint16).reshape(40, 40)
        with patch.object(rn, 'load_cellpose_model', return_value=(worker, {'model': 'nuclei'})):
            _, info = rn.cellpose_labels(image, (.5, .5), self.config['nuclei'])
        sent, parameters = worker.run.call_args.args
        self.assertIs(sent, image)
        self.assertEqual(set(parameters), backend.INFERENCE_KEYS)
        self.assertEqual(parameters['channels'], [0, 0])
        self.assertEqual(parameters['diameter'], 17.0)
        self.assertEqual(info['diameter_um'], 8.5)
        self.assertEqual(info['input'], 'raw DAPI only')
        # Scene metadata carries the echoed keyword set, not a separately typed copy.
        self.assertEqual(info['inference'], parameters)
        self.assertEqual(info['inference']['flow_threshold'], backend.CELLPOSE_INFERENCE['flow_threshold'])

    def test_reference_calibration_never_changes_exact_17_pixel_diameter(self):
        worker = echoing_worker(np.zeros((40, 40), np.int32))
        with patch.object(rn, 'load_cellpose_model', return_value=(worker, {})):
            _, info = rn.cellpose_labels(np.ones((40, 40)),
                                        (.690535481770835, .690535481770835),
                                        self.config['nuclei'])
        self.assertEqual(worker.run.call_args.args[1]['diameter'], 17.0)
        self.assertEqual(info['diameter_px'], 17.0)
        self.assertEqual(info['diameter_um'], 17 * .690535481770835)
        # Avoid the 1023x1023 resize triggered by the former rounded um setting.
        self.assertEqual(int(1024 * (17 / info['diameter_px'])), 1024)

    def fake_subprocess(self, mutate=None, fingerprint=rn.CELLPOSE_WEIGHTS_SHA256):
        """Stand in for the isolated interpreter: honor the CLI contract and write worker outputs."""
        def run(command, **kwargs):
            options = dict(zip(command[2::2], command[3::2]))
            info = {'model': 'nuclei', 'weights_sha256': fingerprint}
            if '--input' in options:
                image = np.load(options['--input'], allow_pickle=False)
                parameters = json.loads(Path(options['--parameters']).read_text())
                np.save(options['--output'], np.zeros(image.shape, np.int32), allow_pickle=False)
                info['inference'] = mutate(parameters) if mutate else parameters
            Path(options['--metadata']).write_text(json.dumps(info))
            return SimpleNamespace(returncode=0, stderr='')
        return run

    def test_worker_requires_faithful_parameter_echo_and_fingerprint(self):
        worker = rn.CellposeWorker('/unused/python', 2)
        image = np.ones((8, 8), np.uint16)
        request = rn.cellpose_inference_parameters(self.config['nuclei'])
        with patch.object(rn.subprocess, 'run', self.fake_subprocess()):
            labels, info = worker.run(image, request)
            self.assertEqual(labels.shape, image.shape)
            self.assertEqual(info['inference'], request)

        def drift(parameters):
            return {**parameters, 'flow_threshold': 0.6}
        with patch.object(rn.subprocess, 'run', self.fake_subprocess(mutate=drift)):
            with self.assertRaisesRegex(RuntimeError, 'did not echo'):
                worker.run(image, request)
        with patch.object(rn.subprocess, 'run', self.fake_subprocess(fingerprint='0' * 64)):
            with self.assertRaisesRegex(RuntimeError, 'fingerprint'):
                worker.run()
        with patch.dict(os.environ, {'CELLPOSE_PYTHON': '/fingerprint/test/python'}), \
             patch.object(rn.subprocess, 'run', self.fake_subprocess(fingerprint='0' * 64)):
            with self.assertRaisesRegex(RuntimeError, 'before analysis'):
                rn.preflight(self.config)

    def test_backend_passes_exactly_the_requested_parameters_and_echoes_them(self):
        model = Mock()
        model.eval.return_value = (np.ones((6, 6), np.int64), None, None)
        request = rn.cellpose_inference_parameters(self.config['nuclei'])
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(backend, 'load_cellpose_model',
                          return_value=(model, {'model': 'nuclei', 'weights_sha256': backend.CELLPOSE_WEIGHTS_SHA256})):
            folder = Path(directory)
            np.save(folder / 'input.npy', np.arange(36, dtype=np.uint16).reshape(6, 6), allow_pickle=False)
            (folder / 'parameters.json').write_text(json.dumps(request))
            argv = ['--threads', '2', '--metadata', str(folder / 'metadata.json'), '--input', str(folder / 'input.npy'),
                    '--output', str(folder / 'labels.npy'), '--parameters', str(folder / 'parameters.json')]
            backend.main(argv)
            self.assertEqual(model.eval.call_args.kwargs, request)
            np.testing.assert_array_equal(model.eval.call_args.args[0], np.arange(36).reshape(6, 6))
            info = json.loads((folder / 'metadata.json').read_text())
            self.assertEqual(info['inference'], request)
            self.assertEqual(info['backend_sha256'], rn.sha256(backend.__file__))
            self.assertEqual(np.load(folder / 'labels.npy').dtype, np.int32)
            # Incomplete or unknown keyword sets are refused rather than filled from Cellpose defaults.
            model.eval.reset_mock()
            for bad in [{**request, 'resize': 1}, {k: v for k, v in request.items() if k != 'min_size'},
                        {**request, 'diameter': 0}]:
                (folder / 'parameters.json').write_text(json.dumps(bad))
                with self.assertRaises(SystemExit):
                    backend.main(argv)
            model.eval.assert_not_called()

    def test_run_manifest_records_worker_script_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / 'bad.czi').write_bytes(b'not a CZI')
            rn.run_batch(root, root / 'out', rn.DEFAULT_CONFIG)
            manifest = json.loads((root / 'out' / 'run_manifest.json').read_text())
        self.assertEqual(manifest['cellpose_backend_sha256'], rn.sha256(backend.__file__))
        self.assertEqual(manifest['script_sha256'], rn.sha256(rn.__file__))

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
        for key, value in [('cellpose_diameter_px', 0), ('cellpose_threads', 1.5),
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
