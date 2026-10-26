"""Behavioral tests with known spot locations and classifications."""
import copy
from contextlib import redirect_stderr
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

import rnascope as rn


def gaussian_spots(shape, positions, sigma=1.2, amplitude=2000, noise=0):
    y, x = np.indices(shape)
    a = np.full(shape, 100.0)
    for py, px in positions:
        a += amplitude * np.exp(-((y - py) ** 2 + (x - px) ** 2) / (2 * sigma ** 2))
    if noise:
        a += np.random.default_rng(42).normal(0, noise, shape)
    return np.clip(a, 0, 65535).astype(np.uint16)


class ConfigurationTests(unittest.TestCase):
    def test_channel_bit_depth_metadata_and_saturation_override(self):
        xml = '''<ImageDocument><Metadata><Information><Image><ComponentBitCount>14</ComponentBitCount>
        <Dimensions><Channels><Channel Name="Cy5"><ComponentBitCount>12</ComponentBitCount></Channel>
        <Channel Name="DAPI"/></Channels></Dimensions></Image></Information></Metadata></ImageDocument>'''
        channels, _ = rn.metadata_from_xml(xml)
        self.assertEqual([c["component_bit_count"] for c in channels], [12, 14])
        raw = np.array([[0, 4095]], np.uint16)
        self.assertEqual(rn.saturation_level(raw, None, 12), 4095)
        self.assertEqual(rn.saturation_level(raw, None, 14), 16383)
        self.assertEqual(rn.saturation_level(raw, None), 65535)
        self.assertEqual(rn.saturation_level(raw, 4000, 12), 4000)
        self.assertEqual(rn.saturation_level(raw, None, 32), 65535)

    def test_scaling_metadata_without_a_usable_value_is_left_unset(self):
        xml = """<ImageDocument><Metadata><Scaling><Items>
        <Distance Id="X"><Value></Value></Distance><Distance Id="Y"><Value>5e-07</Value></Distance>
        </Items></Scaling></Metadata></ImageDocument>"""
        _, sizes = rn.metadata_from_xml(xml)
        self.assertEqual(sizes, {"Y": 0.5})
        _, sizes = rn.metadata_from_xml(xml.replace("<Value></Value>", "<Value>-1</Value>"))
        self.assertEqual(sizes, {"Y": 0.5})
        _, sizes = rn.metadata_from_xml(xml.replace("<Value></Value>", "<Value>2.5e-07</Value>"))
        self.assertEqual(sizes, {"X": 0.25, "Y": 0.5})

    def test_no_config_cannot_silently_select_a_method(self):
        error = io.StringIO()
        with patch.object(rn, "run_batch") as run, redirect_stderr(error), self.assertRaises(SystemExit) as exc:
            rn.main(["--input", "/input", "--output", "/output"])
        self.assertEqual(exc.exception.code, 2)
        self.assertIn("--config is required", error.getvalue())
        run.assert_not_called()

    def test_all_shipped_profiles_validate(self):
        root = Path(rn.__file__).parent
        for path in root.glob("settings*.json"):
            with self.subTest(profile=path.name):
                rn.load_profile(path)
        production = rn.load_profile(root / "settings.json")
        self.assertEqual(production["nuclei"]["method"], "stardist")
        self.assertEqual(production["assignment"]["expand_um"], 2)
        self.assertEqual(rn.load_profile(root / "settings_nucleus_only.json")["assignment"]["expand_um"], 0)

    def test_exported_default_is_the_shipped_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.json"
            self.assertEqual(rn.main(["--write-default-config", str(path)]), 0)
            self.assertEqual(rn.load_profile(path), rn.load_profile(Path(rn.__file__).with_name("settings_cellpose.json")))

    def test_partial_settings_and_typos(self):
        config = rn.merge_config({"spots": {"pnoc": {"positive_min_spots": 4}}})
        rn.validate_config(config)
        self.assertEqual(config["spots"]["cre"]["positive_min_spots"], 3)
        self.assertEqual(rn.DEFAULT_CONFIG["spots"]["pnoc"]["positive_min_spots"], 3)
        with self.assertRaises(ValueError):
            rn.merge_config({"spots": {"pnoc": {"positve_min_spots": 4}}})

    def test_invalid_parameters(self):
        for override in [{"pixel_size_um": [0, 1]}, {"spots": {"cre": {"positive_min_spots": 2.5}}},
                         {"nuclei": {"exclude_border": "false"}}, {"spots": {"pnoc": {"snr_threshold": float("nan")}}}]:
            with self.subTest(override=override), self.assertRaises(ValueError):
                rn.validate_config(rn.merge_config(override))

    def test_metadata_order_takes_precedence_over_fallback(self):
        metadata = [{"name": n, "fluor": ""} for n in ["DAPI-T4", "At520-T3", "Cy3-T2", "Cy5-T1"]]
        resolved, warnings = rn.resolve_channels(metadata, 4, rn.DEFAULT_CONFIG["channels"])
        self.assertEqual(resolved, {"dapi": 0, "pnoc": 1, "cre": 3})
        self.assertFalse(warnings)

    def test_missing_metadata_fallback_and_explicit_override(self):
        resolved, warnings = rn.resolve_channels([], 4, rn.DEFAULT_CONFIG["channels"])
        self.assertEqual(resolved, {"dapi": 3, "pnoc": 2, "cre": 0})
        self.assertEqual(len(warnings), 3)
        config = copy.deepcopy(rn.DEFAULT_CONFIG["channels"])
        for role, index in [("dapi", 0), ("pnoc", 1), ("cre", 2)]:
            config[role]["index"] = index
        self.assertEqual(rn.resolve_channels([{"name": "unknown"}] * 3, 3, config)[0]["cre"], 2)

    def test_unknown_duplicate_and_ambiguous_channels_fail(self):
        with self.assertRaises(ValueError):
            rn.resolve_channels([{"name": "unexpected"}], 4, rn.DEFAULT_CONFIG["channels"])
        config = copy.deepcopy(rn.DEFAULT_CONFIG["channels"])
        config["cre"]["index"] = 3
        with self.assertRaises(ValueError):
            rn.resolve_channels([], 4, config)
        with self.assertRaises(ValueError):
            rn.resolve_channels([{"name": "DAPI"}, {"name": "DAPI"}], 4, config)

    def test_axes_do_not_silently_project_stacks(self):
        a = np.zeros((1, 8, 4, 9))
        self.assertEqual(rn.as_cyx(a, "TYCX").shape, (4, 8, 9))
        with self.assertRaisesRegex(ValueError, "Z has size 2"):
            rn.as_cyx(np.zeros((4, 2, 8, 9)), "CZYX")


class DetectionTests(unittest.TestCase):
    def test_intensity_detector_known_peaks_and_noise(self):
        config = copy.deepcopy(rn.DEFAULT_CONFIG["spots"]["pnoc"])
        config.update(method="intensity", snr_threshold=10, min_distance_um=1.0)
        positions = [(25, 25), (35, 70), (70, 40)]
        image = gaussian_spots((96, 96), positions, noise=6)
        spots, info, _ = rn.detect_spots(image, (0.5, 0.5), config)
        self.assertEqual({(s["y_px"], s["x_px"]) for s in spots}, set(positions))
        self.assertEqual(info["method"], "intensity")
        self.assertLess(len(rn.detect_spots(gaussian_spots((96, 96), [], noise=10), (0.5, 0.5), config)[0]), 3)

    def test_stardist_routes_explicitly_and_never_falls_back(self):
        config = copy.deepcopy(rn.DEFAULT_CONFIG["nuclei"])
        config["method"] = "stardist"
        image = np.zeros((40, 40), np.uint16); image[15:25, 15:25] = 5000
        expected = np.zeros_like(image, np.int32); expected[15:25, 15:25] = 1
        with patch.object(rn, "stardist_labels", return_value=(expected, {"method": "stardist"})) as model:
            labels, records, info = rn.segment_nuclei(image, (0.5, 0.5), config)
            model.assert_called_once(); np.testing.assert_array_equal(labels, expected)
            self.assertEqual(info["method"], "stardist")
            self.assertTrue(records[0]["accepted"])
        with patch.object(rn, "stardist_labels", side_effect=RuntimeError("model unavailable")):
            with self.assertRaisesRegex(RuntimeError, "model unavailable"):
                rn.segment_nuclei(image, (0.5, 0.5), config)

    def test_dapi_contrast_filter_rejects_dim_object(self):
        y, x = np.indices((100, 130))
        labels = np.zeros((100, 130), np.int32)
        labels[(y - 50) ** 2 + (x - 35) ** 2 < 10 ** 2] = 1
        labels[(y - 50) ** 2 + (x - 95) ** 2 < 10 ** 2] = 2
        image = np.full(labels.shape, 100, np.uint16); image[labels == 1] = 5000
        config = copy.deepcopy(rn.DEFAULT_CONFIG["nuclei"]); config["min_dapi_contrast"] = 1.15
        records = rn.nucleus_records(labels, (0.5, 0.5), config, image)
        self.assertTrue(records[0]["accepted"])
        self.assertIn("low_dapi_contrast", records[1]["exclusion_reason"])

    def test_dapi_contrast_smoothing_is_method_neutral(self):
        y, x = np.indices((60, 60))
        labels = np.zeros((60, 60), np.int32); labels[(y - 30) ** 2 + (x - 30) ** 2 < 8 ** 2] = 1
        image = (100 + 3000 * (labels > 0) + 50 * np.sin(x / 2.0)).astype(np.uint16)
        config = copy.deepcopy(rn.DEFAULT_CONFIG["nuclei"]); config["min_dapi_contrast"] = 1.0
        reference = rn.dapi_contrast_ratios(labels, image, (0.5, 0.5), config)[1]
        decoupled = copy.deepcopy(config); decoupled["stardist_smooth_um"] = 4.0
        self.assertEqual(rn.dapi_contrast_ratios(labels, image, (0.5, 0.5), decoupled)[1], reference)
        changed = copy.deepcopy(config); changed["contrast_smooth_um"] = 4.0
        self.assertNotEqual(rn.dapi_contrast_ratios(labels, image, (0.5, 0.5), changed)[1], reference)
        # Saved configurations predating the key still reproduce with their StarDist scale.
        legacy = {k: v for k, v in changed.items() if k != "contrast_smooth_um"}; legacy["stardist_smooth_um"] = 4.0
        self.assertEqual(rn.dapi_contrast_ratios(labels, image, (0.5, 0.5), legacy)[1],
                         rn.dapi_contrast_ratios(labels, image, (0.5, 0.5), changed)[1])
        self.assertEqual(rn.load_profile(Path(rn.__file__).with_name("settings_cellpose.json"))["nuclei"]["contrast_smooth_um"], 1.4)

    def test_log_spots_report_size_scale_and_detection_scale_separately(self):
        narrow, wide = gaussian_spots((80, 80), [(20, 20)], sigma=1.0), gaussian_spots((80, 80), [(60, 60)], sigma=2.4)
        image = (narrow.astype(np.int64) + wide - 100).clip(0, 65535).astype(np.uint16)
        config = rn.DEFAULT_CONFIG["spots"]["pnoc"]
        spots, info, _ = rn.detect_spots(image, (0.5, 0.5), config)
        by_position = {(s["y_px"], s["x_px"]): s for s in spots}
        self.assertIn((20, 20), by_position); self.assertIn((60, 60), by_position)
        # sigma_um stays the scale-selection size estimate (maximal normalized response).
        self.assertLess(by_position[(20, 20)]["sigma_um"], by_position[(60, 60)]["sigma_um"])
        for s in spots:
            self.assertGreater(s["response"], 0)
            self.assertGreaterEqual(s["detection_snr"], info["snr_threshold"])
            self.assertIn(s["detection_sigma_um"], config["sigma_um"])
        self.assertIn("detection_sigma_um", rn.SPOT_FIELDS)
        intensity = copy.deepcopy(config); intensity["method"] = "intensity"
        self.assertIsNone(rn.detect_spots(image, (0.5, 0.5), intensity)[0][0]["detection_sigma_um"])

    def test_blank_image_has_no_nuclei_or_spots(self):
        image = np.zeros((64, 64), np.uint16)
        labels, records, _ = rn.segment_nuclei(image, (0.5, 0.5), rn.DEFAULT_CONFIG["nuclei"])
        self.assertFalse(labels.any())
        self.assertFalse(records)
        spots, _, _ = rn.detect_spots(image, (0.5, 0.5), rn.DEFAULT_CONFIG["spots"]["pnoc"])
        self.assertFalse(spots)

    def test_spots_on_gradient_background_and_no_duplicate_peaks(self):
        positions = [(25, 25), (35, 70), (70, 40)]
        image = gaussian_spots((96, 96), positions, noise=6)
        image += np.arange(96, dtype=np.uint16)[None, :] * 5
        spots, _, _ = rn.detect_spots(image, (0.5, 0.5), rn.DEFAULT_CONFIG["spots"]["pnoc"])
        detected = {(s["y_px"], s["x_px"]) for s in spots}
        self.assertEqual(detected, set(positions))

    def test_noise_only_not_called_positive(self):
        image = gaussian_spots((128, 128), [], noise=10)
        spots, _, _ = rn.detect_spots(image, (0.5, 0.5), rn.DEFAULT_CONFIG["spots"]["pnoc"])
        self.assertLess(len(spots), 3)

    def test_saturated_clusters_are_flagged(self):
        image = gaussian_spots((80, 80), [(40, 40)], sigma=6, amplitude=100000)
        spots, info, mask = rn.detect_spots(image, (0.5, 0.5), rn.DEFAULT_CONFIG["spots"]["pnoc"], 65535)
        self.assertGreater(info["saturated_pixels"], 0)
        self.assertTrue(mask[40, 40])
        self.assertTrue(any("large_bright_cluster" in s["qc_flags"] for s in spots))

    def test_touching_nuclei_split_and_border_excluded(self):
        y, x = np.indices((100, 130))
        mask = (((y - 50) ** 2 + (x - 50) ** 2 < 13 ** 2) |
                ((y - 50) ** 2 + (x - 72) ** 2 < 13 ** 2) |
                ((y - 0) ** 2 + (x - 110) ** 2 < 13 ** 2))
        image = (mask * 5000 + 100).astype(np.uint16)
        labels, records, _ = rn.segment_nuclei(image, (0.5, 0.5), rn.DEFAULT_CONFIG["nuclei"])
        self.assertNotEqual(labels[50, 50], labels[50, 72])
        self.assertEqual(sum(n["accepted"] for n in records), 2)
        self.assertTrue(any("image_border" in n["exclusion_reason"] for n in records))


class CountingTests(unittest.TestCase):
    def test_cre_endpoint_tier_covers_both_cutoffs_and_respects_config(self):
        config = rn.merge_config({"spots": {"pnoc": {"positive_min_spots": 5},
                                             "cre": {"positive_min_spots": 4}}})["spots"]
        labels = np.zeros((12, 84), np.int32)
        nuclei, spots = [], {"pnoc": [], "cre": []}
        cases = [(4, 8, True), (9, 4, True), (9, 3, True), (9, 8, True),
                 (4, 8, False), (0, 8, True), (4, 0, True)]
        for i, (pnoc, cre, accepted) in enumerate(cases, 1):
            x = (i - 1) * 12 + 3
            labels[:, x:x + 5] = i
            nuclei.append({"nucleus_id": i, "accepted": accepted})
            for gene, count in [("pnoc", pnoc), ("cre", cre)]:
                spots[gene].extend({"x_px": x, "y_px": y, "qc_flags": ""} for y in range(count))
        rn.assign_and_classify(labels, nuclei, spots, config, {"pnoc_saturation": labels == 6})
        self.assertEqual([n["review_priority"] for n in nuclei],
                         ["cre_endpoint", "cre_endpoint", "general", "none", "none", "cre_endpoint", "general"])
        self.assertEqual([n["classification"] for n in nuclei],
                         ["cre_only", "double_positive", "pnoc_only", "double_positive", "excluded", "cre_only", "neither"])
        queued = sorted((n for n in nuclei if n["needs_review"]), key=rn.review_sort_key)
        self.assertEqual([n["nucleus_id"] for n in queued], [2, 1, 6, 3, 7])
        summary = rn.summarize(nuclei)
        self.assertEqual(summary["nuclei_needing_review"], 5)
        self.assertEqual(summary["nuclei_cre_endpoint_review"], 3)

    def test_priority_review_separates_information_and_cutoff_risk(self):
        labels = np.zeros((12, 60), np.int32)
        nuclei = []
        spots = {"pnoc": [], "cre": []}
        for i, pnoc_count in enumerate([8, 2, 0, 8, 2], 1):
            x = (i - 1) * 12 + 3
            labels[:, (i - 1) * 12:(i - 1) * 12 + 10] = i
            nuclei.append({"nucleus_id": i, "accepted": i != 5,
                           "qc_flags": "irregular_nuclear_shape" if i == 1 else ""})
            for gene, count in [("pnoc", pnoc_count), ("cre", 8)]:
                spots[gene].extend({"x_px": x, "y_px": y + 1, "qc_flags": ""} for y in range(count))
        masks = {"dapi_saturation": labels == 1, "pnoc_saturation": (labels == 3) | (labels == 4)}
        rn.assign_and_classify(labels, nuclei, spots, rn.DEFAULT_CONFIG["spots"], masks)
        self.assertEqual([n["needs_review"] for n in nuclei], [False, True, True, False, False])
        self.assertTrue(nuclei[0]["flag_dapi_saturation"])
        self.assertTrue(nuclei[0]["flag_irregular_nuclear_shape"])
        self.assertTrue(nuclei[1]["flag_near_pnoc_cutoff"])
        self.assertTrue(nuclei[2]["flag_saturated_pnoc_below_cutoff"])
        self.assertTrue(nuclei[3]["flag_pnoc_saturation"])
        self.assertEqual(rn.summarize(nuclei)["nuclei_needing_review"], 2)

    def test_expansion_preserves_nuclear_counts_and_blocks_rejected_objects(self):
        labels = np.zeros((15, 30), np.int32)
        labels[5:10, 5:10] = 1; labels[5:10, 15:20] = 2; labels[5:10, 23:26] = 3
        nuclei = [{"nucleus_id": i, "accepted": i != 3} for i in [1, 2, 3]]
        assignments = rn.build_assignment_labels(labels, nuclei, (1, 0.5), 2)
        self.assertEqual(assignments[7, 11], 1)
        self.assertEqual(assignments[7, 14], 2)
        self.assertFalse(assignments[labels == 3].any())
        self.assertEqual(assignments[2, 7], 0)  # 3 um away in Y, despite 0.5 um X pixels
        spots = {"pnoc": [{"x_px": x, "y_px": 7, "qc_flags": ""} for x in [6, 10, 11]], "cre": []}
        rn.assign_and_classify(labels, nuclei, spots, rn.DEFAULT_CONFIG["spots"], assignment_labels=assignments)
        self.assertEqual(nuclei[0]["pnoc_nuclear_spots"], 1)
        self.assertEqual(nuclei[0]["pnoc_spots"], 3)
        self.assertEqual(nuclei[0]["classification"], "pnoc_only")
        self.assertEqual(spots["pnoc"][1]["nucleus_id"], 0)
        self.assertEqual(spots["pnoc"][1]["assignment_id"], 1)

    def test_dapi_flags_use_the_nucleus_and_rna_flags_use_the_counting_region(self):
        labels = np.zeros((20, 40), np.int32); labels[5:15, 5:15] = 1
        nuclei = [{"nucleus_id": 1, "accepted": True}]
        assignments = rn.build_assignment_labels(labels, nuclei, (0.5, 0.5), 2.0)
        ring = (assignments == 1) & (labels == 0)
        self.assertTrue(ring.any())
        for masks, dapi_flag, pnoc_flag in [({"dapi_saturation": ring, "pnoc_saturation": ring}, False, True),
                                            ({"dapi_saturation": labels == 1}, True, False)]:
            cells = copy.deepcopy(nuclei)
            rn.assign_and_classify(labels, cells, {"pnoc": [], "cre": []}, rn.DEFAULT_CONFIG["spots"], masks, assignments)
            self.assertEqual((cells[0]["flag_dapi_saturation"], cells[0]["flag_pnoc_saturation"]), (dapi_flag, pnoc_flag))

    def test_thresholds_same_nucleus_outside_and_exclusions(self):
        labels = np.zeros((12, 32), np.int32)
        nuclei = []
        for i in range(1, 6):
            labels[2:10, (i - 1) * 6 + 1:(i - 1) * 6 + 6] = i
            nuclei.append({"nucleus_id": i, "accepted": i != 5})
        spots = {"pnoc": [], "cre": []}
        # both, Pnoc only, Cre only, neither; final nucleus excluded.
        for gene, counts in [("pnoc", [3, 3, 2, 2, 4]), ("cre", [3, 2, 3, 2, 4])]:
            for i, count in enumerate(counts):
                for j in range(count):
                    spots[gene].append({"x_px": i * 6 + 2, "y_px": j + 3, "qc_flags": ""})
            spots[gene].append({"x_px": 0, "y_px": 0, "qc_flags": ""})
        rn.assign_and_classify(labels, nuclei, spots, rn.DEFAULT_CONFIG["spots"])
        summary = rn.summarize(nuclei)
        self.assertEqual([n["classification"] for n in nuclei],
                         ["double_positive", "pnoc_only", "cre_only", "neither", "excluded"])
        self.assertEqual(summary["total_nuclei"], 4)
        self.assertEqual(summary["pnoc_positive"], 2)
        self.assertEqual(summary["cre_positive"], 2)
        self.assertEqual(summary["pct_pnoc_expressing_cre"], 50)
        self.assertEqual(sum(s["counted"] for s in spots["pnoc"]), 10)
        self.assertFalse(spots["pnoc"][-1]["counted"])

    def test_zero_denominator_and_weighted_batch_percentages(self):
        self.assertIsNone(rn.summarize([])["pct_cre_expressing_pnoc"])
        summaries = []
        for p, c, d in [(10, 5, 4), (100, 20, 10)]:
            row = {k: 0 for k in rn.COUNT_FIELDS}
            row.update(status="ok", pnoc_positive=p, cre_positive=c, double_positive=d)
            summaries.append(row)
        summaries.append({"status": "failed"})
        result = rn.aggregate(summaries)
        self.assertAlmostEqual(result["pct_pnoc_expressing_cre"], 100 * 14 / 110)
        self.assertAlmostEqual(result["pct_cre_expressing_pnoc"], 56)


class IntegrationTests(unittest.TestCase):
    def test_imagej_companions_preserve_ids_and_micrometer_calibration(self):
        with tempfile.TemporaryDirectory() as directory:
            for values in [[0, 1, 65535], [0, 65536, 2 ** 24]]:
                labels = np.array([values, values], np.uint32)
                path = Path(directory) / "labels.tif"
                status = rn.write_imagej_labels(path, labels, (0.6, 0.4))
                self.assertTrue(status["written"])
                with rn.tifffile.TiffFile(path) as tif:
                    self.assertTrue(tif.is_imagej)
                    self.assertEqual(tif.imagej_metadata["unit"], "um")
                    numerator, denominator = tif.pages[0].tags["XResolution"].value
                    self.assertAlmostEqual(denominator / numerator, 0.4)
                    np.testing.assert_array_equal(tif.asarray(), labels)
            path = Path(directory) / "too_many.tif"
            with self.assertLogs("rnascope", level="WARNING"):
                status = rn.write_imagej_labels(path, np.array([[2 ** 24 + 1]], np.uint32), (1, 1))
            self.assertFalse(status["written"])
            self.assertFalse(path.exists())

    def test_missing_package_metadata_uses_module_version(self):
        with patch.object(rn.importlib.metadata, "version", side_effect=rn.importlib.metadata.PackageNotFoundError), \
             patch.object(rn.importlib, "import_module", return_value=SimpleNamespace(__version__="2.21.0")):
            self.assertEqual(rn.package_version("tensorflow"), "2.21.0")
        with patch.object(rn.importlib.metadata, "version", side_effect=rn.importlib.metadata.PackageNotFoundError), \
             patch.object(rn.importlib, "import_module", side_effect=ImportError):
            self.assertEqual(rn.package_version("unavailable"), "unknown")

    def test_startup_environment_failure_is_reported_once(self):
        config = rn.load_profile(Path(rn.__file__).with_name("settings.json"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / "one.czi").write_bytes(b"unused")
            with patch.object(rn, "load_stardist_model", side_effect=RuntimeError("missing runtime")) as model, \
                 patch.object(rn, "sha256") as digest:
                with self.assertRaisesRegex(RuntimeError, "before analysis"):
                    rn.run_batch(root, root / "output", config)
                model.assert_called_once(); digest.assert_not_called()
                self.assertFalse((root / "output").exists())

    def test_metadata_failure_does_not_abort_batch_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / "bad.czi").write_bytes(b"invalid CZI")
            with patch.object(rn.importlib.metadata, "version", side_effect=rn.importlib.metadata.PackageNotFoundError):
                self.assertEqual(rn.run_batch(root, root / "output", rn.DEFAULT_CONFIG), 1)
            self.assertTrue((root / "output" / "run_manifest.json").exists())

    def test_shipped_production_profile_exports_with_mocked_nuclei(self):
        config = rn.load_profile(Path(rn.__file__).with_name("settings.json"))
        shape = (100, 100); y, x = np.indices(shape)
        labels = np.zeros(shape, np.int32); labels[(y - 50) ** 2 + (x - 50) ** 2 < 12 ** 2] = 1
        image = np.zeros((4, *shape), np.uint16)
        image[3] = 100 + 5000 * (labels > 0)
        positions = [(46, 46), (46, 54), (54, 50)]
        image[2] = gaussian_spots(shape, positions, noise=3)
        image[0] = gaussian_spots(shape, positions, noise=3)
        metadata = [{"name": n, "fluor": ""} for n in ["Cy5-T1", "Cy3-T2", "At520-T3", "DAPI-T4"]]
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(rn, "stardist_labels", return_value=(labels, {"method": "stardist"})):
            out = Path(directory) / "scene"
            result = rn.analyze_scene(image, (0.5, 0.5), config, metadata, out, "known.czi", 0)
            self.assertEqual(result["total_nuclei"], 1)
            self.assertEqual(result["double_positive"], 1)
            self.assertEqual(result["pct_cre_expressing_pnoc"], 100)
            self.assertIn("stardist/CLAHE", result["analysis_profile"])
            self.assertTrue((out / "assignment_labels.tif").exists())
            self.assertEqual(result["nuclei_cre_endpoint_review"], 1)
            with (out / "review_queue.csv").open() as stream:
                queue = list(csv.DictReader(stream))
            self.assertEqual(len(queue), 1)
            self.assertEqual(queue[0]["review_priority"], "cre_endpoint")
            self.assertEqual(queue[0]["nucleus_id"], "1")

    def test_known_nuclei_to_exported_results(self):
        shape = (150, 150)
        centers = [(40, 40), (40, 105), (105, 75)]
        y, x = np.indices(shape)
        mask = np.zeros(shape, bool)
        for cy, cx in centers:
            mask |= (y - cy) ** 2 + (x - cx) ** 2 <= 12 ** 2
        image = np.zeros((4, *shape), np.uint16)
        image[3] = (100 + 5000 * mask).astype(np.uint16)
        offsets = [(-4, -4), (-4, 4), (4, 0)]
        image[2] = gaussian_spots(shape, [(cy + dy, cx + dx) for cy, cx in centers[:2] for dy, dx in offsets], noise=3)
        image[0] = gaussian_spots(shape, [(centers[0][0] + dy, centers[0][1] + dx) for dy, dx in offsets], noise=3)
        metadata = [{"name": n, "fluor": ""} for n in ["Cy5-T1", "Cy3-T2", "At520-T3", "DAPI-T4"]]
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "scene"
            summary = rn.analyze_scene(image, (0.5, 0.5), rn.DEFAULT_CONFIG, metadata, out, "synthetic.czi", 0)
            self.assertEqual(summary["total_nuclei"], 3)
            self.assertEqual(summary["pnoc_positive"], 2)
            self.assertEqual(summary["cre_positive"], 1)
            self.assertEqual(summary["double_positive"], 1)
            self.assertEqual(summary["pct_pnoc_expressing_cre"], 50)
            self.assertEqual(summary["pct_cre_expressing_pnoc"], 100)
            for name in ["spots.csv", "nuclei.csv", "nuclei_labels.tif", "nuclei_labels_imagej.tif", "metadata.json",
                         "qc_overview.png", "review_crops.png", "annotated_overlay.png", "annotated_overlay.pdf"]:
                self.assertGreater((out / name).stat().st_size, 0)
            from PIL import Image
            with Image.open(out / "annotated_overlay.png") as png:
                self.assertGreaterEqual(png.width, shape[1])
                self.assertGreaterEqual(png.height, shape[0])
            with rn.tifffile.TiffFile(out / "nuclei_labels.tif") as tif:
                self.assertEqual(tif.pages[0].tags["ResolutionUnit"].value, 3)
                numerator, denominator = tif.pages[0].tags["XResolution"].value
                self.assertAlmostEqual(10000 * denominator / numerator, 0.5)

    def test_qc_rendering_helpers_preserve_coordinates_and_pixels(self):
        a = np.arange(5 * 7, dtype=np.uint8).reshape(5, 7)
        self.assertIs(rn.reduce_block(a, 1), a)
        mean = rn.reduce_block(a, 2, "mean")
        self.assertEqual(mean.shape, (2, 3)); self.assertEqual(mean.dtype, np.uint8)
        self.assertEqual(int(mean[0, 0]), round((0 + 1 + 7 + 8) / 4))
        np.testing.assert_array_equal(rn.reduce_block(a, 2, "max")[0], [8, 10, 12])
        np.testing.assert_array_equal(rn.reduce_block(a, 2, "nearest")[0], [0, 2, 4])
        rgba = np.zeros((4, 4, 4), np.uint8); rgba[1, 1] = [255, 0, 0, 255]
        self.assertEqual(tuple(rn.reduce_block(rgba, 2, "max")[0, 0]), (255, 0, 0, 255))
        self.assertEqual(rn.reduced_extent((5, 7), 2), (-0.5, 5.5, 3.5, -0.5))
        self.assertEqual(rn.overview_reduction((1024, 1024)), 1)
        self.assertEqual(rn.overview_reduction((2048, 1024)), 2)
        rgb = np.full((2, 2, 3), 100, np.uint8)
        overlay = np.zeros((2, 2, 4), np.uint8); overlay[0, 0] = [255, 255, 255, 255]; overlay[1, 1] = [0, 255, 255, 128]
        blended = rn.blend_overlay(rgb, overlay)
        self.assertEqual(tuple(blended[0, 0]), (255, 255, 255))
        self.assertEqual(tuple(blended[0, 1]), (100, 100, 100))
        self.assertEqual(tuple(blended[1, 1]), (50, 178, 178))
        self.assertEqual(int(rgb[0, 0, 0]), 100)  # input untouched
        extent, window, window_extent = rn.crop_window(10.4, 3.0, 8.0, (0.5, 0.5), (20, 30))
        np.testing.assert_allclose(extent, [2.4, 18.4, 0, 11.0])
        self.assertEqual((window[0].start, window[0].stop, window[1].start, window[1].stop), (0, 12, 2, 20))
        self.assertEqual(window_extent, (1.5, 19.5, 11.5, -0.5))

    def test_annotated_overlay_is_raster_at_render_scale_with_vector_pdf(self):
        rgb = np.zeros((30, 40, 3), np.uint8)
        overlay = np.zeros((30, 40, 4), np.uint8); overlay[10, :] = [255, 255, 255, 255]
        nuclei = [{"nucleus_id": 7, "accepted": True, "classification": "double_positive", "x_px": 20.0, "y_px": 15.0}]
        spots = {"pnoc": [{"x_px": 5, "y_px": 5, "counted": True}], "cre": [{"x_px": 30, "y_px": 25, "counted": False}]}
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            rn.save_annotated_overlay(out, rgb, overlay, nuclei, spots)
            from PIL import Image
            with Image.open(out / "annotated_overlay.png") as png:
                self.assertEqual((png.width, png.height), (80, 60))  # 2x below 2048 px
                pixels = np.asarray(png.convert("RGB"))
            self.assertEqual(tuple(pixels[20, 2]), (255, 255, 255))  # nuclear boundary row at 2x
            self.assertEqual(tuple(pixels[2, 2]), (0, 0, 0))
            self.assertTrue(np.any(np.all(pixels[25:35, 30:50] == [255, 225, 74], axis=-1)))  # ID text in class colour
            self.assertGreater((out / "annotated_overlay.pdf").stat().st_size, 0)
            legend = json.loads((out / "annotated_overlay_legend.json").read_text())
            self.assertEqual(legend["render_scale"], 2)

    def test_readable_batch_exports_nested_scenes_and_sweeps(self):
        xml = """<ImageDocument><Metadata><Information><Image><ComponentBitCount>16</ComponentBitCount>
        <Dimensions><Channels><Channel Id="0" Name="Cy5-T1"><Fluor>Cy5</Fluor></Channel>
        <Channel Id="1" Name="Cy3-T2"><Fluor>Cy3</Fluor></Channel><Channel Id="2" Name="At520-T3"><Fluor>At520</Fluor></Channel>
        <Channel Id="3" Name="DAPI-T4"><Fluor>DAPI</Fluor></Channel></Channels></Dimensions></Image></Information>
        <Scaling><Items><Distance Id="X"><Value>5e-07</Value></Distance><Distance Id="Y"><Value>5e-07</Value></Distance>
        </Items></Scaling></Metadata></ImageDocument>"""
        shape = (150, 150); centers = [(40, 40), (40, 105), (105, 75)]
        y, x = np.indices(shape); mask = np.zeros(shape, bool)
        for cy, cx in centers:
            mask |= (y - cy) ** 2 + (x - cx) ** 2 <= 12 ** 2
        image = np.zeros((4, *shape), np.uint16); image[3] = (100 + 5000 * mask).astype(np.uint16)
        offsets = [(-4, -4), (-4, 4), (4, 0)]
        image[2] = gaussian_spots(shape, [(cy + dy, cx + dx) for cy, cx in centers[:2] for dy, dx in offsets], noise=3)
        image[0] = gaussian_spots(shape, [(centers[0][0] + dy, centers[0][1] + dx) for dy, dx in offsets], noise=3)

        class FakeScene:
            sizes, dims = {"S": 1, "C": 4, "Y": 150, "X": 150, "0": 1}, "SCYX0"

            def asarray(self, maxworkers=None):
                return image[None, ..., None]

        class FakeCzi:
            def __init__(self, path):
                self.path = path

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def metadata(self):
                return xml

            scenes = {0: FakeScene(), 1: FakeScene()}

        config = rn.load_profile(Path(rn.__file__).with_name("settings_legacy.json"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / "input" / "sub").mkdir(parents=True)
            (root / "input" / "a.czi").write_bytes(b"a"); (root / "input" / "sub" / "b.CZI").write_bytes(b"b")
            with patch.object(rn.czifile, "CziFile", FakeCzi):
                self.assertEqual(rn.run_batch(root / "input", root / "out", config, recursive=True), 0)
            out = root / "out"
            with (out / "summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([(r["file"], r["scene"], r["status"]) for r in rows],
                             [("a.czi", "0", "ok"), ("a.czi", "1", "ok"), ("sub/b.CZI", "0", "ok"), ("sub/b.CZI", "1", "ok")])
            batch = json.loads((out / "batch_summary.json").read_text())
            self.assertEqual((batch["total_nuclei"], batch["cre_positive"], batch["double_positive"], batch["complete"]), (12, 4, 4, True))
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual([i["file"] for i in manifest["inputs"]], ["a.czi", "sub/b.CZI"])
            scene = out / "images" / "sub" / "b.CZI" / "scene_1"
            for name in ["nuclei.csv", "spots.csv", "nuclei_labels.tif", "metadata.json", "annotated_overlay.png", "qc_overview.png"]:
                self.assertGreater((scene / name).stat().st_size, 0)
            self.assertEqual(rn.analysis_root(scene), out.resolve())
            import sweep_parameters as sweep
            with patch.object(sweep, "plot_sweep"):
                totals = sweep.run_sweep(out, root / "sweep", [0, 1], [2, 3])
            self.assertEqual(len(totals), 8)
            baseline = next(r for r in totals if (r["expand_um"], r["pnoc_min_spots"], r["cre_min_spots"]) == (0, 3, 3))
            self.assertEqual((baseline["cre_positive"], baseline["double_positive"]), (4, 4))

    def test_bad_czi_continues_and_returns_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "bad.CZI").write_bytes(b"not a CZI")
            status = rn.run_batch(root, root / "out", rn.DEFAULT_CONFIG)
            self.assertEqual(status, 1)
            result = json.loads((root / "out" / "batch_summary.json").read_text())
            self.assertFalse(result["complete"])
            self.assertEqual(result["failed_items"], 1)
            with self.assertRaisesRegex(ValueError, "not empty"):
                rn.run_batch(root, root / "out", rn.DEFAULT_CONFIG)


if __name__ == "__main__":
    unittest.main()
