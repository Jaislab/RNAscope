"""Behavioral tests with known spot locations and classifications."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

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
    def test_partial_settings_and_typos(self):
        config = rn.merge_config({"spots": {"pnoc": {"positive_min_spots": 4}}})
        rn.validate_config(config)
        self.assertEqual(config["spots"]["cre"]["positive_min_spots"], 3)
        self.assertEqual(rn.DEFAULT_CONFIG["spots"]["pnoc"]["positive_min_spots"], 3)
        with self.assertRaises(ValueError):
            rn.merge_config({"spots": {"pnoc": {"positve_min_spots": 4}}})

    def test_invalid_parameters(self):
        for patch in [{"pixel_size_um": [0, 1]}, {"spots": {"cre": {"positive_min_spots": 2.5}}},
                      {"nuclei": {"exclude_border": "false"}}, {"spots": {"pnoc": {"snr_threshold": float("nan")}}}]:
            with self.subTest(patch=patch), self.assertRaises(ValueError):
                rn.validate_config(rn.merge_config(patch))

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
            for name in ["spots.csv", "nuclei.csv", "nuclei_labels.tif", "metadata.json", "qc_overview.png", "review_crops.png"]:
                self.assertGreater((out / name).stat().st_size, 0)
            with rn.tifffile.TiffFile(out / "nuclei_labels.tif") as tif:
                self.assertEqual(tif.pages[0].tags["ResolutionUnit"].value, 3)
                numerator, denominator = tif.pages[0].tags["XResolution"].value
                self.assertAlmostEqual(10000 * denominator / numerator, 0.5)

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
