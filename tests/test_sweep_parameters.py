"""Sensitivity-analysis behavior with known detections and fixed cohorts."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import tifffile

import rnascope as rn
import sweep_parameters as sweep


class SweepTests(unittest.TestCase):
    def test_complete_batch_uses_summed_denominators_and_preserves_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = {"assignment": {"expand_um": 0},
                      "spots": {"pnoc": {"positive_min_spots": 3}, "cre": {"positive_min_spots": 3}}}
            rn.write_json(root / "config_used.json", config)
            rn.write_json(root / "run_manifest.json", {"config_sha256": rn.sha256(root / "config_used.json")})
            summary = []
            for number, pairs in enumerate([[(6, 4), (0, 0), (4, 3)], [(0, 3), (6, 0)]]):
                scene = root / "images" / "test.czi" / f"scene_{number}"
                scene.mkdir(parents=True)
                labels = np.zeros((10, 10), np.uint32)
                nuclei, spots = [], []
                for i, (p, c) in enumerate(pairs, 1):
                    x = 3 * i - 2
                    labels[1:9, x:x + 2] = i
                    nuclei.append({"nucleus_id": i, "accepted": True, "pnoc_spots": p, "cre_spots": c,
                                   "pnoc_nuclear_spots": p, "cre_nuclear_spots": c})
                    for gene, count in [("pnoc", p), ("cre", c)]:
                        spots.extend({"gene": gene, "y_px": y, "x_px": x} for y in range(1, count + 1))
                tifffile.imwrite(scene / "nuclei_labels.tif", labels)
                rn.write_json(scene / "metadata.json", {"pixel_size_yx_um": [1, 1]})
                rn.write_csv(scene / "nuclei.csv", nuclei, list(nuclei[0]))
                rn.write_csv(scene / "spots.csv", spots, ["gene", "y_px", "x_px"])
                summary.append({"file": "test.czi", "scene": number, "status": "ok", "total_nuclei": len(pairs),
                                "pnoc_positive": sum(p >= 3 for p, c in pairs),
                                "cre_positive": sum(c >= 3 for p, c in pairs),
                                "double_positive": sum(p >= 3 and c >= 3 for p, c in pairs)})
            rn.write_csv(root / "summary.csv", summary, list(summary[0]))
            rn.write_json(root / "batch_summary.json", {"complete": True, "successful_scenes": 2,
                          **{k: sum(s[k] for s in summary) for k in ("total_nuclei", "pnoc_positive", "cre_positive", "double_positive")}})
            before = {p: rn.sha256(p) for p in root.rglob("*") if p.is_file()}
            with patch.object(sweep, "plot_sweep"):
                rows = sweep.run_sweep(root, root / "output", [0], [3])
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["cre_positive"], 3)
            self.assertEqual(rows[0]["double_positive"], 2)
            self.assertAlmostEqual(rows[0]["pct_cre_expressing_pnoc"], 200 / 3)
            self.assertEqual(len(sweep.read_csv(root / "output" / "per_scene_sweep.csv")), 2)
            self.assertEqual(len(sweep.read_csv(root / "output" / "counts_by_nucleus.csv")), 5)
            self.assertEqual(before, {p: rn.sha256(p) for p in before})

    def test_independent_cutoffs_and_fixed_cre_cohort(self):
        counts = {"pnoc": np.array([0, 1, 2, 3, 4, 6]), "cre": np.array([6, 0, 2, 3, 4, 1])}
        baseline = np.array([True, False, False, True, True, False])
        row = sweep.add_percentages(sweep.count_setting(counts, baseline, 3, 3))
        self.assertEqual((row["pnoc_positive"], row["cre_positive"], row["double_positive"]), (3, 3, 2))
        row = sweep.add_percentages(sweep.count_setting(counts, baseline, 1, 4))
        self.assertEqual((row["cre_positive"], row["double_positive"]), (2, 1))
        self.assertEqual(row["pct_cre_expressing_pnoc"], 50)
        self.assertEqual(row["baseline_cre_positive"], 3)
        self.assertAlmostEqual(row["pct_baseline_cre_expressing_pnoc"], 200 / 3)
        self.assertEqual(row["cre_lost_vs_baseline"], 1)
        zero = sweep.add_percentages(sweep.count_setting(counts, baseline, 7, 7))
        self.assertIsNone(zero["pct_cre_expressing_pnoc"])
        self.assertIsNone(zero["pct_pnoc_expressing_cre"])

    def test_reassignment_preserves_nuclei_and_blocks_rejected_pixels(self):
        with tempfile.TemporaryDirectory() as directory:
            scene = Path(directory)
            labels = np.zeros((10, 18), np.uint32)
            labels[4:6, 3:5] = 1
            labels[4:6, 10:12] = 2
            labels[4:6, 6:8] = 3
            tifffile.imwrite(scene / "nuclei_labels.tif", labels)
            rn.write_json(scene / "metadata.json", {"pixel_size_yx_um": [.5, 1]})
            nuclei = [{"nucleus_id": i, "accepted": i != 3, "pnoc_spots": 1, "cre_spots": int(i == 1),
                       "pnoc_nuclear_spots": 1, "cre_nuclear_spots": int(i == 1)} for i in (1, 2, 3)]
            rn.write_csv(scene / "nuclei.csv", nuclei, list(nuclei[0]))
            spots = [{"gene": "pnoc", "y_px": y, "x_px": x} for y, x in [(4, 3), (6, 4), (4, 6), (4, 10)]]
            spots.append({"gene": "cre", "y_px": 4, "x_px": 3})
            rn.write_csv(scene / "spots.csv", spots, list(spots[0]))
            ids, counts, areas = sweep.read_scene(scene, [0, 1, 2], 0)
            np.testing.assert_array_equal(ids, [1, 2])
            np.testing.assert_array_equal(counts[0]["pnoc"], [1, 1])
            np.testing.assert_array_equal(counts[1]["pnoc"], [2, 1])
            np.testing.assert_array_equal(counts[2]["pnoc"], [2, 1])
            self.assertEqual(areas[0], 4)
            self.assertGreater(areas[1], areas[0])
            nuclei[0]["pnoc_spots"] = 9
            rn.write_csv(scene / "nuclei.csv", nuclei, list(nuclei[0]))
            with self.assertRaisesRegex(ValueError, "do not reproduce"):
                sweep.read_scene(scene, [0, 1], 0)

    def test_invalid_grid_and_incomplete_source_fail_before_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "output"
            for radii, cutoffs in [([-1], [3]), ([float("nan")], [3]), ([2], [0]), ([2], [1.5])]:
                with self.assertRaises(ValueError):
                    sweep.run_sweep(root, output, radii, cutoffs)
            (root / "batch_summary.json").write_text(json.dumps({"complete": False}))
            with self.assertRaisesRegex(ValueError, "complete analysis run"):
                sweep.run_sweep(root, output, [2], [3])
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
