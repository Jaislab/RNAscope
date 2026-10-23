"""Known-mask checks for the optional Cellpose second-opinion workflow."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

import cellpose_infer as inference
import compare_cellpose as comparison


class CellposeComparisonTests(unittest.TestCase):
    def test_filtered_assignment_mask_can_omit_high_rejected_ids(self):
        labels = np.array([[0, 2, 2]], np.uint32)
        records = [{"nucleus_id": 2, "accepted": True}, {"nucleus_id": 7, "accepted": False}]
        np.testing.assert_array_equal(comparison.accepted_labels(labels, records), labels)

    def test_cellpose_channels_are_pnoc_then_dapi(self):
        image = np.stack([np.full((4, 5), value, np.uint16) for value in (10, 20, 30, 40)])
        channels = {"cre": 0, "pnoc": 2, "dapi": 3}
        data, specification, axis = inference.model_input(image, channels, "cyto3_pnoc_dapi")
        np.testing.assert_array_equal(data[..., 0], image[2])
        np.testing.assert_array_equal(data[..., 1], image[3])
        self.assertEqual((specification, axis), ([1, 2], -1))
        for variant in ("nuclei_dapi", "cyto3_dapi"):
            data, specification, axis = inference.model_input(image, channels, variant)
            np.testing.assert_array_equal(data, image[3])
            self.assertEqual((specification, axis), ([0, 0], None))
        with self.assertRaises(ValueError):
            inference.model_input(image, channels, "nuclei_dapi_clahe")
        processed = np.ones((4, 5), np.float32)
        data, _, _ = inference.model_input(image, channels, "nuclei_dapi_clahe", processed)
        np.testing.assert_array_equal(data, processed)

    def test_instance_matching_is_one_to_one_for_merges_and_splits(self):
        reference = np.array([[0, 65535, 65535, 99, 99, 0]], np.uint32)
        merged = np.array([[0, 5000, 5000, 5000, 5000, 0]], np.uint32)
        pairs, stats = comparison.match_instances(reference, merged)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["iou"], .5)
        self.assertEqual(stats["reference_objects"], 2)
        self.assertEqual(stats["reference_unmatched"], 1)
        self.assertEqual(stats["query_unmatched"], 0)
        pairs, stats = comparison.match_instances(merged, reference)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(stats["query_unmatched"], 1)
        pairs, stats = comparison.match_instances(reference, reference)
        self.assertEqual({p["reference_id"] for p in pairs}, {99, 65535})
        self.assertTrue(all(p["iou"] == 1 for p in pairs))
        _, stats = comparison.match_instances(reference, np.zeros_like(reference))
        self.assertEqual(stats["matched_at_iou_threshold"], 0)
        self.assertIsNone(stats["median_matched_iou"])

    def test_cell_link_rejects_multiple_nuclei_rejected_nuclei_and_border(self):
        nuclei = np.zeros((12, 30), np.uint32)
        for i, x in enumerate((2, 9, 18, 25), 1):
            nuclei[4:6, x:x + 2] = i
        records = [{"nucleus_id": i, "accepted": i != 3} for i in range(1, 5)]
        cells = np.zeros_like(nuclei)
        cells[3:7, 1:5] = 11
        cells[3:7, 8:21] = 12
        cells[8:10, 2:5] = 13
        cells[:6, 24:28] = 14
        assignments, links = comparison.link_cell_masks(nuclei, records, cells)
        self.assertEqual({r["cell_id"]: r["status"] for r in links},
                         {11: "eligible", 12: "multiple_nuclear_occupants", 13: "no_majority_nucleus", 14: "image_border_cell"})
        np.testing.assert_array_equal(assignments[cells == 11], 1)
        self.assertFalse(assignments[cells != 11].any())
        cells[:] = 0; cells[3:7, 17:21] = 22
        assignments, links = comparison.link_cell_masks(nuclei, records, cells)
        self.assertEqual(links[0]["status"], "rejected_reference_nucleus")
        self.assertFalse(assignments.any())

    def test_ambiguous_half_nucleus_links_are_not_counted_twice(self):
        nuclei = np.zeros((10, 10), np.uint32); nuclei[3:7, 3:7] = 1
        cells = np.zeros_like(nuclei); cells[2:5, 2:8] = 1; cells[5:8, 2:8] = 2
        assignments, links = comparison.link_cell_masks(nuclei, [{"nucleus_id": 1, "accepted": True}], cells)
        self.assertTrue(all(r["status"] == "nucleus_split_across_cells" for r in links))
        self.assertFalse(assignments.any())

    def test_blinded_package_omits_method_key_and_has_empty_templates(self):
        import csv
        import tifffile
        image = np.ones((4, 40, 40), np.uint16)
        labels = np.zeros((40, 40), np.uint32); labels[10:20, 10:20] = 1
        masks = {name: labels for name in ("stardist", "nuclei_dapi", "cyto3_dapi", "stardist_plus2um", "cyto3_pnoc_cell_regions")}
        metadata = {"resolved_channels": {"cre": 0, "pnoc": 2, "dapi": 3}, "pixel_size_yx_um": [.7, .7],
                    "channels": [{"name": name} for name in ("Cy5", "Cy3", "Pnoc", "DAPI")]}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            comparison.save_blinded_review(root, image, metadata, masks)
            folder = root / "blinded_review"
            with (folder / "review_scores.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 5)
            self.assertTrue(all("method" not in row for row in rows))
            self.assertFalse(list(folder.rglob("*key*")))
            self.assertTrue((root / "blinding_key_DO_NOT_SHARE.csv").exists())
            for path in folder.glob("TEMPLATE_*.tif"):
                self.assertFalse(tifffile.imread(path).any())
            with tifffile.TiffFile(folder / "raw_images" / "full_image_raw_channels.tif") as tif:
                self.assertEqual(tif.imagej_metadata["unit"], "um")
                np.testing.assert_array_equal(tif.asarray(), image)


if __name__ == "__main__":
    unittest.main()
