# Changelog

## 0.2.1

- Document Cre-positive → detectable Pnoc association as the primary endpoint,
  with the limitation imposed by predominantly cytoplasmic Pnoc signal.
- Distinguish nuclear overlap from the optional nucleus-plus-2 µm assignment.
- Require an explicit analysis profile, identify it in outputs and export the
  shipped current profile rather than historical defaults.
- Make package-version provenance tolerant of missing distribution metadata,
  including working TensorFlow runtimes installed under alternative names.
- Check runtime capabilities/model initialization before processing images.
- Validate every shipped profile and exercise the production profile in tests.
- Consolidate standalone methods/validation documentation and audit dispositions.
- Preserve the v0.2.0 numerical results; no detection/counting parameters changed.

## 0.2.0

- Add StarDist/CLAHE segmentation, local DAPI-contrast filtering and intensity
  peak detection; retain the earlier watershed/LoG profile separately.
- Add non-overlapping 2 µm assignment, preserve nuclear-only counts and block
  assignment into rejected DAPI objects.
- Add review outputs for every Cre-positive/Pnoc-below-cutoff object.
- Reference results differ by method and region; a higher overlap alone is not
  evidence that one method is more accurate.

## 0.1.0

- Initial CZI batch analysis with watershed nuclei, LoG peaks, independent ≥3
  expression criteria, per-object exports, QC images and provenance.
