# Changelog

## 0.2.3

- Correct attribution: two lab researchers performed the manual counts relayed
  by the auditor; cell, section and animal denominators remain unavailable.
- Prioritize queued Cre-positive objects with `review_priority=cre_endpoint`,
  export a complete sorted `review_queue.csv`, and summarize that tier's size.
- Distinguish the 77 expanded Cre-endpoint review objects from the 22 near the
  Pnoc cutoff; retain general review and every Cre-only review page.
- Document measured low Cre counts and distinguish RNA coexpression from
  functional recombination or reporter-based labeling efficiency.
- Preserve segmentation, detection, assignments and expression counts; add
  tests for configurable-cutoff triage and exported tiers.

## 0.2.2

- Replace the strong localization contrast with measured detected-peak fractions
  and the sparse-Cre / detection-cutoff rationale for the primary endpoint.
- Attribute approximate manual counts to the round-two audit without claiming
  cell-by-cell ground truth or formally validated sensitivity.
- Separate priority review reasons from informational per-reason QC flags.
- Add a native-or-higher-resolution annotated overlay with matching cell IDs.
- Add lossless, µm-calibrated ImageJ companion masks while retaining native IDs.
- Use CZI component bit depth for saturation before falling back to storage type.
- State the Intel-Mac limitation of the pinned TensorFlow 2.21 environment.
- Preserve reference expression counts and add export/QC tests.

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
