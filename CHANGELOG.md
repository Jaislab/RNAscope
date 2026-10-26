## Exact-17-pixel reproducibility fix

- Replace `nuclei.cellpose_diameter_um` with `nuclei.cellpose_diameter_px`, default 17.0. Pass the pixel value directly to inference; retain physical diameter as derived metadata. Existing custom profiles using the old key must migrate explicitly.
- Remove the rounded 11.74 µm conversion that triggered a 1024→1023 pixel internal resize at reference sampling.
- Add regression coverage asserting exactly 17 pixels at the reference calibration. Exact-17-pixel control reproduces historical reference counts: 409 accepted nuclei, 38 Cre-positive, 36 double-positive (94.74%).

## 0.3.0 — Main Cellpose nuclei +2 µm workflow

- Add explicit `settings_cellpose.json` and direct batch routing to legacy Cellpose nuclei on raw DAPI, with the existing ≥3-peak cutoffs, physical 2 µm assignment, QC and exports.
- Isolate inference in a pinned subprocess because the CZI reader and Cellpose have incompatible NumPy requirements; verify weights before loading and record model provenance. No fallback to another method.
- Make exported default configurations use the Cellpose profile; historical shipped profiles remain unchanged.
- Add routing, calibration, failure and export regression tests plus an actual-inference synthetic smoke check. Real-image accuracy remains unvalidated.
- Ensure a zero-noise detector estimate remains JSON serializable.

# Changelog

## Unreleased: optional comparison tools

- Add an isolated, pinned Cellpose 3 workflow for the requested `nuclei` and
  `cyto3` models, with verified source export, model/mask fingerprints and
  DAPI preprocessing controls; preserve the production StarDist environment.
- Compare frozen ≥3-peak calls under nuclear masks/2 µm regions and linked
  Pnoc/DAPI-guided cell estimates. Report missing coverage and ambiguous links.
- Prepare blinded raw-image/outline review files with a separate model key and
  empty annotation templates; no expert-label accuracy claims are made.
- Add six Cellpose-comparison tests; 40 tests pass in total.

- Add `sweep_parameters.py` to reuse complete analysis exports across independent
  gene cutoffs and assignment radii, with baseline reproduction checks and hashes.
- Report full numerators/denominators, fixed-baseline Cre-group comparisons,
  per-object counts, weighted batch totals and PNG/PDF sensitivity plots.
- Document 324 settings on the reference section and a five-object Pnoc ≥2
  validation candidate; retain production profiles and all core analysis behavior.
- Add four sweep tests; the sensitivity update brought the test total to 34.

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
