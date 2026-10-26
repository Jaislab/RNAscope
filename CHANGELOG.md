# Changelog

## Unreleased

### Audit round four: parameter, flag and tooling fixes

- Make the Cellpose worker timeout a profile setting (`nuclei.cellpose_timeout_s`),
  pass thread limits to the worker's environment before it starts, and remove
  partial weight downloads left by a killed worker on the next start.
- Evaluate DAPI review masks inside the nuclear mask rather than the expanded
  counting region; RNA masks still use the counting region.
- Add the method-neutral `nuclei.contrast_smooth_um` (1.4 µm) for the DAPI
  contrast filter, which previously reused the StarDist preprocessing scale.
  Shipped profiles set it explicitly; older saved configurations fall back to
  the scale they were run with.
- Export `detection_sigma_um`, the LoG scale whose SNR triggered each peak,
  alongside the unchanged `sigma_um`/`response` size estimate, so the two
  scales are no longer conflated; detections and existing columns are unchanged.
- Reject `cellpose_diameter_px` outside 3–500 pixels; treat CZI scaling entries
  without a usable value as missing calibration instead of failing with a
  TypeError.
- Add a GitHub Actions workflow running the unit suite, a batch test through
  the reader interface with nested inputs and a sweep of its exports, and
  refresh wording that still described StarDist as the production method.

### Audit round four: QC rendering cost

- Render QC outputs from 8-bit display copies: the annotated overlay PNG is
  drawn with Pillow at the same render scale as before, the PDF embeds the
  native-resolution raster unresampled and keeps IDs/markers as vector
  elements, review crops hand matplotlib only the crop window, overview panels
  are block-reduced to their rendered size, and closed figures are collected
  immediately. On a 2048 x 2048 synthetic scene, rendering fell from about 45 s
  and 2.2 GB peak to 23 s and 0.2 GB; the full scene from 47 s / 2.4 GB RSS to
  27 s / 0.6 GB. Analysis arrays, counts and the exported file set are unchanged.

### Audit round four: inference provenance and comparison tooling

- Define the fixed Cellpose `model.eval` keyword set once in `cellpose_backend.py`.
  The batch runner sends the complete set to the isolated worker as JSON, the
  worker applies exactly that set and echoes it, and the runner refuses masks
  whose echo differs. Scene metadata now records the echoed parameters instead
  of a second hard-coded copy; `cellpose_infer.py` uses the same constant.
- Report the digest of the weight file the worker actually loaded and verify it
  against the reference fingerprint on every worker call; record the worker
  script hash in run manifests and scene metadata.
- Reproduce the csbdeep percentile normalization in NumPy so
  `prepare_cellpose_input.py` runs in the main analysis environment, as documented.
- Locate the analysis run root by walking up from a scene folder, so the
  comparison tools accept scenes exported by `--recursive` runs with subfolders.
- State in the README that the main profile's object count is unvalidated
  against expert labels. No segmentation, detection or counting behavior changed.

### Exact-17-pixel reproducibility fix

- Replace `nuclei.cellpose_diameter_um` with `nuclei.cellpose_diameter_px`, default 17.0. Pass the pixel value directly to inference; retain physical diameter as derived metadata. Existing custom profiles using the old key must migrate explicitly.
- Remove the rounded 11.74 µm conversion that triggered a 1024→1023 pixel internal resize at reference sampling.
- Add regression coverage asserting exactly 17 pixels at the reference calibration. Exact-17-pixel control reproduces historical reference counts: 409 accepted nuclei, 38 Cre-positive, 36 double-positive (94.74%).

## 0.3.0 — Main Cellpose nuclei +2 µm workflow

- Add explicit `settings_cellpose.json` and direct batch routing to legacy Cellpose nuclei on raw DAPI, with the existing ≥3-peak cutoffs, physical 2 µm assignment, QC and exports.
- Isolate inference in a pinned subprocess because the CZI reader and Cellpose have incompatible NumPy requirements; verify weights before loading and record model provenance. No fallback to another method.
- Make exported default configurations use the Cellpose profile; historical shipped profiles remain unchanged.
- Add routing, calibration, failure and export regression tests plus an actual-inference synthetic smoke check. Real-image accuracy remains unvalidated.
- Ensure a zero-noise detector estimate remains JSON serializable.

## 0.3.0 — Optional comparison tools (added before the main-workflow change)

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
