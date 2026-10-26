# Repository audit responses

## Round one: execution fixes (v0.2.1)

The first audit identified a provenance crash on working TensorFlow runtimes
installed under another distribution name, and silent selection of historical
methods without `--config`. Both were fixed and are covered by tests. Analysis
requires an explicit profile; startup, logs, scene/batch summaries and manifests
identify it. Runtime/model checks occur once before image processing.

All shipped profiles are tested, and the current production profile is exercised
with mocked nucleus inference. The round-two auditor reports that both profiles
also ran unpatched on an actual `tensorflow-cpu` installation. That external
result is distinct from the local test simulating missing distribution metadata.

## Round two: evidence and review outputs (v0.2.2)

| Finding | Disposition |
|---|---|
| Localization wording exceeded the evidence | Corrected in README, methods and validation. Recomputed candidate-mask fractions are 66.7% Cre / 65.4% Pnoc. Similar detected-peak distributions do not establish a strong localization contrast. |
| Primary endpoint rationale | Retained Cre-positive → detectable Pnoc, explaining that sparse Cre detection and the ≥3 cutoff can depress the reverse percentage. Conditioning on detected Cre does not remove all selection bias. |
| Manual counting context | Two lab researchers performed approximate manual counts near 90%, relayed by the auditor. Attribution corrected in v0.2.3; no object-level labels/denominators were supplied for a formal accuracy calculation. |
| Most cells flagged for review | Added separate boolean flags and `review_reasons`. `needs_review` now prioritizes near-cutoff counts or below-cutoff genes with saturated RNA. The expanded reference queue falls from 955 to 413; DAPI saturation alone remains informational. |
| Overview too small for IDs | Added a single full-image PNG at native-or-higher resolution and a PDF with zoomable labels, matching accepted IDs in `nuclei.csv`. Existing overview/crops remain. |
| Fiji mask calibration | Added `*_labels_imagej.tif` companions in µm, with lossless uint16/float32 IDs and explicit handling of IDs too large for exact ImageJ floats. Native uint32 masks remain authoritative. |
| 12/14-bit saturation in 16-bit storage | Added CZI `ComponentBitCount` handling, retaining explicit overrides and storage-limit fallback. The 16-bit reference threshold stays 65535. |
| Intel Macs | README/methods now state that the pinned TensorFlow 2.21 environment has no Intel-macOS wheel; only Linux CPU was tested here. |

Segmentation, detection thresholds, regions and expression classes are unchanged.
The review-queue count changes by design; it is not a biological accuracy gain.
All 12 expanded-region Cre-positive/Pnoc-below-cutoff objects remain in dedicated
review pages, regardless of priority. Informational flags stay accessible;
absence from the priority queue does not establish that a mask/class is correct.

## Round three: attribution and endpoint review (v0.2.3)

The audit reports that all 29 v0.2.2 tests passed and both profiles reproduced
the reference counts on Python 3.13 with `tensorflow-cpu`. This is external
execution evidence, not a local cross-platform installation or biological
accuracy benchmark.

| Finding | Disposition |
|---|---|
| Manual counts misattributed | Credit the two lab researchers; the auditor only relayed the approximate 90% values. Cell, section and animal denominators are unavailable, so none are invented. |
| Prioritize endpoint review | Added `review_priority`, a complete sorted `review_queue.csv`, and `nuclei_cre_endpoint_review` in summaries. The expanded first tier contains 77 currently Cre-positive queued objects; the other 336 queued objects remain accessible. |
| Report calls all 77 near Pnoc cutoff | Checked the CSV: 22 are near Pnoc cutoff (17 double-positive, 5 Cre-only). The 77 include near-Cre and saturated-below-cutoff reasons as well. Near-Cre calls affect the denominator; general-tier calls may also affect it after review. |
| Sparse Cre counts | Verified median 1 Cre peak among 371 detected Pnoc-positive objects, with 72.0% ≥1, 48.2% ≥2 and 29.4% ≥3. This describes detections, not the latent expression distribution. |
| Specificity and efficiency interpretation | Documented detected RNA coexpression as the measurement. Functional recombination, reporter specificity and population labeling efficiency require distinct validation. Literature ranges and reporter experimental designs in the report are not software validation evidence. |

All 30 v0.2.3 analysis tests pass, including custom cutoffs and endpoint-tier export.
Both actual CPU profiles preserve their reference masks, peaks and expression
counts. The new tiers only change review order and add output columns/files.

## Round four: code audit of the Cellpose main workflow

| Finding | Disposition |
|---|---|
| Main profile unvalidated | README states that the Cellpose object count (409 vs 1,554 StarDist nuclei; 124 matched at IoU ≥ 0.5) is unvalidated against expert labels and that the blinded review is unscored. Profile choice remains the owner's decision. |
| Recorded Cellpose parameters were a second hard-coded copy | One shared keyword set; the worker applies and echoes it, the runner refuses masks whose echo differs, and manifests record the worker script hash and loaded-weight digest. |
| QC rendering cost scaled with pixel count | 8-bit display copies, Pillow overlay PNG, unresampled PDF raster with vector labels, sliced crops, block-reduced overview panels; 2048² scene rendering 45 s / 2.2 GB → 23 s / 0.2 GB. |
| Comparison prep required csbdeep | NumPy percentile normalization, bit-identical to csbdeep's NumPy path; runs in the main environment. |
| Fixed run-root depth in comparison tools | Walk up to the run manifest; recursive layouts work. |
| Fixed worker timeout, stale partial downloads, thread variables set too late | `nuclei.cellpose_timeout_s`; partial files removed on start; limits passed in the child's environment. |
| DAPI saturation flag used the counting region | Evaluated inside the nuclear mask; RNA masks unchanged. |
| Contrast filter reused a StarDist-named scale | `nuclei.contrast_smooth_um`, with fallback for older saved configurations. |
| LoG scale reported from the raw-response maximum | Reported from the detecting (max-SNR) scale. |
| Validation and messaging gaps | Diameter bounds; scaling metadata guard. |
| Documentation drift, no CI | Wording refreshed; GitHub Actions unit workflow added. License remains the owner's choice. |
| No successful-batch test | Batch through the reader interface with nested inputs, followed by a sweep of its exports. |

## Remaining work and project decisions

- CI across operating systems and a distributable regression dataset remain
  portability/reproducibility work. Cross-platform installs are not claimed tested.
- License and citation metadata require the owner's license/attribution choice;
  this update does not select a license or invent citation details.
- Thread counts remain bounded and fixed, with explicit configurable model tiles.
  Exposing thread settings or selecting tiles automatically is future performance work.
- Optional Z-stack projection, shifted-spot/ring-sweep chance controls, and
  saturated-area count estimates introduce additional analytical methods and
  require separate validation. A projection can combine different depths.
- Matched negative controls and expert object review remain necessary. Reported
  manual agreement and a >90% example are not sensitivity/specificity estimates.

See [VALIDATION.md](../VALIDATION.md) for numerical evidence and denominator
choices. The expanded-region result is not labeled nuclear overlap.
