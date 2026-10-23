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
| Manual counting context | The auditor's two approximate manual counts near 90% are attributed as reported evidence. No object-level labels/denominators were supplied for a formal accuracy calculation. |
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
