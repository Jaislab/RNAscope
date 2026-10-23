# Methods and interpretation

## Intended measurement

The primary endpoint is Pnoc signal associated with Cre-positive nuclei.
The motivating images show nuclear enrichment of Cre RNA and mainly cytoplasmic
Pnoc signal. A DAPI mask directly defines a nuclear region, not an entire cell.
Nucleus-only positivity can therefore miss Pnoc expressed in the cytoplasm.
The reverse conditional percentage applies only to the detected Pnoc subset,
not the complete Pnoc-expressing population. These observations are specific to
the intended images and should be checked when applying the tool to other assays.

The two channel thresholds are applied independently (≥3 peaks in each region).
Double positivity requires both genes in the same region; it does not require
individual peaks to overlap. Counts are detections, not calibrated transcript
abundances. Different gains/laser settings prevent direct abundance comparison
across channels without additional calibration.

## Image reading and calibration

The reader reconstructs individual CZI scenes with C/Y/X axes. Non-singleton
Z, T or other extra axes are rejected; no max projection or first-plane
selection occurs. Input images remain read-only. Multiple scenes are processed
separately, and failed items do not prevent processing later files.

| Role | Metadata match | Fallback 0-based index |
|---|---|---:|
| Cre | Cy5-T1 / Cy5 | 0 |
| Pnoc | At520-T3 / ATTO 520 | 2 |
| DAPI | DAPI-T4 / DAPI | 3 |

Metadata take precedence. Fallback indices are used only when channel names are
missing, with a warning. Ambiguous/unrecognized names require an explicit
`channels.<role>.index` after confirming probe identity. Dye metadata alone do
not establish which biological probe was used. Cy3 is only a QC reference.

X/Y sampling comes from CZI scaling metadata. Missing calibration requires
`pixel_size_um: [y_um, x_um]`; no assumption of square pixels is made silently.
Coordinates are local to the scene's top-left (x=column, y=row), not stage
coordinates. Physical masks and assignment distances account for anisotropy.

## StarDist profiles

1. Smooth DAPI at 1.4 µm, normalize and apply CLAHE with a 25 µm window and
   clip 0.01. StarDist `2D_versatile_fluo` predicts nucleus outlines at target
   sampling 0.69 µm/pixel. Its default probability threshold is used unless
   explicitly overridden. RNA channels do not determine the masks.
2. Accept nuclei of 20–300 µm², excluding image-border objects and objects whose
   smoothed DAPI mean is below 1.15 times the mean in an available 1.5–3.5 µm
   surrounding ring. An unavailable ring is unknown, not grounds for rejection.
   Rejected candidates remain in the CSV and label image with their reasons.
3. Use nuclear masks directly, or expand accepted masks by the configured
   radius with nearest-label assignment and no overlap. Pixels belonging to
   rejected DAPI objects are blocked. The 2 µm expansion is an approximate
   nucleus-associated region and can include another cell's cytoplasm/processes.
4. Independently remove each RNA channel's morphological background at an 8 µm
   radius, smooth at about 0.55 µm, and detect intensity maxima above median plus
   10 × robust noise (`1.4826 × MAD`). Minimum separation uses a 1 µm footprint.
5. Count peak centers strictly inside each selected region and apply the
   per-gene ≥3 criterion. Preserve nuclear counts alongside expanded counts.

StarDist is explicitly selected and never silently replaced by thresholding.
The optional local-watershed backend uses local DAPI thresholding, watershed
and solidity filtering. The historical profile uses Li/Otsu segmentation and
multiscale, scale-normalized LoG responses with per-scale signed noise estimates.
All profile parameters are saved with each run.

## Controls and quality limits

Calibrate masks and spot thresholds with manual annotations, representative
positive fields, and negative-control slides such as dapB under matched imaging
conditions. A near-zero negative-control positive rate supports a threshold but
is not, by itself, a complete accuracy assessment. Fix parameters across
comparable samples; never tune the desired overlap percentage.

There is no automated cell-boundary segmentation, spectral unmixing, chromatic
registration or Cy3 subtraction. Subpixel/merged puncta and saturated clusters
cannot be resolved reliably. Counts include detected peaks, which can include
artifact or cluster maxima. No saturated-area-to-transcript substitution is used.

Saturation uses `qc.saturation_value` when configured, otherwise the integer
storage limit. For 12- or 14-bit data in a 16-bit container, explicitly configure
the acquisition's real ceiling (4095 or 16383); metadata-bit-depth support is an
outstanding audit recommendation. `needs_review` is currently a broad union of
flags, including informational DAPI saturation. Flagged cells are not excluded
automatically. Review all Cre-positive/Pnoc-below-cutoff objects in the dedicated
review CSV/PDF/PNG outputs, plus representative positives and negatives.

TIFF labels are unsigned 32-bit IDs, calibrated in pixels/cm with µm sizes in
metadata; they are not intensity images. Spot CSV `nucleus_id` is the raw
containing nucleus, `assignment_id` is the region owner, and `counted` indicates
contribution to an accepted object. Excluded objects never enter summary totals.

## Runtime and reproducibility

The pinned StarDist environment was tested on Linux CPU and Python 3.12.
Other operating-system/runtime combinations have not been validated here.
Compatible TensorFlow installations may use different distribution names;
version recording falls back to the imported module's version and then
`unknown`, rather than aborting analysis. Missing capabilities/model failures
are reported once during preflight, before reading image pixels.

Model downloads use the publisher's registered archive checksum and normal
TLS verification. `KERAS_HOME` selects the cache; otherwise `.cache/keras/`
beside the script is used. Models, input hashes, source hash, resolved profiles,
dependency versions and numerical detector diagnostics are recorded. Decode
and TensorFlow thread defaults are bounded; model tiling remains configurable
through `nuclei.stardist_tiles` and is not selected automatically.

Corrupt/unsupported images produce failed summary rows and a nonzero batch
exit, while remaining images continue. Verify completeness before using totals.
Batch percentages use summed counts; zero denominators are NA/null. Overlapping
fields are not deduplicated, and biological grouping must come from the study.
