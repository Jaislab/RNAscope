# RNAscope: Pnoc association in Cre-positive nuclei

This script is primarily intended to measure **the fraction of Cre-positive
nuclei with detectable Pnoc RNA signal** in batch CZI images.

Cre is the **sparser detected signal** in the reference image: 1,209 Cre peaks
versus 3,531 Pnoc peaks. Cells expressing Cre weakly can fall below the ≥3-peak
Cre cutoff, increasing the apparent Pnoc-only group and depressing the reverse
percentage. The Cre-conditioned endpoint focuses on the detected Cre-positive
subset; it does not recover missed Cre cells or guarantee an unbiased census.

Detected peak positions do **not** establish a strong localization contrast:
66.7% of Cre peaks and 65.4% of Pnoc peaks lie inside candidate DAPI nuclear
masks. DAPI supports a nuclear-region measurement; assigning whole cytoplasms
to individual cells requires additional boundary information. See the
[measured fractions and interpretation](VALIDATION.md#detected-peak-distribution).

## Install and run

Tested on Linux CPU with Python 3.12; a GPU is not required. The pinned
TensorFlow 2.21 environment is not an Intel-Mac installation target: that
release has no Intel-macOS wheel. Windows and Apple-silicon macOS are not
validated platforms for this repository yet.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-stardist.txt
python rnascope.py --input /path/to/czi_folder --output results_nuclear --config settings_nucleus_only.json
```

First model use downloads official StarDist weights with the registered archive
checksum and normal TLS verification. Output folders must be new or empty.
Add `--recursive` to include subfolders. Input CZIs are read-only.

**`--config` is required.** No analysis method is selected silently. The active
segmentation, detectors and expansion radius are printed at startup and saved
as `analysis_profile` in scene/batch summaries and the run manifest.

## Primary endpoint and biological limits

With the nucleus-only profile, a nucleus is Cre-positive with **≥3 detected Cy5
peaks** and Pnoc-positive with **≥3 detected ATTO 520 peaks** inside that nucleus.
The primary output is:

```text
pct_cre_expressing_pnoc = 100 × double_positive / cre_positive
```

Both genes must meet their thresholds in the same nucleus; their individual
puncta do not need to coincide. This measures nucleus-associated RNA signal,
not molecular colocalization or complete whole-cell expression.

A below-cutoff nucleus can belong to a cell expressing that gene weakly or
outside the mask. The script therefore **does not estimate the total
Pnoc-positive or Cre-positive cell population**. The reverse percentage,
`pct_pnoc_expressing_cre`, is particularly sensitive to missed low-Cre cells
and should not be used as the fraction of all Pnoc cells expressing Cre.
Neither conditional percentage is automatically an unbiased population estimate.

The Cre-conditioned fraction describes **detected RNA coexpression**. It can
inform Cre-driver characterization, but does not directly measure functional
recombination, reporter specificity, or the efficiency of labeling all Pnoc cells.

## Nuclear overlap versus a surrounding region

| Profile | Counting region | Interpretation |
|---|---|---|
| `settings_nucleus_only.json` | DAPI nuclear mask | Primary nuclear-overlap measurement |
| `settings.json` | Nucleus plus a non-overlapping 2 µm region | Additional nucleus-associated cytoplasmic signal; cell ownership is approximate |
| `settings_legacy.json` | Nuclear mask with earlier watershed/LoG methods | Historical method comparison |

The 2 µm profile remains available:

```bash
python rnascope.py --input /path/to/czi_folder --output results_expanded --config settings.json
```

Expansion applies to **both genes**. It can include cytoplasm or neighboring
processes, and it can change which objects meet the Cre cutoff. Its result
must not be labeled nuclear overlap. Nuclear-only counts remain in the expanded
run's cell table for comparison. The profiles and ≥3 criteria are unchanged
from v0.2.0; this documentation distinguishes their biological interpretations.

`--write-default-config my_settings.json` copies the shipped **StarDist + 2 µm**
profile, equivalent to `settings.json`; it does not write the legacy pipeline.
For the primary nuclear endpoint, use `settings_nucleus_only.json` explicitly.

## Example results

The reference single-section `20X_L_CZI.czi` produced:

| Measurement | Nuclear mask | Nuclear mask + 2 µm |
|---|---:|---:|
| Accepted nuclei | 1,554 | 1,554 |
| Cre-positive objects | 89 | 121 |
| Double-positive objects | 78 | 109 |
| Cre-positive with detectable Pnoc | **87.64%** | **90.08%** |

These are provisional detected-peak results. The different percentages use
different Cre-positive denominators; the expanded value is not a more accurate
nuclear-overlap result merely because it exceeds 90%. No target overlap
percentage enters the algorithm. See [methods](docs/METHODS.md) and
[validation](VALIDATION.md) for controls, assumptions and reproducibility.

## Outputs and review

- `summary.csv`, `batch_summary.csv` / `.json`: scene/batch counts, both
  percentages, active profile, assignment radius and failures.
- `images/<CZI-filename>/scene_<id>/nuclei.csv`: accepted and excluded masks,
  nuclear and assigned-region peak counts, classes and QC flags.
- `review_queue.csv`: all queued objects, with the `cre_endpoint` tier first.
  IDs and per-reason flags match `nuclei.csv`.
- `spots.csv`, `nuclei_labels.tif`, and optionally `assignment_labels.tif`:
  individual detections and their exact assignments.
- `annotated_overlay.png` / `.pdf`: a single full-image overlay at at least
  native resolution, with accepted nucleus IDs matching the CSV. The PDF labels
  remain zoomable; `annotated_overlay_legend.json` explains colors and IDs.
- `nuclei_labels_imagej.tif` / `assignment_labels_imagej.tif`: Fiji/ImageJ
  companion masks calibrated directly in µm, with unchanged IDs.
- `qc_overview.png`, `review_crops.png`: nucleus/region boundaries and detections.
- `cre_only_review.csv` / `.pdf` / PNG pages: every accepted Cre-positive object
  below the Pnoc cutoff. Review these for missed signal, background and saturation.
- `config_used.json`, `run_manifest.json`, `metadata.json`, `run.log`: parameters,
  hashes, versions, channel mappings, calibration and diagnostics.

`needs_review` marks accepted objects within ±1 peak
of either cutoff, or below a gene's cutoff with saturated RNA signal for that
gene. `review_reasons` explains the conditions; separate `flag_*` columns retain every
reason, including informational DAPI saturation and irregular mask shape.
Informational flags alone do not add an object to the queue. This is a review
heuristic, not an error probability or exclusion rule; it does not capture
every possible segmentation problem. Every Cre-only object still gets a review
page, whether or not it is queued.

`review_priority` is `cre_endpoint` for queued Cre-positive objects, `general`
for other queued objects, and `none` otherwise. The expanded example has
**77 Cre-endpoint objects out of 413 queued**; 22 of those 77 are near the Pnoc
cutoff. Near-Cre-cutoff calls can affect the denominator too. The general tier
remains relevant: below-cutoff Cre objects could become positive after review.
Summary column `nuclei_cre_endpoint_review` counts the first tier. Review crops
retain class examples and prioritize the Cre tier among the selected objects.

Saturation honors an explicit override, then CZI `ComponentBitCount`, then the
integer storage limit. Saturated clusters are not converted to invented counts.

Zero denominators are `NA` in CSV and `null` in JSON. Batch percentages use
summed counts. Images/nuclei are not automatically independent biological
replicates; group results by animal/sample using the experimental design.

## Verification and audit

```bash
python -m unittest discover -s tests -v
```

The 30 tests cover software behavior, all shipped profiles, missing version
metadata, explicit profile selection, priority/informational QC separation,
Cre-endpoint triage, bit-depth saturation, lossless ImageJ masks and production-profile exports with
mocked nucleus inference. Actual StarDist inference was also executed on
the reference CZI. Passing tests are not a biological sensitivity/specificity
benchmark. See the [audit response](docs/AUDIT_REVIEW.md) for implemented fixes
and outstanding recommendations, and the [changelog](CHANGELOG.md).
