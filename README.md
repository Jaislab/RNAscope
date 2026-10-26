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

The main analysis platform is **Cellpose nuclei +2 µm**, using the pinned
legacy Cellpose 3.1.1.3 `nuclei` model on raw DAPI. Linux CPU / Python 3.12 is
the supported setup. RNA channels do not guide segmentation. Both genes keep
the ≥3-peak threshold. The counting region is a non-overlapping 2 µm expansion
of accepted nuclei; rejected nuclear pixels are blocked. Nuclear-only counts
are also retained. This is nucleus-associated signal, not a validated whole-cell
boundary or a nuclear-overlap measurement.

Cellpose requires NumPy <2.1; the current CZI reader requires NumPy ≥2.1.
Use two environments: the batch runner exchanges numeric NPY arrays with an
isolated Cellpose subprocess. Do not install both dependency sets together.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
python3 -m venv .venv-cellpose
.venv-cellpose/bin/python -m pip install -r requirements-cellpose.txt
.venv/bin/python rnascope.py --input /path/to/czi_folder --output results_cellpose_main --config settings_cellpose.json
```

The worker defaults to `.venv-cellpose/bin/python` beside `rnascope.py`.
Set `CELLPOSE_PYTHON` to an absolute Python path for another environment location.
`CELLPOSE_LOCAL_MODELS_PATH` optionally chooses the weight cache; otherwise it
uses `.cache/cellpose` beside the scripts. The official weight download uses
verified HTTPS and must match the reference SHA-256 recorded in
[the comparison report](docs/CELLPOSE_COMPARISON.md). That fingerprint is a
reproducibility check, not an independently published publisher checksum.
Missing dependencies, weights or worker failures abort without method fallback.

Use a new/empty output folder. Add `--recursive` for subfolders. Input CZIs are
read-only. Only single optical sections and square XY sampling (within 1%) are
supported for Cellpose; anisotropic acquisitions fail explicitly. The nominal
nuclear diameter is 11.74 µm, converted using image calibration (about 17 pixels
at the reference sampling). It can be set with `nuclei.cellpose_diameter_um`;
`nuclei.cellpose_threads` controls CPU parallelism. These are starting settings,
not validated universal parameters.

**`--config` remains required.** Use `settings_cellpose.json` for the main
workflow. `--write-default-config my_settings.json` now copies that profile.
Scene metadata records the model, weight fingerprint, package versions, diameter,
normalization and inference parameters. The active method and expansion are
printed and recorded in summaries and run manifests.

The historical StarDist profiles remain available unchanged: `settings.json`
(+2 µm) and `settings_nucleus_only.json` (no expansion). To run those methods,
install `requirements-stardist.txt` in a separate StarDist environment.
The results below describe historical reference runs, not a new Cellpose
biological validation. Existing reference comparisons show substantially fewer
objects with the legacy Cellpose nuclei model; higher coexpression percentages
do not establish better accuracy. Expert review on representative images remains
necessary before interpreting scientific results.

## Endpoint and biological limits

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

## Historical profiles: nuclear overlap versus a surrounding region

| Profile | Counting region | Interpretation |
|---|---|---|
| `settings_nucleus_only.json` | DAPI nuclear mask | Historical nuclear-overlap measurement |
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


## Historical reference results

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

## Historical Cellpose comparison

The [Cellpose comparison](docs/CELLPOSE_COMPARISON.md) runs the legacy `nuclei`
and `cyto3` models on the same image, preserving ≥3 puncta and frozen detections.
It includes DAPI preprocessing controls, Pnoc/DAPI cell-outline estimates,
per-object disagreements and a blinded annotation template. Cellpose uses a
separate environment; install/run instructions are in the report.

On the reference image, Pnoc-guided cell estimates link only 20 of the original
121 Cre-positive objects. These outlines are not validated whole-cell boundaries
and are not used by the main Cellpose nuclei +2 µm profile. Expert annotations are still needed
for an accuracy comparison.

## Verification and audit

To compare RNA-positivity cutoffs and counting radii using an existing complete
run, without rerunning segmentation or changing settings:

```bash
python sweep_parameters.py --results results_expanded --output results_sweep
```

The reference [324-setting sensitivity analysis](docs/SENSITIVITY.md) shows why a
higher percentage can reflect a smaller Cre-positive denominator. Pnoc ≥2 /
Cre ≥3 at 2 µm gives 114/121 (94.21%); this is a candidate for manual validation,
not an accuracy optimum or a change to the production profile.

```bash
python -m unittest discover -s tests -v
```

The tests cover software behavior, all shipped profiles, missing version
metadata, explicit profile selection, priority/informational QC separation,
Cre-endpoint triage, bit-depth saturation, lossless ImageJ masks and production-profile exports with
mocked nucleus inference. Actual StarDist inference was also executed on
the reference CZI. Passing tests are not a biological sensitivity/specificity
benchmark. See the [audit response](docs/AUDIT_REVIEW.md) for implemented fixes
and outstanding recommendations, and the [changelog](CHANGELOG.md).
Sweep tests additionally cover independent cutoffs, fixed-group comparisons,
physical reassignment, excluded-mask blocking, baseline reproduction and
weighted multi-scene totals.
Cellpose-comparison tests cover channel ordering, unique instance matching,
ambiguous cell links and separation of blinded reviewer files from the model key.

For the main platform, run the suite in the analysis environment, then exercise
actual Cellpose inference and full scene exports on synthetic data:

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python tests/smoke_cellpose.py
```

The smoke check downloads verified weights if needed and checks four synthetic
nuclei, RNA classifications, 2 µm expansion, and CSV/TIFF/PDF exports. It does
not validate real-image segmentation accuracy or CZI end-to-end decoding.
