# RNAscope: Pnoc association in Cre-positive nuclei

This script is primarily intended to measure **the fraction of Cre-positive
nuclei with detectable Pnoc RNA signal** in batch CZI images.

In the images for which it was developed, Cre RNA fluorescence is strongly
localized to the nucleus, whereas Pnoc mRNA appears predominantly cytoplasmic.
DAPI defines nuclei, so **nuclear spatial overlap is the directly supported
measurement**. Reliable assignment of the full cytoplasmic Pnoc population to
individual cells requires additional cell-boundary information. Nuclear
segmentation and puncta detection still require expert review and controls.

## Install and run

Tested on Linux CPU with Python 3.12; a GPU is not required.

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

Because Pnoc is mainly cytoplasmic, a nucleus with insufficient Pnoc signal
can still belong to a Pnoc-expressing cell. The script therefore **does not
estimate the total Pnoc-positive cell population**. The reverse percentage,
`pct_pnoc_expressing_cre`, describes only the detected Pnoc-positive subset;
it should not be used to estimate the fraction of all Pnoc cells expressing Cre.
These two percentages are not equally informative population measurements.

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
- `spots.csv`, `nuclei_labels.tif`, and optionally `assignment_labels.tif`:
  individual detections and their exact assignments.
- `qc_overview.png`, `review_crops.png`: nucleus/region boundaries and detections.
- `cre_only_review.csv` / `.pdf` / PNG pages: every accepted Cre-positive object
  below the Pnoc cutoff. Review these for missed signal, background and saturation.
- `config_used.json`, `run_manifest.json`, `metadata.json`, `run.log`: parameters,
  hashes, versions, channel mappings, calibration and diagnostics.

Review flags are currently broad: `needs_review` includes informational DAPI
saturation as well as possible classification problems. It is not an automatic
exclusion rule or a probability that the class is wrong. Saturated RNA clusters
remain flagged; their area is not converted to invented punctum counts.

Zero denominators are `NA` in CSV and `null` in JSON. Batch percentages use
summed counts. Images/nuclei are not automatically independent biological
replicates; group results by animal/sample using the experimental design.

## Verification and audit

```bash
python -m unittest discover -s tests -v
```

The 26 tests cover software behavior, all shipped profiles, missing version
metadata, explicit profile selection and a production-profile integration test
with mocked nucleus inference. Actual StarDist inference was also executed on
the reference CZI. Passing tests are not a biological sensitivity/specificity
benchmark. See the [audit response](docs/AUDIT_REVIEW.md) for implemented fixes
and outstanding recommendations, and the [changelog](CHANGELOG.md).
