# Response to the v0.2.0 repository audit

The supplied independent audit was reviewed against the repository source and
reference CZI. Findings are evidence to evaluate, not proof of assay accuracy.
This update addresses execution reliability and interpretation without changing
segmentation, peak detection, thresholds, region sizes or example count results.

| Finding | Disposition in v0.2.1 |
|---|---|
| Literal `tensorflow` distribution lookup can crash an otherwise working runtime | Fixed: all package-version recording falls back to the imported module version, then `unknown`. Tests simulate missing distribution metadata. |
| Omitted `--config` silently selects historical methods | Fixed: analysis requires `--config`. Startup, scene/batch summaries and manifests identify the active method and region. |
| Default-config export writes historical parameters | Fixed: it copies the shipped StarDist + 2 µm profile; nuclear-only remains an explicit separate profile. |
| Environment failures repeat for every file | Improved: required CZI API and StarDist/model initialization are checked once before image processing. |
| Shipped profiles are absent from tests | Fixed: every shipped JSON is validated; the current production profile is exercised in an integration test with mocked nucleus inference. |
| Validation notes depend on private paths and absent result files | Rewritten as standalone documentation. Historical v0.1 notes are in `docs/archive/`. Local example outputs and screenshots are not claimed to be repository files. |
| Primary scientific interpretation is unclear | README prioritizes Cre-positive nuclei with detectable Pnoc, explains nuclear/cytoplasmic localization, and distinguishes nuclear overlap from the optional 2 µm region. |

## Remaining work

- Review-flag prioritization: `needs_review` still includes informational flags.
  Per-reason columns and a classification-impact queue need a separately tested
  rule; informational DAPI saturation should not be mistaken for a wrong class.
- Bit-depth-aware saturation, a full-resolution annotated overview and
  ImageJ-native label export remain improvements. Current saturation overrides
  and TIFF calibration are documented in the methods.
- Configurable thread counts, adaptive tiling and additional platform validation
  remain performance/portability work. Current thread limits and tiles are explicit.
- CI, distributable regression data, citation metadata and a license remain
  repository-maintenance work. License and attribution choices belong to the
  project owner; this update does not select a license or invent citation details.
- Optional Z-stack handling, shifted-spot/ring-sweep chance controls and
  area-derived saturated-cluster estimates change or add analytical methods.
  They were not introduced by this audit fix. In particular, a maximum projection
  can superimpose signals from different depths and needs separate validation.
- Matched negative controls and expert review of discrepant/near-cutoff objects
  remain experimental requirements. Neither an independent software audit nor
  a >90% example overlap establishes biological sensitivity/specificity.

The numerical reference example is documented in [VALIDATION.md](../VALIDATION.md).
The higher 2 µm result is not presented as the nuclear-overlap result.
