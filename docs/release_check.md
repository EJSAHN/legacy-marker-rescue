# Release check

Before a public release, confirm that the repository contains:

- reusable implementation code in `src/legacy_marker_rescue/`;
- source-data validation code in `scripts/validate_source_data.py`;
- derived source data in `source_data/Supplementary_Data_S1.xlsx`;
- configuration templates and documentation;
- a provenance map in `docs/provenance/`;
- no raw publisher PDFs or publisher-owned figure crops;
- no manuscript display-item drawing code;
- no local working directories or user-specific paths;
- no exploratory troubleshooting scripts, environment tests, or exploratory run wrappers.

Recommended final checks:

```bash
python -m pip install -e .
python scripts/validate_source_data.py --workbook source_data/Supplementary_Data_S1.xlsx
python -m pytest
```

The repository DOI in `CITATION.cff` can remain as `10.5281/zenodo.TBD` until the archival release is created.
