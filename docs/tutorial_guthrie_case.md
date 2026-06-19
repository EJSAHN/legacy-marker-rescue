# Tutorial: Guthrie et al. schematic RAPD case

This short tutorial runs the reusable implementation steps for a local legacy RAPD case after the user has supplied a genome FASTA inventory, a Mash distance table, and the source PDF locally.

1. Edit `configs/guthrie_1992_cgraminicola.yaml` so that `genomes.fasta_inventory`, `genomes.mash_dist`, and `legacy_figure.pdf` point to local files.
2. Simulate permissive in-silico RAPD:
   ```bash
   python scripts/02_simulate_rapd.py --config configs/guthrie_1992_cgraminicola.yaml
   ```
3. Compare band-space distances with Mash:
   ```bash
   python scripts/03_compare_bandspace_to_mash.py --config configs/guthrie_1992_cgraminicola.yaml
   ```
4. Digitize the configured legacy figure region:
   ```bash
   python scripts/06_digitize_legacy_figure.py --config configs/guthrie_1992_cgraminicola.yaml
   ```
5. Build a simple output manifest:
   ```bash
   python scripts/09_archive_manifest.py --config configs/guthrie_1992_cgraminicola.yaml
   ```

For manuscript-level numerical validation, run:

```bash
python scripts/validate_source_data.py --workbook data/source/Supplementary_Data_S1.xlsx
```

The source PDF itself should not be committed to the repository unless its license explicitly allows redistribution.
