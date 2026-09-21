# Inputs and outputs

Run commands from the repository root. Installation as an importable wheel is not required. `environment.yml` describes the versions used for the local study runs; `requirements.txt` describes compatible dependency ranges. No command changes the active environment automatically.

## Full workflow inputs

Create a tab-separated genome manifest with the columns:

| Column | Meaning |
|---|---|
| `assembly` | Versioned assembly accession, matching `data/genome/inventory.tsv` |
| `fasta` | Path to that genomic FASTA, relative to the manifest or absolute |
| `fasta_sha256` | Expected SHA-256 of the FASTA file; strongly recommended |
| `gff` | Matching GFF annotation path, or an empty cell when unavailable |
| `gff_sha256` | Expected SHA-256 of a supplied GFF |

Missing annotation is not interpreted as absence of genes. Only the supplied annotations are used. Contig names, coordinates, and primer sites are checked against the FASTA.

A local JSON file identifies inputs and caches:

```json
{
  "genome_manifest": "genomes.tsv",
  "scan_cache": "cache/sites",
  "external_cache": "cache/gelgenie"
}
```

Relative paths are resolved against this JSON file. Keep this machine-specific file outside version control. A supplied `mash_table` may replace Mash execution explicitly; it must contain every unique pair as `assembly_a`, `assembly_b`, and `mash_distance`. The output records whether distances were supplied or recalculated. The study release check requires recalculation with Mash.

`genome_data` can select a different prepared input-table directory. That directory must supply the documented inventory, primer sets, comparison matrices, and curated legacy-size inputs with consistent identifiers. The default study uses the bundled tables and does not infer additional samples from a directory listing.

## Stage commands

```bash
python -B run_pipeline.py genome --inputs local_inputs.json --output-dir ../results/genome_run
python -B run_pipeline.py context --inputs context_inputs.json --output-dir ../results/context_run
python -B run_pipeline.py detection --output-dir ../results/detection_run
python -B run_pipeline.py external --inputs image_inputs.json --output-dir ../results/external_run
```

For a separate context run, add `genome_result` pointing to the `genome/` directory of a completed integrated run. `denoyes_data` can point to another compatible observation/input directory; otherwise the supplied source is used.

For candidate extraction from another image with supplied lane rectangles, use the detector directly:

```bash
python -B -m legacy_marker_rescue.gel.detect --image image.png --geometry lanes.json --config config/detection.json --output-dir ../results/candidates
```

The geometry contains the image SHA-256, its native pixel dimensions, and panel-specific lane records. Each record includes `id`, `box: [left, top, right, bottom]`, `label`, `label_source`, `role`, and `assessment`. Coordinates are pixel-edge coordinates with the origin at the top left. `role` distinguishes samples, markers, gaps, and unresolved regions. No band annotations or exclusion regions are accepted by the detector geometry. See `data/denoyes/lanes.example.json` for the supplied source geometry.

## Outputs

`genome/` contains candidate intervals, primer-specific binary matrices, pairwise Mash distances, orientation diagnostics, random-primer comparisons, and internal split results. Candidate products include both site orders so that the inward-only and diagnostic definitions can be compared without a second FASTA scan. The primary matrices are under `product_analysis/band_matrices/inward_min100/`.

`context/` contains selected inward products, interval-level sequence/annotation metrics, background draws, and explicitly scoped summaries. The representative-weighted result and the pooled-interval summaries are different estimands.

`detection/` contains candidate coordinates, lane profiles, all local maxima and rejection reasons, and reader-correspondence tables. No candidate is relabelled as absent solely because it was not marked by the reader.

`external/` contains all image selections, candidates, reference components, match assignments, failed-case records, and per-collection summaries. The two collections and the three analysis methods remain separate. A completed calculation does not imply a particular level of performance.

`STATUS.json`, software versions, code/configuration hashes, and `SHA256SUMS.txt` accompany each run. These records and local inputs can contain machine-specific information and are not added to the public source repository automatically.
