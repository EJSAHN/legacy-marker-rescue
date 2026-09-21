# Analysis definitions

## Genome profiles

The fixed study has 22 assembly records. Exact BioSample or explicit strain/isolate identifiers define groups before resampling. Two distinct accessions labelled M1.001 form one such group; the lexicographically first accession represents that group. This results in 21 metadata representatives. This grouping is not a claim that the two assemblies have identical sequences.

A single 10-mer primer and its reverse complement are searched with at most one mismatch per binding site. Inward-facing products have the forward primer sequence upstream of the reverse-complement sequence on the stored FASTA strand. Primary products have lengths from 100 to 3,000 bp. Opposite-order sites and an 11-bp lower bound are diagnostic calculations, not primary PCR predictions. All fixed random-primer sequences are retained; resource caps stop execution rather than replace a set.

Binning sorts unique integer lengths, joins each length to the current bin if its relative distance from that bin's unrounded mean is within tolerance, and reports a rounded mean as the bin center. Primer identities remain separate. The default family has maximum lengths 2,000, 2,500, or 3,000 bp and relative tolerances 3%, 5%, or 10% in five recorded combinations.

Jaccard distance ignores shared absences. Two empty profiles have undefined between-sample distance. Pearson and Spearman coefficients summarize pairwise distances; pair entries are not treated as independent observations. Matrix-label permutations apply the same permutation to rows and columns. The maximum absolute statistic across the recorded condition family supplies a search-adjusted two-sided comparison. A within-submitter permutation is a separate conditional diagnostic. Neither accounts for unrecorded earlier exploration or establishes universal exchangeability.

Internal splitting occurs at representative level. Bins are fitted using training products only, then applied to held-out representatives; unassigned held-out products are recorded. The condition maximizing training Spearman correlation is selected without using held-out scores. A fixed log-grid representation is a separate diagnostic. Repeated split percentiles describe internal variability, not confidence intervals from independent cohorts.

## Context and size correspondence

For each primary historical product, 25 same-contig, same-length background intervals are drawn with replacement. Seed derivation uses accession and iteration; observed starts remain possible background draws. The primary summary gives each metadata representative equal weight. Gene metrics are missing when annotation coverage is unavailable or inconsistent. Gene-overlap inference is restricted to the documented annotation subset.

Legacy-size matches are compatibility comparisons only. A modern bin present in every assembly is fixed; a bin present in some but not all assemblies is variable. Prevalence-threshold summaries used for random-primer diagnostics are separately named and are not equivalent to this fixed/variable distinction. Results across size tolerances are descriptive. The provided schematic coordinates and published summaries are curated inputs, not new independent annotations of the original experiment.

## Image candidates

The lane detector uses the central fraction of each supplied lane, a vertical opening-based background, smoothing, peak prominence, separation, width, and horizontal support. Saturation and broad or weak signals are flagged for review. All detection parameters are in `config/detection.json`; reader band positions and exclusion masks are not detector inputs.

Reader matching is one-to-one within spatially defined lanes. The original primary localization allowance is 3 native pixels; 1 and 5 pixels are sensitivity analyses. Recorded exclusion masks and uncertain-band neighborhoods are accounted for separately, alongside a no-mask comparison. Because the reader observations are incomplete, unmarked candidates are not labelled false positives and the correspondence fractions are not whole-image accuracy.

For the external benchmark, image-only lane proposals precede access to the reference mask. Reference foreground components use 8-neighbor connectivity with no minimum-area deletion. Candidate centers are matched one-to-one to the component pixel region within a Euclidean native-pixel tolerance, maximizing match count before minimizing distance. All candidates and reference components remain in the denominators. Native resolution is primary. Height-normalized and reference-x-envelope analyses are separate diagnostics. This is a point-detection evaluation, not segmentation Dice or a direct comparison with the GelGenie neural model.
