# Source records

## Genome and legacy marker inputs

Versioned accessions, BioSample/BioProject identifiers, source labels, and expected sequence counts and lengths are in `data/genome/inventory.tsv`. Download the genomic FASTA for those specific versions with NCBI Datasets. Raw FASTA and GFF files are not redistributed here. The local input manifest records file hashes.

The three historical primer sequences are in `data/genome/historical_primers.tsv`; the 25 fixed GC-matched sets are in `random_primer_sequences.tsv`. Guthrie et al. (1992), *Phytopathology* 82:832-835, is the source of the historical primers and schematic records. The legacy size and published-summary tables are derived inputs. The software does not claim to reconstruct original laboratory genotypes from a publisher PDF.

Baseline matrices and comparison tables are retained to diagnose changes in product geometry. They are not the primary inward-only analysis. The auxiliary distance column preserves the earlier comparison layer and is not substituted for Mash in the primary regenerated primer-set analysis.

Ondov, B. D., T. J. Treangen, P. Melsted, et al. 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology* 17:132. https://doi.org/10.1186/s13059-016-0997-x

## Denoyes source and reader observations

Denoyes-Rothan, B., G. Guérin, C. Délye, B. Smith, D. Minz, M. Maymon, and S. Freeman. 2003. Genetic diversity and pathogenic variability among isolates of Colletotrichum species from strawberry. *Phytopathology* 93:219-228. https://doi.org/10.1094/PHYTO.2003.93.2.219

The included image is the original embedded raster of Figure 1, 774 by 674 pixels. The supplied article identifies itself as public domain and permits reprinting with customary source credit. The image is unchanged; no sharpening, reconstruction, or synthetic detail is added. Source-image identity is recorded in the observation and geometry files.

The reader record was supplied by Ezekiel Ahn on 21 September 2026. It retains its draft flags, prior-exposure declaration, uncertain marks, and exclusions. The reader clarified that the marks were approximate rather than pixel-exact or exhaustive. A corrected isolate label is present in the recorded edit history. These observations are not a blinded, independently adjudicated ground-truth dataset.

## External gel benchmark

GelGenie dataset, version 3: https://doi.org/10.5281/zenodo.14641949

Associated article: https://doi.org/10.1038/s41467-025-59189-0

The two named source archives and published MD5 checksums are recorded in `config/benchmark.json`. Every paired image in those archives is evaluated, across the original directory partitions. The expected counts are 35 DNA-ladder images and 25 heterogeneous external gel images. The second collection includes gel types other than RAPD/DNA ladders and is reported separately. Raw external archives are retrieved from the original deposit and remain in the local cache. Refer to the deposit for its data license and attribution requirements.

Some source TIFFs contain malformed ancillary metadata while their image payloads decode consistently. Warnings remain in per-image records; the source files are not rewritten. File decoding, annotation coverage, and analytical performance are distinct checks.
