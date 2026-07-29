# Reproducibility strategy

The project separates three tasks that are often mixed together.

1. **Implementation.** Library functions and example scripts show how band-space simulation, digitization, and comparison are carried out.
2. **Source-data validation.** The source-data workbook contains the derived tables that underlie the manuscript figures and tables. `scripts/validate_source_data.py` recomputes headline results from those tables.
3. **Publisher-source handling.** Raw publisher figures are not redistributed unless their source terms permit it. For copyrighted sources, derived matrices, coordinate tables, and crop metadata are provided instead.

This design keeps the public release clear and inspectable while avoiding redistribution of source figures that may be copyrighted. The validation script is the main public check for numerical agreement between the manuscript and the archived source data.
