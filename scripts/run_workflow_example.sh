#!/usr/bin/env bash
set -euo pipefail
CONFIG=${1:-configs/guthrie_1992_cgraminicola.yaml}
python scripts/01_prepare_genome_inventory.py --config "$CONFIG"
python scripts/02_simulate_rapd.py --config "$CONFIG"
python scripts/03_compare_bandspace_to_mash.py --config "$CONFIG"
python scripts/06_digitize_legacy_figure.py --config "$CONFIG"
python scripts/09_archive_manifest.py --config "$CONFIG"
