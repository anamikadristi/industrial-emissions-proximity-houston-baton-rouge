#!/usr/bin/env bash
# Run the full pipeline from the repository root.
set -euo pipefail
export RAW_DIR="${RAW_DIR:-data_raw}"
python scripts/01_download.py
python scripts/02_prepare.py
python scripts/03_analysis.py
python scripts/04_figures.py
