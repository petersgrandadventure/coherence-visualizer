#!/bin/bash
# Generate the calibration report from the current database (safe while the daemon runs).
cd "$(dirname "$0")" || exit 1
OUT="../data/reg_calibration_report.md"
../monitor/venv/bin/python reg_report.py --db ../data/reg_calibration.db --out "$OUT" "$@" && echo "report written to data/reg_calibration_report.md"
