#!/bin/bash
# Launch the REG calibration daemon detached from this terminal.
#
# RUN THIS FROM Terminal.app (not from an IDE or agent shell): macOS grants
# camera access per launching application, and Terminal is the one you authorised.
#
#   ./run_calibration.sh            # 24-hour calibration (default)
#   ./run_calibration.sh 2          # N hours
#   ./run_calibration.sh stop       # stop a running calibration
#
# Progress: tail -f reg/calibration.log        Report: ./report.sh  (any time)

cd "$(dirname "$0")" || exit 1
PY=../monitor/venv/bin/python
DB=../data/reg_calibration.db
LOG=calibration.log

if [ "$1" = "stop" ]; then
  pkill -INT -f "reg_calibrate.py" && echo "stop signal sent; the daemon closes its database cleanly" || echo "no calibration running"
  exit 0
fi

if pgrep -f "reg_calibrate.py" >/dev/null; then
  echo "a calibration is already running (pid $(pgrep -f reg_calibrate.py | head -1)); use ./run_calibration.sh stop first"
  exit 1
fi

HOURS="${1:-24}"
mkdir -p ../data
nohup "$PY" -u reg_calibrate.py --db "$DB" --hours "$HOURS" >> "$LOG" 2>&1 &
sleep 3
echo "calibration started (pid $!) for $HOURS h"
echo "database: $DB"
echo "log:      reg/$LOG"
echo "--- first log lines:"
tail -n 6 "$LOG"
