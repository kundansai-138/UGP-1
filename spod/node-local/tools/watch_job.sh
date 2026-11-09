#!/bin/bash
# watch_job.sh - while our run_S4_nl.sh run is going, watch the S4 job's Xcompact3d step times (read-only).
# If the job's last 10 steps average more than LIMIT x its baseline, stop OUR run (never the job) and exit.
# Exits when our run ends. One status line per check in $OUT.
# Usage: bash watch_job.sh <pid of run_S4_nl.sh>   (LIMIT=1.20 HOSTS=hosts/coloc8 EVERY=60)
PID=${1:?pid of run_S4_nl.sh}
LIMIT=${LIMIT:-1.20}; EVERY=${EVERY:-60}
HOSTS=${HOSTS:-$HOME/spod-dev/hosts/coloc8}
OUT=${OUT:-$HOME/spod-dev/runs/watch_job.log}

joblog() { ls -t "$HOME"/ugp/runs/S4*/S4_*/run*/x3d_ts.log 2>/dev/null | head -1; }
steps() { awk '/Time for this time step/ {print $NF}' "$1"; }    # one line per step (s)
mean_last() { tail -n "$1" | awk '{s += $1; n++} END {if (n) printf "%.2f", s / n}'; }

L=$(joblog)
BASE=$(steps "$L" | mean_last 20)
[[ -n $BASE ]] || BASE=35.5
echo "$(date +%T) watching $L, baseline $BASE s (last 20 steps before our run), limit x$LIMIT" | tee -a "$OUT"

while kill -0 "$PID" 2>/dev/null; do
  sleep "$EVERY"
  L=$(joblog); n=$(steps "$L" | wc -l)
  age=$(( $(date +%s) - $(stat -c %Y "$L") ))
  if (( n >= 10 && age < 300 )); then
    NOW=$(steps "$L" | mean_last 10)
    echo "$(date +%T) job $(basename "$(dirname "$(dirname "$L")")")/$(basename "$(dirname "$L")") step $n: last 10 avg $NOW s (base $BASE)" >> "$OUT"
    if awk -v a="$NOW" -v b="$BASE" -v l="$LIMIT" 'BEGIN {exit !(a > b * l)}'; then
      echo "$(date +%T) JOB SLOWED: $NOW s > $LIMIT x $BASE s - stopping OUR run (pid $PID)" | tee -a "$OUT"
      kill "$PID" 2>/dev/null
      for h in $(grep -v '^\s*#' "$HOSTS" | cut -d: -f1); do
        ssh -n -o BatchMode=yes -o ConnectTimeout=10 "$h" "pkill -u \$USER -x 'x3d_lab|spod_nl|spod_old'" 2>/dev/null
      done
      exit 2
    fi
  else
    echo "$(date +%T) job not stepping (log age ${age}s, $n steps) - nothing to compare" >> "$OUT"
  fi
done
echo "$(date +%T) our run (pid $PID) has ended" | tee -a "$OUT"
