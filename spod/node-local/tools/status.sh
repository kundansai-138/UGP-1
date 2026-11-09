#!/bin/bash
# status.sh - quick, read-only look at the spod-dev run and the S4 job.
# Usage (on csews180): bash ~/spod-dev/tools/status.sh
D=$HOME/spod-dev
t() { cut -c16- | sed 's/^ *//'; }      # drop the timestamp column

echo "== spod-dev: 301^3 node-local run"
P=$(cat "$D/runs/S4_coloc_nl.pid" 2>/dev/null)
if [[ -n $P ]] && kill -0 "$P" 2>/dev/null; then echo "state   : RUNNING (pid $P)"; else echo "state   : finished / not running"; fi
R=$(ls -td "$D"/runs/S4_coloc_nl/S4_* 2>/dev/null | head -1)
if [[ -n $R ]]; then
  echo "folder  : $R"
  echo "x3d     : $(grep 'Time step =' "$R/run1/x3d_ts.log" | tail -1 | t)   ($(grep -c 'Writing snapshots' "$R/run1/x3d_ts.log") snapshots written)"
  echo "x3d step: $(grep 'Time for this time step' "$R/run1/x3d_ts.log" | tail -1 | t)"
  echo "spod    : $(grep 'Processing Snapshot' "$R/run1/spod_ts.log" | tail -1 | t)"
  [[ -f $R/summary.txt ]] && echo "summary : $R/summary.txt"
fi

echo; echo "== your S4 job (~/ugp)"
L=$(ls -t "$HOME"/ugp/runs/S4*/S4_*/run*/x3d_ts.log 2>/dev/null | head -1)
echo "log     : $L"
echo "x3d     : $(grep 'Time step =' "$L" | tail -1 | t)"
echo "x3d step: $(grep 'Time for this time step' "$L" | tail -1 | t)"

echo; echo "== job step-time monitor (log only)"
tail -2 "$D/runs/watch_job.log" 2>/dev/null
