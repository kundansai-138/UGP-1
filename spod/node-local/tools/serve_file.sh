#!/bin/bash
# serve_file.sh - send FILE once to the first client that connects on PORT, then stop (raw TCP, nc).
# The receiver connects to this machine, the same direction SST uses (reader dials writer).
# Usage: nohup setsid bash serve_file.sh FILE [PORT] [SECONDS] &   -> log in ~/spod-dev/runs/serve_<PORT>.log
F=${1:?file}; P=${2:-26000}; T=${3:-900}; LOG=$HOME/spod-dev/runs/serve_$P.log
echo "serving $F ($(stat -c %s "$F") bytes, md5 $(md5sum < "$F" | cut -d' ' -f1)) on $(hostname):$P from $(date +%T) for ${T}s" > "$LOG"
timeout "$T" nc -N -lv "$P" < "$F" > /dev/null 2>> "$LOG"; rc=$?
if [ "$rc" = 0 ]; then echo "RESULT: SENT at $(date +%T)" >> "$LOG"; else echo "RESULT: nobody connected (rc $rc) at $(date +%T)" >> "$LOG"; fi
