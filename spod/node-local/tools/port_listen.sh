#!/bin/bash
# port_listen.sh - wait for one TCP connection on PORT (firewall test); result in ~/spod-dev/runs/port_<PORT>.result
# Usage: nohup setsid bash port_listen.sh [PORT] [SECONDS] &     (default 26000, 900 s)
P=${1:-26000}; T=${2:-900}; OUT=$HOME/spod-dev/runs/port_$P.result
echo "listening on $(hostname):$P from $(date +%T) for ${T}s" > "$OUT"
timeout "$T" nc -lv "$P" < /dev/null > /dev/null 2>> "$OUT"; rc=$?
if [ "$rc" = 0 ]; then echo "RESULT: CONNECTED (see line above) at $(date +%T)" >> "$OUT"
else echo "RESULT: no connection (rc $rc) at $(date +%T)" >> "$OUT"; fi
