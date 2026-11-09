#!/bin/bash
# rank_map.sh - which MPI ranks of Xcompact3d and SPOD run on each machine (reads PMI_RANK from /proc,
# read-only). Run it while a case is running, with the host file that case uses.
# Usage: bash spod/node-local/tools/rank_map.sh scripts/hostfile
hosts=$(grep -v '^\s*#' "${1:?host file}" | cut -d: -f1)
for h in $hosts; do
  ssh -n -o BatchMode=yes -o ConnectTimeout=5 "$h" '
    for name in xcompact3d streamingSPOD spod_nl spod_old; do
      r=$(for p in $(pgrep -u "$USER" -x $name); do tr "\0" "\n" < /proc/$p/environ 2>/dev/null | sed -n "s/^PMI_RANK=//p"; done | sort -n | paste -sd, -)
      [ -n "$r" ] && printf "%s %-14s ranks %s\n" "$(hostname)" "$name" "$r"
    done' 2>/dev/null
done
