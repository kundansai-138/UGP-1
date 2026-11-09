#!/bin/bash
# rank_map.sh - which MPI ranks of x3d_lab / spod_* run on each host (reads PMI_RANK from /proc, read-only).
# Usage: bash rank_map.sh hosts/coloc8      (run while a test is running)
hosts=$(grep -v '^\s*#' "${1:?host file}" | cut -d: -f1)
for h in $hosts; do
  ssh -n -o BatchMode=yes -o ConnectTimeout=5 "$h" '
    for name in x3d_lab spod_nl spod_old; do
      r=$(for p in $(pgrep -u "$USER" -x $name); do tr "\0" "\n" < /proc/$p/environ 2>/dev/null | sed -n "s/^PMI_RANK=//p"; done | sort -n | paste -sd, -)
      [ -n "$r" ] && printf "%s %-8s ranks %s\n" "$(hostname)" "$name" "$r"
    done' 2>/dev/null
done
