#!/bin/bash
# probe_nodes.sh - read-only look at lab machines: CPU, threads, RAM, load, logged-in users, busiest process.
# Usage: bash probe_nodes.sh [first] [last]   (default 100 200; skips csews160-171 of the S4 job and csews180)
first=${1:-100}; last=${2:-200}
probe() {
  timeout 8 ssh -n -o BatchMode=yes -o ConnectTimeout=4 "$1" '
    cpu=$(lscpu | sed -n "s/^Model name:\s*//p" | sed "s/13th Gen Intel(R) Core(TM) //; s/Intel(R) Core(TM) //")
    mem=$(awk "/MemTotal/{printf \"%.0f\", \$2/1048576}" /proc/meminfo)
    users=$(who | awk "{print \$1}" | sort -u | paste -sd, -)
    top=$(ps -eo pcpu,user,comm --sort=-pcpu --no-headers | head -1 | awk "{print \$1\"%\", \$2, \$3}")
    echo "$(hostname) | $cpu | $(nproc) thr | ${mem}G | load $(cut -d" " -f1-3 /proc/loadavg) | users ${users:--} | top $top"
  ' 2>/dev/null || echo "$1 | unreachable"
}
export -f probe
for i in $(seq "$first" "$last"); do
  case $i in 16[0-9]|17[01]|180) continue ;; esac
  echo "csews$i"
done | xargs -P 25 -I{} bash -c 'probe {}' | sort -V
