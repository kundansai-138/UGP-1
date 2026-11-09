#!/usr/bin/env python3
"""
net_rates.py - network use per machine from the S0/S4 node samplers (run*/node_<host>.csv, one row per ~1 s:
t, cpu_busy_pct, mem_used_gb, net_rx_bytes, net_tx_bytes as cumulative kernel counters of the main card).

For every machine of every run folder given: received and sent MB/s per sample (counter difference / time
difference), then median, 95th percentile, maximum, the share of samples above LIMIT MB/s, and the total GB.
1 GbE carries at most 125 MB/s of payload (about 117 MB/s of real TCP data), so samples near that are
a saturated link.

Usage: python3 net_rates.py <S0/S4 run folder> ...      (LIMIT=100 by default, in MB/s)
"""
import csv
import glob
import os
import sys

LIMIT = float(os.environ.get("LIMIT", "100"))


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    return xs[min(len(xs) - 1, int(q * (len(xs) - 1) + 0.5))]


def rates(path):
    rows = [r for r in csv.DictReader(open(path)) if r.get("net_rx_bytes")]
    rx, tx = [], []
    for a, b in zip(rows, rows[1:]):
        dt = float(b["t"]) - float(a["t"])
        if dt <= 0:
            continue
        rx.append((float(b["net_rx_bytes"]) - float(a["net_rx_bytes"])) / dt / 1e6)
        tx.append((float(b["net_tx_bytes"]) - float(a["net_tx_bytes"])) / dt / 1e6)
    tot_rx = (float(rows[-1]["net_rx_bytes"]) - float(rows[0]["net_rx_bytes"])) / 1e9 if rows else 0
    tot_tx = (float(rows[-1]["net_tx_bytes"]) - float(rows[0]["net_tx_bytes"])) / 1e9 if rows else 0
    return rx, tx, tot_rx, tot_tx


for run_dir in sys.argv[1:]:
    for rd in sorted(glob.glob(os.path.join(run_dir, "run[0-9]*"))):
        files = sorted(glob.glob(os.path.join(rd, "node_*.csv")))
        if not files:
            continue
        print(f"\n== {os.path.relpath(rd, os.path.expanduser('~'))}")
        print(f"   {'host':9s} {'rx med':>7s} {'rx p95':>7s} {'rx max':>7s} {'rx>L%':>6s} {'tx med':>7s} {'tx p95':>7s} "
              f"{'tx max':>7s} {'tx>L%':>6s} {'rx GB':>7s} {'tx GB':>7s}")
        for f in files:
            rx, tx, trx, ttx = rates(f)
            if not rx:
                continue
            host = os.path.basename(f)[5:-4]
            over_rx = 100 * sum(x > LIMIT for x in rx) / len(rx)
            over_tx = 100 * sum(x > LIMIT for x in tx) / len(tx)
            print(f"   {host:9s} {pct(rx, .5):7.1f} {pct(rx, .95):7.1f} {max(rx):7.1f} {over_rx:6.1f} "
                  f"{pct(tx, .5):7.1f} {pct(tx, .95):7.1f} {max(tx):7.1f} {over_tx:6.1f} {trx:7.1f} {ttx:7.1f}")
print(f"\n(MB/s; rx>L% / tx>L% = share of 1-s samples above {LIMIT:.0f} MB/s; 1 GbE limit ~117-125 MB/s)")
