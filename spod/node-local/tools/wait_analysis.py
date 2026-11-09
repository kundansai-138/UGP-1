#!/usr/bin/env python3
"""
wait_analysis.py - how long streaming SPOD (reader) waits for Xcompact3d (writer), the reverse, and the
compute bursts at SPOD block ends, from the per-rank CSV files written by libadios2_wait.so
(adios2_wait_shim.c) plus SPOD's timestamped log (spod_ts.log from run_S4.sh).

Per snapshot and reader rank:
    wait    = time inside adios2_begin_step until it returns the snapshot (reader idle, writer not ready)
    read    = begin_step return -> end_step return (selection + get + data transfer in end_step)
    compute = end_step return -> next begin_step call (SPOD work on that snapshot)
Per snapshot and writer rank:
    end_step = time inside adios2_end_step (marshal + hand-off; grows when the reader is behind and the
               SST queue is full: QueueFullPolicy=Block makes end_step wait)
Derived per run (steady state = all snapshots after the first; the first wait also covers the solver start-up):
    idle %   = sum(wait) / sum(wait + read + compute)
    headroom = median of (wait + read + compute) / (read + compute): how many times more often snapshots
               could arrive before SPOD becomes the bottleneck (1.0 = SPOD is already the bottleneck)
    queue    = reader begin_step return - writer end_step return (time a snapshot waits for the reader;
               needs synced clocks between writer and reader machines, so treat it as approximate)
Bursts (snapshots on which an SPOD block completes):
    found from the 'Processing Block' lines SPOD rank 0 prints (else: compute > 2 x median and > median + 0.5 s)
    mean removal   = 'Processing Block' -> 'DONE : Processing Block'
    frequency loop = 'Processing Each Frequency' -> end of that snapshot's compute
    T              = solver time between snapshots (writer end_step return -> next begin_step call)
    absorbed if  read + burst < (QueueLimit + 1) x T;   needed QueueLimit = max(1, floor((read + burst) / T))
    writer after   = largest writer end_step on the snapshots published during the burst (shows real blocking)

Usage
    python3 scripts/wait_analysis.py runs/S4_coloc/S4_* runs/S4_coloc_nl/S4_*
    python3 wait_analysis.py <S4 folder> --engine data --csv compare.csv
Each S4 folder gets wait_steps.csv, bursts.csv, wait_summary.csv and wait_summary.txt; the comparison
table is printed.
"""
import argparse
import csv
import glob
import math
import os
import re
import statistics as st
import sys
from collections import defaultdict


def med(xs):
    xs = [x for x in xs if x is not None]
    return st.median(xs) if xs else None


def f(x, nd=3):
    return "" if x is None else f"{x:.{nd}f}"


def rank_key(k):
    return int(k) if str(k).isdigit() else 10 ** 9


def load(wait_dir):
    """All records of one run, grouped by tag (reader/writer) and rank."""
    recs = defaultdict(lambda: defaultdict(list))
    for path in sorted(glob.glob(os.path.join(wait_dir, "*_r*.csv"))):
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                try:
                    r["step"], r["status"], r["rc"] = int(r["step"]), int(r["status"]), int(r["rc"])
                    r["t"], r["dur"] = float(r["t_enter"]), float(r["dur"])
                except (KeyError, ValueError, TypeError):
                    continue
                recs[r["tag"]][r["rank"]].append(r)
    for tag in recs:
        for rank in recs[tag]:
            recs[tag][rank].sort(key=lambda r: r["t"])
    return recs


def pick_engine(ranks, want):
    """Engine name to analyse: the one matching --engine, else the one with most successful begin_steps."""
    count = defaultdict(int)
    for rows in ranks.values():
        for r in rows:
            if r["call"] == "begin_step" and r["status"] == 0:
                count[r["engine"]] += 1
    names = [e for e in count if not want or want in e]
    return max(names, key=count.get) if names else None


def reader_steps(rows):
    """One rank's records -> {step: dict(wait, read, compute, t_got, t_end, calls)}, order, eos_wait."""
    out, order = {}, []
    pend, first, ncall, eos = 0.0, None, 0, None

    def close_compute(t_next):
        if order and out[order[-1]]["t_end"] is not None and out[order[-1]]["compute"] is None:
            out[order[-1]]["compute"] = t_next - out[order[-1]]["t_end"]

    for r in rows:
        if r["call"] == "begin_step":
            if first is None:
                first = r["t"]
            pend += r["dur"]; ncall += 1
            if r["status"] == 0:
                close_compute(first)
                s = r["step"] if r["step"] >= 0 else len(order)
                out[s] = dict(wait=pend, read=None, compute=None, t_got=r["t"] + r["dur"], t_end=None, calls=ncall)
                order.append(s)
                pend, first, ncall = 0.0, None, 0
            elif r["status"] == 2:                      # end of stream
                close_compute(first)
                eos = pend
                pend, first, ncall = 0.0, None, 0
        elif r["call"] == "end_step" and order:
            cur = out[order[-1]]
            if cur["t_end"] is None:
                cur["t_end"] = r["t"] + r["dur"]
                cur["read"] = cur["t_end"] - cur["t_got"]
        elif r["call"] == "close" and order:
            close_compute(r["t"])
    return out, order, eos


def writer_steps(rows):
    """One writer rank -> {step: dict(t_begin, end, t_pub)}, close duration."""
    out, cur, close = {}, None, None
    for r in rows:
        if r["call"] == "begin_step" and r["status"] == 0:
            cur = r["step"] if r["step"] >= 0 else len(out)
            out[cur] = dict(t_begin=r["t"], end=None, t_pub=None)
        elif r["call"] == "end_step" and cur is not None and out[cur]["end"] is None:
            out[cur]["end"] = r["dur"]
            out[cur]["t_pub"] = r["t"] + r["dur"]
        elif r["call"] == "close":
            close = r["dur"]
    return out, close


def parse_blocks(log_path):
    """Block-end markers printed by SPOD rank 0, timestamped by run_S4.sh: [dict(t0, t_done, t_freq)]."""
    blocks = []
    if not log_path or not os.path.exists(log_path):
        return blocks
    for line in open(log_path, errors="replace"):
        p = line.split(" ", 1)
        if len(p) < 2:
            continue
        try:
            t = float(p[0])
        except ValueError:
            continue
        txt = p[1].strip()
        if txt == "Processing Block":
            blocks.append(dict(t0=t, t_done=None, t_freq=None))
        elif txt == "DONE : Processing Block" and blocks and blocks[-1]["t_done"] is None:
            blocks[-1]["t_done"] = t
        elif txt == "Processing Each Frequency" and blocks and blocks[-1]["t_freq"] is None:
            blocks[-1]["t_freq"] = t
    return blocks


def queue_limit_of(xml_paths):
    """QueueLimit of the solution-io SST engine (0 = unlimited, also the SST default)."""
    for p in xml_paths:
        if p and os.path.exists(p):
            txt = open(p, errors="replace").read()
            m = re.search(r'<io\s+name="solution-io".*?</io>', txt, re.S)
            q = re.search(r'key="QueueLimit"\s+value="(\d+)"', m.group(0) if m else txt)
            return int(q.group(1)) if q else 0
    return None


def analyse_run(wait_dir, engine_filter, log_path=None, qlimit=None):
    recs = load(wait_dir)
    res = {"steps": [], "summary": {}, "notes": [], "bursts": []}
    if "reader" not in recs:
        res["notes"].append(f"no reader_*.csv in {wait_dir}")
        return res
    rd = recs["reader"]
    eng = pick_engine(rd, engine_filter)
    rd = {k: [r for r in v if r["engine"] == eng] for k, v in rd.items()}
    per_rank = {k: reader_steps(v) for k, v in rd.items()}
    steps = sorted({s for out, order, _ in per_rank.values() for s in order})

    wr, weng, wclose = {}, None, []
    if "writer" in recs:
        weng = pick_engine(recs["writer"], engine_filter)
        for k, v in recs["writer"].items():
            o, c = writer_steps([r for r in v if r["engine"] == weng])
            wr[k] = o
            if c is not None:
                wclose.append(c)

    rows = []
    for i, s in enumerate(steps):
        def col(key):
            return [out[s][key] for out, _, _ in per_rank.values() if s in out and out[s][key] is not None]
        w, rdt, c = col("wait"), col("read"), col("compute")
        got = col("t_got")
        we = [o[s]["end"] for o in wr.values() if s in o and o[s]["end"] is not None]
        pub = [o[s]["t_pub"] for o in wr.values() if s in o and o[s]["t_pub"] is not None]
        rows.append(dict(
            step=s, index=i, ranks=len(w),
            wait_med=med(w), wait_min=min(w) if w else None, wait_max=max(w) if w else None,
            read_med=med(rdt), read_max=max(rdt) if rdt else None,
            compute_med=med(c), compute_max=max(c) if c else None,
            writer_end_med=med(we), writer_end_max=max(we) if we else None,
            queue=(med(got) - max(pub)) if got and pub else None, burst=""))
    res["steps"] = rows
    by_step = {r["step"]: r for r in rows}

    steady = [r for r in rows[1:] if None not in (r["wait_med"], r["read_med"], r["compute_med"])]
    sw = sum(r["wait_med"] for r in steady)
    busy = sum(r["read_med"] + r["compute_med"] for r in steady)
    head = [(r["wait_med"] + r["read_med"] + r["compute_med"]) / (r["read_med"] + r["compute_med"])
            for r in steady if r["read_med"] + r["compute_med"] > 0]

    # snapshot interval T: solver time between snapshots (writer end_step return -> next begin_step call)
    T, T_src = None, ""
    if wr:
        w0 = wr[min(wr, key=rank_key)]
        ws = sorted(s for s in w0 if w0[s]["t_pub"] is not None)
        gaps = [w0[b]["t_begin"] - w0[a]["t_pub"] for a, b in zip(ws, ws[1:]) if b in w0]
        T, T_src = med(gaps), "writer"
    if T is None and steady:
        T, T_src = med([r["wait_med"] + r["read_med"] + r["compute_med"] for r in steady]), "reader cycle"

    # bursts: block ends from SPOD's log, mapped to the snapshot whose compute window holds them (reader rank 0)
    r0_out = per_rank[min(per_rank, key=rank_key)][0]
    windows = sorted((v["t_end"], v["t_end"] + (v["compute"] or 0.0), s) for s, v in r0_out.items() if v["t_end"])
    blocks, source = parse_blocks(log_path), "log"
    found = defaultdict(list)
    for b in blocks:
        cand = [w for w in windows if w[0] - 1.0 <= b["t0"] <= w[1] + 1.0]
        if cand:
            found[cand[-1][2]].append(b)
    if not blocks:
        source = "threshold"
        cm = med([r["compute_med"] for r in steady])
        if cm is not None:
            for r in steady:
                if r["compute_max"] is not None and r["compute_max"] > max(2 * cm, cm + 0.5):
                    found[r["step"]] = []
    order0 = [w[2] for w in windows]
    for s in sorted(found):
        r = by_step.get(s)
        if not r or r["compute_max"] is None:
            continue
        r["burst"] = source
        bl = found[s]
        mean_rm = sum(b["t_done"] - b["t0"] for b in bl if b["t_done"]) if bl else None
        win_end = next((w[1] for w in windows if w[2] == s), None)
        fl = [win_end - b["t_freq"] for b in bl if b["t_freq"] and win_end]
        need = (r["read_med"] or 0.0) + r["compute_max"]
        qn = max(1, math.floor(need / T)) if T else None
        absorbed = None if not T or qlimit is None else (True if qlimit == 0 else need < (qlimit + 1) * T)
        nxt = order0[order0.index(s) + 1: order0.index(s) + 1 + max(1, (qlimit or 0) + 1)] if s in order0 else []
        wafter = [by_step[x]["writer_end_max"] for x in nxt if x in by_step and by_step[x]["writer_end_max"] is not None]
        res["bursts"].append(dict(
            step=s, source=source, compute_med=r["compute_med"], compute_max=r["compute_max"],
            ordinary_med=None, mean_removal=mean_rm, freq_loop=sum(fl) if fl else None, read=r["read_med"],
            need=need, T=T, queue_limit=qlimit, queue_needed=qn, absorbed=absorbed,
            margin=((qlimit + 1) * T - need) if T and qlimit else None,
            writer_after_max=max(wafter) if wafter else None))
    ordinary = [r for r in steady if not r["burst"]]
    om = med([r["compute_med"] for r in ordinary])
    for b in res["bursts"]:
        b["ordinary_med"] = om
    if len(blocks) > len([b for v in found.values() for b in v]):
        res["notes"].append(f"{len(blocks) - sum(len(v) for v in found.values())} block marker(s) in the log "
                            "could not be matched to a snapshot (clock offset > 1 s?)")

    eos = [e for _, _, e in per_rank.values() if e is not None]
    wend = [r["writer_end_med"] for r in rows if r["writer_end_med"] is not None]
    bs = res["bursts"]
    res["summary"] = dict(
        engine=eng, writer_engine=weng or "", reader_ranks=len(per_rank), writer_ranks=len(wr),
        snapshots=len(rows),
        first_wait_s=rows[0]["wait_med"] if rows else None,
        wait_med_s=med([r["wait_med"] for r in steady]),
        wait_max_s=max([r["wait_max"] for r in steady], default=None),
        wait_total_s=sw if steady else None,
        wait_spread_med_s=med([r["wait_max"] - r["wait_min"] for r in steady]),
        read_med_s=med([r["read_med"] for r in steady]),
        read_max_s=max([r["read_max"] for r in steady], default=None),
        compute_med_s=med([r["compute_med"] for r in steady]),
        compute_ordinary_s=om,
        compute_max_s=max([r["compute_max"] for r in steady], default=None),
        last_compute_s=rows[-1]["compute_med"] if rows else None,
        idle_pct=100 * sw / (sw + busy) if steady and sw + busy > 0 else None,
        headroom_x=med(head),
        eos_wait_s=med(eos),
        writer_end_med_s=med(wend),
        writer_end_max_s=max([r["writer_end_max"] for r in rows if r["writer_end_max"] is not None], default=None),
        writer_end_total_s=sum(wend) if wend else None,
        writer_close_med_s=med(wclose),
        queue_med_s=med([r["queue"] for r in rows[1:]]),
        snapshot_interval_s=T,
        bursts=len(bs),
        burst_max_s=max([b["compute_max"] for b in bs], default=None),
        burst_mean_removal_max_s=max([b["mean_removal"] for b in bs if b["mean_removal"] is not None], default=None),
        burst_freq_loop_max_s=max([b["freq_loop"] for b in bs if b["freq_loop"] is not None], default=None),
        queue_limit=qlimit,
        queue_needed=max([b["queue_needed"] for b in bs if b["queue_needed"] is not None], default=None),
        burst_margin_min_s=min([b["margin"] for b in bs if b["margin"] is not None], default=None),
        writer_after_burst_max_s=max([b["writer_after_max"] for b in bs if b["writer_after_max"] is not None],
                                     default=None),
    )
    res["T_src"] = T_src
    calls = [out[s]["calls"] for out, _, _ in per_rank.values() for s in out]
    if calls and max(calls) > 1:
        res["notes"].append(f"begin_step was retried (timeout / not ready) up to {max(calls)} times per snapshot; "
                            "wait includes all tries")
    if not wr:
        res["notes"].append("no writer_*.csv: writer side not measured (T taken from the reader cycle)")
    if log_path and os.path.exists(log_path) and not blocks:
        res["notes"].append("no 'Processing Block' lines in spod_ts.log: bursts found by compute threshold, "
                            "no mean-removal / frequency-loop split")
    return res


def layout_of(folder):
    m = re.search(r"S4_r(\d+x\d+)", os.path.basename(folder.rstrip("/")))
    if m:
        return m.group(1)
    sc = os.path.join(folder, "summary.csv")
    if os.path.exists(sc):
        with open(sc, newline="") as fh:
            for r in csv.DictReader(fh):
                if r.get("readers"):
                    return r["readers"]
    return os.path.basename(folder.rstrip("/"))


def runs_in(path):
    """(label, wait_dir, run_dir) for an S4 folder (run*/wait), a run folder (wait/) or a wait folder."""
    if glob.glob(os.path.join(path, "run*", "wait")):
        return [(os.path.basename(os.path.dirname(w)), w, os.path.dirname(w))
                for w in sorted(glob.glob(os.path.join(path, "run*", "wait")))]
    if os.path.isdir(os.path.join(path, "wait")):
        return [(os.path.basename(path.rstrip("/")), os.path.join(path, "wait"), path)]
    if glob.glob(os.path.join(path, "*_r*.csv")):
        return [(os.path.basename(path.rstrip("/")), path, os.path.dirname(path.rstrip("/")))]
    return []


STEP_KEYS = ["step", "index", "ranks", "wait_med", "wait_min", "wait_max", "read_med", "read_max",
             "compute_med", "compute_max", "writer_end_med", "writer_end_max", "queue", "burst"]
BURST_KEYS = ["step", "source", "compute_med", "compute_max", "ordinary_med", "mean_removal", "freq_loop", "read",
              "need", "T", "queue_limit", "queue_needed", "absorbed", "margin", "writer_after_max"]
SUM_KEYS = ["snapshots", "reader_ranks", "writer_ranks", "first_wait_s", "wait_med_s", "wait_max_s", "wait_total_s",
            "wait_spread_med_s", "read_med_s", "read_max_s", "compute_med_s", "compute_ordinary_s", "compute_max_s",
            "last_compute_s", "idle_pct", "headroom_x", "eos_wait_s", "writer_end_med_s", "writer_end_max_s",
            "writer_end_total_s", "writer_close_med_s", "queue_med_s", "snapshot_interval_s", "bursts",
            "burst_max_s", "burst_mean_removal_max_s", "burst_freq_loop_max_s", "queue_limit", "queue_needed",
            "burst_margin_min_s", "writer_after_burst_max_s", "engine", "writer_engine"]


def cell(v, nd=6):
    return f(v, nd) if isinstance(v, float) else ("" if v is None else v)


def burst_lines(res):
    s, bs = res["summary"], res["bursts"]
    if not bs:
        return ["  bursts       none (no SPOD block completed, or none stood out)"]
    q = s["queue_limit"]
    qtxt = "unknown" if q is None else ("0 (unlimited)" if q == 0 else str(q))
    T = s["snapshot_interval_s"]
    cap = "" if not T or not q else f" -> absorbs read + burst up to {f((q + 1) * T, 1)} s"
    out = [f"  bursts       {len(bs)} block end(s) (found from {bs[0]['source']}): max {f(s['burst_max_s'], 2)} s"
           f" vs ordinary {f(s['compute_ordinary_s'])} s;"
           f" mean removal max {f(s['burst_mean_removal_max_s'], 2)} s, frequency loop max {f(s['burst_freq_loop_max_s'], 2)} s",
           f"               snapshot interval T {f(T, 1)} s ({res.get('T_src', '')}), QueueLimit {qtxt}{cap};"
           f" QueueLimit needed >= {s['queue_needed']}; smallest margin {f(s['burst_margin_min_s'], 1)} s;"
           f" writer end_step after a burst max {f(s['writer_after_burst_max_s'])} s"]
    for b in bs:
        out.append(f"               step {b['step']:>4}: compute {f(b['compute_max'], 2)} s"
                   f" (mean removal {f(b['mean_removal'], 2)}, frequency loop {f(b['freq_loop'], 2)}),"
                   f" read {f(b['read'], 2)} s, needs QueueLimit >= {b['queue_needed']},"
                   f" {'absorbed' if b['absorbed'] else ('NOT absorbed' if b['absorbed'] is False else '?')},"
                   f" writer end_step after {f(b['writer_after_max'])} s")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", help="S4 folders (with run*/wait), run folders or wait folders")
    ap.add_argument("--engine", default="", help="only engines whose name contains this (default: busiest engine)")
    ap.add_argument("--csv", help="write the comparison table (one row per layout) to this file")
    a = ap.parse_args()

    table = []
    for folder in a.folders:
        runs = runs_in(folder)
        if not runs:
            print(f"WARNING: no wait CSV files under {folder}", file=sys.stderr)
            continue
        results = []
        for name, w, rdir in runs:
            ql = queue_limit_of([os.path.join(folder, "adios2_config.xml"), os.path.join(rdir, "adios2_config.xml"),
                                 os.path.join(os.path.dirname(rdir.rstrip("/")), "adios2_config.xml")])
            results.append((name, analyse_run(w, a.engine, os.path.join(rdir, "spod_ts.log"), ql)))
        lines = [f"Reader wait analysis: {folder}  (layout {layout_of(folder)})", ""]
        with open(os.path.join(folder, "wait_steps.csv"), "w", newline="") as fh:
            wcsv = csv.writer(fh)
            wcsv.writerow(["run"] + STEP_KEYS)
            for name, res in results:
                for r in res["steps"]:
                    wcsv.writerow([name] + [cell(r[k]) for k in STEP_KEYS])
        with open(os.path.join(folder, "bursts.csv"), "w", newline="") as fh:
            wcsv = csv.writer(fh)
            wcsv.writerow(["run"] + BURST_KEYS)
            for name, res in results:
                for b in res["bursts"]:
                    wcsv.writerow([name] + [cell(b[k]) for k in BURST_KEYS])
        with open(os.path.join(folder, "wait_summary.csv"), "w", newline="") as fh:
            wcsv = csv.writer(fh)
            wcsv.writerow(["run"] + SUM_KEYS)
            for name, res in results:
                wcsv.writerow([name] + [cell(res["summary"].get(k)) for k in SUM_KEYS])
        for name, res in results:
            s = res["summary"]
            if not s:
                lines += [f"{name}: " + "; ".join(res["notes"]), ""]
                continue
            lines += [
                f"{name}: {s['snapshots']} snapshots, {s['reader_ranks']} reader ranks, {s['writer_ranks']} writer ranks"
                f" (engine {s['engine']})",
                f"  reader wait  first {f(s['first_wait_s'], 1)} s (includes solver start-up),"
                f" then per snapshot median {f(s['wait_med_s'], 2)} s, max {f(s['wait_max_s'], 2)} s,"
                f" total {f(s['wait_total_s'], 1)} s",
                f"  reader busy  read {f(s['read_med_s'])} s (max {f(s['read_max_s'])}),"
                f" compute {f(s['compute_med_s'])} s (max {f(s['compute_max_s'])}), last snapshot {f(s['last_compute_s'], 2)} s",
                f"  reader idle  {f(s['idle_pct'], 1)} %   headroom {f(s['headroom_x'], 1)}x"
                f"   rank spread of wait {f(s['wait_spread_med_s'])} s",
                f"  writer       end_step median {f(s['writer_end_med_s'])} s, max {f(s['writer_end_max_s'])} s,"
                f" total {f(s['writer_end_total_s'], 2)} s; close {f(s['writer_close_med_s'], 2)} s;"
                f" queue time {f(s['queue_med_s'])} s",
            ] + burst_lines(res) + [f"  note: {n}" for n in res["notes"]] + [""]
        ok = [res["summary"] for _, res in results if res["summary"].get("wait_med_s") is not None]

        def agg(k):
            return med([s.get(k) for s in ok])
        row = dict(name=os.path.basename(folder.rstrip("/")), layout=layout_of(folder), folder=folder, runs=len(ok))
        for k in SUM_KEYS[:-2]:
            row[k] = agg(k)
        table.append(row)
        txt = "\n".join(lines)
        print(txt)
        with open(os.path.join(folder, "wait_summary.txt"), "w") as fh:
            fh.write(txt + "\n")

    if not table:
        sys.exit("no data")
    cols = [("name", "S4 folder", None), ("layout", "layout", None), ("runs", "runs", None),
            ("reader_ranks", "ranks", 0), ("snapshots", "snaps", 0), ("wait_med_s", "wait/snap s", 2),
            ("read_med_s", "read s", 3), ("compute_ordinary_s", "compute s", 3), ("idle_pct", "idle %", 1),
            ("headroom_x", "headroom x", 1), ("bursts", "bursts", 0), ("burst_max_s", "burst max s", 2),
            ("queue_needed", "Q needed", 0), ("writer_end_med_s", "W end_step s", 3),
            ("writer_end_max_s", "W end_step max", 3)]
    print("COMPARISON (median over runs; wait/read = median per snapshot after the first; "
          "compute = ordinary snapshots; burst = block-end snapshots)")
    widths = [max([len(h), 8] + [len(str(r.get(k, ''))) for r in table if nd is None]) for k, h, nd in cols]
    print("  ".join(h.rjust(w) for (_, h, _), w in zip(cols, widths)))
    for row in table:
        cells = []
        for (k, _, nd), w in zip(cols, widths):
            v = row.get(k)
            cells.append((f(v, nd) if isinstance(v, (float, int)) and not isinstance(v, bool) and nd is not None
                          else str(v if v is not None else "")).rjust(w))
        print("  ".join(cells))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            wcsv = csv.DictWriter(fh, list(table[0].keys()))
            wcsv.writeheader()
            wcsv.writerows(table)
        print(f"\nwritten {a.csv}")


if __name__ == "__main__":
    main()
