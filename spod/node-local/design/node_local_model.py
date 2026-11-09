"""
Reference model for the node-local SPOD input layout.

It answers three questions without needing MPI, PETSc or ADIOS2:

  1. Mapping: does every SPOD rank read only boxes that lie inside a single
     Xcompact3d (X3D) writer's pencil, with that writer on the same node,
     and do all ranks together read every grid point exactly once?
  2. Maths: does streaming SPOD give the same energies and modes when the
     rows are stored in the node-local order instead of the natural order?
  3. Traffic: how many snapshot bytes cross nodes with the old reader
     (contiguous PETSc rows -> whole x-y planes) vs the new reader?

The algorithms here are ported 1:1 to node_local_map.F90 and
streamingSPODModule.F90, so a pass here is a pass for the logic there.

Run:  python tools/node_local_model.py      (needs numpy)
"""
import itertools
import numpy as np


# --------------------------------------------------------------------------
# 2decomp-fft v2.0.4 (the version Xcompact3d pins) pencil boxes
# --------------------------------------------------------------------------
def distribute(n, p, extra_first=False):
    """1-D split used by 2decomp. Returns 0-based starts and sizes.
    v2.0.4 always gives the extra points to the LAST ranks."""
    sz = [n // p] * p
    r = n % p
    rng = range(r) if extra_first else range(p - r, p)
    for i in rng:
        sz[i] += 1
    st = [0] * p
    for i in range(1, p):
        st[i] = st[i - 1] + sz[i - 1]
    return st, sz


def writer_box(b, nx, ny, nz, prow, pcol, pencil=1, extra_first=False):
    """Box (0-based start, count) written by X3D rank b.
    MPI_Cart_create(dims=[prow,pcol], reorder=.false.) is row-major:
    coord1 = b // pcol, coord2 = b % pcol."""
    c1, c2 = b // pcol, b % pcol
    if pencil == 1:      # x local, y over prow, z over pcol
        (sy, ny_), (sz, nz_) = distribute(ny, prow, extra_first), distribute(nz, pcol, extra_first)
        return (0, sy[c1], sz[c2]), (nx, ny_[c1], nz_[c2])
    if pencil == 2:      # x over prow, y local, z over pcol
        (sx, nx_), (sz, nz_) = distribute(nx, prow, extra_first), distribute(nz, pcol, extra_first)
        return (sx[c1], 0, sz[c2]), (nx_[c1], ny, nz_[c2])
    if pencil == 3:      # x over prow, y over pcol, z local
        (sx, nx_), (sy, ny_) = distribute(nx, prow, extra_first), distribute(ny, pcol, extra_first)
        return (sx[c1], sy[c2], 0), (nx_[c1], ny_[c2], nz)
    raise ValueError(pencil)


def writer_node(b, n_writers, n_nodes, per_node=0):
    """Node of X3D rank b when ranks are placed sequentially (block)."""
    if per_node > 0:
        return (b // per_node) % n_nodes
    base, extra = divmod(n_writers, n_nodes)
    if b < extra * (base + 1):
        return b // (base + 1)
    return extra + (b - extra * (base + 1)) // base


def split_range(total, nparts, part):
    """Balanced contiguous split, same rule as PETSC_DECIDE."""
    base, extra = divmod(total, nparts)
    lo = part * base + min(part, extra)
    return lo, lo + base + (1 if part < extra else 0)


def line_range_to_boxes(wstart, wcount, a, e):
    """Lines [a,e) of a writer box -> at most 3 boxes, in line order.
    A line is one run along dim 1 (wcount[0] points).
    Line l sits at (y = l % wcount[1], z = l // wcount[1]) inside the box."""
    c1, c2 = wcount[0], wcount[1]
    za, ya = divmod(a, c2)
    zb, yb = divmod(e - 1, c2)
    out = []
    if za == zb:
        out.append(((0, ya, za), (c1, yb - ya + 1, 1)))
    else:
        if ya > 0:                       # partial first plane
            out.append(((0, ya, za), (c1, c2 - ya, 1)))
            za += 1
        last = None
        if yb < c2 - 1:                  # partial last plane
            last = ((0, 0, zb), (c1, yb + 1, 1))
            zb -= 1
        if za <= zb:                     # full planes in between
            out.append(((0, 0, za), (c1, c2, zb - za + 1)))
        if last:
            out.append(last)
    return [(tuple(ws + s for ws, s in zip(wstart, st)), ct) for st, ct in out]


def node_local_boxes(my_node, node_rank, node_size, n_nodes, nx, ny, nz,
                     prow, pcol, pencil=1, extra_first=False, per_node=0):
    """Boxes this SPOD rank reads: the lines of all writers on its node,
    split evenly over the SPOD ranks of that node."""
    nw = prow * pcol
    writers = [b for b in range(nw) if writer_node(b, nw, n_nodes, per_node) == my_node]
    wb = [writer_box(b, nx, ny, nz, prow, pcol, pencil, extra_first) for b in writers]
    nlines = [c[1] * c[2] for _, c in wb]
    lo, hi = split_range(sum(nlines), node_size, node_rank)
    boxes, off = [], 0
    for b, (ws, wc), nl in zip(writers, wb, nlines):
        a, e = max(lo, off), min(hi, off + nl)
        if a < e:
            boxes += [(b, s, c) for s, c in line_range_to_boxes(ws, wc, a - off, e - off)]
        off += nl
    return boxes


def natural_indices(boxes, nvars, nx, ny, nz):
    """Natural row index (q = [var1; var2; ...], Fortran order inside a var)
    of every local row, in local order: var-major, then box, then k, j, i."""
    npts = nx * ny * nz
    out = []
    for v in range(nvars):
        for _, (sx, sy, sz), (cx, cy, cz) in boxes:
            k, j, i = np.meshgrid(np.arange(cz), np.arange(cy), np.arange(cx), indexing="ij")
            g = v * npts + (sx + i) + nx * ((sy + j) + ny * (sz + k))
            out.append(g.ravel())          # C order of (k,j,i) == Fortran order of (i,j,k)
    return np.concatenate(out) if out else np.zeros(0, dtype=np.int64)


def spod_node_ids(l, n, placement):
    """Physical node of each SPOD rank."""
    if placement == "block":
        return [writer_node(r, l, n) for r in range(l)]
    if placement == "cyclic":
        return [r % n for r in range(l)]
    raise ValueError(placement)


def node_info(phys):
    """What MPI_Comm_split_type + a leaders communicator give each rank:
    node index ordered by the lowest rank on that node, rank in node, node size."""
    order = []
    for p in phys:
        if p not in order:
            order.append(p)
    idx = {p: i for i, p in enumerate(order)}
    out, seen = [], {}
    for p in phys:
        q = seen.get(p, 0)
        seen[p] = q + 1
        out.append((idx[p], q))
    sizes = {idx[p]: phys.count(p) for p in order}
    return [(m, q, sizes[m]) for m, q in out], len(order)


# --------------------------------------------------------------------------
# 1. Mapping check
# --------------------------------------------------------------------------
def check_mapping(nx, ny, nz, prow, pcol, n, l, nvars=2, pencil=1,
                  extra_first=False, placement="block", per_node=0):
    nw = prow * pcol
    phys = spod_node_ids(l, n, placement)
    info, n_nodes = node_info(phys)
    # node index -> physical node (only equal to identity if launch order agrees)
    idx_to_phys = {}
    for (m, _, _), p in zip(info, phys):
        idx_to_phys[m] = p
    perm, rows_per_rank = [], []
    for r in range(l):
        m, q, s = info[r]
        boxes = node_local_boxes(m, q, s, n_nodes, nx, ny, nz, prow, pcol,
                                 pencil, extra_first, per_node)
        for b, st, ct in boxes:
            ws, wc = writer_box(b, nx, ny, nz, prow, pcol, pencil, extra_first)
            inside = all(ws[d] <= st[d] and st[d] + ct[d] <= ws[d] + wc[d] for d in range(3))
            assert inside, ("box outside writer", r, b, st, ct, ws, wc)
            wphys = writer_node(b, nw, n, per_node)
            assert wphys == idx_to_phys[m], ("remote read", r, b)
        g = natural_indices(boxes, nvars, nx, ny, nz)
        rows_per_rank.append(len(g))
        perm.append(g)
    perm = np.concatenate(perm)
    N = nvars * nx * ny * nz
    assert len(perm) == N, (len(perm), N)
    assert np.array_equal(np.sort(perm), np.arange(N)), "not a permutation"
    return perm, rows_per_rank


# --------------------------------------------------------------------------
# 2. Streaming SPOD, same steps as streamingSPODModule.F90
# --------------------------------------------------------------------------
def streaming_spod(snapshots, sqrtw, nDFT, ovlp, k, hamming=False):
    N = snapshots.shape[1]
    nOvlp = int(ovlp * nDFT)
    nFreq = nDFT // 2 + 1
    dn = nDFT - nOvlp
    t_idx = [1, 1 - dn]
    i = np.arange(nDFT)
    F = np.exp(-2j * np.pi * np.outer(i, np.arange(nFreq)) / nDFT)   # fft of identity rows
    F[:, 1:] *= 2
    w = (0.54 - 0.46 * np.cos(2 * np.pi * i / (nDFT - 1))) if hamming else np.ones(nDFT)
    winWeight = 1.0 / (w.sum() / nDFT)
    mu = np.zeros(N, complex)
    X = [np.zeros((N, nFreq), complex), np.zeros((N, nFreq), complex)]
    U_hat = np.zeros((N, k * nFreq), complex)
    S_hat = np.zeros((k, nFreq))
    state = {"block": 1}

    def handle_block(XH):
        state["block"] += 1
        cb = state["block"]
        XH = XH - np.outer(mu, w @ F)          # mean removal, row-wise
        XH = XH * (winWeight / nDFT)
        for fi in range(nFreq):
            U = U_hat[:, fi * k:(fi + 1) * k]
            x = XH[:, fi] * sqrtw
            Ux = U.conj().T @ x                 # global reduction
            up = x - U @ Ux
            abs_up = np.linalg.norm(up)         # global reduction
            unew = up / abs_up
            c = np.sqrt(cb / (cb + 1) ** 2)
            K = np.zeros((k + 1, k + 1), complex)
            K[np.arange(k), np.arange(k)] = c * np.sqrt(cb + 1) * S_hat[:, fi]
            K[:k, k] = c * Ux
            K[k, k] = c * abs_up
            Um, s, _ = np.linalg.svd(K)
            S_hat[:, fi] = s[:k]
            U_hat[:, fi * k:(fi + 1) * k] = (np.column_stack([U, unew]) @ Um)[:, :k]

    for x in snapshots:
        x = x.astype(complex)
        t = state.setdefault("t", 0) + 1
        state["t"] = t
        mu[:] = (1.0 / t) * ((t - 1.0) * mu + x)
        for b in range(2 if nOvlp > 0 else 1):
            if t_idx[b] > 0:
                X[b] += np.outer(x * w[t_idx[b] - 1], F[t_idx[b] - 1])
        for b in range(2 if nOvlp > 0 else 1):
            if t_idx[b] == nDFT:
                t_idx[b] = min(t_idx) - dn
                handle_block(X[b])
                X[b][:] = 0
            else:
                t_idx[b] += 1
    return S_hat, U_hat / sqrtw[:, None]


def check_spod_equivalence(nx=10, ny=9, nz=7, nvars=2, prow=2, pcol=3, n=2, l=5,
                           nsnap=31, nDFT=8, k=3, seed=0):
    perm, _ = check_mapping(nx, ny, nz, prow, pcol, n, l, nvars)
    N = nvars * nx * ny * nz
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((nsnap, N)) + np.sin(np.arange(nsnap))[:, None] * rng.standard_normal(N)
    wts = rng.uniform(0.5, 2.0, N)
    S_nat, U_nat = streaming_spod(q, np.sqrt(wts), nDFT, 0.5, k)
    # node-local run: row r of the PETSc vector holds natural point perm[r]
    S_loc, U_loc = streaming_spod(q[:, perm], np.sqrt(wts[perm]), nDFT, 0.5, k)
    U_back = np.empty_like(U_loc)
    U_back[perm] = U_loc                      # the reverse VecScatter
    ds = np.abs(S_nat - S_loc).max() / np.abs(S_nat).max()
    # A SPOD mode is only defined up to a factor e^{i*theta} (SVD phase), and
    # LAPACK may pick a different one after round-off-level changes (the old
    # code does the same when you change the number of ranks). Align it first.
    ip = np.einsum("ij,ij->j", U_nat.conj(), U_back)
    phase = np.where(np.abs(ip) > 0, ip / np.where(ip == 0, 1, np.abs(ip)), 1)
    du = np.abs(U_nat * phase - U_back).max() / np.abs(U_nat).max()
    nflip = int(np.sum(np.abs(phase - 1) > 1e-6))
    return ds, du, nflip


# --------------------------------------------------------------------------
# 3. Inter-node traffic: old reader vs node-local reader
# --------------------------------------------------------------------------
def box_overlap(s1, c1, s2, c2):
    v = 1
    for d in range(3):
        lo, hi = max(s1[d], s2[d]), min(s1[d] + c1[d], s2[d] + c2[d])
        if hi <= lo:
            return 0
        v *= hi - lo
    return v


def old_reader_boxes(r, l, nvars, nx, ny, nz):
    """spod-adios2.txt: PETSC_DECIDE rows, whole x-y planes per variable."""
    npts, plane = nx * ny * nz, nx * ny
    rstart, rend = split_range(nvars * npts, l, r)
    out = []
    for v in range(nvars):
        lo, hi = max(rstart, v * npts), min(rend - 1, (v + 1) * npts - 1)
        if lo <= hi:
            k0 = (lo - v * npts) // plane
            k1 = (hi - v * npts) // plane
            out.append(((0, 0, k0), (nx, ny, k1 - k0 + 1)))
    return out


def traffic(nx, ny, nz, nvars, prow, pcol, n, l):
    nw = prow * pcol
    wboxes = [writer_box(b, nx, ny, nz, prow, pcol) for b in range(nw)]
    wnode = [writer_node(b, nw, n) for b in range(nw)]
    phys = spod_node_ids(l, n, "block")
    info, n_nodes = node_info(phys)
    res = {}
    for name in ("old", "new"):
        same = cross = 0
        writers_per_reader = []
        for r in range(l):
            if name == "old":
                boxes = [(None, s, c) for s, c in old_reader_boxes(r, l, nvars, nx, ny, nz)]
                mult = 1
            else:
                m, q, s = info[r]
                boxes = node_local_boxes(m, q, s, n_nodes, nx, ny, nz, prow, pcol)
                mult = nvars
            touched = set()
            for _, s, c in boxes:
                for b in range(nw):
                    ov = box_overlap(s, c, *wboxes[b])
                    if ov:
                        touched.add(b)
                        if wnode[b] == phys[r]:
                            same += ov * mult
                        else:
                            cross += ov * mult
            writers_per_reader.append(len(touched))
        res[name] = (same * 8, cross * 8, max(writers_per_reader))
    return res


if __name__ == "__main__":
    print("== 1. mapping checks (tiling, inside one writer box, same node) ==")
    cases = []
    for (nx, ny, nz), (prow, pcol), (n, l), pencil, xf, plc in itertools.product(
            [(10, 9, 7), (17, 13, 11), (8, 8, 8)],
            [(1, 1), (2, 3), (3, 2), (4, 4)],
            [(1, 1), (2, 4), (2, 5), (3, 7), (4, 32)],
            [1, 2, 3], [False, True], ["block", "cyclic"]):
        if prow > min(nx, ny) or pcol > min(ny, nz):
            continue
        check_mapping(nx, ny, nz, prow, pcol, n, l, 2, pencil, xf, plc)
        cases.append(1)
    # explicit X3D ranks per node (x3dRanksPerNode > 0)
    for (nx, ny, nz), (prow, pcol), (n, l) in [((17, 13, 11), (4, 4), (4, 8)), ((10, 9, 7), (2, 3), (3, 6)),
                                               ((16, 16, 16), (4, 8), (2, 6))]:
        check_mapping(nx, ny, nz, prow, pcol, n, l, 2, per_node=prow * pcol // n)
        cases.append(1)
    # more SPOD ranks than x-lines on a node: some ranks own 0 rows, still exact
    _, rows = check_mapping(3, 2, 2, 1, 1, 1, 7, 2)
    assert min(rows) == 0
    cases.append(1)
    print(f"   {len(cases)} configurations: all PASS (incl. explicit ranks-per-node and 0-row ranks)")

    print("== 2. SPOD with node-local row order vs natural row order ==")
    for cfg in [dict(), dict(l=4, n=2, prow=3, pcol=2), dict(l=7, n=3, prow=2, pcol=2, nsnap=40, nDFT=10, k=2)]:
        ds, du, nflip = check_spod_equivalence(**cfg)
        print(f"   {cfg or 'default'}: max rel diff energies {ds:.1e}, "
              f"modes {du:.1e} (after phase alignment; {nflip} mode(s) had a different phase)")
        assert ds < 1e-10 and du < 1e-9

    print("== 3. snapshot bytes per step that cross nodes ==")
    for (nx, ny, nz, nv, prow, pcol, n, l) in [(301, 301, 301, 4, 8, 16, 4, 64),
                                               (301, 301, 301, 4, 16, 16, 8, 128),
                                               (256, 128, 256, 4, 4, 16, 2, 32)]:
        t = traffic(nx, ny, nz, nv, prow, pcol, n, l)
        tot = sum(t["old"][:2])
        print(f"   grid {nx}x{ny}x{nz} x{nv} vars, X3D {prow}x{pcol}={prow*pcol} ranks, "
              f"{n} nodes, SPOD {l} ranks  (step = {tot/2**20:.0f} MiB)")
        for name in ("old", "new"):
            same, cross, mw = t[name]
            print(f"      {name}: cross-node {cross/2**20:8.1f} MiB ({100*cross/tot:5.1f}%), "
                  f"same-node {same/2**20:8.1f} MiB, max writers per SPOD rank {mw}")
