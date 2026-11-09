"""
Compare the outputs of two SPOD runs on the same data, e.g. the current
reader (spod-adios2.txt) and the node-local reader.

    python tools/compare_spod_outputs.py RUN_DIR_A RUN_DIR_B

Compares
  * EnegySpectra.txt   (singular values per mode and frequency)
  * SPODMODE_freq_<i>_sigma_*.dat  (PETSc binary vectors, natural order)

A SPOD mode is only defined up to a factor e^{i*theta}, so each mode pair
is compared after removing that phase. Expect ~1e-12 or smaller; the energy
file is printed by PETSc with 6 digits, so expect ~1e-6 there.
"""
import glob
import os
import re
import sys

import numpy as np

VEC_FILE_CLASSID = 1211214


def read_petsc_vec(path):
    raw = open(path, "rb").read()
    for itype in (">i4", ">i8"):
        isz = np.dtype(itype).itemsize
        cid, n = np.frombuffer(raw[:2 * isz], dtype=itype)
        if cid == VEC_FILE_CLASSID and len(raw) == 2 * isz + 16 * n:
            return np.frombuffer(raw[2 * isz:], dtype=">c16").astype(complex)
    raise ValueError(f"{path}: not a complex PETSc binary vector")


def read_spectra(path):
    """MatView ASCII of S_hat: 'row r: (c, value) ...' -> dict[(r, c)] = real part."""
    out = {}
    for line in open(path):
        m = re.match(r"\s*row\s+(\d+):(.*)", line)
        if not m:
            continue
        r = int(m.group(1))
        for c, val in re.findall(r"\((\d+),\s*([^)]*)\)", m.group(2)):
            num = re.match(r"\s*([-+]?[0-9.]+(?:[eE][-+]?\d+)?)", val)
            out[(r, int(c))] = float(num.group(1))
    return out


def modes_by_index(run_dir):
    files = {}
    for f in glob.glob(os.path.join(run_dir, "SPODMODE_freq_*_sigma_*.dat")):
        i = int(re.search(r"_freq_(\d+)_", os.path.basename(f)).group(1))
        files[i] = f
    return files


def main(a, b):
    ok = True
    sa, sb = read_spectra(os.path.join(a, "EnegySpectra.txt")), read_spectra(os.path.join(b, "EnegySpectra.txt"))
    keys = sorted(set(sa) | set(sb))
    va = np.array([sa.get(key, np.nan) for key in keys])
    vb = np.array([sb.get(key, np.nan) for key in keys])
    ds = np.nanmax(np.abs(va - vb)) / np.nanmax(np.abs(va))
    print(f"EnegySpectra.txt: {len(keys)} values, max rel diff {ds:.2e}")
    ok &= bool(ds < 1e-5) and not np.isnan(va - vb).any()

    ma, mb = modes_by_index(a), modes_by_index(b)
    for i in sorted(set(ma) | set(mb)):
        if i not in ma or i not in mb:
            print(f"mode {i}: missing in one run")
            ok = False
            continue
        x, y = read_petsc_vec(ma[i]), read_petsc_vec(mb[i])
        if x.shape != y.shape:
            print(f"mode {i}: lengths differ {x.shape} vs {y.shape}")
            ok = False
            continue
        ip = np.vdot(x, y)
        cosang = abs(ip) / (np.linalg.norm(x) * np.linalg.norm(y))
        phase = ip / abs(ip) if abs(ip) > 0 else 1.0
        d = np.abs(x * phase - y).max() / np.abs(x).max()
        print(f"mode {i}: |cos angle| = {cosang:.15f}, max rel diff after phase alignment {d:.2e}"
              f"  ({os.path.basename(ma[i])} vs {os.path.basename(mb[i])})")
        ok &= bool(d < 1e-8)
    print("SAME RESULTS" if ok else "RESULTS DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
