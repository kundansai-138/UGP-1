#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
PY=/usr/bin/python3
B=$U/build/spod
W=$U/build/spod-env
N=20000; NSNAP=96; NDFT=32; NP=4
D=$U/data/synth_n${N}_t${NSNAP}_d${NDFT}

[ -x $B/streamingSPOD ] || { echo "STOP: run step4.sh first"; exit 1; }

# ---- 1. patched copy: inputs overridable by environment variables
rm -rf $W; mkdir -p $W $U/tools
cp $B/streamingSPODModule.F90 $B/Makefile $W/
cat > $U/tools/patch_userinput.py <<'PY1'
import re, sys
path = sys.argv[1]
with open(path, newline="") as fh:
    src = fh.read()
nl = "\r\n" if "\r\n" in src else "\n"
if "UGP: run-time overrides" in src:
    print("already patched:", path); sys.exit(0)
decl_re = re.compile(r"(Subroutine\s+UserInput[ \t]*\r?\n[ \t]*implicit[ \t]+none[ \t]*\r?\n)", re.I)
end_re = re.compile(r"(\r?\n[ \t]*End[ \t]+Subroutine[ \t]+UserInput)", re.I)
if len(decl_re.findall(src)) != 1 or len(end_re.findall(src)) != 1:
    sys.exit("STOP: UserInput layout not recognised; patch not applied")
pad = " " * 32
decl = pad + "character(len=500) :: ugp_tmp" + nl + pad + "integer :: ugp_st" + nl
lines = ["", pad + "! ---- UGP: run-time overrides from environment variables ----"]
for var, target in [("SPOD_NPOINTS", "numberOfDataInSnapShot"), ("SPOD_NSNAP", "numberOfSnapShot"),
                    ("SPOD_NDFT", "nDFT"), ("SPOD_K", "k"), ("SPOD_NFREQSAVE", "nFreqToSave"),
                    ("SPOD_OVERLAP", "percentageOfOverall")]:
    lines.append(pad + "call get_environment_variable('%s', ugp_tmp, status=ugp_st)" % var)
    lines.append(pad + "if (ugp_st == 0) read(ugp_tmp, *) %s" % target)
lines += [
    pad + "call get_environment_variable('SPOD_DATADIR', ugp_tmp, status=ugp_st)",
    pad + "if (ugp_st == 0) filenameFormat = '(\"' // trim(ugp_tmp) // '/\", I0, \".bin\")'",
    pad + "call get_environment_variable('SPOD_WEIGHTS', ugp_tmp, status=ugp_st)",
    pad + "if (ugp_st == 0) weightMatrixFilePath = trim(ugp_tmp)",
    pad + "if (trim(weightMatrixFilePath) == 'none') weightMatrixFilePath = ''",
    ""]
block = nl.join(lines)
src = decl_re.sub(lambda m: m.group(1) + decl, src, count=1)
src = end_re.sub(lambda m: block + m.group(1), src, count=1)
with open(path, "w", newline="") as fh:
    fh.write(src)
print("patched:", path)
PY1
$PY $U/tools/patch_userinput.py $W/streamingSPODModule.F90
(cd $W && make > build.log 2>&1) || { tail -30 $W/build.log; echo "STOP: patched build failed"; exit 1; }
[ -x $W/streamingSPOD ] || { echo "STOP: no patched executable"; exit 1; }
echo "patched build ok"

# ---- 2. synthetic snapshots with a known answer (bins 4 10 14), raw float64 like Example/*.bin
if [ ! -f $D/done ]; then
  mkdir -p $D
  $PY - $D $N $NSNAP $NDFT <<'PY2'
import sys, math, random, array
out, N, T, ndft = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
bins = [int(round(ndft * f)) for f in (0.125, 0.3125, 0.4375)]
amps = [3.0, 2.0, 1.2]; kap = [1, 2, 3]
random.seed(12345)
cb = [[math.cos(2*math.pi*k*i/N) for i in range(N)] for k in kap]
sb = [[math.sin(2*math.pi*k*i/N) for i in range(N)] for k in kap]
for t in range(1, T + 1):
    f = [0.05 * random.gauss(0, 1) for _ in range(N)]
    for m, a, c, s in zip(bins, amps, cb, sb):
        ang = 2*math.pi*(m/ndft)*(t-1)
        ca, sa = math.cos(ang), math.sin(ang)
        for i in range(N):
            f[i] += a * (ca*c[i] + sa*s[i])
    with open("%s/%d.bin" % (out, t), "wb") as fh:
        array.array("d", f).tofile(fh)
print("wrote %d snapshots, expected bins %s" % (T, bins))
PY2
  touch $D/done
fi

# ---- 3. run and check
R=$U/runs/smoke-$(date +%Y%m%d-%H%M%S); mkdir -p $R; cd $R
export SPOD_NPOINTS=$N SPOD_NSNAP=$NSNAP SPOD_NDFT=$NDFT SPOD_K=3 SPOD_OVERLAP=0.5 SPOD_NFREQSAVE=3
export SPOD_DATADIR=$D SPOD_WEIGHTS=none
env -u DISPLAY -u XAUTHORITY mpiexec -n $NP $W/streamingSPOD > spod.log 2>&1 \
  || { echo "RUN FAILED"; tail -30 spod.log; exit 1; }
got=$(grep "Eigen values at index" spod.log | awk '{print $5}' | tr '\n' ' ')
echo "run dir : $R"
echo "expected: 4 10 14"
echo "got     : $got"
ls SPODMODE_freq_*.dat EnegySpectra.txt
if [ "$got" = "4 10 14 " ]; then echo "SMOKE TEST: PASS"; else echo "SMOKE TEST: DIFFERENT (paste: tail -40 $R/spod.log)"; fi
