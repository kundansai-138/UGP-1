#!/usr/bin/env bash
# =====================================================================
# setup_adios_case.sh - prepare the ADIOS2 version of your case, end to end.
#
#   A. Build  ~/ugp/src/x3d_lab_adios  (Xcompact3d with IO=adios2)
#        1. int_to_str fix in src/visu.f90   (the "End of record" crash)
#        2. snapshot-0 guard in src/case.f90 (ivisu, same as the MPI build)
#        3. SST support: decomp2d/io.f90 knows only BP4/HDF5 names -> add SST/BP5;
#           snapshot writes in visu.f90 -> sync (deferred writes read the array only at end_step)
#        4. make clean; make IO=adios2 ADIOS2DIR=~/ugp/install; check the binary
#   B. Case folder  ~/ugp/runs/testcase_adios  from the master case (inputs only)
#        same input.i3d as everywhere (ilast 150, ioutput 5), plus the settings
#        S0 used: ilist = 1 (step timing), icheckpoint = ilast+1 (no restart writes)
#   C. ADIOS2 configs in the case folder
#        adios2_config.bp4.xml  all IOs -> BP4 files
#        adios2_config.sst.xml  solution-io -> SST stream "data" (for SPOD), rest BP4
#        adios2_config.xml      the active one (ENGINE=bp4 or sst)
#   D. Checks: binary, IO names used by the code, open mode of solution-io
#
# Usage:   bash setup_adios_case.sh
# Options: ENGINE=sst|bp4 (active config, default sst)   SKIP_BUILD=1 (only B-D)
#          SRC=~/ugp/src/x3d_lab_adios  MASTER=~/ugp/runs/testcase_mpi  CASE=~/ugp/runs/testcase_adios
#          ADIOS2DIR=~/ugp/install  MPIFC=/usr/local/bin/mpif90
# =====================================================================
set -euo pipefail

SRC=${SRC:-$HOME/ugp/src/x3d_lab_adios}
MASTER=${MASTER:-$HOME/ugp/runs/testcase_mpi}
CASE=${CASE:-$HOME/ugp/runs/testcase_adios}
ADIOS2DIR=${ADIOS2DIR:-$HOME/ugp/install}
# MPI Fortran wrapper. Not taken from $FC: ~/ugp/env.sh exports FC=gfortran-11 (plain compiler, no mpi.mod)
MPIFC=${MPIFC:-/usr/local/bin/mpif90}
ENGINE=${ENGINE:-sst}
SKIP_BUILD=${SKIP_BUILD:-0}
TOOLS=$(cd "$(dirname "$0")" && pwd)

say() { printf '\n==> %s\n' "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
[[ -d $SRC/src && -f $SRC/Makefile ]]   || die "no Xcompact3d source tree in $SRC"
[[ -f $MASTER/input.i3d ]]              || die "no input.i3d in $MASTER"
[[ -x $ADIOS2DIR/bin/adios2-config ]]   || die "ADIOS2 not found in $ADIOS2DIR"
[[ -x $MPIFC ]]                         || die "MPI Fortran wrapper not found: $MPIFC"
[[ $ENGINE == sst || $ENGINE == bp4 ]]  || die "ENGINE must be sst or bp4"

# ===================================================================== A. build
if [[ $SKIP_BUILD != 1 ]]; then
  say "A1. int_to_str fix in $SRC/src/visu.f90"
  python3 - "$SRC/src/visu.f90" <<'PY'
import sys, shutil
p = sys.argv[1]; s = open(p).read()
old = '''    character(len=(1 + int(log10(real(max(i, 1)))))) :: int_to_str

    write(int_to_str, "(I0)") i'''
new = '''    character(len=:), allocatable :: int_to_str
    character(len=32) :: buf

    write(buf, "(I0)") i
    int_to_str = trim(buf)'''
if "character(len=32) :: buf" in s:
    print("already fixed")
elif s.count(old) == 1:
    shutil.copy(p, p + ".orig"); open(p, "w").write(s.replace(old, new)); print("fixed (backup visu.f90.orig)")
else:
    sys.exit("int_to_str block not found in visu.f90 - send me: grep -n -A8 'function int_to_str' src/visu.f90")
PY

  say "A2. snapshot-0 guard in $SRC/src/case.f90"
  [[ -f $TOOLS/patch_snapshot.py ]] || die "patch_snapshot.py must be next to this script"
  python3 "$TOOLS/patch_snapshot.py" "$SRC/src/case.f90"

  say "A3. SST support in decomp2d/io.f90 + sync snapshot writes in src/visu.f90"
  python3 - "$SRC" <<'PY'
import re, shutil, sys, os
src = sys.argv[1]
# (a) old 2decomp only knows BP4/HDF5 output names -> add BP5 and SST (stream name = io_dir, no extension)
p = os.path.join(src, 'decomp2d/io.f90'); s = open(p).read()
if '"SST"' in s:
    print('io.f90: SST already supported')
else:
    old = """    else if (io%engine_type .eq. "HDF5") then
       ext = ".hdf5"
    else"""
    new = """    else if (io%engine_type .eq. "HDF5") then
       ext = ".hdf5"
    else if (io%engine_type .eq. "BP5") then
       ext = ".bp5"
    else if (io%engine_type .eq. "SST" .or. io%engine_type .eq. "Sst" .or. io%engine_type .eq. "sst") then
       ext = ""        ! SST stream: name = io_dir ("data"), contact file data.sst
    else"""
    if s.count(old) != 1:
        sys.exit("io.f90: engine-type block not found - send me: grep -n -B10 -A5 'Unkown engine type' decomp2d/io.f90")
    shutil.copy(p, p + '.pre_sst'); open(p, 'w').write(s.replace(old, new)); print('io.f90: added BP5 and SST')
# (b) full-resolution snapshot writes were deferred (array read only at end_step) -> make them sync
p = os.path.join(src, 'src/visu.f90'); L = open(p).read().split('\n')
idx = [k for k, l in enumerate(L) if 'decomp_2d_write_one(1,f1,"data"' in l and not l.strip().startswith('!')]
if not idx:
    sys.exit('visu.f90: write_one(1,f1,"data"...) line not found')
changed = 0
for i in idx:
    if 'opt_deferred_writes' not in L[i]:
        L[i] = re.sub(r',io_name\)\s*$', ',io_name,opt_deferred_writes=.false.)', L[i])
        if 'opt_deferred_writes' not in L[i]:
            sys.exit('visu.f90: could not change line %d: %s' % (i + 1, L[i]))
        changed += 1
if changed:
    shutil.copy(p, p + '.pre_sync'); open(p, 'w').write('\n'.join(L)); print('visu.f90: %d write(s) -> sync' % changed)
else:
    print('visu.f90: snapshot writes already sync')
PY
  n_ready=$(cat "$SRC"/src/*.f90 | grep -c 'call visu_ready' || true)
  [[ $n_ready == 1 ]] || echo "WARNING: 'call visu_ready' found $n_ready times - solution-io must be opened once for SST"

  say "A4. build with IO=adios2 (serial make, ~5 min)"
  ( cd "$SRC" && make clean >/dev/null 2>&1 || true
    echo "compiler: $MPIFC -> $("$MPIFC" -show 2>/dev/null | cut -d' ' -f1)" > build.log
    make CMP=gcc FC="$MPIFC" IO=adios2 ADIOS2DIR="$ADIOS2DIR" >> build.log 2>&1 ) \
    || die "build failed - see $SRC/build.log (grep -n -i error $SRC/build.log | head)"
fi
BIN=$SRC/xcompact3d
[[ -x $BIN ]] || die "no binary $BIN"
grep -q -- '-DADIOS2' "$SRC/build.log" 2>/dev/null || echo "WARNING: -DADIOS2 not seen in $SRC/build.log"
n=$(nm "$BIN" | grep -ci adios2_init || true)
[[ $n -gt 0 ]] || die "binary has no ADIOS2 (adios2_init symbols: 0)"
ldd "$BIN" | grep -q 'not found' && die "missing libraries: $(ldd "$BIN" | grep 'not found')"
echo "binary ok: $BIN  (adios2_init symbols: $n; MPI: $(ldd "$BIN" | awk '/libmpi.so/{print $3}'))"

# ===================================================================== B. case folder
say "B. case folder $CASE"
if [[ -e $CASE ]]; then
  OLD=$CASE.old.$(date +%Y%m%d_%H%M%S); mv "$CASE" "$OLD"; echo "existing folder moved to $OLD"
fi
mkdir -p "$CASE"
tar -C "$MASTER" --exclude='./data' --exclude='./probes' --exclude='./statistics' \
    --exclude='./streaklines' --exclude='./out' --exclude='./prof_*' --exclude='*.log' \
    --exclude='*.log1' --exclude='*.bin' --exclude='*.bp4' --exclude='*.bp5' --exclude='*.sst' \
    --exclude='restart*' --exclude='checkpoint*' --exclude='adios2_config*.xml' -cf - . | tar -C "$CASE" -xf -
cd "$CASE"
for d in data probes statistics streaklines; do mkdir -p "$d"; done
cp input.i3d input.i3d.master

getv() { grep -ioE "(^|[ ,])$1\s*=\s*[-0-9]+" "$2" | head -1 | sed -E 's/.*=\s*//'; }
setv() {
  grep -qiE "(^|[ ,])$1\s*=" "$3" || die "$1 not found in $3"
  sed -i -E "s/(^|[ ,])($1\s*=\s*)[-0-9]+/\1\2$2/I" "$3"
  [[ $(getv "$1" "$3") == "$2" ]] || die "could not set $1=$2 in $3"
}
IL=$(getv ilast input.i3d); IO=$(getv ioutput input.i3d)
setv icheckpoint $((IL + 1)) input.i3d     # same as S0: no restart writes
setv ilist 1 input.i3d                     # same as S0: one step line per step
echo "ilast = $IL, ioutput = $IO  ->  snapshots at steps $IO, $((2*IO)), ..., plus snapshot 0 at start-up"
echo "changes to input.i3d (vs master):"; diff input.i3d.master input.i3d || true

# ===================================================================== C. ADIOS2 configs
say "C. ADIOS2 configs"
IONAMES=$(grep -h -o -i '"[a-z0-9_-]*-io"' "$SRC"/src/*.f90 "$SRC"/decomp2d/*.f90 2>/dev/null | tr -d '"' | sort -u)
[[ -n $IONAMES ]] || IONAMES="solution-io restart-io restart-forces-io in-outflow-io statistics-io turb-io BC-lock-exchange-io"
echo "IO names in the code: $(echo $IONAMES)"

bp4_io() { printf '  <io name="%s">\n    <engine type="BP4"/>\n  </io>\n' "$1"; }
{
  echo '<?xml version="1.0"?>'
  echo '<!-- all output as BP4 files (written by setup_adios_case.sh) -->'
  echo '<adios-config>'
  for n in $IONAMES; do bp4_io "$n"; done
  echo '</adios-config>'
} > adios2_config.bp4.xml
{
  echo '<?xml version="1.0"?>'
  echo '<!-- snapshots (solution-io) as SST stream "data" for streaming SPOD; everything else BP4 -->'
  echo '<adios-config>'
  cat <<'EOF'
  <io name="solution-io">
    <engine type="SST">
      <!-- the solver waits at start until 1 reader (SPOD) has connected -->
      <parameter key="RendezvousReaderCount" value="1"/>
      <!-- at most 2 snapshots buffered on the writer; then the solver waits for SPOD -->
      <parameter key="QueueLimit" value="2"/>
      <parameter key="QueueFullPolicy" value="Block"/>
      <!-- TCP sockets over the normal 1 GbE network -->
      <parameter key="DataTransport" value="WAN"/>
      <!-- every csews machine also has docker0: advertise and use eno1 only -->
      <parameter key="NetworkInterface" value="eno1"/>
      <!-- reader side: wait up to 1 h for the writer to appear -->
      <parameter key="OpenTimeoutSecs" value="3600"/>
    </engine>
  </io>
EOF
  for n in $IONAMES; do [[ $n == solution-io ]] || bp4_io "$n"; done
  echo '</adios-config>'
} > adios2_config.sst.xml
cp "adios2_config.$ENGINE.xml" adios2_config.xml
echo "written: adios2_config.bp4.xml, adios2_config.sst.xml; active: adios2_config.xml = $ENGINE"

# ===================================================================== D. checks
say "D. checks"
echo "how solution-io is opened (SST cannot use append mode):"
grep -n -i -E 'open_io|adios2_mode|decomp_2d_(write|append)_mode' "$SRC/src/visu.f90" | head -8 || true
IR=$(getv irestart input.i3d || true)
[[ ${IR:-0} == 0 ]] && echo "irestart = 0: fresh start, write mode (fine for SST)" \
                    || echo "WARNING: irestart = $IR - a restarted run may open snapshots in append mode, which SST does not support"
cat <<EOF

Done. Case folder: $CASE
  binary : $BIN
  input  : ilast $IL, ioutput $IO, ilist 1, no restart writes
  ADIOS2 : $(readlink -f adios2_config.xml) (active = $ENGINE)
Switch engine:  cp adios2_config.bp4.xml adios2_config.xml   (or .sst.xml)
SST stream name: data  -> contact file $CASE/data.sst while the solver runs
EOF
