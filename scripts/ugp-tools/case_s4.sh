#!/usr/bin/env bash
# =====================================================================
# run_S4.sh - S4: Xcompact3d (ADIOS2 build) streams snapshots over SST straight into streaming SPOD.
# Same layout as basecase.sh (S0), so the two summaries compare line by line.
#
#   writer : x3d_lab_adios/xcompact3d   hostfile       (csews160..167, 8 x 24 = 192 ranks, under APS)
#   profiling: Intel APS on both sides (APS=1, APS_SPOD=1) -> aps_report.txt, aps_spod_report.txt
#   reader : streamingSPOD (spod-adios2) hostfile_spod  (csews168..171, 4 x RPPN ranks: 12 -> 48, 8 -> 32, 4 -> 16)
#   input  : ilast = 150, ioutput = 5   -> 31 snapshots (step 0, 5, ..., 150), ilist = 1, no restart writes
#   SPOD   : 31 snapshots, nDFT 8, overlap 4 -> 6 blocks (real values are hard-coded in streamingSPODModule.F90)
#
# Steps
#   1. settings and checks
#   2. fresh run folder from testcase_adios (inputs + adios2_config.sst.xml; SPOD uses adios2_config.xml too)
#   3. machine information of all 12 machines (meta.txt)
#   4. REPEATS runs: writer starts -> data.sst appears -> SPOD starts -> both finish
#      every run: x3d_ts.log, spod_ts.log, node_<host>.csv, APS raw + reports (SPOD outputs deleted, KEEP_SPOD=1 keeps)
#   5. summary over all runs (summary.csv / summary.txt)
#   6. reader/writer wait per snapshot (libadios2_wait.so times every begin_step/end_step per rank)
#      -> run<N>/wait/*.csv, wait_steps.csv, wait_summary.csv, wait_summary.txt
#
# Usage:   bash run_S4.sh
# Options: ILAST=150 IOUTPUT=5 REPEATS=1  WPPN=24 RPPN=12  APS=1 (writer under APS)
#          reader sweep: RPPN=8 bash run_S4.sh ; RPPN=4 bash run_S4.sh   (-ppn overrides the :12 in hostfile_spod)
#          SPOD_NDFT=16 SPOD_OVERLAP=8 SPOD_K=3 SPOD_NSAVE=2 SPOD_VARS="ux uy uz pp"
#          S0_STEP=<S0 median step s>  (prints the S4 overhead per step)
#          short test: ILAST=10 SPOD_NDFT=2 SPOD_OVERLAP=1 REPEATS=1
#          WAIT=0 (no reader wait timing)  WAIT_WRITER=0 (no writer-side timing)
#          APS_IMBALANCE_TYPE=2 (default: MPI imbalance via a barrier before each collective; works with MPICH)
#          APS_IMBALANCE_TYPE=0 (off)
#          APS_STAT_LEVEL=2 (default: MPI time + data moved per MPI function; 1 = time only, 3-5 = more detail)
#          Disk: after each run only APS (raw folders + text reports) and the analysis are kept (summaries, logs,
#          wait CSVs, node CSVs). Removed: SPOD outputs run<N>/spod/ (spectrum + mode files, ~1.75 GB per mode)
#          and the Xcompact3d output folders (after the last run).
#          KEEP_SPOD=1 keeps the spectrum (+ KEEP_MODES=N leading modes of run1); KEEP_X3D_OUT=1 keeps the
#          Xcompact3d outputs; KEEP_APS_RAW=0 removes the raw APS folders once the reports are written
#          WPIN=0-19 RPIN=20-23  CPU lists for every writer / reader rank (taskset -c), e.g. SPOD on the writer
#          machines: WPPN=20 RPPN=4 RHOSTS=~/ugp/tools/hostfile WPIN=0-19 RPIN=20-23 (check lscpu -e first)
# Compare layouts afterwards: python3 wait_analysis.py ~/ugp/runs/S4/S4_r4x12_* ~/ugp/runs/S4/S4_r4x8_* ...
# =====================================================================
set -euo pipefail

# ---------------------------------------------------------------- 1. settings
CASE=${CASE:-$HOME/ugp/runs/testcase_adios}           # from setup_adios_case.sh
BIN=${BIN:-$HOME/ugp/src/x3d_lab_adios/xcompact3d}    # writer (IO=adios2 build with SST support)
SBIN=${SBIN:-$HOME/ugp/build/spod-adios2/streamingSPOD}   # reader (streaming SPOD, ADIOS2 input)
WHOSTS=${WHOSTS:-$HOME/ugp/tools/hostfile}            # csews160..167
RHOSTS=${RHOSTS:-$HOME/ugp/tools/hostfile_spod}       # csews168..171
WPPN=${WPPN:-24}; RPPN=${RPPN:-12}
WPIN=${WPIN:-}; RPIN=${RPIN:-}                         # CPU lists (taskset -c); empty = no pinning
REPEATS=${REPEATS:-1}                                  # runs per case (REPEATS=2 for a repeat)
ILAST=${ILAST:-150}; IOUTPUT=${IOUTPUT:-5}
APS=${APS:-1}                                         # Xcompact3d under Intel APS
APS_SPOD=${APS_SPOD:-$APS}                            # SPOD under Intel APS
APS_IMB=${APS_IMBALANCE_TYPE:-2}                      # 2 = MPI_Barrier before each collective (type 1 needs Intel MPI)
APS_LEVEL=${APS_STAT_LEVEL:-2}                        # 2 = also data moved per MPI function (APS default 1 = time only)
SPOD_NDFT=${SPOD_NDFT:-8}; SPOD_OVERLAP=${SPOD_OVERLAP:-4}
SPOD_K=${SPOD_K:-3};        SPOD_NSAVE=${SPOD_NSAVE:-2}
SPOD_VARS=${SPOD_VARS:-ux uy uz pp}
S0_STEP=${S0_STEP:-}
MPIRUN=${MPIRUN:-/usr/local/bin/mpirun}
APS_ENV=${APS_ENV:-/opt/intel/oneapi/vtune/2025.7/env/vars.sh}
RUNROOT=${RUNROOT:-$HOME/ugp/runs/S4}
SPOD_XML=${SPOD_XML:-adios2_config.xml}               # SPOD reads the same config as the solver (solution-io = SST)
SST_WAIT=${SST_WAIT:-600}                             # s to wait for data.sst after writer start
NET=(-iface eno1 -genv FI_PROVIDER tcp -genv FI_TCP_IFACE eno1 -genv UCX_LOG_LEVEL error)
SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WAIT=${WAIT:-1}                                       # time SPOD's begin_step/end_step (reader wait)
WAIT_WRITER=${WAIT_WRITER:-$WAIT}                     # same on the Xcompact3d side (writer blocked by SPOD)
WAIT_SHIM=${WAIT_SHIM:-$SELF_DIR/libadios2_wait.so}   # built from adios2_wait_shim.c if missing
WAIT_PY=${WAIT_PY:-$SELF_DIR/wait_analysis.py}
KEEP_SPOD=${KEEP_SPOD:-0}                             # 1 = keep run<N>/spod (spectrum; modes per KEEP_MODES)
KEEP_MODES=${KEEP_MODES:-0}                           # with KEEP_SPOD=1: leading mode files kept in run1
KEEP_X3D_OUT=${KEEP_X3D_OUT:-0}                       # 1 = keep Xcompact3d output folders after the last run
KEEP_APS_RAW=${KEEP_APS_RAW:-1}                       # 0 = remove raw APS folders once the reports exist

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -f $CASE/input.i3d ]]             || die "no input.i3d in $CASE (run setup_adios_case.sh)"
[[ -f $CASE/adios2_config.sst.xml ]] || die "no adios2_config.sst.xml in $CASE (run setup_adios_case.sh)"
[[ -x $BIN ]]   || die "writer binary not found: $BIN"
[[ -x $SBIN ]]  || die "SPOD binary not found: $SBIN"
[[ -f $WHOSTS ]] || die "hostfile not found: $WHOSTS"
[[ -f $RHOSTS ]] || die "hostfile not found: $RHOSTS"
[[ -x $MPIRUN ]] || die "mpirun not found: $MPIRUN"
command -v perl >/dev/null    || die "perl is needed for timestamps"
command -v python3 >/dev/null || die "python3 is needed for the summary"
(( ILAST % IOUTPUT == 0 )) || echo "NOTE: ilast is not a multiple of ioutput"
if [[ $APS == 1 || $APS_SPOD == 1 ]]; then
  [[ -f $APS_ENV ]] || die "APS environment not found: $APS_ENV"
  # shellcheck disable=SC1090
  source "$APS_ENV"; command -v aps >/dev/null || die "aps not found after sourcing $APS_ENV"
fi
if [[ -n $WPIN || -n $RPIN ]]; then command -v taskset >/dev/null || die "taskset is needed for WPIN/RPIN"; fi
if [[ $WAIT == 1 || $WAIT_WRITER == 1 ]]; then
  if [[ ! -f $WAIT_SHIM && -f $SELF_DIR/adios2_wait_shim.c ]]; then
    gcc -O2 -fPIC -shared -o "$WAIT_SHIM" "$SELF_DIR/adios2_wait_shim.c" -ldl || die "could not build $WAIT_SHIM"
  fi
  [[ -f $WAIT_SHIM ]] || die "wait timing library not found: $WAIT_SHIM (put adios2_wait_shim.c next to this script, or WAIT=0)"
  [[ -f $WAIT_PY ]]   || echo "NOTE: $WAIT_PY not found - wait CSVs are written, analysis skipped"
fi

mapfile -t WH < <(grep -v '^\s*#' "$WHOSTS" | grep -v '^\s*$' | cut -d: -f1)
mapfile -t RH < <(grep -v '^\s*#' "$RHOSTS" | grep -v '^\s*$' | cut -d: -f1)
mapfile -t ALLH < <(printf '%s\n' "${WH[@]}" "${RH[@]}" | awk '!seen[$0]++')   # each machine once (SPOD may share writer machines)
WNP=$(( ${#WH[@]} * WPPN )); RNP=$(( ${#RH[@]} * RPPN ))
# CPU lists must exist on every machine (taskset fails otherwise and the run would hang at the SST rendezvous)
check_pin() {  # check_pin <cpu list> <hosts...>
  local pin=$1 h rc; shift
  [[ -n $pin ]] || return 0
  for h in "$@"; do
    rc=0; ssh -o BatchMode=yes -o ConnectTimeout=10 "$h" "taskset -c '$pin' true" 2>/dev/null || rc=$?
    (( rc == 0 || rc == 255 || rc == 127 )) || die "CPU list '$pin' is not valid on $h (see: ssh $h lscpu -e)"
  done
}
check_pin "$WPIN" "${WH[@]}"; check_pin "$RPIN" "${RH[@]}"
OUTPUT_DIRS="data probes statistics streaklines out"

getv() { grep -ioE "(^|[ ,])$1\s*=\s*[-0-9]+" "$2" | head -1 | sed -E 's/.*=\s*//'; }
setv() {
  grep -qiE "(^|[ ,])$1\s*=" "$3" || die "$1 not found in $3"
  sed -i -E "s/(^|[ ,])($1\s*=\s*)[-0-9]+/\1\2$2/I" "$3"
  [[ $(getv "$1" "$3") == "$2" ]] || die "could not set $1=$2 in $3"
}

# ---------------------------------------------------------------- 2. run folder
RUN=$RUNROOT/S4_r${#RH[@]}x${RPPN}_$(date +%Y%m%d_%H%M)          # e.g. S4_r4x8_20261108_1400
mkdir -p "$RUN"
tar -C "$CASE" --exclude='./data' --exclude='./probes' --exclude='./statistics' \
    --exclude='./streaklines' --exclude='./out' --exclude='./prof_*' --exclude='*.log' \
    --exclude='*.log1' --exclude='*.bin' --exclude='*.bp4' --exclude='*.bp5' --exclude='*.sst' \
    --exclude='*.bp' --exclude='restart*' --exclude='checkpoint*' --exclude='*.old.*' -cf - . | tar -C "$RUN" -xf -
cd "$RUN"
for d in $OUTPUT_DIRS; do mkdir -p "$d"; done
cp input.i3d input.i3d.orig
cp adios2_config.sst.xml adios2_config.xml           # solution-io -> SST stream "data"

setv ilast "$ILAST" input.i3d
setv ioutput "$IOUTPUT" input.i3d
grep -qiE '(^|[ ,])ivisu\s*=' input.i3d && setv ivisu 1 input.i3d   # snapshots ON (incl. snapshot 0)
setv icheckpoint $((ILAST + 1)) input.i3d                           # no restart writes (as S0)
setv ilist 1 input.i3d                                              # one step line per step (as S0)
NX=$(getv nx input.i3d); NY=$(getv ny input.i3d); NZ=$(getv nz input.i3d)
NSNAP=$(( ILAST / IOUTPUT + 1 ))                                    # + snapshot 0 at start-up
NBLK=$(( SPOD_NDFT > NSNAP ? 0 : (NSNAP - SPOD_NDFT) / (SPOD_NDFT - SPOD_OVERLAP) + 1 ))
(( NBLK >= 1 )) || die "nDFT $SPOD_NDFT > $NSNAP snapshots: no SPOD block would finish"

# SPOD settings: exported (MPICH passes the environment to all ranks) and also given with -genv.
# Both name sets are set, so either SPOD build finds its own (unused names are ignored).
export SPOD_INPUT=adios2                                  # backend: ADIOS2 instead of .bin files
export SPOD_ADIOS_CONFIG="$RUN/adios2_config.xml"         # same config as the solver
export SPOD_CONFIG="$SPOD_ADIOS_CONFIG"
export SPOD_IO=solution-io                                # IO name -> SST engine block in adios2_config.xml
export SPOD_STREAM="$RUN/data"                            # stream name (contact file data.sst)
export SPOD_VAR=${SPOD_VAR:-ux}                           # one field per snapshot
export SPOD_VARS SPOD_NX=$NX SPOD_NY=$NY SPOD_NZ=$NZ
export SPOD_NPOINTS=$(( NX * NY * NZ ))
export SPOD_NSNAP=$NSNAP SPOD_NDFT SPOD_OVERLAP SPOD_K SPOD_NSAVE
export SPOD_WEIGHTS=""                                    # equal weights
SPOD_ENV=()
for v in SPOD_INPUT SPOD_ADIOS_CONFIG SPOD_CONFIG SPOD_IO SPOD_STREAM SPOD_VAR SPOD_VARS SPOD_NX SPOD_NY SPOD_NZ \
         SPOD_NPOINTS SPOD_NSNAP SPOD_NDFT SPOD_OVERLAP SPOD_K SPOD_NSAVE; do
  SPOD_ENV+=(-genv "$v" "${!v}")
done

echo "run folder : $RUN"
echo "writer     : $WNP ranks on ${WH[*]}$([[ $APS == 1 ]] && echo ' (APS)')"
echo "reader     : $RNP ranks on ${RH[*]}$([[ $APS_SPOD == 1 ]] && echo ' (APS)')"
echo "CPU pins   : writer ${WPIN:-none}, reader ${RPIN:-none}"
echo "case       : ${NX}x${NY}x${NZ}, ilast $ILAST, ioutput $IOUTPUT -> $NSNAP snapshots"
echo "SPOD       : $SBIN, var $SPOD_VAR, points $SPOD_NPOINTS, IO $SPOD_IO, nDFT $SPOD_NDFT, overlap $SPOD_OVERLAP -> $NBLK block(s), k $SPOD_K"
echo "changes to input.i3d:"; diff input.i3d.orig input.i3d || true

# ---------------------------------------------------------------- 3. machine information
{
  echo "date    : $(date -Is)"
  echo "writer  : $BIN (md5 $(md5sum "$BIN" | cut -d' ' -f1))"
  echo "reader  : $SBIN (md5 $(md5sum "$SBIN" | cut -d' ' -f1))"
  echo "launch  : W $MPIRUN -f $WHOSTS -np $WNP -ppn $WPPN | R $MPIRUN -f $RHOSTS -np $RNP -ppn $RPPN"
  echo "SPOD    : NSNAP $NSNAP NDFT $SPOD_NDFT OVERLAP $SPOD_OVERLAP K $SPOD_K NSAVE $SPOD_NSAVE VARS [$SPOD_VARS]"
  echo "pins    : writer ${WPIN:-none} reader ${RPIN:-none} (taskset -c)"
  echo "aps     : writer $APS reader $APS_SPOD, APS_IMBALANCE_TYPE $APS_IMB, APS_STAT_LEVEL ${APS_LEVEL:-default}"
  echo "wait    : reader $WAIT writer $WAIT_WRITER $WAIT_SHIM$([[ -f $WAIT_SHIM ]] && echo " (md5 $(md5sum "$WAIT_SHIM" | cut -d' ' -f1))")"
  echo "case fs : $(findmnt -n -o FSTYPE,SOURCE -T "$RUN" 2>/dev/null || echo unknown)"
  echo "--- writer libraries"; ldd "$BIN"  | grep -iE 'mpi|adios|fftw|not found' || true
  echo "--- reader libraries"; ldd "$SBIN" | grep -iE 'mpi|adios|petsc|slepc|fftw|not found' || true
  echo "--- input (S4)"
  grep -iE '(^|[ ,])(nx|ny|nz|p_row|p_col|ilast|ilist|ioutput|ivisu|icheckpoint|iibm|nraf)\s*=' input.i3d || true
  echo "--- solution-io engine"; grep -A2 'solution-io' adios2_config.xml | grep -i engine || true
  echo "--- machines (name, cpu, logical cpus, RAM GB, network card speed Mb/s, load average)"
  for h in "${ALLH[@]}"; do
    ssh -o BatchMode=yes "$h" 'IF=$(ip -o route show default | awk "{print \$5; exit}");
      echo "$(hostname) | $(lscpu | sed -n "s/^Model name:\s*//p") | $(nproc) |" \
           "$(awk "/MemTotal/{printf \"%.1f\", \$2/1048576}" /proc/meminfo) | $IF $(cat /sys/class/net/$IF/speed 2>/dev/null) |" \
           "$(cut -d" " -f1-3 /proc/loadavg)"' 2>/dev/null || echo "$h | ssh failed"
  done
} > meta.txt
cat meta.txt

# ---------------------------------------------------------------- helpers for one run
SAMPLER='
STOP=$1
IF=$(ip -o route show default 2>/dev/null | awk "{print \$5; exit}")
echo "t,cpu_busy_pct,mem_used_gb,net_rx_bytes,net_tx_bytes"
read -r _ a b c d e f g h _ < /proc/stat; pt=$((a+b+c+d+e+f+g+h)); pi=$((d+e))
end=$(( $(date +%s) + 86400 ))
while [ ! -e "$STOP" ] && [ "$(date +%s)" -lt "$end" ]; do
  sleep 1
  read -r _ a b c d e f g h _ < /proc/stat; tt=$((a+b+c+d+e+f+g+h)); ti=$((d+e))
  busy=$(awk -v t=$((tt-pt)) -v i=$((ti-pi)) "BEGIN{ if (t>0) printf \"%.1f\", 100*(t-i)/t; else print 0 }")
  pt=$tt; pi=$ti
  mem=$(awk "/MemTotal/{t=\$2} /MemAvailable/{a=\$2} END{printf \"%.3f\", (t-a)/1048576}" /proc/meminfo)
  net=$(awk -v IF="$IF" "NR>2{sub(/:/,\" \"); if (\$1!=\"lo\" && (IF==\"\" || \$1==IF)) {r+=\$2; x+=\$10}} END{print r+0\",\"x+0}" /proc/net/dev)
  echo "$(date +%s.%N | cut -c1-14),$busy,$mem,$net"
done'

ts() { perl -MTime::HiRes=time -ne '$|=1; printf "%.3f %s", time, $_'; }
now() { date +%s.%N; }
kill_stale() {
  for h in "${ALLH[@]}"; do
    ssh -o BatchMode=yes -o ConnectTimeout=10 "$h" "pkill -u \$USER -f '$(basename "$BIN")|$(basename "$SBIN")'" 2>/dev/null || true
  done
}
clean_out() { for d in $OUTPUT_DIRS; do rm -rf "${d:?}"/*; done; rm -f restart.info checkpoint* data.sst; }

# run_one <name> -> <name>/ with x3d_ts.log, spod_ts.log, node_<host>.csv, aps/, spod/
run_one() {
  local name=$1 dir=$RUN/$1 h pids=() wpid i wrc rrc
  mkdir -p "$dir/spod"; kill_stale; clean_out
  touch "$dir/.marker"
  for h in "${ALLH[@]}"; do
    ssh -o BatchMode=yes "$h" "bash -s -- '$dir/.stop'" <<< "$SAMPLER" > "$dir/node_$h.csv" 2>/dev/null &
    pids+=($!)
  done
  # wait timing: LD_PRELOAD the shim (keeps an existing LD_PRELOAD, e.g. the perfstubs fix), one CSV per rank
  local wenv=() renv=() pre="$WAIT_SHIM${LD_PRELOAD:+:$LD_PRELOAD}" nr nw
  mkdir -p "$dir/wait"
  if [[ $WAIT_WRITER == 1 ]]; then
    wenv=(-genv LD_PRELOAD "$pre" -genv ADIOS_WAIT_DIR "$dir/wait" -genv ADIOS_WAIT_TAG writer)
  fi
  if [[ $WAIT == 1 ]]; then
    renv=(-genv LD_PRELOAD "$pre" -genv ADIOS_WAIT_DIR "$dir/wait" -genv ADIOS_WAIT_TAG reader)
  fi
  # APS settings go to every rank (aps runs inside each rank)
  local apsenv=(-genv APS_IMBALANCE_TYPE "$APS_IMB")
  [[ -n $APS_LEVEL ]] && apsenv+=(-genv APS_STAT_LEVEL "$APS_LEVEL")
  [[ $APS == 1 ]] && wenv+=("${apsenv[@]}")
  [[ $APS_SPOD == 1 ]] && renv+=("${apsenv[@]}")
  local wcmd=("$MPIRUN" "${NET[@]}" ${wenv[@]+"${wenv[@]}"} -f "$WHOSTS" -np "$WNP" -ppn "$WPPN")
  [[ -n $WPIN ]] && wcmd+=(taskset -c "$WPIN")
  [[ $APS == 1 ]] && wcmd+=(aps -r "$dir/aps")
  wcmd+=("$BIN" input.i3d)
  local rcmd=("$MPIRUN" "${NET[@]}" "${SPOD_ENV[@]}" ${renv[@]+"${renv[@]}"} -f "$RHOSTS" -np "$RNP" -ppn "$RPPN")
  [[ -n $RPIN ]] && rcmd+=(taskset -c "$RPIN")
  [[ $APS_SPOD == 1 ]] && rcmd+=(aps -r "$dir/aps_spod")
  rcmd+=("$SBIN")

  # writer in the background (blocks at stream open until SPOD connects: RendezvousReaderCount=1)
  echo "$(now) # START writer" > "$dir/x3d_ts.log"
  set +e
  { GFORTRAN_UNBUFFERED_PRECONNECTED=y "${wcmd[@]}" 2>&1 | ts >> "$dir/x3d_ts.log"
    echo "$(now) # END rc=${PIPESTATUS[0]}" >> "$dir/x3d_ts.log"; } &
  wpid=$!

  # wait for the SST contact file
  for ((i = 0; i < SST_WAIT; i += 2)); do
    [[ -s data.sst ]] && break; kill -0 "$wpid" 2>/dev/null || break; sleep 2
  done
  if [[ ! -s data.sst ]]; then
    echo "$name: writer created no data.sst within ${SST_WAIT}s:"
    grep -n -E 'ERROR|At line|EXCEPTION|Fortran runtime' "$dir/x3d_ts.log" | head -5
    kill_stale; wait "$wpid" 2>/dev/null
    touch "$dir/.stop"; wait "${pids[@]}" 2>/dev/null; set -e; return 1
  fi
  echo "$(now) # SST data.sst ready" >> "$dir/x3d_ts.log"

  # reader in the foreground, in the run folder
  echo "$(now) # START spod" > "$dir/spod_ts.log"
  # unbuffered Fortran output, so every SPOD line gets its real time (block / frequency markers for the burst report)
  GFORTRAN_UNBUFFERED_PRECONNECTED=y "${rcmd[@]}" 2>&1 | ts >> "$dir/spod_ts.log"
  rrc=${PIPESTATUS[0]}
  echo "$(now) # END rc=$rrc" >> "$dir/spod_ts.log"
  [[ $rrc == 0 ]] || { echo "$name: SPOD failed (rc $rrc) - stopping the writer"; kill_stale; }
  wait "$wpid"
  set -e
  wrc=$(sed -n 's/.*# END rc=\([0-9]*\).*/\1/p' "$dir/x3d_ts.log" | tail -1)
  touch "$dir/.stop"; wait "${pids[@]}" 2>/dev/null || true

  # SPOD outputs written into the run folder during this run -> <name>/spod/
  find "$RUN" -mindepth 1 -maxdepth 1 -newer "$dir/.marker" \
       ! -name 'run[0-9]*' ! -name 'data' ! -name 'probes' ! -name 'statistics' \
       ! -name 'streaklines' ! -name 'out' ! -name 'data.sst' ! -name 'input.i3d*' \
       -exec mv -t "$dir/spod/" {} + 2>/dev/null || true
  rm -f data.sst

  # disk space: SPOD outputs (mode files ~1.75 GB each); by default the whole run<N>/spod/ goes
  local keep=0 f m nrm=0 brm=0
  if [[ $KEEP_SPOD != 1 ]]; then
    echo "$name: removed SPOD outputs ($(du -sh "$dir/spod" 2>/dev/null | cut -f1)) - KEEP_SPOD=1 keeps them"
    rm -rf "${dir:?}/spod"
  fi
  [[ $name == run1 ]] && keep=$KEEP_MODES
  for f in "$dir"/spod/SPODMODE_freq_*; do
    [[ -e $f ]] || continue
    m=$(basename "$f" | sed -n 's/^SPODMODE_freq_\([0-9]*\)_.*/\1/p')
    [[ $m =~ ^[0-9]+$ ]] || continue
    if (( m > keep )); then brm=$(( brm + $(stat -c %s "$f") )); rm -f "$f"; nrm=$(( nrm + 1 )); fi
  done
  (( nrm == 0 )) || echo "$name: removed $nrm SPOD mode file(s), $(awk -v b=$brm 'BEGIN{printf "%.1f", b/1e9}') GB (kept the $keep leading mode(s); spectrum kept)"

  # APS reports: aps_report.txt (Xcompact3d), aps_spod_report.txt (SPOD)
  local a aok
  for a in aps aps_spod; do
    [[ -d $dir/$a ]] || continue
    aok=1
    aps --report "$dir/$a" > "$dir/${a}_report.txt" 2>&1 \
      || aps-report "$dir/$a" > "$dir/${a}_report.txt" 2>&1 \
      || { echo "WARNING: aps report failed for $name/$a - try: aps --report $dir/$a"; aok=0; }
    aps-report -t "$dir/$a" > "$dir/${a}_mpi_per_rank.txt" 2>&1 \
      || echo "NOTE: aps-report -t (MPI time per rank) failed for $name/$a"
    aps-report -f "$dir/$a" > "$dir/${a}_mpi_functions.txt" 2>&1 \
      || echo "NOTE: aps-report -f (time and data per MPI function) failed for $name/$a"
    if [[ $KEEP_APS_RAW != 1 && $aok == 1 ]]; then         # raw data only kept if the report failed
      echo "$name: removed raw APS data $a ($(du -sh "$dir/$a" 2>/dev/null | cut -f1); reports kept)"
      rm -rf "${dir:?}/$a"
    fi
  done
  nr=$(find "$dir/wait" -name 'reader_r*.csv' | wc -l); nw=$(find "$dir/wait" -name 'writer_r*.csv' | wc -l)
  if [[ $WAIT == 1 ]] && (( nr != RNP )); then
    echo "WARNING: $name: $nr of $RNP reader wait files - if 0 with APS_SPOD=1, APS may drop LD_PRELOAD: retry with APS_SPOD=0"
  fi
  if [[ $WAIT_WRITER == 1 ]] && (( nw != WNP )); then
    echo "WARNING: $name: $nw of $WNP writer wait files - if 0 with APS=1, retry with APS=0"
  fi
  echo "$name: writer rc ${wrc:-?}, SPOD rc $rrc, $(awk '/# START/{s=$1} /# END/{e=$1} END{printf "%.1f", e-s}' "$dir/x3d_ts.log") s"
  [[ ${wrc:-1} == 0 && $rrc == 0 ]]
}

# ---------------------------------------------------------------- 4. runs
RUNS=()
for ((r = 1; r <= REPEATS; r++)); do
  echo; echo "=== run $r of $REPEATS (S4: SST -> SPOD) ==="
  run_one "run$r" || echo "WARNING: run$r failed - see $RUN/run$r/x3d_ts.log and spod_ts.log"
  RUNS+=("run$r")
done

# ---------------------------------------------------------------- 5. summary
echo; echo "=== summary ==="
python3 - "$RUN" "$NSNAP" "${S0_STEP:-none}" "${#RH[@]}x${RPPN}" "${RUNS[@]}" <<'PY'
import csv, glob, os, re, statistics as st, sys
run_dir, nsnap_exp, s0, rlay, names = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5:]
s0 = None if s0 == "none" else float(s0)

def marks(log):
    out = {}
    for line in open(log, errors="replace"):
        p = line.split(" ", 1)
        if len(p) > 1 and p[1].startswith("# "):
            try: out.setdefault(p[1].split()[1], float(p[0]))
            except ValueError: pass
            m = re.search(r"rc=(\d+)", p[1])
            if m: out["rc"] = int(m.group(1))
    return out

rows, spod_lines = [], {}
for n in names:
    d = os.path.join(run_dir, n); log = os.path.join(d, "x3d_ts.log")
    if not os.path.exists(log): continue
    mk = marks(log); steps, snaps = [], []
    for line in open(log, errors="replace"):
        p = line.split(" ", 1)
        try: t = float(p[0])
        except ValueError: continue
        txt = p[1] if len(p) > 1 else ""
        m = re.search(r"Time step\s*=\s*(\d+)\s*/", txt)
        if m: steps.append((int(m.group(1)), t))
        if "Writing snapshots" in txt: snaps.append(t)
    plain, withsnap = [], []
    for (n1, t1), (n2, t2) in zip(steps, steps[1:]):
        if n2 > n1:
            (withsnap if any(t1 <= s <= t2 for s in snaps) else plain).append((t2 - t1) / (n2 - n1))
    dur = plain + withsnap
    rlog = os.path.join(d, "spod_ts.log"); rk = marks(rlog) if os.path.exists(rlog) else {}
    if os.path.exists(rlog):
        spod_lines[n] = [l.split(" ", 1)[1].rstrip() for l in open(rlog, errors="replace")
                         if re.search(r"wait|read|compute|block|total|time", l, re.I) and " # " not in l][-12:]
    mem, cpu = {}, []
    for f in glob.glob(os.path.join(d, "node_*.csv")):
        h = os.path.basename(f)[5:-4]
        rr = [r for r in csv.DictReader(open(f)) if r.get("mem_used_gb")]
        if rr:
            mem[h] = max(float(r["mem_used_gb"]) for r in rr)
    start, end = mk.get("START"), mk.get("END")
    first_t = steps[0][1] if steps else None; last_t = steps[-1][1] if steps else None
    med = st.median(plain) if plain else (st.median(dur) if dur else None)
    rows.append({
        "run": n, "readers": rlay, "writer_rc": mk.get("rc", ""), "spod_rc": rk.get("rc", ""),
        "wall_s": round(end - start, 2) if start and end else "",
        "sst_ready_s": round(mk["SST"] - start, 2) if "SST" in mk and start else "",
        "init_s": round(first_t - start, 2) if first_t and start else "",
        "loop_s": round(last_t - first_t, 2) if steps else "",
        "end_s": round(end - last_t, 2) if steps and end else "",
        "spod_wall_s": round(rk["END"] - rk["START"], 2) if "START" in rk and "END" in rk else "",
        "spod_after_writer_s": round(rk["END"] - end, 2) if "END" in rk and end else "",
        "steps": steps[-1][0] if steps else 0, "snapshots": len(snaps),
        "step_median_s": round(med, 3) if med is not None else "",
        "snap_step_extra_s": round(st.mean(withsnap) - med, 3) if withsnap and med is not None else "",
        "step_max_s": round(max(dur), 3) if dur else "",
        "mem_peak_writer_gb": "", "mem_peak_reader_gb": "",
        "mem_peak_host": max(mem, key=mem.get) if mem else "",
    })
    meta = open(os.path.join(run_dir, "meta.txt")).read()
    w = re.search(r"-f (\S+) -np", meta); rr_ = re.search(r"R \S+ -f (\S+)", meta)
    wh = [l.split(":")[0].strip() for l in open(w.group(1)) if l.strip() and not l.startswith("#")] if w else []
    rh = [l.split(":")[0].strip() for l in open(rr_.group(1)) if l.strip() and not l.startswith("#")] if rr_ else []
    rows[-1]["mem_peak_writer_gb"] = round(max([mem[h] for h in wh if h in mem] or [0]), 2)
    rows[-1]["mem_peak_reader_gb"] = round(max([mem[h] for h in rh if h in mem] or [0]), 2)

keys = list(rows[0].keys()) if rows else []
with open(os.path.join(run_dir, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(rows)
lines = ["  ".join(f"{k:>19s}" for k in keys)] + ["  ".join(f"{str(r[k]):>19s}" for k in keys) for r in rows]
ok = [r for r in rows if r["step_median_s"] != "" and r["wall_s"] != ""]
if ok:
    ws, sm = st.median([r["wall_s"] for r in ok]), st.median([r["step_median_s"] for r in ok])
    ex = [r["snap_step_extra_s"] for r in ok if r["snap_step_extra_s"] != ""]
    lines += ["", f"S4 result, reader {rlay} (median of {len(ok)} runs): wall {ws:.1f} s, plain step {sm:.4f} s, "
              f"extra per snapshot step {st.median(ex):.3f} s" if ex else
              f"S4 result, reader {rlay} (median of {len(ok)} runs): wall {ws:.1f} s, step {sm:.4f} s"]
    if s0:
        lines += [f"vs S0 step {s0:.4f} s: plain-step overhead {100*(sm-s0)/s0:+.1f} %"]
for r in rows:
    if r["snapshots"] != nsnap_exp:
        lines += [f"NOTE {r['run']}: {r['snapshots']} snapshot writes seen, expected {nsnap_exp}"]
for n, ls in spod_lines.items():
    if ls: lines += ["", f"--- SPOD timing lines ({n}, last {len(ls)})"] + ["  " + l for l in ls]
txt = "\n".join(lines); print(txt)
open(os.path.join(run_dir, "summary.txt"), "w").write(txt + "\n")
PY
if [[ ($WAIT == 1 || $WAIT_WRITER == 1) && -f $WAIT_PY ]]; then
  echo; echo "=== reader / writer wait (per snapshot) ==="
  python3 "$WAIT_PY" "$RUN" || echo "WARNING: wait analysis failed - try: python3 $WAIT_PY $RUN"
fi
if [[ $KEEP_X3D_OUT != 1 ]]; then
  clean_out; echo; echo "removed Xcompact3d output folders ($OUTPUT_DIRS) - KEEP_X3D_OUT=1 keeps them"
fi
echo
echo "Wait times  : $RUN/wait_summary.txt, wait_steps.csv (per snapshot), run*/wait/*.csv (per rank)"
echo "APS reports : $RUN/run*/aps_report.txt (Xcompact3d), $RUN/run*/aps_spod_report.txt (SPOD), *_mpi_per_rank.txt, *_mpi_functions.txt"
[[ $KEEP_SPOD == 1 ]] && echo "SPOD output : $RUN/run*/spod/"
echo "All results : $RUN ($(du -sh "$RUN" 2>/dev/null | cut -f1))"
