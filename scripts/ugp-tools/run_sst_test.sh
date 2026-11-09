#!/usr/bin/env bash
# =====================================================================
# run_sst_test.sh - direct SST test: Xcompact3d streams snapshots over ADIOS2 SST,
# a reader on other machines receives them (no SPOD yet).
#
#   writer : ~/ugp/src/x3d_lab_adios/xcompact3d on csews160..167 (hostfile, 192 ranks)
#   reader : adios2_reorganize_mpi on csews168..171 (hostfile_spod, 1 rank per machine)
#            KEEP=1 (default): saves what it receives to sst_copy.bp  -> check with bpls
#            KEEP=0          : receives and discards (Null engine), no disk use
#   steps  : ILAST=10 (snapshots at 0, 5, 10 with ioutput 5); input.i3d is restored afterwards
#
# Every mpirun uses eno1 (all csews machines also have docker0).
# Usage:  bash run_sst_test.sh            Options: ILAST=10 KEEP=1 RNP=4 CASE=~/ugp/runs/testcase_adios
# =====================================================================
set -uo pipefail

CASE=${CASE:-$HOME/ugp/runs/testcase_adios}
BIN=${BIN:-$HOME/ugp/src/x3d_lab_adios/xcompact3d}
WHOSTS=${WHOSTS:-$HOME/ugp/tools/hostfile}
RHOSTS=${RHOSTS:-$HOME/ugp/tools/hostfile_spod}
WPPN=${WPPN:-24}
RNP=${RNP:-4}
ILAST=${ILAST:-10}
KEEP=${KEEP:-1}
MPIRUN=/usr/local/bin/mpirun
A2=$HOME/ugp/install/bin
F=(-iface eno1 -genv FI_PROVIDER tcp -genv FI_TCP_IFACE eno1)

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -x $BIN ]] || die "no binary $BIN (run setup_adios_case.sh first)"
[[ -f $CASE/adios2_config.sst.xml ]] || die "no adios2_config.sst.xml in $CASE (run setup_adios_case.sh first)"
grep -q '"SST"' "$(dirname "$BIN")/decomp2d/io.f90" || die "x3d_lab_adios has no SST support in io.f90 (rerun setup_adios_case.sh)"
grep -q NetworkInterface "$CASE/adios2_config.sst.xml" || die "adios2_config.sst.xml has no NetworkInterface=eno1 (rerun setup_adios_case.sh SKIP_BUILD=1)"
mapfile -t WH < <(grep -v '^\s*$' "$WHOSTS" | cut -d: -f1)
mapfile -t RH < <(grep -v '^\s*$' "$RHOSTS" | cut -d: -f1)
WNP=$(( ${#WH[@]} * WPPN ))

cd "$CASE"
# ---------------------------------------------------------------- 1. clean start
echo "== 1. clean up old processes and files"
for h in "${WH[@]}" "${RH[@]}"; do
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$h" 'pkill -u $USER -f "xcompact3d|adios2_reorganize"' 2>/dev/null
done
rm -rf data.sst sst_copy.bp x3d_sst_test.log reader_test.log
cp adios2_config.sst.xml adios2_config.xml
cp input.i3d input.i3d.before_sst_test
sed -i -E "s/^(\s*ilast\s*=\s*)[0-9]+/\1$ILAST/" input.i3d
trap 'cp input.i3d.before_sst_test input.i3d' EXIT          # always restore the input
echo "   ilast=$ILAST, ioutput=$(grep -oiE 'ioutput\s*=\s*[0-9]+' input.i3d | grep -oE '[0-9]+$')"

# ---------------------------------------------------------------- 2. writer
echo "== 2. start Xcompact3d (writer): $WNP ranks on ${WH[*]}"
T0=$(date +%s)
"$MPIRUN" "${F[@]}" -f "$WHOSTS" -np "$WNP" -ppn "$WPPN" "$BIN" input.i3d > x3d_sst_test.log 2>&1 &
WPID=$!

echo "== 3. wait for data.sst (writer opens the stream after start-up and geometry)"
for i in $(seq 120); do
  [[ -s data.sst ]] && break
  kill -0 $WPID 2>/dev/null || break
  sleep 5
done
if [[ ! -s data.sst ]]; then
  echo "writer did not create data.sst - it stopped or is stuck:"
  grep -n -E 'ERROR|At line|EXCEPTION|Fortran runtime|Launch|ssh:' x3d_sst_test.log | head
  kill $WPID 2>/dev/null; exit 1
fi
echo "   data.sst after $(( $(date +%s) - T0 )) s:"; sed 's/^/   | /' data.sst | head -5

# ---------------------------------------------------------------- 4. reader
if [[ $KEEP == 1 ]]; then OUT=sst_copy.bp; ENG=BP4; else OUT=/dev/null; ENG=Null; fi
echo "== 4. start reader: $RNP ranks on ${RH[*]} -> $OUT ($ENG)"
"$MPIRUN" "${F[@]}" -f "$RHOSTS" -np "$RNP" -ppn 1 \
  "$A2/adios2_reorganize_mpi" data "$OUT" SST "OpenTimeoutSecs=600,NetworkInterface=eno1" "$ENG" "" \
  > reader_test.log 2>&1
RRC=$?
wait $WPID; WRC=$?
echo "   writer exit $WRC, reader exit $RRC, total $(( $(date +%s) - T0 )) s"

# ---------------------------------------------------------------- 5. results
echo "== 5. writer: snapshots"
grep -E 'Writing snapshots|Time for writing' x3d_sst_test.log
grep -E 'Time step\s*=\s*'"$ILAST"'/' x3d_sst_test.log | tail -1
grep -n -E 'ERROR|EXCEPTION|At line' x3d_sst_test.log | head -5
echo "== 5. reader"
grep -v -E 'Authorization|^\s*$' reader_test.log | tail -8
if [[ $KEEP == 1 && -d sst_copy.bp ]]; then
  echo "== 5. what arrived (bpls -l sst_copy.bp): expect each variable with $(( ILAST / 5 + 1 ))*{301, 301, 301} and different min/max"
  "$A2/bpls" -l sst_copy.bp
  du -sh sst_copy.bp
  echo "delete when done: rm -rf $CASE/sst_copy.bp"
fi
[[ $WRC == 0 && $RRC == 0 ]] && echo "SST TEST: both sides finished without error" || echo "SST TEST: FAILED (see x3d_sst_test.log / reader_test.log)"
