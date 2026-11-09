#!/bin/bash
# =====================================================================
# setup.sh - one-time set-up of ~/spod-dev.
#
# Copies from ~/ugp (read only, nothing there is changed):
#   src/spod_old.F90   the ADIOS2 reader S4 runs (~/ugp/src/SPOD-adios2, working copy)
#   bin/x3d_lab        the Xcompact3d ADIOS2 build (~/ugp/src/x3d_lab_adios/xcompact3d)
#   cases/base/       the S4 case inputs (~/ugp/runs/testcase_adios, no outputs)
# src/spod_nl.F90 (node-local reader) is uploaded separately.
# The binary gets a new name, so the pkill in run_S4.sh can never match it.
# Then: bash build.sh ; bash run_compare.sh
# =====================================================================
set -euo pipefail
D=$HOME/spod-dev
U=$HOME/ugp
cd "$D"
mkdir -p src build bin cases runs tools

OLD_SRC="$U/src/SPOD-adios2/Streaming SPOD/streamingSPODModule.F90"
X3D_BIN=$U/src/x3d_lab_adios/xcompact3d
CASE=$U/runs/testcase_adios
[[ -f $OLD_SRC ]] || { echo "STOP: $OLD_SRC missing"; exit 1; }
[[ -x $X3D_BIN ]] || { echo "STOP: $X3D_BIN missing"; exit 1; }
[[ -f $CASE/input.i3d && -f $CASE/adios2_config.sst.xml ]] || { echo "STOP: $CASE incomplete"; exit 1; }
[[ -f src/spod_nl.F90 ]] || { echo "STOP: upload src/spod_nl.F90 first"; exit 1; }

sed 's/\r$//' "$OLD_SRC" > src/spod_old.F90
cp "$X3D_BIN" bin/x3d_lab
rm -rf cases/base && mkdir -p cases/base
tar -C "$CASE" --exclude='./data' --exclude='./probes' --exclude='./statistics' \
    --exclude='./streaklines' --exclude='./out' --exclude='./prof_*' --exclude='*.log' \
    --exclude='*.log1' --exclude='*.bin' --exclude='*.bp4' --exclude='*.bp5' --exclude='*.sst' \
    --exclude='*.bp' --exclude='restart*' --exclude='checkpoint*' --exclude='*.old.*' -cf - . \
  | tar -C cases/base -xf -
cp cases/base/adios2_config.sst.xml cases/base/adios2_config.xml

{
  echo "set up : $(date -Is) on $(hostname)"
  echo "old src: $OLD_SRC (md5 $(md5sum "$OLD_SRC" | cut -d' ' -f1))"
  echo "x3d bin: $X3D_BIN (md5 $(md5sum "$X3D_BIN" | cut -d' ' -f1))"
  echo "case   : $CASE"
  echo "nl src : src/spod_nl.F90 (md5 $(md5sum src/spod_nl.F90 | cut -d' ' -f1))"
} | tee origin.txt
echo "set-up done: $D"
