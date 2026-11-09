#!/usr/bin/env python3
"""Let the SPOD settings in Subroutine UserInput be set by environment variables.

    python3 patch_env_overrides.py FILE.F90 [FILE.F90 ...]

The same patch goes into the old and the node-local reader, so a test gives both
the same settings without editing code. A variable that is unset or empty keeps
the value written in UserInput. Running it twice changes nothing.

  SPOD_NX SPOD_NY SPOD_NZ          grid of each field          (gridNx, gridNy, gridNz)
  SPOD_STREAM                      SST stream name or path     (streamName)
  SPOD_CONFIG                      ADIOS2 XML for 'spod-input' (adiosConfigFile)
  SPOD_NSNAP SPOD_NDFT SPOD_K SPOD_OVERLAP SPOD_NSAVE
  SPOD_WEIGHTS                     weight file, 'none' = equal weights
  X3D_PROW X3D_PCOL X3D_RANKS_PER_NODE   node-local reader only
"""
import re
import sys

MARK = "spod-dev: settings from environment variables"
NUM = [("SPOD_NX", "gridNx"), ("SPOD_NY", "gridNy"), ("SPOD_NZ", "gridNz"),
       ("SPOD_NSNAP", "numberOfSnapShot"), ("SPOD_NDFT", "nDFT"), ("SPOD_K", "k"),
       ("SPOD_OVERLAP", "percentageOfOverall"), ("SPOD_NSAVE", "nFreqToSave")]
NL_ONLY = [("X3D_PROW", "x3dProw"), ("X3D_PCOL", "x3dPcol"), ("X3D_RANKS_PER_NODE", "x3dRanksPerNode")]
TXT = [("SPOD_STREAM", "streamName"), ("SPOD_CONFIG", "adiosConfigFile"),
       ("SPOD_WEIGHTS", "weightMatrixFilePath")]


def patch(path):
    with open(path, newline="") as fh:
        src = fh.read().replace("\r\n", "\n")
    if MARK in src:
        print("already patched:", path)
        return
    head = re.compile(r"^[ \t]*Subroutine[ \t]+UserInput[ \t]*\n[ \t]*implicit[ \t]+none[ \t]*\n", re.I | re.M)
    tail = re.compile(r"^[ \t]*End[ \t]+Subroutine[ \t]+UserInput\b", re.I | re.M)
    if len(head.findall(src)) != 1 or len(tail.findall(src)) != 1:
        sys.exit(f"STOP: Subroutine UserInput not recognised in {path}")
    if re.search(r"\benv(Txt|St)\b", src):
        sys.exit(f"STOP: {path} already uses the names envTxt/envSt")
    nodelocal = re.search(r"::\s*x3dProw\b", src) is not None

    pad = " " * 32
    decl = f"{pad}character(len=500) :: envTxt\n{pad}integer :: envSt\n"
    lines = ["", f"{pad}! ---- {MARK} (unset or empty = value above) ----"]
    for var, target in NUM + (NL_ONLY if nodelocal else []):
        lines.append(f"{pad}call get_environment_variable('{var}', envTxt, status=envSt)")
        lines.append(f"{pad}if (envSt == 0 .and. len_trim(envTxt) > 0) read(envTxt, *) {target}")
    for var, target in TXT:
        lines.append(f"{pad}call get_environment_variable('{var}', envTxt, status=envSt)")
        lines.append(f"{pad}if (envSt == 0 .and. len_trim(envTxt) > 0) {target} = trim(envTxt)")
    lines.append(f"{pad}if (trim(weightMatrixFilePath) == 'none') weightMatrixFilePath = ''")
    lines.append(f"{pad}numberOfDataInSnapShot = int(nStreamVars * gridNx * gridNy * gridNz)")
    block = "\n".join(lines) + "\n\n"

    src = head.sub(lambda m: m.group(0) + decl, src, count=1)
    m = tail.search(src)
    src = src[:m.start()] + block + src[m.start():]
    with open(path, "w", newline="\n") as fh:
        fh.write(src)
    print(f"patched: {path} ({'node-local' if nodelocal else 'old'} reader)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        patch(p)
