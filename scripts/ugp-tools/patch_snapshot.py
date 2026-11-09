#!/usr/bin/env python3
"""Guard the start-up snapshot (snapshot 0) in Xcompact3d's preprocessing with ivisu.

Usage:  python3 patch_snapshot0.py src/case.f90

In subroutine preprocessing, each call to write_snapshot / end_snapshot (with its
'&' continuation lines) is wrapped in
    if (ivisu.ne.0) then
       ...
    endif
With the default ivisu = 1 nothing changes; ivisu = 0 in &BasicParam gives zero snapshots.
A backup <file>.orig is written. Running it twice does nothing the second time.
"""
import re
import shutil
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "src/case.f90"
lines = open(path).read().split("\n")

# locate subroutine preprocessing
start = next((i for i, l in enumerate(lines)
              if re.match(r"\s*subroutine\s+preprocessing\b", l, re.I)), None)
if start is None:
    sys.exit("subroutine preprocessing not found")
end = next((i for i in range(start + 1, len(lines))
            if re.match(r"\s*end\s+subroutine\s+preprocessing\b", lines[i], re.I)), None)
if end is None:
    sys.exit("end of subroutine preprocessing not found")


def statement_end(i):
    """last line of the statement starting at line i (follows '&' continuations)."""
    while lines[i].split("!")[0].rstrip().endswith("&"):
        i += 1
    return i


calls = [i for i in range(start, end)
         if re.match(r"\s*call\s+(write_snapshot|end_snapshot)\b", lines[i], re.I)]
if not calls:
    sys.exit("no call write_snapshot / end_snapshot inside preprocessing - nothing to patch")

if any(re.search(r"ivisu", lines[c - 1], re.I) for c in calls):
    print("already guarded by ivisu - no change")
    sys.exit(0)

# wrap each call (with its continuation lines) on its own, so other lines in between stay untouched
shutil.copy(path, path + ".orig")
for c in reversed(calls):
    e = statement_end(c)
    indent = re.match(r"\s*", lines[c]).group(0)
    block = ["   " + l if l.strip() else l for l in lines[c:e + 1]]
    lines[c:e + 1] = [f"{indent}if (ivisu.ne.0) then   ! zero-snapshot runs: ivisu = 0"] + block + [f"{indent}endif"]
open(path, "w").write("\n".join(lines))
print(f"patched {path}: {len(calls)} call(s) in preprocessing now inside 'if (ivisu.ne.0)'; backup {path}.orig")
