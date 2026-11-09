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
