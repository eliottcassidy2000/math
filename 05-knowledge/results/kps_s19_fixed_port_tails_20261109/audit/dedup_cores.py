"""S19 audit: distinct core states (h0, R0, Delta0, hu) per base and maximal core order, from cores.exe classes output.
Also cross-checks against the labelled (no-isomorphism) enumeration files if given.
Usage: python dedup_cores.py CLASS_STATS OUTDIR [LABELLED_FILE_n ...]
Writes OUTDIR/cores_b{2,3}_le{n}.txt (sorted) and prints the counts.
"""
import sys
src, outdir = sys.argv[1], sys.argv[2]
rows = []
for line in open(src):
    if line.startswith('#'): continue
    f = list(map(int, line.split()))
    rows.append(f)
summary = []
for x in (2, 3):
    for nmax in range(1, 8):
        s = set()
        for f in rows:
            if f[0] > nmax: continue
            st = tuple(f[3:7]) if x == 2 else tuple(f[7:11])
            s.add(st)
        with open(f"{outdir}/cores_b{x}_le{nmax}.txt", "w") as fh:
            for t in sorted(s): fh.write(" ".join(map(str, t)) + "\n")
        summary.append((x, nmax, len(s)))
for x, nmax, c in summary:
    print(f"base {x}, cores of order <= {nmax}: {c} distinct states (h0, R0, Delta0, hu)")
# labelled cross-check: each file holds the distinct states of all labelled tournaments of one order n with port 0
for path in sys.argv[3:]:
    n = int(path.rsplit('_', 1)[1].split('.')[0])
    lab = {2: set(), 3: set()}
    for line in open(path):
        x, *t = map(int, line.split()); lab[x].add(tuple(t))
    for x in (2, 3):
        cls = set((tuple(f[3:7]) if x == 2 else tuple(f[7:11])) for f in rows if f[0] == n)
        print(f"order {n}, base {x}: labelled distinct {len(lab[x])}, class-based distinct {len(cls)}, equal: {lab[x] == cls}")
