# Over complete circulant class tables (odd n <= 19): which class maximises each OCF packing number alpha_k?
# Test of the packing reading of the two peaks (results note section 13).  mac-mini-2026-10-21-fanopsp.
# usage: python3 alpha_all.py ALPHA_DECOMP_BINARY [NMIN NMAX WORKERS]  (defaults 9 19 1)
import re, glob, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from math import gcd
from collections import defaultdict
BIN = sys.argv[1]; NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 9; NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 19; WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        mm = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if mm: rows[int(mm.group(1))][tuple(sorted(int(x) for x in mm.group(2).split(',')))] = int(mm.group(4))
def canon(n, S): return min(tuple(sorted((a * s) % n for s in S)) for a in range(1, n) if gcd(a, n) == 1)
def E(n, S):
    Ss = set(S); return sum(1 for a in S for b in S for c in S if (a + b - c) % n in Ss)
for n in sorted(rows):
    if n < NMIN or n > NMAX: continue
    m = (n - 1) // 2; R = canon(n, range(1, m + 1))
    reps = {canon(n, S): (S, h) for S, h in rows[n].items()}
    def one(item):
        c, (S, h) = item
        out = subprocess.run([BIN, str(n)] + [str(s) for s in S], capture_output=True, text=True, check=True).stdout
        a = {int(k): int(v) for k, v in re.findall(r'alpha_(\d+)=(\d+)', out)}
        Hs = int(re.search(r'H=sum=(\d+)', out).group(1)); assert Hs == h, (n, S, Hs, h)
        return c, a
    with ThreadPoolExecutor(WORKERS) as ex: al = dict(ex.map(one, reps.items()))
    emin = min(reps, key=lambda c: E(n, c))
    K = max(max(a) for a in al.values())
    desc = []
    for k in range(1, K + 1):
        best = max(al[c].get(k, 0) for c in al); arg = [c for c in al if al[c].get(k, 0) == best]
        tag = "R" if arg == [R] else ("Emin" if arg == [emin] else f"{len(arg)} class(es) {[list(a) for a in arg][:1]}")
        desc.append(f"a{k}:{tag}")
    print(f"n={n:2d} ({len(al)} classes; H checked): argmax " + "  ".join(desc), flush=True)
