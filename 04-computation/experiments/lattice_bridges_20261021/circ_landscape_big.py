# Complete circulant H-landscape at one odd order n (default 25): every connection set S (one of k, n-k for each k)
# up to multipliers, H by inclusion-exclusion (../paley_extremal_20261021/ie_ham.c, exact in Z/2^128), then the global
# maximum, all strict local maxima of the flip graph, and their flip distance to R_n.  Pre-registered test of the
# one-flip conjecture HYP-9362(a).  mac-mini-2026-10-21-fanopsp.
# usage: python3 circ_landscape_big.py N IE_HAM_BINARY OUTFILE [WORKERS]
import sys, subprocess, re
from math import gcd
from itertools import product
from concurrent.futures import ThreadPoolExecutor
n = int(sys.argv[1]); BIN = sys.argv[2]; OUT = sys.argv[3]; W = int(sys.argv[4]) if len(sys.argv) > 4 else 2
m = (n - 1) // 2
units = [a for a in range(1, n) if gcd(a, n) == 1]
def canon(S): return min(tuple(sorted((a * s) % n for s in S)) for a in units)
classes = {}
for choice in product((0, 1), repeat=m):
    S = [k if c == 0 else n - k for k, c in zip(range(1, m + 1), choice)]
    c = canon(S); classes[c] = classes.get(c, 0) + 1
reps = sorted(classes)
print(f"n={n}: {2**m} connection sets, {len(reps)} multiplier classes", flush=True)
done = {}
try:
    for line in open(OUT):
        mm = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
        if mm: done[tuple(int(x) for x in mm.group(1).split(','))] = int(mm.group(2))
except FileNotFoundError: pass
todo = [S for S in reps if S not in done]
print(f"{len(done)} already computed, {len(todo)} to go", flush=True)
def run(S):
    r = subprocess.run([BIN, 'circ', str(n)] + [str(s) for s in S], capture_output=True, text=True, check=True)
    return S, int(re.search(r'H=(\d+)', r.stdout).group(1))
with open(OUT, 'a') as f, ThreadPoolExecutor(W) as ex:
    for S, h in ex.map(run, todo):
        done[S] = h; f.write(f"S={{{','.join(map(str, S))}}} H={h}\n"); f.flush()
H = {S: done[S] for S in reps}
R = canon(range(1, m + 1))
def flip(S, k):
    S2 = set(S)
    if k in S2: S2.discard(k); S2.add(n - k)
    else: S2.discard(n - k); S2.add(k)
    return canon(S2)
def dist(A, B):
    A = set(A); return min(len(A - {(a * s) % n for s in B}) for a in units)
best = max(H, key=H.get)
print(f"global max H={H[best]} S={list(best)} flip distance to R_{n} = {dist(best, R)}; H(R_{n}) = {H[R]}")
loc = [(H[S], S) for S in reps if all(H[S] > H[flip(S, k)] for k in range(1, m + 1))]
loc.sort(reverse=True)
print(f"{len(loc)} strict local maxima:")
for h, S in loc: print(f"  H={h} dR={dist(S, R)} S={list(S)}" + ("  [R]" if S == R else ""))
nb = sorted(((H[flip(R, k)], k) for k in range(1, m + 1)), reverse=True)
print(f"R_{n} is a strict local max: {all(H[R] > h for h, _ in nb)}; best neighbour ratio = {nb[0][0] / H[R]:.6f} (flip k={nb[0][1]})")
