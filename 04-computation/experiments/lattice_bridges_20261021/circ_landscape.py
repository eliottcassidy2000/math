# The H-landscape on circulant tournaments: vertices = multiplier classes of connection sets S (one of k, n-k for each k),
# edges = reversing one difference class {k, n-k}.  Lists all local maxima (strict) with their flip distance to R_n and,
# for prime n = 3 mod 4, to the Paley set QR_n.  Data from circ_counts.c (complete class tables).  mac-mini-2026-10-21-fanopsp.
import re, glob
from math import gcd
from collections import defaultdict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m_ = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m_: rows[int(m_.group(1))][frozenset(int(x) for x in m_.group(2).split(','))] = int(m_.group(4))
def canon(n, S):
    best = None
    for a in range(1, n):
        if gcd(a, n) != 1: continue
        t = tuple(sorted((a * s) % n for s in S))
        if best is None or t < best: best = t
    return best
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))
for n in sorted(rows):
    m = (n - 1) // 2
    H = {canon(n, S): h for S, h in rows[n].items()}
    R = canon(n, range(1, m + 1))
    P = canon(n, [x * x % n for x in range(1, n)]) if isprime(n) and n % 4 == 3 else None
    def dist(A, B):
        A = set(A); best = m
        for a in range(1, n):
            if gcd(a, n) != 1: continue
            aB = {(a * s) % n for s in B}
            best = min(best, len(A - aB))
        return best
    locmax = []
    for S, h in H.items():
        nbrs = []
        for k in range(1, m + 1):
            S2 = set(S)
            if k in S2: S2.discard(k); S2.add(n - k)
            else: S2.discard(n - k); S2.add(k)
            nbrs.append(H[canon(n, S2)])
        if all(h > x for x in nbrs): locmax.append((h, S))
    locmax.sort(reverse=True)
    desc = []
    for h, S in locmax:
        tag = "R" if S == R else ("Paley" if S == P else "")
        desc.append(f"H={h} dR={dist(S, R)}" + (f" dP={dist(S, P)}" if P else "") + (f" [{tag}]" if tag else f" S={list(S)}"))
    print(f"n={n:2d}: {len(H)} classes, {len(locmax)} strict local maxima: " + " | ".join(desc))
