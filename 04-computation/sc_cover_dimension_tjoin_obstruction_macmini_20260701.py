#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S85 -- does the T-join boundary/parity obstruct low-dimensional covers of the SC classes?

Two threads joined:
 - S84: the BLUE line subgraph on the SC (self-complementary) merged nodes is ALL-ODD-degree = a T-JOIN with
   T = every SC node (boundary = sum of all SC nodes in GF(2)).
 - S81: kappa(n) = min dimension of an axis-aligned subcube (choose k free arcs, fix the rest) realizing a
   target set of iso classes.

DEFINE kappa_SC(n) = min subcube dimension realizing ALL SC iso classes. Compare to:
   info floor ceil(log2 #SC),  #SC=A051337,  half-tiling dim floor((n-1)^2/4) (the grid-sym subcube, a natural
   complement-symmetric SC cover),  kappa_all=1+C(n-2,2).
Also analyze the blue T-join: degrees, is it bipartite, does it have a perfect matching (min T-join = #SC/2),
and the parity/boundary. OBSTRUCTION = kappa_SC exceeding the floor, and whether the gap is parity-forced.
"""
import itertools
from math import comb, log2, ceil
from collections import defaultdict

def build(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}; P = len(pairs)
    perms = list(itertools.permutations(range(n)))
    def relabel(mask, perm):
        out = 0
        for (i, j), k in idx.items():
            bit = (mask >> k) & 1; pi, pj = perm[i], perm[j]
            a, b = (pi, pj) if pi < pj else (pj, pi)
            out |= (bit if pi < pj else 1 - bit) << idx[(a, b)]
        return out
    # transpose (reverse all arcs) = complement: for tournament mask, complement flips every arc bit
    comp = lambda mask: mask ^ ((1 << P) - 1)
    canon_of = [0]*(1 << P); cid = {}; nid = 0
    for mask in range(1 << P):
        c = min(relabel(mask, p) for p in perms)
        if c not in cid: cid[c] = nid; nid += 1
        canon_of[mask] = cid[c]
    # SC classes: class c is self-complementary iff canon(complement of a rep) == c
    rep = {}
    for mask in range(1 << P):
        rep.setdefault(canon_of[mask], mask)
    sc_class = set()
    for c, mask in rep.items():
        if canon_of[comp(mask)] == c: sc_class.add(c)
    return pairs, idx, P, canon_of, nid, sc_class

def kappa_cover(P, canon_of, targets, kmax=None):
    """min subcube dim (free arcs) whose 2^k tournaments realize every class in `targets`."""
    lo = ceil(log2(len(targets)))
    for k in range(lo, (kmax or P) + 1):
        for F in itertools.combinations(range(P), k):
            Fb = list(F); rest = [b for b in range(P) if b not in F]
            subs = [sum((1 << Fb[t]) for t in range(k) if (s >> t) & 1) for s in range(1 << k)]
            for fx in range(1 << len(rest)):
                fm = sum((1 << rest[t]) for t in range(len(rest)) if (fx >> t) & 1)
                seen = set()
                for mm in subs:
                    seen.add(canon_of[mm | fm])
                if targets <= seen:
                    return k, F, fm
    return None

print(f"{'n':>2} {'#SC':>4} {'floor':>5} {'kappa_SC':>8} {'half-dim':>8} {'kappa_all':>9}  obstruction?")
for n in [4, 5, 6]:
    pairs, idx, P, canon_of, nclasses, sc = build(n)
    floor = ceil(log2(len(sc)))
    res = kappa_cover(P, canon_of, sc, kmax=P)
    kSC = res[0] if res else None
    half = ((n-1)**2)//4; kall = 1 + comb(n-2, 2)
    obstr = "YES (>floor)" if kSC is not None and kSC > floor else ("no (=floor)" if kSC == floor else "?")
    print(f"{n:>2} {len(sc):>4} {floor:>5} {str(kSC):>8} {half:>8} {kall:>9}  {obstr}")
    if res:
        print(f"      cover: free arcs {[pairs[b] for b in res[1]]} (dim {kSC})")
