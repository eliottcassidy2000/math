#!/usr/bin/env python3
"""Which tournaments have VC dimension 1?  (mac-mini-2026-07-20-S131)
The S131 run found EXACTLY TWO iso classes with VC = 1 at every n = 4,5,6,7.
This identifies them and tests the structural characterisation."""
import numpy as np
from itertools import permutations, combinations

def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return pairs, {p: k for k, p in enumerate(pairs)}, len(pairs)

def out_nbrs(code, pairs, n):
    N = [0] * n
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: N[j] |= 1 << i
        else:             N[i] |= 1 << j
    return N

def vc(code, pairs, n):
    N = out_nbrs(code, pairs, n); best = 0
    for k in range(1, n + 1):
        if (1 << k) > n: break
        if any(len({tuple((N[v] >> s) & 1 for s in S) for v in range(n)}) == (1 << k)
               for S in combinations(range(n), k)): best = k
        else: break
    return best

def canon_codes(n):
    reps = {0}
    for k in range(2, n + 1):
        pk, ik, Ek = scaffold(k)
        op, _, _ = scaffold(k - 1)
        cand = []
        for r in reps:
            base = 0
            for e, (i, j) in enumerate(op):
                if r >> e & 1: base |= 1 << ik[(i, j)]
            for mask in range(1 << (k - 1)):
                v = base
                for b in range(k - 1):
                    if mask >> b & 1: v |= 1 << ik[(b, k - 1)]
                cand.append(v)
        p2 = (1 << np.arange(Ek, dtype=np.int64))
        A = ((np.array(cand, dtype=np.int64)[:, None] >> np.arange(Ek)) & 1).astype(np.uint8)
        best = None
        for p in permutations(range(k)):
            src = np.empty(Ek, dtype=np.int64); fl = np.zeros(Ek, dtype=np.uint8)
            for e, (i, j) in enumerate(pk):
                a, b = p[i], p[j]
                t = ik[(min(a, b), max(a, b))]
                src[t] = e; fl[t] = 1 if a > b else 0
            c = (A[:, src] ^ fl) @ p2
            best = c if best is None else np.minimum(best, c)
        reps = set(int(x) for x in best.tolist())
    return sorted(reps)

for n in range(4, 8):
    pairs, idx, E = scaffold(n)
    low = [r for r in canon_codes(n) if vc(r, pairs, n) == 1]
    print(f"n={n}: {len(low)} classes with VC=1")
    for r in low:
        N = out_nbrs(r, pairs, n)
        sc = sorted(bin(x).count('1') for x in N)
        # structure probes
        transitive = (sorted(sc) == list(range(n)))
        # 3-cycle count
        c3 = n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in sc)
        # is it "transitive plus one reversed arc"?  count arcs to flip to reach transitive
        print(f"    code {r:>7}  scores {sc}  3-cycles {c3}  transitive={transitive}")
    # the structural test: VC>=2 iff some arc a->b has BOTH a vertex beating a not b
    # AND a vertex beating both.
    def has_pattern(r):
        N = out_nbrs(r, pairs, n)
        for a in range(n):
            for b in range(n):
                if a == b or not (N[a] >> b & 1): continue      # need a->b
                up  = any((N[u] >> a & 1) and not (N[u] >> b & 1) for u in range(n))
                bot = any((N[w] >> a & 1) and (N[w] >> b & 1) for w in range(n))
                if up and bot: return True
        return False
    ok = all((vc(r, pairs, n) >= 2) == has_pattern(r) for r in canon_codes(n))
    print(f"    structural criterion (exists a->b with a beater-of-a-only and a "
          f"beater-of-both) matches VC>=2 on all classes: {ok}")
