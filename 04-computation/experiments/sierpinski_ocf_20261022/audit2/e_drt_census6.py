# Audit2 (E), extra: 6-vertex subtournament censuses of the order-23 doubly regular tournaments of drt_results.out
# (D_23.2 and D_23.3 have identical 5-vertex censuses; does the difference enter at 6 vertices?).  Canonical form of a
# 6-vertex subtournament = min over the 720 vertex orders of its 15-bit upper-triangle code (numpy, chunked).
import itertools, sys
import numpy as np
sys.path.insert(0, '.')
from common import code_to_A
rows = [l.split() for l in open('../drt_results.out') if l.strip() and not l.startswith('#')]
k = 6
perms = np.array(list(itertools.permutations(range(k))), dtype=np.int64)          # 720 x 6
pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
PI = perms[:, [i for i, j in pairs]]; PJ = perms[:, [j for i, j in pairs]]          # 720 x 15
W = (1 << np.arange(len(pairs), dtype=np.int64))
res = {}
for r in rows:
    n, lab, code = int(r[0]), r[1], r[-1]
    if n != 23: continue
    A = code_to_A(code).astype(np.int8)
    subs = np.array(list(itertools.combinations(range(n), k)), dtype=np.int64)     # C(23,6) x 6
    cnt = {}
    for s in range(0, len(subs), 4000):
        X = subs[s:s + 4000]                                                         # B x 6
        U = X[:, PI]; V = X[:, PJ]                                                   # B x 720 x 15
        bits = A[U, V].astype(np.int64)
        codes = (bits * W).sum(axis=2).min(axis=1)
        u, c = np.unique(codes, return_counts=True)
        for a, b in zip(u.tolist(), c.tolist()): cnt[a] = cnt.get(a, 0) + b
    res[lab] = cnt
    print(f"{lab}: {len(cnt)} isomorphism types of 6-subtournaments present", flush=True)
labs = list(res)
keys = sorted(set().union(*[set(c) for c in res.values()]))
for i in range(len(labs)):
    for j in range(i + 1, len(labs)):
        diff = sum(1 for kk in keys if res[labs[i]].get(kk, 0) != res[labs[j]].get(kk, 0))
        print(f"  {labs[i]} vs {labs[j]}: 6-vertex census {'identical' if diff == 0 else f'differs in {diff} of {len(keys)} types'}")
