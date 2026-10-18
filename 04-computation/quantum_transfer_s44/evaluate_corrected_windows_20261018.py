#!/usr/bin/env python3
"""Post hoc check (NOT pre-registered) of the audit-corrected stability windows on all S44 transfer-order tables.
m = L + 2i (relative level).  (U) universality: T_(p,q)(L,i) = T_(p',q')(L,i) for m < min(2p-1, 2p'-1).
(Q) q-stability: same p, m < 2 min(q,q') - 1.  (S) saturation: orders k <= min(p,p') - 2 agree for m < 2 min(q,q') - 1.
Also reports the first disagreement beyond each window (to show the windows are sharp).  Session opus-2026-10-18-S44.
"""
import re, os, ast, itertools, glob
HERE = os.path.dirname(os.path.abspath(__file__))
def load(fn):
    T, N, Dmax, pq = {}, None, None, None
    for line in open(fn):
        m = re.match(r'=== Vir_\((\d+),(\d+)\).*N = (\d+)', line)
        if m: pq = (int(m.group(1)), int(m.group(2))); N = int(m.group(3))
        m = re.match(r'D=\s*(\d+) (.*)', line)
        if m:
            Dmax = int(m.group(1))
            for part in m.group(2).split(';'):
                mm = re.search(r's=(\d+)\(Li (\d+)\) ker (\d+) (\{[^}]*\})', part)
                if mm: T[(int(mm.group(2)), int(mm.group(1)) - N)] = ast.literal_eval(mm.group(4))
    return pq, N, Dmax, T
best = {}
for fn in glob.glob(os.path.join(HERE, 'orders_*.out')):
    pq, N, Dmax, T = load(fn)
    if pq and (pq not in best or Dmax > best[pq][1]):
        best[pq] = (N, Dmax, T)
def val(pq, k): return best[pq][2].get(k, {})
def keys(pq):
    N, Dmax, _ = best[pq]
    return {(L, i) for i in range(0, 40) for L in range(0, 60) if 2 * (N + i) + L <= Dmax}
res = {'U': [0, 0], 'Q': [0, 0], 'S': [0, 0]}
beyond = []
for a, b in itertools.combinations(sorted(best), 2):
    common = keys(a) & keys(b)
    mU = min(2 * a[0] - 1, 2 * b[0] - 1)
    mQ = 2 * min(a[1], b[1]) - 1
    for (L, i) in sorted(common, key=lambda k: (k[0] + 2 * k[1], k)):
        m = L + 2 * i
        va, vb = val(a, (L, i)), val(b, (L, i))
        if not (va or vb):
            continue
        if m < mU:
            res['U'][0] += 1; res['U'][1] += (va != vb)
        if a[0] == b[0] and m < mQ:
            res['Q'][0] += 1; res['Q'][1] += (va != vb)
        if m < mQ:
            kmax = min(a[0], b[0]) - 2
            for k in range(1, kmax + 1):
                res['S'][0] += 1; res['S'][1] += (va.get(k, 0) != vb.get(k, 0))
    # first disagreement for same-p pairs
    if a[0] == b[0]:
        diffs = sorted(((L + 2 * i, (L, i)) for (L, i) in common if val(a, (L, i)) != val(b, (L, i))))
        beyond.append((a, b, mQ, diffs[0] if diffs else None))
print('models:', sorted(best))
print(f"(U) universality, m < min(2p-1,2p'-1): {res['U'][0]} nonempty bidegree comparisons, {res['U'][1]} disagreements")
print(f"(Q) q-stability, m < 2min(q,q')-1:     {res['Q'][0]} comparisons, {res['Q'][1]} disagreements")
print(f"(S) saturation k <= p_min-2, m < 2q_min-1: {res['S'][0]} order-count comparisons, {res['S'][1]} disagreements")
print('same-p pairs: (pair, window bound 2q_min-1, first disagreement (m,(L,i)) if any):')
for a, b, mQ, d in beyond:
    print('  ', a, b, mQ, d)
