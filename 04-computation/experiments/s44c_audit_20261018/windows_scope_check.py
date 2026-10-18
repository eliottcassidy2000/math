#!/usr/bin/env python3
"""Scope check of the 'windows are sharp' claims (li note section 3, HYP-9315(f)) on the session's own order tables
(orders_*.out in 04-computation/quantum_transfer_s44; read only).  For each pair of models: the largest relative level
m covered by BOTH tables, the window bounds mU = min(2p-1,2p'-1), mQ = 2 min(q,q') - 1, and the first disagreement."""
import re, os, ast, itertools, glob, json
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'quantum_transfer_s44')
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
for fn in glob.glob(os.path.join(SRC, 'orders_*.out')):
    pq, N, Dmax, T = load(fn)
    if pq and (pq not in best or Dmax > best[pq][1]):
        best[pq] = (N, Dmax, T)
def keys(pq):
    N, Dmax, _ = best[pq]
    return {(L, i) for i in range(0, 40) for L in range(0, 80) if 2 * (N + i) + L <= Dmax}
rows = []
for a, b in itertools.combinations(sorted(best), 2):
    common = keys(a) & keys(b)
    mmax = max(L + 2 * i for (L, i) in common)
    # largest m such that ALL bidegrees with that m are in both tables (complete level)
    def complete(m):
        need = {(m - 2 * i, i) for i in range(0, m // 2 + 1)}
        return need <= common
    mcomplete = max([m for m in range(0, mmax + 1) if all(complete(k) for k in range(0, m + 1))], default=-1)
    diffs = sorted((L + 2 * i, (L, i)) for (L, i) in common if best[a][2].get((L, i), {}) != best[b][2].get((L, i), {}))
    mU = min(2 * a[0] - 1, 2 * b[0] - 1)
    mQ = 2 * min(a[1], b[1]) - 1
    rows.append({'a': a, 'b': b, 'same_p': a[0] == b[0], 'mU': mU, 'mQ': mQ, 'max_common_m': mmax,
                 'complete_through_m': mcomplete, 'first_disagreement': diffs[0] if diffs else None})
same = [r for r in rows if r['same_p']]
reach = [r for r in same if r['max_common_m'] >= r['mQ']]
print(f'same-p pairs: {len(same)}; pairs whose common data reach m = 2q_min-1: {len(reach)}')
print('  all of those first disagree exactly at 2q_min-1:', all(r['first_disagreement'] and r['first_disagreement'][0] == r['mQ'] for r in reach))
for r in same:
    if r['max_common_m'] < r['mQ']:
        print(f"  untested: {r['a']} {r['b']}: common data reach m <= {r['max_common_m']} (complete through {r['complete_through_m']}) < 2q_min-1 = {r['mQ']}")
diffp = [r for r in rows if not r['same_p']]
reachU = [r for r in diffp if r['max_common_m'] >= r['mU']]
atU = [r for r in reachU if r['first_disagreement'] and r['first_disagreement'][0] == r['mU']]
print(f'different-p pairs: {len(diffp)}; reaching m = min(2p-1,2p\'-1): {len(reachU)}; first disagreement exactly at that m: {len(atU)}')
for r in reachU:
    if r not in atU:
        print(f"  {r['a']} {r['b']}: mU = {r['mU']}, first disagreement {r['first_disagreement']}, max common m {r['max_common_m']}")
json.dump(rows, open(os.path.join(HERE, 'windows_scope_check.json'), 'w'), indent=1, default=str)
