#!/usr/bin/env python3
"""Mechanical evaluation of the pre-registered predictions P4-P8 (PREREGISTRATION_stability_20261018.md).
Reads orders_*.out (transfer-order tables) and checks q-stability, the length-N row, saturation and monotonicity in p
on the stable windows W(N) = {(L, i): L <= N + i - 1}.  Session opus-2026-10-18-S44.
"""
import re, os, ast, itertools, json
HERE = os.path.dirname(os.path.abspath(__file__))
FILES = {(3, 8): 'orders_38_D30.out', (3, 10): 'orders_310_D30.out', (3, 11): 'orders_311_D32.out',
         (3, 13): 'orders_313_D36.out', (4, 7): 'orders_47_D30.out', (4, 9): 'orders_49_D36.out',
         (4, 11): 'orders_411_D42.out', (5, 6): 'orders_56_D30.out', (5, 7): 'orders_57_D36.out',
         (5, 8): 'orders_58_D40.out', (5, 9): 'orders_59_D42.out', (6, 7): 'orders_67_D41.out'}
NEW = {(3, 11), (3, 13), (4, 9), (4, 11), (5, 7), (5, 8), (5, 9), (6, 7)}

def load(fn):
    T, N, Dmax, mism = {}, None, None, 0
    for line in open(os.path.join(HERE, fn)):
        m = re.match(r'=== Vir_\((\d+),(\d+)\).*N = (\d+)', line)
        if m:
            N = int(m.group(3))
        if 'MISMATCH' in line:
            mism += 1
        m = re.match(r'D=(\d+) (.*)', line)
        if m:
            D = int(m.group(1)); Dmax = D
            for part in m.group(2).split(';'):
                mm = re.search(r's=(\d+)\(Li (\d+)\) ker (\d+) (\{[^}]*\})', part)
                if mm:
                    s, L, k = int(mm.group(1)), int(mm.group(2)), int(mm.group(3))
                    T[(L, s - N)] = ast.literal_eval(mm.group(4))
                    assert sum(T[(L, s - N)].values()) == k
    return T, N, Dmax, mism

data = {pq: load(fn) for pq, fn in FILES.items()}
def computed(pq):
    T, N, Dmax, _ = data[pq]
    # every (L,i) with 2(N+i)+L <= Dmax counts as computed (absent entry = empty multiset)
    return {(L, i) for i in range(0, Dmax) for L in range(0, Dmax) if 2 * (N + i) + L <= Dmax and L >= 0}
def val(pq, key):
    return data[pq][0].get(key, {})
def stable(pq, key):
    L, i = key
    return L <= data[pq][1] + i - 1

res = {'P4': [], 'P5': [], 'P6': [], 'P7': [], 'P8': []}
# P8
for pq in FILES:
    res['P8'].append((pq, data[pq][3] == 0))
# P5 (new models)
for pq in sorted(NEW):
    p = pq[0]; T, N, Dmax, _ = data[pq]
    ok = True
    for L in range(3, Dmax - 2 * N + 1):
        m = min(p - 2, (L - 1) // 2)
        want = {k: 1 for k in range(1, m + 1)}
        if val(pq, (L, 0)) != want:
            ok = False
            res['P5'].append((pq, L, val(pq, (L, 0)), want))
    res['P5'].append((pq, ok))
# P4
fam = {}
for pq in FILES:
    fam.setdefault(pq[0], []).append(pq)
for p, lst in fam.items():
    for a, b in itertools.combinations(sorted(lst), 2):
        common = [k for k in computed(a) & computed(b) if stable(a, k) and stable(b, k) and (k[0] >= 3 or k[1] > 0)]
        bad = [(k, val(a, k), val(b, k)) for k in sorted(common) if val(a, k) != val(b, k)]
        if a in NEW or b in NEW:
            res['P4'].append((a, b, len(common), bad))
# P6, P7 across p
pqs = sorted(FILES)
for a, b in itertools.combinations(pqs, 2):
    if a[0] == b[0]:
        continue
    lo, hi = (a, b) if a[0] < b[0] else (b, a)
    if not (lo in NEW or hi in NEW):
        continue
    common = [k for k in computed(lo) & computed(hi) if stable(lo, k) and stable(hi, k)]
    bad6, bad7 = [], []
    for k in sorted(common):
        vl, vh = val(lo, k), val(hi, k)
        for kk in set(vl) | set(vh):
            if kk <= lo[0] - 2 and vl.get(kk, 0) != vh.get(kk, 0):
                bad6.append((k, kk, vl, vh))
            if vl.get(kk, 0) > vh.get(kk, 0):
                bad7.append((k, kk, vl, vh))
    res['P6'].append((lo, hi, len(common), bad6))
    res['P7'].append((lo, hi, len(common), bad7))

def summary(name, items, badidx):
    fails = [it for it in items if it[badidx]]
    return fails
print('P8 no mismatches:', all(ok for _, ok in res['P8']))
print('P5 length-N row:', [x for x in res['P5'] if len(x) == 2])
p5fails = [x for x in res['P5'] if len(x) == 4]
if p5fails: print('  P5 FAILURES:', p5fails)
print('P4 q-stability (pair, #common stable bidegrees, failures):')
for a, b, n, bad in res['P4']:
    print(f'   {a} vs {b}: {n} common, {len(bad)} failures', bad[:4])
print('P6 saturation (k <= p-2):')
for lo, hi, n, bad in res['P6']:
    print(f'   {lo} vs {hi}: {n} common, {len(bad)} failures', bad[:3])
print('P7 monotonicity:')
for lo, hi, n, bad in res['P7']:
    print(f'   {lo} vs {hi}: {n} common, {len(bad)} failures', bad[:3])
json.dump({k: [str(x) for x in v] for k, v in res.items()}, open(os.path.join(HERE, 'evaluate_stability_prereg.json'), 'w'), indent=1)

# ---- informative counts: bidegrees where at least one side is nonempty
print('\nNonempty common stable bidegrees per P4 pair (where the test has content):')
for a, b, n, bad in res['P4']:
    common = [k for k in computed(a) & computed(b) if stable(a, k) and stable(b, k)]
    ne = [k for k in common if val(a, k) or val(b, k)]
    hi = [k for k in ne if max(list(val(a, k)) + list(val(b, k)) + [0]) >= 2]
    print(f'   {a} vs {b}: {len(ne)} nonempty, {len(hi)} with an order >= 2 class')
print('Nonempty common stable bidegrees with an order-k class, k <= p_lo - 2, k >= 2 (P6 content):')
for lo, hi_, n, bad in res['P6']:
    common = [k for k in computed(lo) & computed(hi_) if stable(lo, k) and stable(hi_, k)]
    cont = [k for k in common if any(kk >= 2 and kk <= lo[0] - 2 for kk in list(val(lo, k)) + list(val(hi_, k)))]
    if cont:
        print(f'   {lo} vs {hi_}: {len(cont)} bidegrees with order 2..{lo[0]-2} content')
