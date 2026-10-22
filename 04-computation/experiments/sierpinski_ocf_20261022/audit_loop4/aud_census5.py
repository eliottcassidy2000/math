# AUDIT (C1) / Section 2: own 5-vertex census by TRUE isomorphism class (lookup table over all 1024 labelled
# 5-tournaments, canonical form = min over the 120 relabellings), for the 13 DRT(23) classes and the n = 19 pair.
#  * Are the 12 isomorphism classes of 5-tournaments separated by the scripts' key (sorted scores | #Ham paths)?
#  * Rank of the census differences, the per-TT5 direction at n = 23 and n = 19 (claim: identical, (C1) table);
#  * converse symmetry of the census (c_tau = c_(tau converse));
#  * 4-census against the closed form for doubly regular tournaments (HYP-9365 (d)).
import json, itertools
from fractions import Fraction
from math import comb
from collections import Counter
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
def adj_from_code10(c):
    A = [[0] * 5 for _ in range(5)]
    for b, (i, j) in enumerate(PAIRS):
        if c >> b & 1: A[i][j] = 1
        else: A[j][i] = 1
    return A
def code10(A):
    return sum(1 << b for b, (i, j) in enumerate(PAIRS) if A[i][j])
perms = list(itertools.permutations(range(5)))
canon = {}
for c in range(1024):
    A = adj_from_code10(c)
    canon[c] = min(code10([[A[p[i]][p[j]] for j in range(5)] for i in range(5)]) for p in perms)
classes = sorted(set(canon.values()))
def key_of(cc):
    A = adj_from_code10(cc)
    sc = tuple(sorted(sum(A[i]) for i in range(5)))
    h = sum(1 for p in perms if all(A[p[i]][p[i + 1]] for i in range(4)))
    return f"{sc}|{h}"
def conv(cc):
    A = adj_from_code10(cc)
    return canon[code10([[A[j][i] for j in range(5)] for i in range(5)])]
keys = {cc: key_of(cc) for cc in classes}
print(f"isomorphism classes of 5-tournaments: {len(classes)}; distinct (scores|H) keys: {len(set(keys.values()))}")
for cc in classes:
    print(f"   class {cc:4d}: key {keys[cc]:22s} converse class {conv(cc):4d} (key {keys[conv(cc)]}){'  self-converse' if conv(cc) == cc else ''}")
dup = [k for k, v in Counter(keys.values()).items() if v > 1]
print("   keys shared by >1 class:", dup if dup else 'none')

def census(code):
    m = len(code); n = int(round((1 + (1 + 8 * m) ** 0.5) / 2))
    A = np.zeros((n, n), dtype=int); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1': A[i, j] = 1
            else: A[j, i] = 1
            k += 1
    cnt = Counter()
    for X in itertools.combinations(range(n), 5):
        c = 0
        for b, (i, j) in enumerate(PAIRS):
            if A[X[i], X[j]]: c |= 1 << b
        cnt[canon[c]] += 1
    c4 = Counter()
    for X in itertools.combinations(range(n), 4):
        c4[tuple(sorted(int(A[np.ix_(X, X)][r].sum()) for r in range(4)))] += 1
    return n, cnt, c4

cls = json.load(open(f'{HERE}/drt23_classes.json'))
codes19 = {}
for line in open(f'{HERE}/drt_results.out'):
    t = line.split()
    if len(t) > 6 and t[0] == '19': codes19[t[1]] = t[6]
data = {}
for c in cls: data[c['label']] = census(c['code'])
for lab in ('P_19', 'D_19.1'): data[lab] = census(codes19[lab])
TT = [cc for cc in classes if keys[cc] == '(0, 1, 2, 3, 4)|1'][0]
print("\nconverse symmetry c_tau = c_(tau^T) for every class:", all(all(cnt[cc] == cnt[conv(cc)] for cc in classes) for (_, cnt, _) in data.values()))
# census agrees with the JSON (by key)
def by_key(cnt):
    out = Counter()
    for cc, v in cnt.items(): out[keys[cc]] += v
    return out
ok = all(by_key(data[c['label']][1]) == Counter({k: v for k, v in c['census5'].items() if v}) for c in cls)
print("own census (by key) equals drt23_classes.json census5 for all 13 classes:", ok)
for n_, ref, others in ((23, 'P_23', [c['label'] for c in cls if c['label'] != 'P_23']), (19, 'P_19', ['D_19.1'])):
    cP = data[ref][1]
    M = np.array([[data[o][1][cc] - cP[cc] for cc in classes] for o in others])
    dirs = set()
    for o in others:
        dt = data[o][1][TT] - cP[TT]
        dirs.add(tuple(Fraction(data[o][1][cc] - cP[cc], dt) for cc in classes))
    print(f"n = {n_}: rank of census differences {np.linalg.matrix_rank(M)}; distinct per-TT5 directions {len(dirs)}: "
          + ', '.join(f"{keys[cc]}:{d}" for cc, d in zip(classes, list(dirs)[0])))
    n = n_
    m = (n - 1) // 2; l = (n - 3) // 4
    tt4 = n * m * comb(l, 2); src = n * comb(m, 3)
    pred = {(0, 1, 2, 3): tt4, (1, 1, 1, 3): src - tt4, (0, 2, 2, 2): src - tt4, (1, 1, 2, 2): comb(n, 4) - 2 * src + tt4}
    c4s = {lab: dict(data[lab][2]) for lab in data if (lab.startswith('P_') and int(lab[2:].split('.')[0]) == n) or lab.startswith(f'D{n}') or lab.startswith(f'D_{n}')}
    print(f"   4-census closed form (n C((n-1)/2) ... ): {pred}; all classes of order {n} equal it: {all(v == pred for v in c4s.values())} ({len(c4s)} classes)")

print("\nthe two isomorphism classes with key (1, 1, 2, 3, 3)|9:")
for cc in classes:
    if keys[cc] != '(1, 1, 2, 3, 3)|9': continue
    A = adj_from_code10(cc)
    sc = [sum(A[i]) for i in range(5)]
    lo = [i for i in range(5) if sc[i] == 1]; hi = [i for i in range(5) if sc[i] == 3]; mid = [i for i in range(5) if sc[i] == 2][0]
    sub4 = Counter(tuple(sorted(sum(A[i][j] for j in X if j != i) for i in X)) for X in itertools.combinations(range(5), 4))
    print(f"   class {cc}: arcs {[(i, j) for i in range(5) for j in range(5) if A[i][j]]}; score-2 vertex beats "
          f"{['score-1' if sc[j] == 1 else 'score-3' for j in range(5) if A[mid][j]]}; 4-subtournament score profile {dict(sub4)}")
