# Direction of the 5-vertex census differences of doubly regular tournaments per transitive 5-set: n = 23 (12
# non-Paley classes against P_23, drt23_classes.json) and n = 19 (D_19 against P_19, drt_results.out).  The census is
# by isomorphism type of the 5-subtournament (complete invariant: canonical form = minimum 10-bit code over the 120
# relabellings, via a lookup table over all 1024 labelled 5-tournaments).  The older key (sorted scores | number of
# Hamiltonian paths) merges two of the 12 types (audit_loop4: both have key (1,1,2,3,3)|9).
# Output: drt_census5_direction.out, drt_census5_direction.json (the per-TT5 direction, used by drt_c1_test.py).
# mac-mini-2026-10-22-loop4 (rewritten after audit_loop4/REPORT.md).
import sys, json, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from drt_switch import code_to_S

PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
PERMS = list(itertools.permutations(range(5)))
def _adj(c):
    A = [[0] * 5 for _ in range(5)]
    for b, (i, j) in enumerate(PAIRS):
        if c >> b & 1: A[i][j] = 1
        else: A[j][i] = 1
    return A
def _code(A): return sum(1 << b for b, (i, j) in enumerate(PAIRS) if A[i][j])
CANON = [min(_code([[A[p[i]][p[j]] for j in range(5)] for i in range(5)]) for p in PERMS) for A in map(_adj, range(1024))]
def _label(c):
    A = _adj(c); sc = tuple(sorted(sum(r) for r in A))
    h = sum(1 for p in PERMS if all(A[p[i]][p[i + 1]] for i in range(4)))
    lab = f"{sc}|{h}"
    if sc == (1, 1, 2, 3, 3):                      # the two types sharing (1,1,2,3,3)|9
        v2 = [v for v in range(5) if sum(A[v]) == 2][0]
        beats = sorted(sum(A[w]) for w in range(5) if A[v2][w])
        lab += '|score-2 vertex beats scores ' + ','.join(map(str, beats))
    return lab
LABEL = {c: _label(c) for c in set(CANON)}
TT5_LABEL = LABEL[CANON[_code([[1 if j > i else 0 for j in range(5)] for i in range(5)])]]

def census5_full(S):
    """counts of 5-subtournaments by isomorphism type (labels as above); S skew +-1"""
    n = S.shape[0]; cnt = {}
    A = (S == 1)
    for X in itertools.combinations(range(n), 5):
        c = 0
        for b, (i, j) in enumerate(PAIRS):
            if A[X[i], X[j]]: c |= 1 << b
        lab = LABEL[CANON[c]]; cnt[lab] = cnt.get(lab, 0) + 1
    return cnt

if __name__ == '__main__':
    assert len(LABEL) == 12
    cls = json.load(open('drt23_classes.json'))
    P = [c for c in cls if c['label'] == 'P_23'][0]
    cP = census5_full(code_to_S(P['code']))
    keys = sorted(LABEL.values())
    dirs = set()
    rows = []
    for c in cls:
        if c is P: continue
        cD = census5_full(code_to_S(c['code']))
        dt = cD[TT5_LABEL] - cP[TT5_LABEL]
        assert dt == c['TT5'] - P['TT5']
        rows.append([cD.get(k, 0) - cP.get(k, 0) for k in keys])
        dirs.add(tuple(Fraction(cD.get(k, 0) - cP.get(k, 0), dt) for k in keys))
    print(f"n = 23: 12 isomorphism types; rank of the census differences {np.linalg.matrix_rank(np.array(rows, float))}; "
          f"{len(dirs)} distinct direction vector(s) per TT5 over 12 classes")
    d23 = dict(zip(keys, list(dirs)[0]))
    codes19 = {}
    for line in open('drt_results.out'):
        t = line.split()
        if len(t) > 6 and t[0] == '19': codes19[t[1]] = t[6]
    c19P, c19D = census5_full(code_to_S(codes19['P_19'])), census5_full(code_to_S(codes19['D_19.1']))
    dt19 = c19D[TT5_LABEL] - c19P[TT5_LABEL]
    d19 = {k: Fraction(c19D.get(k, 0) - c19P.get(k, 0), dt19) for k in keys}
    print(f"n = 19: Delta TT5 = {dt19}")
    print(f"{'isomorphism type (scores|Hamiltonian paths)':62s} {'n=23 per TT5':>13s} {'n=19 per TT5':>13s}")
    for k in keys:
        print(f"{k:62s} {str(d23[k]):>13s} {str(d19[k]):>13s}")
    same = all(d23[k] == d19[k] for k in keys)
    print('identical directions at n = 19 and n = 23:', same)
    json.dump({k: str(v) for k, v in d23.items()}, open('drt_census5_direction.json', 'w'), indent=1)
