# Direction of the 5-vertex census differences of doubly regular tournaments per transitive 5-set: n = 23 (12
# non-Paley classes against P_23, drt23_classes.json) and n = 19 (D_19 against P_19, drt_results.out).  If the vector
# Delta census / Delta TT5 is the same at both orders, the 5-census of a DRT is (on this evidence) affine in TT5 with
# an n-independent direction.  Census key: sorted score sequence | number of Hamiltonian paths.  mac-mini-2026-10-22-loop4.
import sys, json, itertools
from collections import Counter
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from drt_switch import code_to_S
perms5 = list(itertools.permutations(range(5)))
def census5(S):
    n = S.shape[0]; cnt = Counter()
    for X in itertools.combinations(range(n), 5):
        sub = S[np.ix_(X, X)]
        sc = tuple(sorted(int((sub[i] == 1).sum()) for i in range(5)))
        h = sum(1 for p in perms5 if all(sub[p[i], p[i + 1]] == 1 for i in range(4)))
        cnt[f'{sc}|{h}'] += 1
    return cnt
cls = json.load(open('drt23_classes.json'))
P23 = [c for c in cls if c['label'] == 'P_23'][0]
keys = sorted(set(k for c in cls for k in c['census5']))
dirs23 = set()
for c in cls:
    if c is P23: continue
    dt = c['TT5'] - P23['TT5']
    dirs23.add(tuple(Fraction(c['census5'].get(k, 0) - P23['census5'].get(k, 0), dt) for k in keys))
print(f"n = 23: {len(dirs23)} distinct direction vector(s) Delta census / Delta TT5 over 12 classes")
d23 = list(dirs23)[0]
codes19 = {}
for line in open('drt_results.out'):
    t = line.split()
    if len(t) > 6 and t[0] == '19': codes19[t[1]] = t[6]
cP, cD = census5(code_to_S(codes19['P_19'])), census5(code_to_S(codes19['D_19.1']))
keys19 = sorted(set(cP) | set(cD) | set(keys))
dt19 = cD['((0, 1, 2, 3, 4), 1)'.replace("'", '')] if False else None
tt = '(0, 1, 2, 3, 4)|1'
dt19 = cD[tt] - cP[tt]
d19 = {k: Fraction(cD.get(k, 0) - cP.get(k, 0), dt19) for k in keys19}
print(f"n = 19: Delta TT5 = {dt19}")
print(f"{'census type (scores|H)':28s} {'n=23 per TT5':>14s} {'n=19 per TT5':>14s}")
same = True
for k in keys19:
    a = dict(zip(keys, d23)).get(k, Fraction(0)); b = d19.get(k, Fraction(0))
    same &= a == b
    print(f"{k:28s} {str(a):>14s} {str(b):>14s}")
print('identical directions at n = 19 and n = 23:', same)
