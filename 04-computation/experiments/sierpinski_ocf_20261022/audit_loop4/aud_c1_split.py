# AUDIT: explicit counts of the two isomorphism classes hidden behind the census key (1, 1, 2, 3, 3)|9, for P_23 and
# D23_0 (and the per-TT5 change), to state the corrected (C1) row.
import sys, json
sys.path.insert(0, '.')
from aud_newclaims_lib import census5, adj10, CLS
from fractions import Fraction
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
cls = json.load(open(f'{HERE}/drt23_classes.json'))
def S_of(code):
    n = 23; S = np.zeros((n, n), dtype=int); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = 1 if code[k] == '1' else -1; S[i, j] = s; S[j, i] = -s; k += 1
    return S
cP = census5(S_of(cls[12]['code'])); c0 = census5(S_of(cls[0]['code']))
for cc in (8, 10):
    A = adj10(cc); sc = [sum(A[i]) for i in range(5)]
    mid = sc.index(2)
    desc = 'score-2 vertex beats both score-1 vertices' if all(sc[j] == 1 for j in range(5) if A[mid][j]) else 'score-2 vertex beats one score-1 and one score-3 vertex'
    print(f"class {cc} ({desc}): P_23 {cP[cc]}, D23_0 {c0[cc]}, change per TT5 {Fraction(c0[cc] - cP[cc], c0[0] - cP[0])}")
print("merged key total: P_23", cP[8] + cP[10], "D23_0", c0[8] + c0[10], "(json key (1, 1, 2, 3, 3)|9:", cls[12]['census5'].get('(1, 1, 2, 3, 3)|9'), cls[0]['census5'].get('(1, 1, 2, 3, 3)|9'), ")")
