# AUDIT: "The Paley switching families of P_27 and P_31 are single classes" (HYP-9371 (C1) evidence).  Own GF(27)
# (GF(3)[x]/(x^3 + 2x + 1), an irreducible cubic chosen independently of drt_c1_test.py) and P_31; every row
# normalisation of the bordered skew-Hadamard matrix, nauty classes.
import sys, itertools
import numpy as np
sys.path.insert(0, '.')
from aud_newclaims_lib import paley_S, skew_hadamard_from_drt, drts_from, canon
def gf27_S():
    # elements a + b x + c x^2, modulus x^3 + 2x + 1 (irreducible over GF(3): no roots 0,1,2)
    assert all((t ** 3 + 2 * t + 1) % 3 for t in range(3))
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    def mul(u, v):
        p = [0] * 5
        for i in range(3):
            for j in range(3): p[i + j] += u[i] * v[j]
        for d in (4, 3):                       # x^3 = -2x - 1 = x + 2
            c = p[d]; p[d] = 0; p[d - 3] += 2 * c; p[d - 2] += c
        return tuple(x % 3 for x in p[:3])
    sq = {mul(u, u) for u in els if u != (0, 0, 0)}
    assert len(sq) == 13 and (2, 0, 0) not in sq          # -1 is a non-square in GF(27): tournament
    S = np.array([[0 if u == v else (1 if tuple((b - a) % 3 for a, b in zip(u, v)) in sq else -1) for v in els] for u in els])
    return S
for name, S in (('P_27', gf27_S()), ('P_31', paley_S(31))):
    n = S.shape[0]
    assert (S @ S.T == n * np.eye(n, dtype=int) - 1).all()
    fam = drts_from(skew_hadamard_from_drt(S))
    c = canon([S] + fam)
    print(f"{name}: doubly regular; {len(fam)} row normalisations; nauty classes in family incl. Paley: {len(set(c))}; all isomorphic to Paley: {set(c) == {c[0]}}")
