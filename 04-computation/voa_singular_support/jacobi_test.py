import random, itertools
from fractions import Fraction
from small_n4 import SmallN4, ODD
from lie_voa_ss import *
P = 2147483647
alg = SmallN4(-9, P)
basis = alg.pbw_basis(4)
states = [mo for key, ms in basis.items() for mo in ms]
random.seed(1)
modes = []
for t in range(8):
    h = alg.h[t]
    for w in [h - 3, h - 2, h - 1, h, h + 1, h - Fraction(1,1)*0]:
        modes.append((w if not isinstance(w, Fraction) or w.denominator != 1 else int(w), t))
modes = sorted(set(modes))
bad = 0; tot = 0
for _ in range(1500):
    X = random.choice(modes); Y = random.choice(modes); v = random.choice(states)
    s = -1 if (alg.par[X[1]] and alg.par[Y[1]]) else 1
    lhs = alg.apply_vec(X, alg.apply_vec(Y, {v: 1}))
    rhs2 = alg.apply_vec(Y, alg.apply_vec(X, {v: 1}))
    diff = dict(lhs)
    for kk, val in rhs2.items():
        diff[kk] = (diff.get(kk, 0) - s * val) % P
    terms, cen = alg.bracket(X, Y)
    br = {}
    for (m, co) in terms:
        for kk, val in alg.apply(m, v).items():
            br[kk] = (br.get(kk, 0) + co * val) % P
    if cen:
        br[v] = (br.get(v, 0) + cen) % P
    keys = set(diff) | set(br)
    ok = all((diff.get(kk, 0) - br.get(kk, 0)) % P == 0 for kk in keys)
    tot += 1
    if not ok:
        bad += 1
        if bad <= 5: print("FAIL", X, Y, v)
print("Jacobi/consistency checks:", tot, "failures:", bad)
