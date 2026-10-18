import sys, random
from fractions import Fraction
from n2_sca import N2SCA, L, J, GP, GM
from lie_voa_ss import *
P = 2147483647
kk = int(sys.argv[1]); Dmax = Fraction(sys.argv[2])
c = Fraction(3 * kk, kk + 2)
alg = N2SCA(c, P)
# consistency test
basis = alg.pbw_basis(4)
states = [mo for key, ms in basis.items() for mo in ms]
modes = sorted({(w if not isinstance(w, Fraction) or w.denominator != 1 else int(w), t) for t in range(4) for w in [alg.h[t] - 3, alg.h[t] - 2, alg.h[t] - 1, alg.h[t], alg.h[t] + 1]})
random.seed(2); bad = 0
for _ in range(600):
    X, Y, v = random.choice(modes), random.choice(modes), random.choice(states)
    s = -1 if (alg.par[X[1]] and alg.par[Y[1]]) else 1
    lhs = alg.apply_vec(X, alg.apply_vec(Y, {v: 1})); r2 = alg.apply_vec(Y, alg.apply_vec(X, {v: 1}))
    diff = dict(lhs)
    for kq, val in r2.items(): diff[kq] = (diff.get(kq, 0) - s * val) % P
    terms, cen = alg.bracket(X, Y); br = {}
    for (m, co) in terms:
        for kq, val in alg.apply(m, v).items(): br[kq] = (br.get(kq, 0) + co * val) % P
    if cen: br[v] = (br.get(v, 0) + cen) % P
    if any((diff.get(q, 0) - br.get(q, 0)) % P for q in set(diff) | set(br)): bad += 1
print(f"N=2 c={c}: bracket consistency failures {bad}/600")
basis = alg.pbw_basis(int(Dmax) + 1)
for key in sorted(k for k in basis if 0 < k[0] <= Dmax):
    vecs = find_singular(alg, key[0], key[1])
    if vecs:
        print(f"  singular: weight {key[0]}, charge {key[1]}, dim {len(vecs)}, block {len(basis[key])}")
