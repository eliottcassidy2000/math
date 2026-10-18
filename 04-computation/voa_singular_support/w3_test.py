import random, sys
from fractions import Fraction
from w3 import W3
from lie_voa_ss import *
from max_submodule import simple_quotient_dims, per_weight
P = 2147483647
p_, q_ = int(sys.argv[1]), int(sys.argv[2])
c = 2 * (1 - Fraction(12 * (p_ - q_) ** 2, p_ * q_))
alg = W3(c, P)
basis = alg.pbw_basis(6)
states = [mo for key, ms in basis.items() for mo in ms]
modes = [(w, t) for t in (0, 1) for w in range(alg.h[t] - 4, alg.h[t] + 2)]
random.seed(3); bad = 0
for _ in range(400):
    X, Y, v = random.choice(modes), random.choice(modes), random.choice(states)
    lhs = alg.apply_vec(X, alg.apply_vec(Y, {v: 1})); r2 = alg.apply_vec(Y, alg.apply_vec(X, {v: 1}))
    diff = dict(lhs)
    for kq, val in r2.items(): diff[kq] = (diff.get(kq, 0) - val) % P
    terms, cen = alg.bracket(X, Y); br = {}
    for (m, co) in terms:
        for kq, val in alg.apply(m, v).items(): br[kq] = (br.get(kq, 0) + co * val) % P
    if cen: br[v] = (br.get(v, 0) + cen) % P
    for kq, val in alg.bracket_extra(X, Y, v).items(): br[kq] = (br.get(kq, 0) + val) % P
    if any((diff.get(q, 0) - br.get(q, 0)) % P for q in set(diff) | set(br)): bad += 1
print(f"W3({p_},{q_}) c={c}: bracket consistency failures {bad}/400")
d = per_weight(simple_quotient_dims(alg, 10, [(-1, 0), (-2, 0), (-1, 1), (-2, 1)]))
print("max-submodule dims by weight:", d)
