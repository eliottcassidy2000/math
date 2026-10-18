#!/usr/bin/env python3
"""python run_km.py {sl N | so2n n} kp kq Dsing 'charge' Dmax [prime]
   e.g. run_km.py sl 3 -3 2 2 '1,1' 6"""
import sys, ast
from fractions import Fraction
from affine_km import AffineKM, sl_basis, so_even_basis
from lie_voa_ss import SingularSupport, find_singular
fam, rk, kp, kq, Dsing = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
charges = [tuple(int(x) for x in c.split(',')) for c in sys.argv[6].split(';')]
Dmax = int(sys.argv[7]); P = int(sys.argv[8]) if len(sys.argv) > 8 else 2147483647
bd = sl_basis(rk) if fam == 'sl' else so_even_basis(rk)
alg = AffineKM(bd, Fraction(kp, kq), P)
sv = []
for ch in charges:
    vecs = find_singular(alg, Dsing, ch)
    print(f"{fam}{rk} k={Fraction(kp,kq)}: singular vectors at weight {Dsing}, charge {ch}: {len(vecs)}")
    sv += vecs
ss = SingularSupport(alg, Dmax, sv)
gm = [(1, t) for t in range(len(alg.gens))]
ex = [(0, t) for t in alg.simple_neg] + [(0, t) for t in alg.simple_pos]
ss.run(gen_modes=gm, extra_ops=ex)
print("generators (D, Li, charge, #):", [g for g in ss.gens_found if g[1] > 0] or "none at positive Li")
print("Li-0 generators (D, #):", sorted({(g[0], sum(h[3] for h in ss.gens_found if h[0]==g[0] and h[1]==0)) for g in ss.gens_found if g[1]==0}))
