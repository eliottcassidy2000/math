#!/usr/bin/env python3
"""Cone-restricted run: python run_km_cone.py {sl N|so2n n} kp kq Dsing 'ch1;ch2' Dtop 'topcharge' [prime]"""
import sys, time
from fractions import Fraction
from affine_km import AffineKM, sl_basis, so_even_basis
from lie_voa_ss import SingularSupport, find_singular
fam, rk, kp, kq, Dsing = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
charges = [tuple(int(x) for x in c.split(',')) for c in sys.argv[6].split(';')]
Dtop = int(sys.argv[7]); top = tuple(int(x) for x in sys.argv[8].split(','))
P = int(sys.argv[9]) if len(sys.argv) > 9 else 2147483647
t0 = time.time()
bd = sl_basis(rk) if fam == 'sl' else so_even_basis(rk)
alg = AffineKM(bd, Fraction(kp, kq), P)
sv = []
for ch in charges:
    vecs = find_singular(alg, Dsing, ch)
    print(f"{fam}{rk} k={Fraction(kp,kq)}: singular vectors at weight {Dsing}, charge {ch}: {len(vecs)}", flush=True)
    sv += vecs
ss = SingularSupport(alg, Dtop, sv)
print(f"PBW basis built ({time.time()-t0:.1f}s)", flush=True)
gm = [(1, t) for t in range(len(alg.gens))]
ex = [(0, t) for t in alg.simple_neg] + [(0, t) for t in alg.simple_pos]
ss.run(gen_modes=gm, extra_ops=ex, cone_top=(Dtop, top))
# report top block
D = Dtop
cols = ss.cols[(D, top)]
dI = len(ss.I[(D, top)])
kers = {}
for p in range(0, D + 1):
    g = ss.grI.get((D, p, top), []); Jb = ss.J.get((D, p, top), [])
    if len(g) != len(Jb):
        kers[p] = len(g) - len(Jb)
print(f"TOP block weight {D} charge {top}: dim V={len(cols)} dim I={dI} dim L={len(cols)-dI}; ker pi by Li: {kers or 'none'}")
for D2 in range(Dsing, Dtop + 1):
    tot = 0
    for (DD, p, c), g in ss.grI.items():
        if DD == D2 and c == top:
            tot += len(g) - len(ss.J.get((DD, p, c), []))
    print(f"  weight {D2}, charge {top}: dim ker pi = {tot}")
