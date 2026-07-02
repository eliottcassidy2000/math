#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3954: THE TORUS-BAND THEOREM. kind-pasteur-2026-07-01-S31.

The c-averaged ledger A(U) = E_c[ meas{x : ||u x - c|| >= h for all u in U} ]  (h = 1/14, w = 2h = 1/7)
equals 1 - vol( Union_i B_i ) on the (x,c)-torus, B_i = {(x,c) : ||c - u_i x|| < h} the band of integer
slope u_i and width w.

THEOREM (elementary, this session):
 (i)   vol(B_i) = w.
 (ii)  PAIR INDEPENDENCE: (x,c) -> (c - u_i x, c - u_j x) has determinant u_j - u_i != 0, hence is a
       measure-preserving surjection T^2 -> T^2:  vol(B_i ∩ B_j) = w^2 EXACTLY. (The THM-594-B
       fluctuations live on the x-circle at fixed c; the c-average removes them.)
 (iii) d-FOLD = SUBTORUS BOX VOLUME: the joint law of (c - u_i x)_{i in S} is Haar on the 2-dim subtorus
       cut out by the SATURATED sum-zero relation lattice Λ_S = {m : Σm_i = 0, Σ m_i u_i = 0} (rank d-2);
       by Poisson/orthogonality
           vol(∩_S B_i) = Σ_{m ∈ Λ_S} Π_i  ŵ_h(m_i),   ŵ_h(0) = 2h, ŵ_h(s) = sin(2π s h)/(π s).
       For a triple with primitive relation (1,-2,1) (AP): closed form 2h^2 (computed by hand; the
       series must agree). All terms are exact rationals at h = 1/14 in principle (sin at π-rationals).
 (iv)  A(U) is translation- and dilation-invariant (a difference-pattern functional); small-coefficient
       relations (AP) carry the max overlap weight and minimize A.
Consequence: HYP-3953's (⋆)-census is SYMBOLIC — finite inclusion-exclusion, no sampling.

TESTS: V1 pair independence (grid). V2 triple Fourier vs grid vs the 2h^2 hand value. V3 the full
inclusion-exclusion vs the S30 sampled ledger (36/49 and 61/98 checks). V4 Bonferroni-3 floors for the
S30 argmin patterns k=4..13 vs witnessMP (all triples via the Fourier engine).
"""
import sys, itertools, random
from math import gcd, floor, sin, pi
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
rng = random.Random(3954)
h = 1.0/14.0; w = 2*h

# ------------------------------------------------------------------ Fourier engine (d-fold overlaps)
def what(s):
    return 2*h if s == 0 else sin(2*pi*s*h)/(pi*s)

_tricache = {}
def vol3_fourier(U, T=4000):
    """triple overlap via the rank-1 saturated relation: m = ((u3-u2), (u1-u3), (u2-u1))/g."""
    u1, u2, u3 = sorted(U)
    al, be, ga = u3-u2, u1-u3, u2-u1
    g = gcd(gcd(abs(al), abs(be)), abs(ga))
    key = tuple(sorted((abs(al)//g, abs(be)//g, abs(ga)//g)))
    if key in _tricache: return _tricache[key]
    a_, b_, c_ = al//g, be//g, ga//g
    tot = (2*h)**3
    for t in range(1, T+1):
        tot += 2 * what(a_*t) * what(b_*t) * what(c_*t)
    _tricache[key] = tot
    return tot

def dfold_grid(U, ngrid=1000):
    """d-fold overlap by exact-pushforward sheets + fine grid (verification engine)."""
    U = sorted(U); D = U[1] - U[0]
    tot = 0.0
    for r in range(D):
        cnt = 0
        for ia in range(ngrid):
            a = -h + 2*h*(ia+0.5)/ngrid
            for ib in range(ngrid):
                b = -h + 2*h*(ib+0.5)/ngrid
                x = (a - b + r)/D
                c = a + U[0]*x
                ok = True
                for ui in U[2:]:
                    y = c - ui*x
                    y = y - floor(y + 0.5)
                    if abs(y) >= h: ok = False; break
                if ok: cnt += 1
        tot += (2*h)**2 * cnt/(ngrid*ngrid)
    return tot/D

print("="*100)
print(" V1: PAIR INDEPENDENCE  vol(B_i ∩ B_j) = w^2 = 1/49 = 0.02040816  (theorem; grid check)")
print("="*100)
for (a_, b_) in [(1,2), (1,14), (3,10), (5,37), (7,21), (100,101)]:
    v = dfold_grid((a_, b_), ngrid=1000)
    print(f"  slopes ({a_:3d},{b_:3d}): grid vol = {v:.6f}")

print("\n" + "="*100)
print(" V2: TRIPLES — Fourier vs grid; AP relation (1,-2,1): hand value 2h^2 = 1/98 = 0.01020408")
print("="*100)
for U in [(0,1,2), (0,7,14), (10,24,38), (0,1,3), (0,2,5), (0,1,7), (0,3,11)]:
    vf = vol3_fourier(U)
    vg = dfold_grid(U, ngrid=900)
    al, be, ga = U[2]-U[1], U[0]-U[2], U[1]-U[0]
    g = gcd(gcd(abs(al), abs(be)), abs(ga))
    rel = (al//g, be//g, ga//g)
    print(f"  U={str(U):>12} rel {str(rel):>12}:  fourier = {vf:.6f}   grid = {vg:.6f}"
          + ("   [2h^2 = 0.010204]" if abs(vf - 2*h*h) < 1e-6 else ""))

print("\n" + "="*100)
print(" V3: FULL INCLUSION-EXCLUSION vs the S30 sampled ledger")
print("="*100)
def A_ie(U):
    """A = 1 - vol(union): singles w, pairs w^2 (theorem), triples Fourier, >=4-fold grid."""
    U = sorted(set(U)); k = len(U)
    tot = 1.0 - k*w + (k*(k-1)/2)*w*w
    for S in itertools.combinations(U, 3):
        tot -= vol3_fourier(S)
    for d in range(4, k+1):
        for S in itertools.combinations(U, d):
            tot += (-1)**d * dfold_grid(S, ngrid=420)
    return tot
def A_sampled(U, NC=1200):
    BAND = h
    def clip(arcs, v, c):
        out = []
        for (lo, hi) in arcs:
            k0 = int(floor(lo*v - c)) - 1; k1 = int(floor(hi*v - c)) + 1
            for kk in range(k0, k1+1):
                aa = (kk + c + BAND)/v; bb = (kk + 1 + c - BAND)/v
                l = lo if lo > aa else aa
                hh = hi if hi < bb else bb
                if l < hh - 1e-15: out.append((l, hh))
        return out
    tot = 0.0
    for ic in range(NC):
        c = (ic+0.5)/NC
        arcs = [(0.0, 1.0)]
        for u in sorted(U): arcs = clip(arcs, u, c)
        tot += sum(hi-lo for lo, hi in arcs)/NC
    return tot
for U in [(5,37), (1,2,3), (10,24,38), (1,2,4), (5,6,7,8), (1,2,3,4)]:
    Ae = A_ie(U); As = A_sampled(U)
    note = ""
    if len(U) == 2: note = f"  [36/49 = {36/49:.6f}]"
    if U in [(1,2,3),(10,24,38)]: note = f"  [61/98 = {61/98:.6f}]"
    print(f"  U={str(U):>12}: A_IE = {Ae:.6f}   A_sampled = {As:.6f}{note}")

print("\n" + "="*100)
print(" V4: BONFERRONI FLOORS (torus, exact terms) for the S30 argmin patterns, k = 4..13")
print("="*100)
argmins = { 4: (5,6,7,8), 5: (3,5,6,7,9), 6: (2,4,6,7,8,10), 7: (3,5,6,7,8,9,11),
            8: (1,3,5,6,7,8,9,11), 9: (1,3,4,5,6,7,8,9,11), 10: (1,3,4,5,6,7,8,9,11,13),
            11: (2,3,4,5,6,7,8,9,10,12,14), 12: (4,8,9,10,12,13,14,15,16,17,22,24),
            13: (3,7,8,9,11,13,15,16,17,19,23,24,27) }
MP = 14249/252252
print(f"  {'k':>3} {'Bonf-2':>8} {'Bonf-3':>8} {'sampled A':>10}  vs witnessMP={MP:.4f}")
for k in range(4, 14):
    U = argmins[k]
    b2 = 1 - k*w + (k*(k-1)/2)*w*w
    tri = sum(vol3_fourier(S) for S in itertools.combinations(U, 3))
    b3 = b2 - tri
    As = A_sampled(U, NC=600)
    print(f"  {k:>3} {b2:>8.4f} {b3:>8.4f} {As:>10.4f}   Bonf-3 {'>= MP' if b3 >= MP else '< MP'}"
          f"   (triple sum {tri:.4f}, {k*(k-1)*(k-2)//6} triples, {len(_tricache)} distinct relations cached)")
print("""
 READING: every inclusion-exclusion term is a subtorus box volume — exact, rational at h = 1/14, with
 NO arithmetic in singles and pairs, and the first arithmetic entering at triples through the primitive
 sum-zero relation (AP (1,-2,1) carries the max weight 2h^2). The (⋆)-census is symbolic. Bonferroni-3
 is a certified signed truncation when the 4-fold tail is controlled — for the floor at witnessMP the
 table shows the working margins.""")
print("DONE.")
