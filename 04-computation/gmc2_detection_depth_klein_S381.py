#!/usr/bin/env python3
"""
klein-2026-07-20-S381 -- IS THE GMC(2) DETECTION DEPTH BOUNDED BY THE SPAN ALONE, uniform in
radial degree?  (HYP-8540, the radial analog of THM-1710's toral 'depth = M+N'.)

Owner: "now work to prove the GMC(2)."

FRAME (opus/kp THM-1740, THM-1710, undisputed Wiener-Hopf half of the bridge): GMC(2) reduces
to NC2 = 'E[P^m]=0 for all m => P one-sided'. E[P^m] is P-recursive in m. For the TORAL nullcone
the detection depth is EXACTLY the span D=M+N (THM-1710): a_1..a_D = 0 forces all a_m=0. The
open uniform bound (HYP-8540) is the RADIAL analog: is the depth for E[P^m] bounded by the SPAN
alone, uniformly over the radial degree d?  If YES, GMC(2) holds for ALL bounded-span P (any d).

MEASURED HERE: for genuine (charge-radius-LOCKED) two-sided P on a fixed charge span, at radial
degrees d = 0,1,2,..., the DETECTION DEPTH = smallest D such that {E[P^1]=..=E[P^D]=0} already
forces one-sidedness (Rabinowitsch-empty on the first D moments).  Compared to the span.
"""
import sympy as sp
from math import factorial

Z, Zb = sp.symbols('Z Zb')
def Ewick(poly):
    poly = sp.expand(poly)
    if poly == 0: return sp.Integer(0)
    tot = sp.Integer(0)
    for (a, b), c in sp.Poly(poly, Z, Zb).as_dict().items():
        if a == b: tot += c * factorial(a)
    return tot

def build_P(charges, dmax):
    """genuine locked P: for each charge q in `charges`, and each radius h with h ≡ q (mod 2),
       |q| <= h <= |q| + 2*dmax, add coeff * Z^{(h+q)/2} Zb^{(h-q)/2}."""
    P = sp.Integer(0); coeffs = {}
    for q in charges:
        cs = []
        for k in range(dmax + 1):           # radial depth k -> radius |q| + 2k
            h = abs(q) + 2 * k
            a = (h + q) // 2; b = (h - q) // 2
            name = sp.symbols(f'g_{"m" if q<0 else ""}{abs(q)}_{k}')
            cs.append((q, name)); coeffs.setdefault(q, []).append(name)
            P += name * Z**a * Zb**b
    return P, coeffs

def detection_depth(charges, dmax, mmax=12):
    P, coeffs = build_P(charges, dmax)
    allc = [c for cl in coeffs.values() for c in cl]
    moms = []
    for m in range(1, mmax + 1):
        moms.append(sp.nsimplify(sp.expand(Ewick(sp.expand(P**m)))))
    # one-sided test: for every pos-charge coeff p and neg-charge coeff n, (p*n) nilpotent mod I;
    # and every charge-0 coeff nilpotent.  Detection depth = smallest D s.t. this holds with E[P^1..D].
    pos = [c for q in coeffs for c in coeffs[q] if q > 0]
    neg = [c for q in coeffs for c in coeffs[q] if q < 0]
    zer = [c for q in coeffs for c in coeffs[q] if q == 0]
    targets = [p*n for p in pos for n in neg] + zer
    from sympy import groebner
    for D in range(1, mmax + 1):
        gens = [g for g in moms[:D] if g != 0]
        if not gens: continue
        try:
            G = groebner(gens, *allc, order='grevlex')
        except Exception:
            continue
        ok = all(any(sp.simplify(G.reduce((t)**k)[1]) == 0 for k in range(1, 8)) for t in targets)
        if ok: return D, len(allc), len(targets)
    return None, len(allc), len(targets)

print("=" * 90)
print("DETECTION DEPTH vs (charge span, radial degree d)")
print("=" * 90)
print(f"{'charges':>18} {'span':>5} {'d':>3} {'#coeffs':>8} {'#targets':>9} {'DETECTION DEPTH':>16}")
tests = [
    ([-1,0,1], 0), ([-1,0,1], 1), ([-1,0,1], 2),
    ([-1,1], 0), ([-1,1], 1), ([-1,1], 2),
    ([-2,-1,1,2], 0), ([-2,-1,1,2], 1),
    ([-2,-1,0,1,2], 0),
    ([-3,3], 0), ([-3,3], 1),
]
for charges, d in tests:
    span = max(charges) - min(charges)
    D, nc, nt = detection_depth(charges, d)
    print(f"{str(charges):>18} {span:>5} {d:>3} {nc:>8} {nt:>9} {str(D):>16}")
print("""
 READING: if DETECTION DEPTH depends only on the SPAN (constant down each fixed-charges block
 as d grows), that is the uniform bound HYP-8540 -- and GMC(2) holds for all bounded-span P at
 EVERY radial degree.  If it GROWS with d, the radial layer genuinely needs more than the toral
 depth and HYP-8540 as stated is false (the bound must include d).
""")
