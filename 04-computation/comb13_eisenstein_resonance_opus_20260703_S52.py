#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE 13-SPACED COMB LEVER as the EISENSTEIN RESONANCE, and the rational/irrational + odd/even duality of
covering times (OPEN-Q-108, the wide-far/slow-base residual of CoveringFarLonely / the two-sided architecture).

opus-2026-07-03-S52. Owner: work the 13-comb lever creatively; think Vitali, rational/irrational & odd/even
dualities, and how TYPES OF NUMBERS relate to tournaments.

CLAIMS TO VERIFY (all exact):
 1. 183 = Phi_6(14) = 14^2 - 14 + 1 = 3*61.  14 is a PRIMITIVE 6th ROOT OF UNITY mod 183 (order 6):
    powers of 14 mod 183 = {14, 13, -1(=182), -14(=169), -13(=170), 1} = the 6 Eisenstein units.
    In particular 13 = 14^2 mod 183 and 14^3 = -1 mod 183 (so 14*13 = 14^3 = -1: the comb resonance).
 2. t* = 14/183.  ||13 t*|| = 1/183 (since 13*14 = 182 = -1 mod 183).  So a 13-SPACED comb
    {w, w+13, ..., w+13(r-1)} has phases at t* forming an AP with step 1/183, spanning (r-1)/183.
    For r <= 13 the span <= 12/183 = 0.0656 < 1/14 = 0.0714: the whole comb fits in ONE danger-arc width.
 3. The base {1..12} at t*: min_v ||v t*|| = 14/183 = 1/14 + 13/2562 (slack 13/2562 > 0).
 4. t* is NEAR the obvious 13-resonance 1/13 (=14/182): at 1/13 the comb COLLAPSES to one point (all
    equal mod 1), but a covering comb (13|w) then sits AT the origin (dangerous). t* shifts off 1/13 to
    lift the comb into the safe region while keeping the base safe -- the "core-slack absorbs the comb".
 5. THE CERTIFICATE IS RATIONAL AT DENOMINATOR 183: a covering family with a 13-spaced far comb is lonely
    at t* = 14/183 (denominator 183 = Phi_6(14)), which the small-q census (q<=~45) CANNOT see -- the
    rational/irrational duality: resonant families are certified by a SPECIFIC Eisenstein rational, not the
    free irrational sweep.  (Complement: klein HYP-3791 -- NON-13 combs equidistribute = the irrational side.)
"""
import sys
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):  # ||x|| = distance to nearest integer, x a Fraction
    f = x - (x.numerator // x.denominator)   # fractional part in [0,1)
    return min(f, 1 - f)

print("="*98); print(" (1) THE EISENSTEIN STRUCTURE OF 183 = Phi_6(14)"); print("="*98)
N = 183
print(f"  Phi_6(14) = 14^2-14+1 = {14*14-14+1};  factor 183 = 3 * 61 (both primes = 3 mod 4 ... 3, 61=1 mod4)")
pw = []; x = 1
for j in range(1, 13):
    x = (x * 14) % N; pw.append(x)
    if x == 1: order = j; break
signed = [p if p <= N//2 else p - N for p in pw]
print(f"  powers of 14 mod 183: {pw}")
print(f"  signed:               {signed}   <- the 6 Eisenstein units {{+-1,+-13,+-14}}; ORDER of 14 = {order}")
print(f"  => 14^2 = {pow(14,2,N)} (=13), 14^3 = {pow(14,3,N)} (=182=-1), so 14*13 = 14^3 = -1 mod 183.")
print(f"  ODD/EVEN: order 6 = 2 * 3 = (antipode Z/2: 14^3=-1) x (Eisenstein Z/3: 14^2=13). The '2' is the")
print(f"           complement/antipode (14 even =2*7), the '3' is the Eisenstein 6th-root (3 | 183).")

print("\n"+"="*98); print(" (2) THE COMB RESONANCE at t* = 14/183"); print("="*98)
tstar = Fr(14, 183)
print(f"  ||13 t*|| = ||{13*14}/183|| = ||182/183|| = {norm(Fr(13*14,183))}  (= 1/183: the tiny step)")
print(f"  {'r':>3} {'comb span (r-1)/183':>20} {'< 1/14?':>8}")
for r in [2,3,5,8,13]:
    span = Fr(r-1,183); print(f"  {r:>3} {str(span):>20} {'YES' if span < Fr(1,14) else 'no':>8}  ({float(span):.4f} vs 0.0714)")
# comb phases at t* for a concrete far w (covering: 13|w)
for w in [169, 13*20, 13*44]:
    ph = sorted(set(norm(Fr((w+13*k)*14, 183)) for k in range(5)))
    print(f"    w={w} (=13*{w//13}): comb phases ||.|| at t* (k=0..4) = {[str(p) for p in ph]}")

print("\n"+"="*98); print(" (3) THE BASE {1..12} SLACK AT t*, and (4) the 1/13 collapse"); print("="*98)
base12 = list(range(1,13))
mins = min(norm(Fr(v*14,183)) for v in base12)
print(f"  min_(v=1..12) ||v*14/183|| = {mins} = {float(mins):.5f};  1/14 = {float(Fr(1,14)):.5f};  slack = {mins - Fr(1,14)} = 13/2562? {mins-Fr(1,14) == Fr(13,2562)}")
base13 = list(range(1,14))
print(f"  min_(v=1..13) ||v*14/183|| = {min(norm(Fr(v*14,183)) for v in base13)}  (the covering-min value M=14/183 is attained)")
# the 1/13 collapse
print(f"  at t=1/13: 13-spaced comb collapses -- ||(w+13k)/13|| = ||w/13|| for all k (13k/13 in Z).")
print(f"     covering comb has 13|w => ||w/13||=0 = AT THE ORIGIN (dangerous). t*=14/183 shifts off 1/13")
print(f"     (1/13={float(Fr(1,13)):.5f}, t*={float(tstar):.5f}, gap {float(Fr(1,13)-tstar):.5f}) to lift the comb into safe.")

print("\n"+"="*98); print(" (5) THE RATIONAL CERTIFICATE AT DENOMINATOR 183 (what the small-q census misses)"); print("="*98)
def M_and_witness(S, Dmax):
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w)
            if v!=w: dens.add(abs(v-w))
    best=Fr(0); bt=None
    for b in sorted(d for d in dens if 0<d<=Dmax):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            md=min(norm(Fr(a*v,b)) for v in S)
            if md>best: best=md; bt=Fr(a,b)
    return best, bt
# COVERING 13-family with a 13-spaced far TRIPLE {156,169,182}=13*{12,13,14}:
#   156=2^2*3*13 (covers 2,3,4,6,12,13), 169=13^2, 182=2*7*13 (covers 2,7,13,14). Triple covers {2,3,4,6,7,12,13,14}.
#   Missing {5,8,9,10,11} => base must cover them. base {1,2,3,5,6,7,8,9,10,11} (10) covers 5(5,10),8,9,10,11 + rest.
def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2,15))
Sgood = sorted(set(range(1,13)) | {182})   # THE DEEP WELL {1..12, 182}, 182=13*14 (kps-S23), 13 runners
print(f"  S = {Sgood}  (THE DEEP WELL; 182 = 13*14 the resonant far runner); covering? {covering(Sgood)}, |S|={len(Sgood)}")
M, wt = M_and_witness(Sgood, 2*max(Sgood)+2)
print(f"    M(S) = {M} = {float(M):.5f} at t = {wt} (DENOMINATOR {wt.denominator if wt else '-'}); lonely(>=1/14)? {M >= Fr(1,14)}")
print(f"    witness denom {wt.denominator} vs small-q census cutoff ~45: {'>>45 => MISSED by the census (needs the Eisenstein certificate)' if wt and wt.denominator>45 else 'within census'}")
# show the phases of S at the witness
ph = sorted((norm(Fr(wt.numerator*v, wt.denominator)), v) for v in Sgood)
print(f"    binding runners at the witness (smallest ||v t||): {[(str(p),v) for p,v in ph[:4]]}")
# and the census-scale value: best small-q (q<=45) witness -- to show the gap
Msmall, wsmall = M_and_witness(Sgood, 45)
print(f"    best q<=45 witness: min-dist = {float(Msmall):.5f} at {wsmall} (q={wsmall.denominator if wsmall else '-'}) -- {'LONELY already' if Msmall>=Fr(1,14) else 'NOT lonely at small q => the census fails, the Eisenstein 183 certificate is needed'}")

print("\n"+"="*98); print(" SYNTHESIS: types of numbers <-> tournaments <-> covering times"); print("="*98)
print("  RATIONAL t=a/q (roots of unity)  <-> TRANSITIVE tournament (Cayley->roots of unity) <-> RESONANT/pinned")
print("     census; the 13-comb resonates here at the EISENSTEIN t*=14/183=14/Phi_6(14).")
print("  QUADRATIC IRRATIONAL (Gauss sum i*sqrt p) <-> PALEY tournament (skew spec +-i sqrt p) <-> mixing.")
print("  GENERIC IRRATIONAL (equidistributed)      <-> regular/mixing tournament <-> FREE sweep (THM-608).")
print("  The lever: the hardest configs sit at the RESONANT rational (Eisenstein 183); certified either")
print("  by that exact rational (denominator 183, missed by small-q) OR by a sweep off 1/13 (the two-sided")
print("  architecture). ODD/EVEN: 14=2*7, order-6 = 2(antipode)*3(Eisenstein); RATIONAL/IRRATIONAL = census/sweep.")
print("DONE.")
