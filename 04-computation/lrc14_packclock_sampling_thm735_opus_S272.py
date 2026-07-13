#!/usr/bin/env python3
"""
lrc14_packclock_sampling_thm735_opus_S272.py
============================================
opus-2026-07-13-S272.  THM-735: the PACK-CLOCK SAMPLING LEMMA (measure form of the detuned
dispatch).  Exact-rational verification of every claim.

THEOREM (THM-735).  Let v = c*R u D with c >= 2, R a finite set of nonzero integers,
D = {u_1..u_d}, and let G = G_R^(1/14) = {s : ||w s|| >= 1/14 for all w in R} (pack good set).
Then the 1/14-safe measure of v satisfies

    L(v)  >=  |G| * ( 1 - Sum_i (2*(1/14)*c + gcd(u_i,c)) / c )
          =   |G| * ( 1 - d/7 - Sum_i gcd(u_i,c)/c ),
and in the coprime case (gcd(u_i,c)=1 for all i), with the sharp integer count:
    L(v)  >=  |G| * ( 1 - d * (floor(c/7)+1) / c )        [ floor(c/7) -> c/7 when 7|c ].

PROOF.  Partition [0,1) into the c branches t = (j+s)/c, j = 0..c-1, s in [0,1).  The pack
condition ||c w t|| >= 1/14 reads ||w s|| >= 1/14: s in G, THE SAME for every branch (the pack
shares one clock).  For s fixed, branch j is u_i-dangerous iff u_i (j+s) in (-c/14, c/14) mod c.
The points {u_i j mod c} form a gcd(u_i,c)-spaced lattice on R/cZ; a g-spaced lattice meets an
open arc of length c/7 in at most floor((c/7)/g)+1 points, so at most g*(floor((c/7)/g)+1)
<= c/7 + g of the c branches are u_i-dangerous (coprime case: floor(c/7)+1 exactly; 7|c, g=1:
at most c/7 since an open arc of integer length m meets a unit lattice in at most m points).
Union bound over D and integrate over s in G.  QED.

Contrast THM-668 (monad-S3, the point-witness detuned dispatch g*H u {delta}, M >= 1/13):
 - THM-668 gives M >= 1/13 via LRC(13)-citation + branch pigeonhole; its Lipschitz fattening
   yields only L >~ 1/(1092 c) -- DECAYING in c.  THM-735 gives a UNIFORM measure floor
   L >= |G|/2 (c=2) rising to (6/7)|G| (7|c): the covering-route residue object (L) itself.
 - d >= 2 is OPEN in THM-668 (degenerate branch subgroup: detuned elements congruent mod g
   move diagonally).  THM-735's union bound needs NO joint condition on the u_i mod c:
   the d >= 2 case closes for all c with d*(floor(c/7)+1) < c (e.g. d=2: all c >= 3).
 - lambda-form: at threshold lambda, L_lambda >= |G^lambda|(1 - d(floor(2 lambda c)+1)/c);
   for the tower this is positive for every lambda < 1/13, recovering M >= 1/13 (THM-668's
   constant) in the limit; at lambda = 1/14 the pack input |G^(1/14)| > 0 is a FINITE exact
   computation (no LRC(13) citation needed).
 - SHARP at the binding case: L(2{1..12} u {13}) = |G_{1..12}|/2 EXACTLY (S270 exact fractions;
   equality forces N_s = 1 for a.e. s in G: on the pack clock one cannot be 13-safe on both
   branches -- the c=2 rigidity behind M = 1/13 exact, MISTAKE-141).

This script verifies (all Fractions, no floats in logic):
  A. the tower bound vs exact L for c = 2,4,6,8 (and the c=2 EQUALITY);
  B. the lambda-form floor at lambda = 5/66 < 1/13 (tower, c=2): still positive;
  C. d=2 coprime instances: 3{1..11} u {13,14} and 4{1..11} u {13,21} (both covering,
     primitive) -- bound vs exact L;
  D. a THM-668-DEGENERATE d=2 instance (both detuned CONGRUENT mod c -- the diagonal case
     open in THM-668): 5{1..11} u {13,18} (13 = 18 = 3 mod 5) -- bound vs exact L;
  E. a gcd>1 instance: 6{1..11} u {13,21} (gcd(21,6)=3): the general-gcd bound vs exact L.
"""
from fractions import Fraction as F

def normalize(ivs):
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1] = (out[-1][0], b)
        else: out.append((a, b))
    return out

def intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def measure(A): return sum(b - a for a, b in A)

def safe_set(w, lam):
    return [((k + lam) / w, (k + 1 - lam) / w) for k in range(w)]

def good_intervals(speeds, lam=F(1, 14)):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w, lam))
    return cur

def gcd(a, b):
    while b: a, b = b, a % b
    return a

def bound_thm735(G, c, us, lam=F(1, 14)):
    """the coprime-sharp / general-gcd THM-735 lower bound"""
    tot = F(0)
    for u in us:
        g = gcd(u, c)
        arc = 2 * lam * c          # window length on R/cZ
        cnt = g * ((arc / g).__floor__() + 1)   # g-spaced lattice in open arc
        if g == 1 and (arc % 1 == 0):           # integer-length open arc, unit lattice
            cnt = arc
        tot += F(cnt, c)
    return G * (1 - tot)

print("=" * 100)
print("THM-735 (pack-clock sampling) -- exact verification.  All quantities Fractions.")
print("=" * 100)

G12 = measure(good_intervals(range(1, 13)))
G11 = measure(good_intervals(range(1, 12)))
print("\npack good sets:  |G_{1..12}| = %s = %.6f ;  |G_{1..11}| = %s = %.6f"
      % (G12, float(G12), G11, float(G11)))

print("\nA. THE COMPRESSED TOWER  c*{1..12} u {13}  (d=1, coprime):")
print("   %3s %26s %12s %26s %12s %8s" % ("c", "bound (exact)", "bound", "true L (exact)", "true L", "OK"))
for c in [2, 4, 6, 8]:
    body = [c * w for w in range(1, 13)] + [13]
    Lex = measure(good_intervals(body))
    b = bound_thm735(G12, c, [13])
    ok = "EQUAL!" if Lex == b else ("OK" if Lex >= b else "**VIOLATED**")
    print("   %3d %26s %12.6f %26s %12.6f %8s" % (c, str(b), float(b), str(Lex), float(Lex), ok))

print("\nB. LAMBDA-FORM (tower c=2, lambda = 5/66 < 1/13 = 0.0769):")
lam = F(5, 66)
G12l = measure(good_intervals(range(1, 13), lam))
body2 = [2 * w for w in range(1, 13)] + [13]
Ll = measure(good_intervals(body2, lam))
bl = bound_thm735(G12l, 2, [13], lam)
print("   |G^lam_{1..12}| = %s = %.6f ; bound = %s = %.6f ; true L_lam = %s = %.6f ; OK=%s"
      % (G12l, float(G12l), bl, float(bl), Ll, float(Ll), Ll >= bl))
print("   (positivity for every lam < 1/13 gives M(tower) >= 1/13 = THM-668's constant, measure-style)")

print("\nC. d=2 COPRIME (open in THM-668 for small g):")
for c, us in [(3, [13, 14]), (4, [13, 21])]:
    body = [c * w for w in range(1, 12)] + us
    cov = all(any(w % q == 0 for w in body) for q in range(2, 15))
    Lex = measure(good_intervals(body))
    b = bound_thm735(G11, c, us)
    print("   %d*{1..11} u %-9s covering=%s  bound=%s=%.6f  true L=%.6f  OK=%s"
          % (c, us, cov, b, float(b), float(Lex), Lex >= b > 0))

print("\nD. d=2 DEGENERATE-DIAGONAL (13=18=3 mod 5 -- THM-668's open obstruction):")
c, us = 5, [13, 18]
body = [c * w for w in range(1, 12)] + us
Lex = measure(good_intervals(body))
b = bound_thm735(G11, c, us)
print("   5*{1..11} u {13,18}:  bound=%s=%.6f  true L=%.6f  OK=%s  (diagonal case CLOSED)"
      % (b, float(b), float(Lex), Lex >= b > 0))

print("\nE. gcd>1 detuned (6*{1..11} u {13,21}, gcd(21,6)=3):")
c, us = 6, [13, 21]
body = [c * w for w in range(1, 12)] + us
Lex = measure(good_intervals(body))
b = bound_thm735(G11, c, us)
print("   bound=%s=%.6f  true L=%.6f  OK=%s" % (b, float(b), float(Lex), Lex >= b > 0))

print("\ndone.")
