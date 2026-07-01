#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S68
=======================
A DIVERSE METRICS TOOLBOX for the LRC(2p) family (apex prime p): beyond the genus, which niche
invariants fire at the p=7 frontier, which detect other structure (p mod 4 / QR), which are scale?

Organized by WHAT each metric detects.  All formulas noted.  (Extends S65/S66/S67.)
"""
from math import gcd, cos, pi, isqrt
from fractions import Fraction as F

def phi(n): return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)
def primes(N): return [p for p in range(2, N+1) if N % p == 0 and all(p % d for d in range(2, p))]
def psi(N):
    r = N
    for p in set(primes(N)): r = r*(p+1)//p
    return r
def nu2(N):
    if N % 4 == 0: return 0
    r = 1
    for p in set(primes(N)): r *= 1 if p == 2 else (1+(1 if p % 4 == 1 else -1))
    return r
def nu3(N):
    if N % 9 == 0: return 0
    r = 1
    for p in set(primes(N)): r *= 1 if p == 3 else (1+(1 if p % 3 == 1 else -1))
    return r
def cusps(N): return sum(phi(gcd(d, N//d)) for d in range(1, N+1) if N % d == 0)
def genus(N): return 1+F(psi(N), 12)-F(nu2(N), 4)-F(nu3(N), 3)-F(cusps(N), 2)
def hclass(D):  # class number of discriminant -D (D>0) via reduced primitive forms
    if D <= 0 or D % 4 in (1, 2): return 0
    h = 0
    for a in range(1, isqrt(D//3)+1):
        for b in range(-a+1, a+1):
            if (b*b+D) % (4*a): continue
            c = (b*b+D)//(4*a)
            if c < a or gcd(gcd(abs(a), abs(b)), c) != 1: continue
            if b < 0 and (a == c or -b == a): continue
            h += 1
    return h
def tri_finite(p):
    d = F(1, 2)+F(1, 3)+F(1, p)-1
    return int(4/d) if d > 0 else None

print("="*94)
print("LRC(2p) INVARIANT TOOLBOX -- p=3,5,7,11,13,17,19,23")
print("="*94)
print(f"{'p':>2}|{'genus':>5} {'triGrp':>7} {'gon':>3}|{'nu2':>3} {'h(-p)':>5}|{'psi':>3} {'covol/pi':>8} {'phi':>3}|{'apexgap':>8} {'2p-1':>4}")
for p in [3, 5, 7, 11, 13, 17, 19, 23]:
    N = 2*p; g = genus(N); tf = tri_finite(p)
    gon = 1 if g == 0 else (2 if g <= 2 else "?")   # genus0=rational(1), genus1,2=hyperelliptic(2)
    q = 2*p-1; qp = q > 1 and all(q % d for d in range(2, isqrt(q)+1))
    qtag = ("Paley-graph" if q % 4 == 1 else "tournament") if qp else "composite"
    print(f"{p:>2}|{str(g):>5} {str(tf) if tf else 'INF':>7} {str(gon):>3}|{nu2(N):>3} {hclass(p):>5}|"
          f"{psi(N):>3} {str(F(psi(N),3)):>8} {phi(N):>3}|{2-2*cos(pi/p):.5f} {q:>4}={qtag}")

print()
print("="*94)
print("SORTED BY WHAT EACH METRIC DETECTS (the diversity):")
print("="*94)
print("""  A. FRONTIER-DETECTORS (spherical -> hyperbolic at p=7, = genus jump 0->1):
       genus g(X0(2p)) = 1+psi/12-nu2/4-nu3/3-cusps/2  = dim S2^new = # cusp-form obstructions (f14)
       (2,3,p) triangle group order 4/(1/2+1/3+1/p-1): FINITE p<=5 (24,120), INFINITE p>=7
       gonality of X0(2p): 1 (P^1) for p<=5, 2 (elliptic/hyperelliptic) for p>=7
       HURWITZ THRESHOLD: (2,3,7) is the minimal hyperbolic triangle group; Klein quartic genus 3,
         Aut=PSL(2,7)=168=84*(3-1) attains the Hurwitz bound 84(g-1). p=7 = the Hurwitz prime.
    B. QR / p-mod-4 DETECTORS (a DIFFERENT axis -- the Paley/tournament dichotomy on 2p-1):
       nu2(2p)=1+(-1/p): 2 if p=1mod4, 0 if p=3mod4  (p=7 -> 0)
       h(-p): nonzero (class number of Q(sqrt-p)) for p=3mod4; the QR structure
       2p-1: p=7 -> 13 = 1mod4 -> PALEY GRAPH (Ramanujan, the covering-min circulant); else tournament
    C. SCALE (monotone, no transition): psi, phi(2p), J2, hyperbolic covolume (=psi*pi/3), apex gap ~1/p^2
    D. THE UNIQUE-ELLIPTIC METRIC at p=7:  X0(14) has genus EXACTLY 1 => it is an ELLIPTIC CURVE (14a),
       the ONLY apex prime where X0(2p) is elliptic (p<=5 rational, p>=11 genus>=2). And:
         14a: conductor 14, RANK 0, TORSION = Z/6, |torsion| = 6 = phi(14) = #(Z/14)^* = the lonely set!
       So the modular curve's torsion order recovers the triple-6 (S66).""")
print("="*94)
print("NEW FRAMING: genus = dim(un-regularizable residual) (S67). The 'metric like genus' IS the genus,")
print("read as the RESIDUAL DIMENSION: 0 (fully zeta-regularizable, easy) for p=3,5; 1 (one cusp form,")
print("first hard) at p=7. The frontier-detectors (triGrp/gonality/Hurwitz) are its geometric shadows;")
print("the QR-detectors (nu2/class number/Paley) are an ORTHOGONAL axis (the covering-min circulant type).")
print("\nDONE.")
