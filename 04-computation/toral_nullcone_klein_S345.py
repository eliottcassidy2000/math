#!/usr/bin/env python3
"""
klein-2026-07-20-S345 -- THE TORAL NULLCONE: an exact Lagrange-inversion proof that settles
the extreme-weight case, and the reduction of NC2 to it.

EXTENDS (does not duplicate):
  mac-mini THM-1500  -- master theorem + Lagrange uniqueness of f = f(0)/(1+s).  CEDED.
  boxeph HYP-8320    -- GMC(2) radial-family descent blocked (bounded search).
  death-star HYP-8330-- GMC(2) evidence TRUE (532 kernel, 0 counterexamples).
  klein THM-1510     -- EMP + NC2 + the two-weight theorem.  THIS FILE CONTINUES THAT LINE.

SETUP (n = 2).  One complex Gaussian Z, E[Z^a Zb^b] = delta_ab a!.  Polar: Z = sqrt(r) e^{i th}
with r ~ Exp(1) INDEPENDENT of th ~ Unif.  Any P in C[Z,Zb] is a function P(r,th), and
E[P^m] = E_{r,th}[P^m] -- the th-average is exactly the weight-0 projection.

Writing u = e^{i th} and P(r,u) = sum_q g_q(r) u^q, and letting D = max_q deg g_q with
LEADING SYMBOL  Lam(u) = sum_{deg g_q = D} c_q u^q  (a Laurent polynomial), the large-m
behaviour is governed by

        E[P^m]  ~  Gamma(Dm+1) * CT_u( Lam(u)^m ),

so NC2 forces the TORAL NULLCONE condition CT(Lam^m) = 0 for all m >= 1.

TORAL NULLCONE CONJECTURE (TNC).  For a Laurent polynomial Lam, CT(Lam^m) = 0 for every
m >= 1  iff  all exponents of Lam have the same STRICT sign.
"""
from fractions import Fraction as Fr
import itertools, cmath, math

def ct_pow(lam, m):
    """lam: dict exponent -> coeff.  Returns CT(lam^m) exactly."""
    cur = {0: 1}
    for _ in range(m):
        nxt = {}
        for e1, c1 in cur.items():
            for e2, c2 in lam.items():
                nxt[e1 + e2] = nxt.get(e1 + e2, 0) + c1 * c2
        cur = {e: c for e, c in nxt.items() if c != 0}
    return cur.get(0, 0)

print("=" * 80)
print("PART 1 -- TNC IS PROVED WHEN THE MINIMUM WEIGHT IS -1 (or the maximum is +1)")
print("=" * 80)
print("""
 Let Lam have minimum exponent exactly -1: Lam(u) = u^{-1} R(u) with R a polynomial,
 R(0) != 0.  Then Lam^m = u^{-m} R^m, so

        CT(Lam^m) = [u^m] R(u)^m.

 LAGRANGE-BUERMANN.  For the equation w = t R(w) (which has a unique small root w(t),
 analytic, with w'(0) = R(0) != 0):

        sum_{m>=0} t^m [u^m] R(u)^m  =  1 / (1 - t R'(w(t))).

 Demanding CT(Lam^m) = 0 for every m >= 1 makes the left side identically 1, so

        t R'(w(t)) = 0  for all small t   =>   R'(w(t)) = 0 for t != 0.

 But w is NONCONSTANT (w'(0) = R(0) != 0), so w(t) sweeps a neighbourhood of 0 and R'
 vanishes on an open set.  A polynomial vanishing on an open set is 0, so R' = 0, i.e.
 R is a CONSTANT and Lam = c u^{-1}: a single weight.  []

 By u -> 1/u the same argument settles maximum exponent exactly +1.

 NOTE THE CONTRAST with min exponent -M, M >= 2: there Lam = u^{-M} R and CT(Lam^m) =
 [u^{Mm}] R^m, and the same Lagrange computation (branching over the M-th roots of unity,
 with H = R^{1/M}) gives the condition [u^m] H(u)^m = 0 only for m divisible by M -- a
 STRICTLY WEAKER demand.  So M >= 2 is genuinely a different problem, not a harder version
 of the same one.  Tested below rather than assumed.
""")
# machine check of the M = 1 proof: no counterexample should exist
print(" machine check (M=1, i.e. min exponent -1), exhaustive over small R:")
bad = []
for deg in (1, 2, 3):
    for coeffs in itertools.product([0, 1, -1, 2, -2, 1j, -1j], repeat=deg + 1):
        if coeffs[0] == 0 or coeffs[deg] == 0: continue      # R(0) != 0, deg exact
        lam = {j - 1: c for j, c in enumerate(coeffs) if c != 0}
        if len(lam) == 1: continue                            # single weight: allowed
        if all(ct_pow(lam, m) == 0 for m in range(1, 9)):
            bad.append(lam)
print(f"   Laurent polys with min exponent -1, >1 term, CT(Lam^m)=0 for m=1..8: {len(bad)}"
      f"   {'<-- would REFUTE the proof' if bad else '(none, as proved)'}")

print("\n" + "=" * 80)
print("PART 2 -- THE OPEN CASE: min exponent <= -2 AND max exponent >= 2")
print("=" * 80)
found = []
tested = 0
for M in (2, 3):
    for N in (2, 3):
        for supp in itertools.combinations(range(-M, N + 1), 3):
            if min(supp) != -M or max(supp) != N: continue
            for cs in itertools.product([1, -1, 2, -2, 1j, -1j], repeat=3):
                lam = {e: c for e, c in zip(supp, cs)}
                tested += 1
                if all(ct_pow(lam, m) == 0 for m in range(1, 11)):
                    found.append(lam)
print(f"   3-term Laurent polys with min <= -2 and max >= 2, tested: {tested}")
print(f"   with CT(Lam^m) = 0 for m = 1..10: {len(found)}")
if found:
    print(f"   EXAMPLES (would be TNC counterexamples): {found[:6]}")
else:
    print("   NONE -- consistent with TNC in the hard range too (bounded, not a proof).")

print("\n" + "=" * 80)
print("PART 3 -- WHY THE M>=2 HITS ARE (OR ARE NOT) REAL: the gcd escape")
print("=" * 80)
print(" If every exponent of Lam is divisible by g > 1 then Lam(u) = Mu(u^g) and")
print(" CT(Lam^m) = CT(Mu^m), so nothing new -- such hits are RESCALINGS, not counterexamples.")
if found:
    real = [l for l in found if math.gcd(*[abs(e) for e in l if e != 0]) == 1]
    print(f"   of the {len(found)} hits, those with gcd(exponents) = 1: {len(real)}")
    if real: print(f"   GENUINE candidates: {real[:6]}")
    else: print("   ALL hits are u -> u^g rescalings.  TNC survives.")
else:
    print("   (no hits to classify)")

print("\n" + "=" * 80)
print("PART 4 -- CONSEQUENCE FOR GMC(2): the {-1,0,1} stratum's leading symbol is SETTLED")
print("=" * 80)
print("""
 The {-1,0,1} stratum -- the shape of BOTH n>=3 witnesses, and the case THM-1510 left open --
 has leading symbol with min exponent -1 and max exponent +1.  Part 1 applies at BOTH ends,
 so its leading symbol must be a single weight.  Concretely, for that stratum the exact
 ordinary generating function is (derived by residues, and checked below)

        sum_m L_m(alpha,beta) t^m  =  ((1 - beta t)^2 - 4 alpha t^2)^{-1/2},
        alpha = r a(r) c(r),  beta = b(r),

 so   sum_m E[P^m] t^m = E_r[ ((1-beta t)^2 - 4 alpha t^2)^{-1/2} ],  and setting it to 1
 with CONSTANT alpha, beta forces (1-beta t)^2 - 4 alpha t^2 = 1, i.e. beta = 0 and
 alpha = 0.  That is TNC at M = N = 1, by hand.
""")
def Lm(alpha, beta, m):
    tot = 0
    for k in range(m // 2 + 1):
        tot += (math.factorial(m) // (math.factorial(k) ** 2 * math.factorial(m - 2 * k))) \
               * alpha ** k * beta ** (m - 2 * k)
    return tot
print(" check of the generating function  sum L_m t^m = ((1-bt)^2-4at^2)^{-1/2}:")
for (a, b) in ((2, 3), (1, 0), (0, 5), (-1, 2)):
    ser = [Lm(a, b, m) for m in range(0, 7)]
    # expand ((1-bt)^2-4at^2)^{-1/2} as a power series
    import numpy as np
    N = 7
    q = np.zeros(N); q[0] = 1.0
    poly = np.zeros(N); poly[0] = 1 - 0.0
    base = np.zeros(N); base[0] = 1.0; base[1] = -2 * b; 
    if N > 2: base[2] = b * b - 4 * a
    # (base)^{-1/2} by Newton iteration on power series
    out = np.zeros(N); out[0] = 1.0
    for n in range(1, N):
        s = 0.0
        for k in range(1, n + 1):
            s += (k * (-0.5) - (n - k)) * base[k] * out[n - k]
        out[n] = s / n
    print(f"   alpha={a:>3} beta={b:>3}:  L_m = {ser}")
    print(f"                  series = {[round(float(x),6) for x in out]}   match: "
          f"{all(abs(float(out[m])-ser[m])<1e-9 for m in range(N))}")
