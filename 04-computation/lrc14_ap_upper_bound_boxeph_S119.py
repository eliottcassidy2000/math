#!/usr/bin/env python3
"""
boxeph-2026-07-18-S119 — proof of the S118 UPPER bound M(AP) <= (q-11)/(2q).

Completes HYP-7710: the centering witness is the GLOBAL maximizer, so
    M({a, a+d, ..., a+11d}) = (q-11)/(2q) = 1/2 - 11/(2q)   (d odd, q=2a+11d)  EXACTLY.

THE ARGUMENT (d odd, gcd(a,d)=1, q=2a+11d, mu=(q-11)/(2q)).
Suppose some t has ||v_k t|| > mu for all k (v_k=a+dk); derive a contradiction.
 - All y_k := v_k t mod 1 lie in the safe arc (mu, 1-mu), of length 1-2mu = 11/q.
 - beta := <d t>, alpha := a t (mod 1); y_k = alpha + k beta, so
   confinement of y_0,y_11 gives ||11 beta|| < 11/q.
 - [THRESHOLD q>132]  b:=|<beta>| <= ||beta|| < 11/q < 1/11, and the wrap zone
   (some k<=11 with k*b near 1) needs 11b >= 1 - 11/q, impossible once
   121/q < 1 - 11/q i.e. q>132.  So no wrap: 11b = ||11 beta|| < 11/q => b < 1/q.
 - ||a beta|| <= a*||beta|| = a*b < a/q  (subadditivity of ||.||).
 - COUPLING: a*beta == d*alpha (mod 1)  [a(dt) - d(at) = 0], so ||d alpha|| < a/q.
 - Hence alpha is within a/(dq) of some i/d; and alpha in (mu,1-mu) => |alpha-1/2| < 11/(2q).
 - So |i/d - 1/2| < a/(dq) + 11/(2q); but d ODD => |i/d-1/2| >= 1/(2d).
 - Multiply 1/(2d) < a/(dq)+11/(2q) by 2dq:  q < 2a+11d = q.  CONTRADICTION.
So M <= mu for q>132; the finitely many APs with q<=132 are verified here exhaustively.

This script:
  (A) exhaustively verifies M(AP) = (q-11)/(2q) for ALL primitive odd-d APs with q<=132
      (the finite residue of the proof), via exact maximization;
  (B) traces the confinement-coupling quantities (b, ||a beta||, ||d alpha||, the q<q
      collapse) at the witness and at perturbed times, to confirm the mechanism;
  (C) spot-checks a few LARGE-q spread APs (q>132) where the clean argument applies.
"""
from fractions import Fraction as F
from math import gcd


def frac_dist(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(speeds, t: F) -> F:
    return min(frac_dist(F(v) * t) for v in speeds)


def true_M_exact(speeds):
    """Exact M = max_t min_i ||v_i t||.  Candidate maximizers are t=m/D with
    D = v_i+v_j or |v_i-v_j| (the tent-pole denominators); this set provably
    contains every local max of the piecewise-linear min-of-triangle-waves."""
    Ds = set()
    n = len(speeds)
    for i in range(n):
        for j in range(i + 1, n):
            Ds.add(speeds[i] + speeds[j])
            Ds.add(abs(speeds[i] - speeds[j]))
    best = F(0)
    bestt = F(0)
    for D in Ds:
        if D == 0:
            continue
        for m in range(1, D):
            t = F(m, D)
            if t >= 1:
                continue
            val = fmin(speeds, t)
            if val > best:
                best, bestt = val, t
    return best, bestt


print("=" * 74, flush=True)
print("(A) EXHAUSTIVE finite-residue check: M(AP)=(q-11)/(2q) for all q<=132", flush=True)
print("    (primitive, d odd) — the cases the q>132 argument does not cover", flush=True)
print("=" * 74, flush=True)
bad = []
count = 0
maxq = 132
for d in range(1, maxq, 2):            # d odd
    a = 1
    while 2 * a + 11 * d <= maxq:
        if gcd(a, d) == 1:
            q = 2 * a + 11 * d
            speeds = [a + d * k for k in range(12)]
            M, tstar = true_M_exact(speeds)
            formula = F(q - 11, 2 * q)
            count += 1
            if M != formula:
                bad.append((a, d, q, M, formula, tstar))
        a += 1
print(f"  checked {count} primitive odd-d APs with q<=132", flush=True)
print(f"  mismatches (M != (q-11)/(2q)): {len(bad)}", flush=True)
for b in bad[:20]:
    print("   ", b, flush=True)
print(f"  ALL q<=132 verified: {len(bad) == 0}", flush=True)

print("", flush=True)
print("=" * 74, flush=True)
print("(B) MECHANISM TRACE: the confinement-coupling quantities", flush=True)
print("=" * 74, flush=True)
print(f"{'a':>3}{'d':>3}{'q':>5} | {'t':>9} {'min_k':>8} {'b=<beta>':>9} "
      f"{'||a b||':>8} {'||d al||':>9} {'q<q? (2a+11d vs bound)':>10}", flush=True)
for (a, d) in [(2, 41), (3, 17), (5, 19), (1, 23)]:
    q = 2 * a + 11 * d
    speeds = [a + d * k for k in range(12)]
    p = pow(d, -1, q)
    tw = F(p, q)                       # the witness
    # a few times: the witness and small perturbations
    cand = [tw, tw + F(1, 3 * q), tw - F(1, 5 * q), F(p + 1, q)]
    for t in cand:
        beta = frac_dist(F(d) * t)     # |<d t>|
        alpha = (F(a) * t) - (F(a) * t).numerator // (F(a) * t).denominator  # frac(a t)
        ab = frac_dist(F(a) * (F(d) * t))     # ||a*(d t)|| = ||a beta||
        da = frac_dist(F(d) * (F(a) * t))     # ||d*(a t)|| = ||d alpha||
        mn = fmin(speeds, t)
        note = "witness" if t == tw else ""
        print(f"{a:>3}{d:>3}{q:>5} | {str(t):>9} {str(mn):>8} {str(beta):>9} "
              f"{str(ab):>8} {str(da):>9}  {note}", flush=True)
    # confirm the coupling identity a*beta == d*alpha (mod 1) exactly at several t
    ok = all(frac_dist(F(a) * F(d) * t - F(d) * F(a) * t) == 0 for t in cand)  # trivially 0
    coup = all(frac_dist(F(a) * (F(d) * t)) == frac_dist(F(d) * (F(a) * t)) for t in cand)
    print(f"      coupling ||a*(dt)|| == ||d*(at)|| for all traced t: {coup}", flush=True)

print("", flush=True)
print("=" * 74, flush=True)
print("(C) LARGE-q spread APs (q>132, clean-argument regime): M = (q-11)/(2q)", flush=True)
print("=" * 74, flush=True)
for (a, d) in [(1, 31), (2, 41), (3, 37), (7, 45), (5, 51)]:
    if gcd(a, d) != 1 or d % 2 == 0:
        continue
    q = 2 * a + 11 * d
    speeds = [a + d * k for k in range(12)]
    M, tstar = true_M_exact(speeds)
    formula = F(q - 11, 2 * q)
    print(f"  a={a:>2} d={d:>2} q={q:>4}: M={str(M):>9}  (q-11)/(2q)={str(formula):>9}  "
          f"match={M == formula}  >1/13={M > F(1,13)}", flush=True)

print("\nDONE.", flush=True)
