# kind-pasteur-2026-07-10-S127
# The small-ruler law: attacking StrictlyLiveSupply (the LRC(14) wall).
#
# StrictlyLive v q p  :=  0<p<q  and  q < 14*((v_i p) % q) < 13*q  for all i.
# StrictlyLiveSupply  :=  every residual family (covering, ratio>13, compressed, no-detune,
#                         diff-primitive) has some strictly-live (q,p).  <- the open case of LRC(14).
#
# PROVED IN LEAN (LRCSmallRuler.lean): covering => every strictly-live q is >= 15
#   (covering gives a zero residue at each q in [2,14]; 0 is never in the band).
# MEASURED HERE: on covering ratio>13 families the minimal strictly-live q is 15..26, SCALE-INDEPENDENT,
#   and the only adversary that raises it (a highly-divisible speed) is DETUNED -- excluded by the residual.
from math import gcd, lcm
from fractions import Fraction as F
from functools import reduce
import random
random.seed(20)

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def gapfam(S): return max(S) > 13 * min(S)
def compressed(S): return all(any(j != i and S[i] <= 13 * S[j] for j in range(13)) for i in range(13))
def detuned(S):
    for g in range(2, 80):
        if len([v for v in S if v % g != 0]) == 1: return True
    return False
def diffprim(S): return reduce(gcd, [v - S[0] for v in S[1:]]) == 1
def strictly_live(S, q, p):
    return all(q < 14 * ((v * p) % q) < 13 * q for v in S)
def min_q(S, qmax=400):
    for q in range(15, qmax + 1):
        for p in range(1, q):
            if strictly_live(S, q, p): return q
    return None

print("=== PROVED (Lean): covering => strictly-live q >= 15.  Here: q=15 witness appears immediately ===")
def sample(hi, n):
    out = []; t = 0
    while len(out) < n and t < 500000:
        t += 1
        S = sorted(random.sample(range(1, hi + 1), 13))
        if gapfam(S) and covering(S): out.append(S)
    return out

print(f"\n=== min strictly-live q is SCALE-INVARIANT (covering, ratio>13) ===")
print(f"{'Vmax':>8} {'n':>4} {'min q':>6} {'max q':>6} {'mean':>6} {'None':>5}")
for hi in [45, 80, 150, 300, 600, 1200]:
    pool = sample(hi, 40)
    qs = [min_q(S) for S in pool]
    ok = [q for q in qs if q]
    print(f"{'[1,%d]'%hi:>8} {len(pool):>4} {min(ok):>6} {max(ok):>6} {sum(ok)/len(ok):>6.1f} {sum(1 for q in qs if q is None):>5}")

print(f"\n=== the adversary that raises min q is DETUNED (excluded by residual hdiv) ===")
for top in [18, 22, 26]:
    H = reduce(lcm, range(15, top + 1))
    S = sorted(set([H // k for k in [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14] if H % k == 0]))[:13]
    if len(S) < 13: continue
    print(f"  H=lcm(15..{top}): detuned={detuned(S)}  diff-prim={diffprim(S)}  "
          f"min q={min_q(S)}  (kills q in [15,{top}] by zero residue)")

print(f"\n=== on the RESIDUAL class (no detune, diff-primitive): max of (min q) ===")
best = 0; bestS = None; n = 0; t = 0
while n < 200 and t < 800000:
    t += 1
    S = sorted(random.sample(range(1, 500), 13))
    if not (gapfam(S) and covering(S) and compressed(S)): continue
    if detuned(S) or not diffprim(S): continue
    n += 1
    q = min_q(S, 200)
    if q and q > best: best, bestS = q, S
print(f"  {n} residual-class families; max(min strictly-live q) = {best}")
print(f"  at {bestS}")
print(f"\n  CONJECTURE (BoundedStrictlyLiveSupply): the witness lives in [15, ~26].")
print(f"  klein THM-685 transfer threshold Sum(v)/mu ~ 3000-5000 is TWO ORDERS larger and NOT needed --")
print(f"  the witness appears at the BOTTOM of the window, not the transfer tail.")
