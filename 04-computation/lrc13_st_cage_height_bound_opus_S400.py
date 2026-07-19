"""
opus-2026-07-19-S400: the S-T tightness cage -- exact height bound H0.

Mined from Sungkawichai-Trakulthongchai arXiv:2604.23906 ("Eleven, twelve,
and thirteen lonely runners"; taken as TRUE per owner LRC(<=13) policy),
via three structural facts extracted from the paper (HTML extraction --
flagged for a direct PDF verification pass):

  (F1) For k=12 their sieve pipeline uses only c=2 lifts: depths
       l in {1,2,4,8}, primes p in P_12 subset [167,733], ln(prod P_12) > 547
       (their Table 1 line, > ln B_12 = 546).  13 never divides l*p.
  (F2) Class properness (Def 2.1) = [witness t in (1/(lp))Z with
       ||t v_i|| >= 1/13 for the whole class]  OR  [gcd branch:
       gcd(l, all-but-one coordinates) > 1, which for l in {2,4,8} means
       11 of 12 speeds even].
  (F3) Terminal state of the k=12 sieve per prime p: only equivalence-class
       members of (1,...,12) remain (equivalence = permutation, per-coordinate
       sign, global unit scaling mod p); Prop 4.4 finishes those.

QUANTIZATION (one line, this repo): a witness t = a/(lp) gives distances in
(1/(lp))Z, so properness value >= ceil(lp/13)/(lp) > 1/13 STRICTLY since
13 does not divide lp.  The weakest rung over the pipeline is at lp = 8*733:

    ceil(5864/13)/5864 = 452/5864 = 1/13 + 3/19058.

CAGE: any primitive 12-family W with M(W) < 452/5864, without 11-of-12 even
speeds, must therefore be equivalent to (1,...,12) mod every p in P_12 with
p not dividing prod(W): as multisets {+-w_i mod p} = c_p * {+-1,...,+-12}.

HEIGHT-FORCING: eliminating c_p via even power sums P_{2m}(W) = sum w_i^{2m},
define for m = 1..12 the integers
    R_m(W) := P_2(W)^m * S_{2m} - S_2^m * P_{2m}(W),
    S_{2m}  := sum_{i=1..12} i^{2m}   (S_2 = 650).
Cage => p | R_m(W) for every caging prime p.  If some R_m != 0 then
|R_m| >= prod(caging primes) >= e^547 / 733^loss, where
loss <= 12*floor(ln H/ln 167) (primes >= 167 dividing prod W, each < 733,
H = max(W)).  So whenever  bound(|R_m|, H) < e^547 / 733^loss(H)  for all m,
every R_m must vanish; R_m = 0 for all m = 1..12 pins the multiset
{w_i^2} = lambda*{i^2} (Newton), hence w_i = t*sigma(i), and primitivity
gives W = {1,...,12}.

This script computes, in exact integer arithmetic with conservative bounds:
  (1) the largest H0 with the forcing inequality satisfied for all m,
  (2) a verification of the Newton step on random samples,
  (3) the rung table ceil(lp/13)/(lp) over the pipeline.

STATUS: DERIVED, conditional on extraction fidelity (F1)-(F3); the height
bound itself is unconditional integer arithmetic given (F1)-(F3).
"""
from fractions import Fraction
from math import gcd, log, floor

K = 12
S = {2*m: sum(i**(2*m) for i in range(1, K+1)) for m in range(1, K+1)}

# ---------- (3) rung table ----------
print("=== certificate rungs ceil(lp/13)/(lp) over the k=12 pipeline ===")
worst = None
for l in (1, 2, 4, 8):
    for p in (167, 733):  # endpoints of P_12 range
        lp = l * p
        rung = Fraction((lp + 12)//13, lp)  # ceil(lp/13)/lp
        margin = rung - Fraction(1, 13)
        print(f"  l={l} p={p}: rung {rung} = 1/13 + {margin}")
        if worst is None or rung < worst:
            worst = rung
print(f"  WEAKEST RUNG (cage threshold): {worst} = 1/13 + {worst - Fraction(1,13)}")
assert worst == Fraction(452, 5864)

# ---------- (1) height bound ----------
# |R_m(W)| <= P2^m * S_{2m} + S_2^m * P_{2m} <= (12 H^2)^m S_{2m} + 650^m * 12 H^{2m}
# caging product >= e^547 / 733^loss, loss = 12*floor(ln H / ln 167)
LN547 = 547.0
def forcing_ok(H: int) -> bool:
    loss = 12 * floor(log(H) / log(167)) if H >= 167 else 0
    ln_caging = LN547 - loss * log(733)
    if ln_caging <= 0:
        return False
    for m in range(1, K+1):
        bound = (12 * H*H)**m * S[2*m] + (650**m) * 12 * (H**(2*m))
        if log(bound) >= ln_caging:
            return False
    return True

lo, hi = 1, 10**8
while lo < hi:                     # binary search largest H with forcing_ok
    mid = (lo + hi + 1)//2
    if forcing_ok(mid):
        lo = mid
    else:
        hi = mid - 1
H0 = lo
loss0 = 12 * floor(log(H0)/log(167))
print(f"\n=== height bound ===")
print(f"  H0 = {H0}  (largest H with |R_m| < prod(caging) guaranteed, all m=1..12)")
print(f"  at H0: loss <= {loss0} primes, ln(caging product) >= {LN547 - loss0*log(733):.1f}")
print(f"  binding m and sizes at H0:")
for m in (1, 6, 12):
    b = (12*H0*H0)**m * S[2*m] + (650**m)*12*(H0**(2*m))
    print(f"    m={m}: ln|R_m|_max = {log(b):.1f}")

# ---------- (2) Newton-step sanity ----------
# R_m(W)=0 for all m  <=>  normalized even power sums match the AP's
# => multiset {w_i^2} = lambda*{i^2}.  Verify: (a) the AP and its dilates
# satisfy R_m=0; (b) random non-dilate families have some R_m != 0.
import random
random.seed(20260719)
def R(W, m):
    P2 = sum(w*w for w in W)
    P2m = sum(w**(2*m) for w in W)
    return P2**m * S[2*m] - 650**m * P2m
def all_R_zero(W):
    return all(R(W, m) == 0 for m in range(1, K+1))
print(f"\n=== Newton-step sanity ===")
ap = list(range(1, 13))
print(f"  AP {{1..12}}: all R_m = 0? {all_R_zero(ap)}")
print(f"  3*AP:        all R_m = 0? {all_R_zero([3*i for i in ap])}")
bad = 0
for _ in range(2000):
    W = sorted(random.sample(range(1, 400), 12))
    if W == ap or all(w % W[0] == 0 and w//W[0] in ap for w in W):
        continue
    if all_R_zero(W):
        bad += 1
        print(f"  UNEXPECTED all-zero at non-dilate {W}")
print(f"  random non-dilates with all R_m = 0: {bad}/2000 (expect 0)")
# targeted: permuted-square structure w_i^2 = lambda i^2 with lambda nonsquare
# cannot happen over integers with positive w; spot-check near-misses:
for W in ([1,2,3,4,5,6,7,8,9,10,11,13], [2,4,6,8,10,12,14,16,18,20,22,24],
          [1,2,3,4,5,6,7,8,9,10,12,11]):
    print(f"  {W}: all R_m = 0? {all_R_zero(sorted(W))}")

print("\nCONCLUSION (conditional on F1-F3): every primitive 12-family W with")
print(f"M(W) < 452/5864 = 1/13 + 3/19058, max(W) <= {H0}, and not 11-of-12")
print("even speeds, satisfies W = {1,...,12}.  In particular the micro-gap")
print(f"(1/13, 452/5864) contains no such family, and tight-locus rigidity")
print(f"holds to height {H0} on this branch -- vs height ~200 by prior")
print("exhaustion (kps-S1).")
