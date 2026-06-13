#!/usr/bin/env python3
"""collatz_defect_rapidity.py — The EXACT rapidity conservation law for Collatz.

Builds on collatz_rapidity_s116n.py, which observed that in rapidity space
rho(n) = (1/2)ln(n) = arctanh((n-1)/(n+1)) = log_F(Cayley address of n),
the Collatz map is a near-additive walk with drift ln(3/4)/2, "plus corrections".

THIS SESSION nails the corrections. They are not a fudge term: they are an
EXACT, telescoping harmonic defect. We prove the conservation law

        ln n  =  K*ln2  -  L*ln3  -  D(n)                              (EQ*)

for the Syracuse (odd-to-odd) map, where along the trajectory
    n = a_0 -> a_1 -> ... -> a_L = 1,   a_i = (3 a_{i-1}+1)/2^{k_i},
    L = number of odd steps,  K = sum k_i (total halvings),
    D(n) = sum_{i=0}^{L-1} ln(1 + 1/(3 a_i))  >  0   (the HARMONIC DEFECT).

Consequences proved/verified here:
  (1) EQ* is exact to machine precision (and its multiplicative form is an
      exact integer identity, verified with Fraction).
  (2) The +1 is EXACTLY the multiplicative fudge prod 3a_i/(3a_i+1) = e^{-D(n)}
      that makes the otherwise-impossible ratio 2^K / 3^L land on the integer n.
  (3) UNCONDITIONAL inequality: every n>1 that reaches 1 satisfies
      K/L > log2(3) = 1.58496...,  with  K - L*log2(3) = log2(n) + D(n)/ln2.
      The mean halving exponent strictly exceeds log2(3) -- this is forced, not
      statistical. (The conjecture is only that L,K are finite for all n.)
  (4) D(n) is dominated by the TERMINAL small-odd visits, not by large n.
      It is the arithmetic residue of the formal-group obstruction: Collatz
      would be a pure rapidity translation (linear in log_F) but for D(n).

Framework tie-in (repo thesis "arctanh is the universal linearizer"):
  log_F = arctanh = rapidity linearizes multiplication-by-3 and halving into
  translations; D(n) is the EXACT measure of how far the +1 is from being a
  formal-group operation. cf. 07-reflections/grand-trichotomy.md (Collatz row),
  collatz_rapidity_s116n.py, rapidity_numbertheory_s116b.py.

Session: (new) collatz-defect-rapidity, continuing the S116 rapidity thread.
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from math import log, log2
from fractions import Fraction
from collections import Counter
import statistics

LN2, LN3 = log(2), log(3)
LOG2_3 = log2(3)

print("\n  THE EXACT RAPIDITY CONSERVATION LAW FOR COLLATZ\n")
print("=" * 70)

# ----------------------------------------------------------------------
def syracuse_trajectory(n):
    """Odd n -> list of odd values a_0..a_L=1 and exponents k_1..k_L.
    Returns (odds, ks) where 2^{k_i} a_i = 3 a_{i-1} + 1."""
    assert n % 2 == 1
    odds = [n]
    ks = []
    a = n
    while a != 1:
        val = 3 * a + 1
        k = 0
        while val % 2 == 0:
            val //= 2
            k += 1
        ks.append(k)
        a = val
        odds.append(a)
    return odds, ks

# ============================================================
print("\n  I. THE IDENTITY, PROVED BY TELESCOPING")
print("  " + "-" * 50)
print("""
  Each Syracuse step is, by definition,   2^{k_i} a_i = 3 a_{i-1} + 1.
  Take logs:
     k_i ln2 + ln a_i = ln(3 a_{i-1} + 1)
                      = ln3 + ln a_{i-1} + ln(1 + 1/(3 a_{i-1})).
  Sum i = 1..L.  The sum of (ln a_i - ln a_{i-1}) telescopes to
     ln a_L - ln a_0 = ln 1 - ln n = -ln n.
  Hence
     K ln2 - ln n = L ln3 + sum_{i=0}^{L-1} ln(1 + 1/(3 a_i)),
  i.e.
     ln n = K ln2 - L ln3 - D(n),     D(n) = sum ln(1 + 1/(3 a_i)) > 0.   (EQ*)

  Dividing by 2 gives the RAPIDITY form (rho = (1/2)ln = log_F = arctanh):
     rho(n) = (K ln2 - L ln3)/2 - D(n)/2.
  Since rho(1)=0, the orbit spends exactly  (L ln3 - K ln2)/2 + D(n)/2  of
  rapidity descending to the origin. The linear part is a formal-group
  translation; D(n)/2 is the ENTIRE non-formal-group content of Collatz.
""")

# ============================================================
print("  II. NUMERICAL VERIFICATION OF EQ* (machine precision)")
print("  " + "-" * 50)
print(f"  {'n':>7} {'L':>4} {'K':>5} {'K/L':>9} {'D(n)':>9} "
      f"{'K ln2 - L ln3 - D':>18} {'ln n':>10} {'err':>9}")
worst = 0.0
for n in [3, 7, 27, 97, 255, 871, 6171, 9663, 77031, 837799]:
    odds, ks = syracuse_trajectory(n)
    L = len(ks)
    K = sum(ks)
    D = sum(log(1 + 1.0/(3*a)) for a in odds[:-1])  # i=0..L-1
    lhs = K*LN2 - L*LN3 - D
    err = abs(lhs - log(n))
    worst = max(worst, err)
    print(f"  {n:>7} {L:>4} {K:>5} {K/L:>9.5f} {D:>9.5f} "
          f"{lhs:>18.10f} {log(n):>10.5f} {err:>9.1e}")
print(f"\n  Worst floating error over sample: {worst:.2e}  (== rounding)\n")

# ============================================================
print("  III. EXACT (RATIONAL) FORM:  the +1 is the integer-making fudge")
print("  " + "-" * 50)
print("""
  Exponentiate EQ*:   n = (2^K / 3^L) * e^{-D(n)},   and
     e^{-D(n)} = prod_{i} 1/(1 + 1/(3 a_i)) = prod_i  3 a_i / (3 a_i + 1).
  So, EXACTLY (rational identity, no logs):
     n * 3^L * prod_i (3 a_i + 1)/(3 a_i)  =  2^K.
  The pure map x->3x (no +1) would demand 2^K = 3^L * n, impossible for n>1
  since 3 never divides 2^K. The product prod (3a_i+1)/(3a_i) = e^{D(n)} > 1 is
  the EXACT rational correction that repairs the divisibility -- the arithmetic
  shadow of the formal-group obstruction.
""")
for n in [3, 27, 97, 871]:
    odds, ks = syracuse_trajectory(n)
    L, K = len(ks), sum(ks)
    # exact rational e^{D} = prod (3a_i+1)/(3a_i)
    eD = Fraction(1)
    for a in odds[:-1]:
        eD *= Fraction(3*a + 1, 3*a)
    lhs = Fraction(n) * (3**L) * eD          # should equal 2^K exactly
    ok = (lhs == 2**K)
    print(f"  n={n:>4}: L={L:>2} K={K:>2}  n*3^L*prod(3a+1)/(3a) = 2^K ? {ok}  "
          f"(e^D = {float(eD):.6f})")
print()

# ============================================================
print("  IV. THE FORCED INEQUALITY  K/L > log2(3)  (unconditional)")
print("  " + "-" * 50)
print("""
  From EQ*: K ln2 = L ln3 + ln n + D(n), and D(n) > 0, ln n >= 0 for n>=1.
  Divide by L ln2:
       K/L = log2(3) + (ln n + D(n))/(L ln2)  >  log2(3) = 1.5849625...
  EVERY n>1 reaching 1 has mean halving exponent strictly above log2(3).
  Equivalently  K - L*log2(3) = log2(n) + D(n)/ln2  > 0  exactly.
  This is not a statistical 'on average' -- it holds for each individual orbit.
""")
print(f"  {'n':>7} {'K/L':>10} {'-log2(3)=excess':>16} "
      f"{'(ln n+D)/(L ln2)':>18}")
for n in [3, 27, 255, 6171, 77031, 837799, 8400511]:
    odds, ks = syracuse_trajectory(n)
    L, K = len(ks), sum(ks)
    D = sum(log(1 + 1.0/(3*a)) for a in odds[:-1])
    excess = K/L - LOG2_3
    pred = (log(n) + D)/(L*LN2)
    print(f"  {n:>7} {K/L:>10.6f} {excess:>16.6f} {pred:>18.6f}")
print()

# ============================================================
print("  V. WHAT DOMINATES D(n)?  (the terminal small-odd visits)")
print("  " + "-" * 50)
print("""
  ln(1+1/(3a)) ~ 1/(3a): a single visit to a is worth ~1/(3a) in D.
  Large values along the orbit contribute negligibly; the defect is
  accumulated almost entirely in the last few odd values before 1.
  Single-visit weights:""")
for a in [1, 3, 5, 7, 11, 21, 85, 341]:
    print(f"     a={a:>4}: ln(1+1/(3a)) = {log(1+1/(3*a)):.6f}"
          + ("   <- terminal a=1 dominates" if a == 1 else ""))
print()

# Decompose D(27) by contribution
odds, ks = syracuse_trajectory(27)
contribs = sorted(((log(1+1/(3*a)), a) for a in odds[:-1]), reverse=True)
Dtot = sum(c for c, _ in contribs)
print(f"  D(27) = {Dtot:.6f}; top contributors (value a : share of D):")
cum = 0.0
for c, a in contribs[:6]:
    cum += c
    print(f"     a={a:>4}: {c:.6f}  ({100*c/Dtot:5.1f}%)  cum {100*cum/Dtot:5.1f}%")
print()

# ============================================================
print("  VI. STATISTICS OF D(n) OVER A LARGE SAMPLE")
print("  " + "-" * 50)
N = 200001
Ds, KoverL, ratios_drift = [], [], []
maxD = (0.0, 0)
sumL = sumK = 0
for n in range(3, N, 2):
    odds, ks = syracuse_trajectory(n)
    L, K = len(ks), sum(ks)
    D = sum(log(1 + 1.0/(3*a)) for a in odds[:-1])
    Ds.append(D); KoverL.append(K/L)
    sumL += L; sumK += K
    if D > maxD[0]:
        maxD = (D, n)
print(f"  Sample: odd n in [3, {N}).")
print(f"  D(n):  mean={statistics.mean(Ds):.5f}  median={statistics.median(Ds):.5f}"
      f"  stdev={statistics.pstdev(Ds):.5f}")
print(f"         min={min(Ds):.5f}  max={maxD[0]:.5f} (at n={maxD[1]})")
print(f"  K/L :  mean={statistics.mean(KoverL):.6f}  (cf. log2(3)={LOG2_3:.6f})")
print(f"  Aggregate K/L over all steps = {sumK/sumL:.6f}  (-> log2(3) from above)")
print(f"  Min K/L in sample = {min(KoverL):.6f}  (still > log2(3)? "
      f"{min(KoverL) > LOG2_3})")
print()
# histogram of D
print("  Histogram of D(n):")
buckets = Counter()
for D in Ds:
    buckets[round(D, 1)] += 1
for key in sorted(buckets):
    if buckets[key] >= len(Ds)//200:
        bar = "#" * (buckets[key] * 60 // max(buckets.values()))
        print(f"    D~{key:>4.1f}: {bar} {buckets[key]}")
print()

# ============================================================
print("  VII. THE DEFECT IS BOUNDED  (CONJ-defect-bound)")
print("  " + "-" * 50)
print("""
  CONJ-defect-bound:  there is an absolute constant  D* ~ 0.2257  with
        D(n) < D*   for ALL n >= 1,
  equivalently the exact rational fudge satisfies
        1 <=  prod_i (3 a_i + 1)/(3 a_i)  <  e^{D*} = 1.2531...   for all n.
  The supremum is approached (champion) at  n = 993  (D = 0.225654), which is
  the record holder for ALL odd n up to 5,000,000 (separate search, see
  results/collatz_defect_rapidity.out notes). The record PROGRESSION converges:
        n=3:.1699  n=9:.2220  n=559:.2249  n=745:.2253  n=993:.2257
  so the worst case is NOT growing -- D(n) is (conjecturally) uniformly bounded.

  WHY this is believable, not just observed:  ln(1+1/(3a)) ~ 1/(3a), so D(n) is
  dominated by the few small odd values in the TERMINAL cascade to 1. The
  Syracuse predecessor tree of 1 is {1,5,21,85,341,...} = (4^j-1)/3; any one
  trajectory threads a single root-path and can collect only a bounded budget
  of small odds. Visiting them all (impossible on one path) would give the
  divergent ~(1/6)ln M, but the dynamics forbid it -- hence a finite sup.

  This is the falsifiable edge:  any n with D(n) >= 0.2257 would (a) beat the
  champion and (b) reveal a trajectory lingering anomalously among small odds.
  Tested over the sample below; champion n=993 unbeaten.""")
for hi in [10**3, 10**4, 10**5, 2*10**5]:
    best = max(((sum(log(1+1/(3*a)) for a in syracuse_trajectory(n)[0][:-1]), n)
                for n in range(3, hi, 2)))
    print(f"     odd n<{hi:>7}: max D(n) = {best[0]:.6f} at n={best[1]:>6}  "
          f"(bound D* = 0.2257)")
print()

# ============================================================
print("  VIII. OEIS / FRAMEWORK HOOKS")
print("  " + "-" * 50)
print("  Per-n data the framework cares about (L=odd-step count, K=total")
print("  halvings, and the exact defect numerator). L(n) is A006667-related;")
print("  K(n)+L(n) is the total stopping time A006577. The pair (L,K) are the")
print("  two coordinates of the rapidity translation; D(n) is the third,")
print("  non-formal-group coordinate. Sample table:")
print(f"  {'n':>5} {'L':>4} {'K':>4} {'K+L (A006577?)':>16} {'D(n)':>9}")
for n in [3, 5, 7, 9, 27, 97]:
    odds, ks = syracuse_trajectory(n)
    L, K = len(ks), sum(ks)
    # total stopping time = full Collatz steps to 1 = K + L (each odd step is
    # 1 tripling-step then k_i halvings)
    D = sum(log(1 + 1.0/(3*a)) for a in odds[:-1])
    print(f"  {n:>5} {L:>4} {K:>4} {K+L:>16} {D:>9.5f}")
print()

print("=" * 70)
print("""
  SUMMARY
  -------
  * Collatz in rapidity space obeys the EXACT law  ln n = K ln2 - L ln3 - D(n).
  * The previously-hand-waved 'corrections' = the harmonic defect D(n) > 0,
    equal to log of the exact rational fudge prod (3a_i+1)/(3a_i) that repairs
    the impossible divisibility 2^K = 3^L n.
  * Forced, unconditional:  K/L > log2(3)  for every orbit reaching 1.
  * D(n) is small, terminal-dominated, and (conjecturally) BOUNDED by an
    absolute constant D* ~ 0.2257 (champion n=993, unbeaten to 5,000,000).
    It is the precise, quantified obstruction to Collatz being a formal-group
    (arctanh-linear) translation. This sharpens the repo thesis: arctanh
    linearizes the *3 and /2 into rapidity translations; D(n) is exactly --
    and boundedly -- what the +1 costs.
""")
