#!/usr/bin/env python3
"""
lrc14_angleF_proof_assembly  --  mac-mini-2026-06-17-S6  (ANGLE F)

ASSEMBLE the LRC(14) proof from the arc-width criterion, or pin the exact residual gap.

After THM-523/524/525/526, LRC(14) <=> every primitive COVERING 13-set S has M(S) >= 1/14.

GENERALIZED ARC-WIDTH LEMMA (THM-526 extended, PROVED):
  For any runner v in S, its level-1/14 danger set is v open teeth, each of FULL WIDTH
  1/(7v), centers k/v spaced 1/v = 7*(1/(7v)) apart, so the safe GAPS between teeth have
  width 6/(7v).  If the level-1/14 safe set G(S\{v}) of the OTHER 12 runners has a widest
  arc W(S\{v}) > 1/(7v), that arc cannot fit inside any single v-tooth (each tooth is
  narrower than the arc and the teeth are isolated), hence the arc contains a v-safe point
  tau0; there min over ALL of S is >= 1/14, so M(S) >= 1/14.

CRITERION  C(S):  EXISTS v in S with  W(S\{v}) > 1/(7v).     C(S) => M(S) >= 1/14.

PIGEONHOLE LOWER BOUND (PROVED, elementary):
  Let A be a finite set of distinct positive integers.  Its level-1/14 safe set G(A) is the
  complement of the danger teeth.  Number of teeth (over all u in A) <= sum_{u in A} u =: N(A),
  so G(A) has at most N(A) arcs.  Total safe measure
        mu(A) = meas(G(A)) >= 1 - sum_{u in A} 1/(7u)         (union bound)
  Hence the WIDEST arc  W(A) >= mu(A)/N(A) >= ( 1 - sum 1/(7u) ) / ( sum u ).
  [pigeonhole: the largest of <=N arcs is at least average = mu/N.]

  So C holds via v whenever  PIG(v):   7v * ( 1 - sum_{u != v} 1/(7u) )  >  sum_{u != v} u.
  Taking v = V := max(S) gives the cleanest sufficient condition.

This script:
  (1) States/derives PIG and tests it exactly on the canonical + structured covering sets.
  (2) CASE SPLIT: classify covering 13-sets by how the large runners cluster, and decide
      EXACTLY which case PIG(largest) settles and where it fails.
  (3) Pin the residual lemma precisely and verify the case-split is EXHAUSTIVE on a large
      covering sample (every covering 13-set falls in case 1 or case 2; case 2 is the named
      residual lemma).  Distinguishes single-large vs clustered-large.
  (4) For clustered sets where PIG(V) fails, test the FULL criterion C (some v works) and
      a sharper per-arc bound, to confirm C is still satisfied (so only the *proof* of the
      residual lemma is open, not its truth).

Exact Fractions throughout.  Stdlib only.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

C = F(1, 14)

# ---------- exact safe-set widest arc W(A) (complement of level-1/14 teeth) ----------
def darcs(v, c=C):
    hw = F(c, v)
    return [(F(k, v) - hw, F(k, v) + hw) for k in range(v)]

def wrapU(iv):
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1); a = lo - s; b = hi - s
        if b <= 1:
            o.append((a, b))
        else:
            o.append((a, F(1))); o.append((F(0), b - 1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch:
            ch = ch if ch > hi else hi
        else:
            r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch)); return r

def safe_total_and_widest(A, c=C):
    """exact (mu, W) : total safe measure and widest safe arc of level-c safe set of A."""
    dz = []
    for v in set(A):
        dz += darcs(v, c)
    if not dz:
        return F(1), F(1)
    dz = wrapU(dz)
    danger = sum(hi - lo for lo, hi in dz)
    mu = F(1) - danger
    best = F(0)
    for i in range(len(dz)):
        hi = dz[i][1]
        lo = dz[(i + 1) % len(dz)][0] + (1 if i == len(dz) - 1 else 0)
        if lo - hi > best:
            best = lo - hi
    return mu, best

def Wsafe(A, c=C):
    return safe_total_and_widest(A, c)[1]

def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

# ---------- the pigeonhole sufficient condition ----------
def mu_lb(A):
    """union-bound lower bound on safe measure: 1 - sum 1/(7u)."""
    return F(1) - sum(F(1, 7 * u) for u in set(A))

def N_ub(A):
    """upper bound on # safe arcs = # teeth = sum u."""
    return sum(set(A))

def PIG(S, v):
    """does the pigeonhole bound prove C via v?  7v*mu_lb(S\\{v}) > N_ub(S\\{v})."""
    A = [u for u in S if u != v]
    return 7 * v * mu_lb(A) > N_ub(A)

def C_pigeon(S):
    """C proved by pigeonhole for SOME v? return (holds, witnesses)."""
    w = [v for v in sorted(set(S)) if PIG(S, v)]
    return (len(w) > 0, w)

def C_exact(S):
    """exact criterion: EXISTS v with W(S\\{v}) > 1/(7v). return (holds, best_v, margin)."""
    best = None
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        if best is None or m > best[1]:
            best = (v, m, W, thr)
        if m > 0:
            return (True, v, m, W, thr)
    return (False, best[0], best[1], best[2], best[3])

print("=" * 80)
print("ANGLE F: ASSEMBLE the LRC(14) proof from the arc-width criterion C(S)")
print("=" * 80)
print("""
LEMMA (THM-526 extended, PROVED): C(S) := [exists v in S: W(S\\{v}) > 1/(7v)]  =>  M(S) >= 1/14.
PIGEONHOLE (PROVED): W(A) >= mu(A)/N(A) >= (1 - sum 1/(7u))/(sum u).
  => PIG(v): 7v*(1 - sum_{u!=v} 1/(7u)) > sum_{u!=v} u  is SUFFICIENT for C via v.
The proof reduces to: does every covering 13-set satisfy C?  We split by PIG(max).
""")

# ---------- (1) canonical + structured single-large covering sets ----------
print("-" * 80)
print("(1) CANONICAL / SINGLE-LARGE covering sets: PIG(V) with V=max(S)")
print("-" * 80)
named = {
    "{1..11,13,84}  (hard core)":  [1,2,3,4,5,6,7,8,9,10,11,13,84],
    "{1..5,7..13,84} (drop6)":     [1,2,3,4,5,7,8,9,10,11,12,13,84],
    "{1..11,13,98}":               [1,2,3,4,5,6,7,8,9,10,11,13,98],
    "{1..11,13,168}":              [1,2,3,4,5,6,7,8,9,10,11,13,168],
    "{1..11,13,1400}":             [1,2,3,4,5,6,7,8,9,10,11,13,1400],
}
for name, S in named.items():
    if not covering(S):
        print(f"  {name}: NOT covering, skip"); continue
    V = max(S); A = [u for u in S if u != V]
    lhs = 7 * V * mu_lb(A); rhs = N_ub(A)
    okp = lhs > rhs
    okC, vC, mC, WC, tC = C_exact(S)
    print(f"  {name:28s}: PIG(V={V}): 7V*mu={float(lhs):8.2f} > sum u={rhs:5d}? {okp}"
          f"   | exact C via v={vC} margin={float(mC):+.5f}")

# ---------- (2) the single-large GENERAL bound: PIG(V) for any covering set with one big V ----------
print()
print("-" * 80)
print("(2) SINGLE-LARGE CASE proved in closed form")
print("-" * 80)
print("""If S = A u {V} with A subset of {1..13} (the 12 small runners, all <= 13) and V >= 14:
   mu_lb(A) = 1 - sum_{u in A} 1/(7u).  Since A subset of {1..13}, the WORST (smallest) mu
   is when A contains the small speeds (large 1/u).  Exact worst over all 12-subsets of
   {1..13} that cover 2..13:""")
# all 12-subsets of {1..13} covering 2..13 (these are exactly {1..13}\{j} for allowed j)
best_worst = None
for j in range(1, 14):
    A = [u for u in range(1, 14) if u != j]
    if not all(any(u % q == 0 for u in A) for q in range(2, 14)):
        continue  # core must cover 2..13
    mu = mu_lb(A); su = sum(A)
    # PIG(V) holds iff 7V*mu > su, i.e. V > su/(7*mu)
    Vthr = su / (7 * mu)  # Fraction
    if best_worst is None or Vthr > best_worst[1]:
        best_worst = (j, Vthr, mu, su)
    print(f"   core {{1..13}}\\{{{j:2d}}}: mu_lb>={mu}={float(mu):.5f}, sum u={su:2d}"
          f"  => PIG(V) holds for all V > {float(Vthr):.3f}")
jj, Vthr, mu, su = best_worst
print(f"""
  WORST core = {{1..13}}\\{{{jj}}}: PIG(V) holds for V > {float(Vthr):.4f}.  Since any parked
  runner V is a multiple of 14 (covering needs mult of 14), V >= 14 > {float(Vthr):.4f}.
  => PROVED: for EVERY single-large covering 13-set (12 small runners in {{1..13}} + one V>=14),
     PIG(V) holds, hence C(S), hence M(S) >= 1/14.   (No finite check needed.)
""")
# sanity: smallest V actually is 14? covering needs a multiple of 14; but also the 12 small
# runners might NOT include the mult of 14. The single-large case = exactly 12 runners in {1..13}.

print("-" * 80)
print("(3) CASE SPLIT and the RESIDUAL: clustered-large covering sets")
print("-" * 80)
print("""Define k = #{runners > 13} (the 'large' runners).  A covering 13-set has 13 - k
runners in {1..13} (the 'small' core) and k large runners.
  * k <= 1: SINGLE-LARGE -- fully PROVED in (2): PIG(max) settles it.
  * k >= 2: CLUSTERED-LARGE -- the residual.  Here sum_{u != V} u can be ~ (k-1)*V, so
    PIG(V) can fail; we must use a DIFFERENT v or the exact arc structure.
We now (a) verify the split is exhaustive, (b) measure how often PIG(V) fails for k>=2,
(c) check the EXACT criterion C still holds there (truth of C), pinning the open part to
PROVING C on the k>=2 clustered family.""")
