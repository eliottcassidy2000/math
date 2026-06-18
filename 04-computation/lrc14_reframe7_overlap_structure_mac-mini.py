#!/usr/bin/env python3
"""
REFRAME 7, part 2 — probing for NON-tautological leverage in the band system.

Two leads from part 1:
 (A) The clearing-num count is ALWAYS EVEN (2,4,6,...). Why?  -> num <-> D-num pairing:
     if num clears at level k/D then so does D-num (||v(D-num)/D|| = ||v num/D||).
     So clearing nums come in pairs {num, D-num}; the count is 2*(#pairs).  We CONFIRM this
     and ask: does the GUARANTEE "count is even & nonneg" + "covering" force count >= 2?
     I.e., can we LOWER-BOUND the number of clearing pairs by a covering/CRT argument?

 (B) Inclusion-exclusion / overlap: union bound fails because forbidden caps overlap.
     Compute the SECOND-order Bonferroni term (sum of pairwise intersections) -- does
     I-E with 2 terms give a positive lower bound on the clearing count?  And: is the
     overlap structure DETERMINED by the covering condition (i.e., does "S has a mult of
     each q" pin down |F_u cap F_v|)?  If overlaps are forced large, clearance follows.

 (C) The KEY structural test: forbidden set F_v at modulus D depends only on gcd(v,D) and
     the residue v mod D.  The covering condition says: for each q in 2..14, some v ≡ 0
     mod q.  At the M-attaining D, which runners are "small forbidden" (gcd(v,D) small)?
     Is there a runner whose forbidden set is SMALL enough that, intersected with the
     others' allowed sets, leaves a clearing num by a pigeonhole on a SINGLE modulus?
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
from collections import Counter
import random

N = 14

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

def gcd_all(S): return reduce(gcd, S)
def is_covering(S, Nval=N):
    return all(any(v % q == 0 for v in S) for q in range(2, Nval + 1))
def make_covering_13set(maxspeed=120, rng=random):
    for _ in range(4000):
        S = set(rng.sample(range(1, 30), rng.randint(6, 11)))
        while len(S) < 13: S.add(rng.randint(1, maxspeed))
        S = set(list(S)[:13])
        if len(S) == 13 and gcd_all(S) == 1 and is_covering(S):
            return S
    return None

def clearing_nums(S, D, k):
    return [num for num in range(1, D) if all(k <= (v*num) % D <= D-k for v in S)]

# ============================================================
print("="*78); print("REFRAME 7 part 2 — overlap structure / non-tautological leverage"); print("="*78)

# ----- (A) clearing nums come in pairs {num, D-num}; count is even -----
print("\n[A] num <-> D-num pairing (explains EVEN clearing counts)")
random.seed(3)
allpaired = True; checked = 0
for _ in range(300):
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N): continue
    D = at.denominator; k = ceil(F(1, N) * D)
    cl = set(clearing_nums(S, D, k))
    # check num in cl  <=> D-num in cl
    paired = all((D - num) in cl for num in cl)
    checked += 1
    if not paired: allpaired = False; print("  UNPAIRED:", sorted(S), D, sorted(cl)); break
print(f"  checked {checked} covering sets: clearing nums always closed under num->D-num: {allpaired}")
print("  => clearing count = 2*(#unordered pairs); count>=1 forces count>=2. (structural, proved by hand:")
print("     ||v(D-num)/D|| = ||(-v num)/D|| = ||v num/D||, so D-num clears iff num clears.)")

# ----- (B) Bonferroni 2nd-order: does inclusion-exclusion (2 terms) lower-bound the count? -----
print("\n[B] Inclusion-Exclusion lower bound on clearing count (Bonferroni)")
print("    #clearing >= (D-1) - sum|F_v| + (correction).  Bonferroni-2 lower bound:")
print("    #good = #{num: avoid all F_v} >= (D-1) - sum|F_v|  (union bound, term 1)")
print("    A *better* bound uses the structure of F_v as unions of gcd(v,D) caps.")
def forbidden_set(v, D, k):
    return set(num for num in range(1, D) if (v*num) % D < k or (v*num) % D > D-k)
random.seed(5)
# For each set: union-bound value, true count, and the 'deficit' that I-E must overcome.
deficits = []; nsets = 0
for _ in range(40000):
    if nsets >= 120: break
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N): continue
    D = at.denominator; k = ceil(F(1, N) * D)
    if 2*k-1 >= D: continue
    Sl = sorted(S)
    forb = [forbidden_set(v, D, k) for v in Sl]
    union = set().union(*forb)
    truecount = (D-1) - len(union)
    sum_single = sum(len(fv) for fv in forb)
    ub = (D-1) - sum_single  # union bound (always negative)
    # pairwise intersections (positive correction)
    pair_corr = 0
    for i in range(len(forb)):
        for j in range(i+1, len(forb)):
            pair_corr += len(forb[i] & forb[j])
    bonf2 = (D-1) - sum_single + pair_corr  # Bonferroni: overcounts, this is an UPPER bound on #good...
    # Actually inclusion-exclusion: |union| <= sum_single - pair_corr + ...
    # |union| >= sum_single - pair_corr (Bonferroni lower bound on union -> upper bound on #good)
    # #good = (D-1) - |union|. Lower bound on #good needs UPPER bound on |union| = sum_single (union bd).
    # So I-E pair term does NOT directly help lower-bound #good. Record the gap.
    nsets += 1
    deficits.append((-ub, truecount, pair_corr))  # (how negative UB is, actual count, overlap mass)
import statistics
neg = [d[0] for d in deficits]; tc = [d[1] for d in deficits]; pc=[d[2] for d in deficits]
print(f"  over {nsets} sets at M-attaining D:")
print(f"    union-bound deficit  (sum|F_v| - (D-1)) : min={min(neg)} max={max(neg)} mean={statistics.mean(neg):.1f}")
print(f"    true clearing count                      : min={min(tc)} max={max(tc)} mean={statistics.mean(tc):.1f}")
print(f"    total pairwise overlap mass sum|F_i&F_j| : min={min(pc)} max={max(pc)} mean={statistics.mean(pc):.0f}")
print("  Bonferroni inequality direction: |union| >= sum|F| - sum|F_i&F_j|.  That gives an")
print("  UPPER bound on #good, not a lower bound.  I-E does NOT yield a clearance proof.")
print("  The overlap is large precisely BECAUSE clearance is hard; counting can't see it.")

# ----- (C) single-modulus pigeonhole via covering: smallest forbidden runner -----
print("\n[C] Covering -> single-modulus structure of F_v (does a SMALL F_v save us?)")
print("    |F_v| at modulus D depends on gcd(v,D). The cap [0,k-1]U[D-k+1,D-1] pulled back")
print("    by mult-by-v has size ~ (2k-1) when gcd(v,D)=1, larger when gcd>1.")
print("    Covering forces some v=14*t (mult of 14) and some v=q*t for each q<=14.")
print("    Q: at the M-attaining D, does any runner give |F_v| << (2k-1)? (would shrink union)")
random.seed(8)
samples = []
for _ in range(40000):
    if len(samples) >= 8: break
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N): continue
    D = at.denominator; k = ceil(F(1, N)*D)
    if 2*k-1 >= D: continue
    Sl = sorted(S)
    rows = []
    for v in Sl:
        gg = gcd(v, D)
        fv = len(forbidden_set(v, D, k))
        rows.append((v, gg, fv))
    samples.append((Sl, D, k, str(m), rows))
for Sl, D, k, mstr, rows in samples[:5]:
    print(f"  S={Sl}")
    print(f"    D={D} k={k} M={mstr}  band size = {D-2*k+1} residues")
    print(f"    (v, gcd(v,D), |F_v|): {[(v,gg,fv) for v,gg,fv in rows]}")
    minfv = min(r[2] for r in rows)
    print(f"    smallest |F_v|={minfv}; base cap size (gcd=1) = {2*k-1}; "
          f"sum|F_v|={sum(r[2] for r in rows)} vs (D-1)={D-1}")

# ----- (D) THE decisive non-tautology check: is there ALWAYS a single modulus D in a
#           BOUNDED family (e.g. D = lcm-like, or D <= poly(N)) where covering forces clearance? -----
print("\n[D] Is clearance forced at a SMALL/structured D (vs needing the full pair-search)?")
print("    Test: for covering minimizers, what is the M-attaining D? bounded by what?")
random.seed(21)
Ds_seen = Counter(); maxD = 0; maxD_case=None; nsets=0
for _ in range(60000):
    if nsets >= 400: break
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N): continue
    nsets += 1
    D = at.denominator
    Ds_seen[D] += 1
    if D > maxD: maxD = D; maxD_case = (sorted(S), str(m), D)
print(f"  over {nsets} covering sets: M-attaining D ranges; max D = {maxD}")
print(f"    max-D case: S={maxD_case[0]} M={maxD_case[1]} D={maxD_case[2]}")
print(f"    most common attaining D: {Ds_seen.most_common(10)}")
print("  Conclusion line printed at end.")
print("\nDONE.")
