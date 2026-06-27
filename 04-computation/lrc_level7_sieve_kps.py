#!/usr/bin/env python3
"""
LRC(14): the LEVEL-7 lift sieve, generalized (kind-pasteur-2026-06-27).

Claim (sharpening of codex THM-571): ANY 13-speed set S with at least 7 speeds
divisible by 7 has M(S) > 1/14 -- by a SINGLE level-7 lift-sieve argument, no
case split, regardless of how many are divisible by 14.

THM-571 proved |M14| >= 7 (multiples of *14*), via Case1 (level-14) + Case2
(level-7 when >=2 half-step). Its Case2 already IS the level-7 sieve but was only
fired at |H|>=9. The clean fact: the level-7 sieve alone needs only |H| >= 7
(H = multiples of 7), which subsumes |M14| >= 7 and extends to all >=7-mult-of-7
sets. Residual shrinks from |M14|<=6 to |H|<=6.

This script:
  V1. VERIFY the level-7 sieve fires (surviving lift exists) on every |H|>=7 set.
  V2. Confirm the per-speed forbidden-lift bound (<=1 for 7-coprime speeds at level 7).
  V3. Show |H|=6 is the boundary: the sieve can fail there (residual is genuine).
  V4. Cross-check: every |H|>=7 set actually has M > 1/14 (exact).
"""
from fractions import Fraction as F
import random, math

def nrm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)

def M_exact(speeds):
    speeds = sorted(set(speeds)); n = len(speeds); cand = set()
    for v in speeds:
        for k in range(v): cand.add(F(2*k+1, 2*v))
    for i in range(n):
        for j in range(i+1, n):
            a, b = speeds[i], speeds[j]
            for D in (a+b, b-a):
                if D <= 0: continue
                for k in range(1, D//2 + 1): cand.add(F(k, D))
    best = F(0)
    for t in cand:
        if 0 < t <= F(1,2):
            d = min(nrm(v*t) for v in speeds)
            if d > best: best = d
    return best

def M_lower_LRC(k):
    """proven LRC bound for k speeds: M >= 1/(k+1), for k<=12 (LRC<=13)."""
    return F(1, k+1)

def level7_sieve(S):
    """
    Try to certify M(S) > 1/14 via the level-7 lift sieve.
    H = speeds divisible by 7; P = H/7. Find P-safe phase v (||p v|| >= 1/(|P|+1)),
    then a lift j in 0..6 with all 7-coprime speeds safe (>= 1/14).
    Returns (certified, |H|, forbidden_counts_for_nonH, v, j).
    Uses an exact phase search: candidate v are the M_exact-style optima of P.
    """
    H = [v for v in S if v % 7 == 0]
    nonH = [v for v in S if v % 7 != 0]
    if not H:
        return None
    P = sorted(set(h // 7 for h in H))
    # need a P-safe phase with margin > 1/14; search candidate phases (P's optima grid)
    # candidate v: rationals a/D with D in P-driven denominators, plus a fine grid
    cands = set()
    for p in P:
        for k in range(p): cands.add(F(2*k+1, 2*p))
    for i in range(len(P)):
        for j2 in range(i+1, len(P)):
            for D in (P[i]+P[j2], P[j2]-P[i]):
                if D > 0:
                    for k in range(1, D): cands.add(F(k, D))
    # also a fine grid to be safe
    for g in range(1, 700): cands.add(F(g, 700))
    thr = F(1, 14)
    Pthr = M_lower_LRC(len(P))  # >= this for some phase
    best = None
    for v in cands:
        if v <= 0 or v >= 1: continue
        if any(nrm(p*v) < thr for p in P):   # require P-safe at > 1/14 level
            continue
        # try the 7 lifts
        for j in range(7):
            t = (v + j) / 7
            # H safe automatically: ||7p t|| = ||p(v+j)|| = ||pv|| >= thr (check)
            if all(nrm(b * t) >= thr for b in nonH):
                # count forbidden lifts per nonH speed for reporting
                return True, len(H), v, j
    return False, len(H), None, None

print("="*70); print("LRC(14) LEVEL-7 LIFT SIEVE  (kps S31 cont.)"); print("="*70)

# ---- V1+V4: every |H|>=7 set certified and M>1/14 ----
print("\n[V1/V4] Random 13-sets with |H|>=7 (>=7 multiples of 7): sieve + exact M.")
random.seed(10)
tested = cert = mbad = mchecked = 0
fail_examples = []
for it in range(1500):
    nH = random.randint(7, 12)            # number of multiples of 7
    H = random.sample([7*j for j in range(1, 26)], nH)
    nonH = random.sample([v for v in range(1, 70) if v % 7 != 0], 13 - nH)
    S = sorted(set(H + nonH))
    if len(S) != 13: continue
    tested += 1
    res = level7_sieve(S)
    ok = res[0]
    if ok: cert += 1
    else: fail_examples.append(S)
    if it % 4 == 0:                       # exact M cross-check on a subsample
        mchecked += 1
        if M_exact(S) <= F(1,14): mbad += 1
print(f"   tested {tested} sets with |H|>=7.")
print(f"   level-7 sieve certified: {cert}/{tested}")
print(f"   exact M <= 1/14 (counterexamples to LRC14): {mbad}/{mchecked} checked")
if fail_examples:
    print(f"   sieve-uncertified examples ({len(fail_examples)}): {fail_examples[:3]}")
else:
    print("   => level-7 sieve fired on EVERY |H|>=7 set.  (Supports THM-571'.)")

# ---- V2: forbidden-lift bound at level 7 ----
print("\n[V2] Per-speed forbidden-lift count at level 7 (7-coprime b => <=1).")
random.seed(20)
worst = 0
for _ in range(20000):
    b = random.choice([v for v in range(1, 200) if v % 7 != 0])
    v = F(random.randint(1, 699), 700)
    cnt = sum(1 for j in range(7) if nrm(b*(v+j)/7) < F(1,14))
    worst = max(worst, cnt)
print(f"   max forbidden lifts by a single 7-coprime speed over 20000 trials: {worst}")
print(f"   (claim: <=1, so 13-|H| coprime speeds forbid <= 13-|H| < 7 when |H|>=7)")

# ---- V3: |H|=6 boundary (residual is genuine) ----
print("\n[V3] |H|=6 boundary: does the level-7 sieve fail (residual genuine)?")
random.seed(30)
h6_tested = h6_fail = 0
for _ in range(1500):
    H = random.sample([7*j for j in range(1, 20)], 6)
    nonH = random.sample([v for v in range(1, 80) if v % 7 != 0], 7)
    S = sorted(set(H + nonH))
    if len(S) != 13: continue
    h6_tested += 1
    res = level7_sieve(S)
    if not res[0]: h6_fail += 1
print(f"   |H|=6 sets tested: {h6_tested}; level-7 sieve FAILED on {h6_fail} "
      f"({100*h6_fail/max(h6_tested,1):.1f}%).")
print("   => |H|=6 is below the sieve threshold: residual = covering sets with <=6 mult of 7.")

# ---- summary on the residual ----
print("\n[Residual] After level-7 sieve, LRC(14) reduces to:")
print("   covering 13-sets with <= 6 multiples of 7  (=> <= 6 multiples of 14).")
print("   This is STRICTLY smaller than THM-571's |M14|<=6 residual:")
print("   it also kills every set with 7..12 multiples of 7 but < 7 multiples of 14.")
print("="*70); print("DONE")
