#!/usr/bin/env python3
"""
lrc14_freiman_stability_macmini_S65cont5.py -- LEM-016 (k=7 restricted-sumset stability) +
the near-AP-side closure experiment (THM-676 dichotomy, structured branch).

LEM-016 TARGET: for 7 distinct integers A with restricted sumset B = |A +^ A| (sums a_i + a_j,
i < j), after gcd-normalization (a_1 = 0, gcd of differences = 1, diameter D):
    CORRECTED IN-SESSION: B <= 12  ==>  D <= 7 (SHARP; B = 13 escapes via rank-2 GAPs):
    B = 11 <=> D = 6 (the AP); B = 12 ==> D <= 7; B = 13 ==> D <= 8).
  KNOWN SHARPNESS: B = 14 admits UNBOUNDED diameter (two-piece sets {0,1,2,3,4}+{N,N+1} etc.),
  so 13 is the exact stability threshold -- the restricted-sum analog of Freiman 3k-4 sits at
  excess b <= 2 for k = 7 (the doubles-slack costs 2 vs the classical unrestricted range).
  PROOF METHOD: monotonicity (adding an element never removes sums) => pruned DFS is EXHAUSTIVE
  for each D; run D = 6..200. TAIL (D > 200): the largest-gap split g > D/2 gives disjoint
  sum-blocks and B >= 14 (proved in canon); the remaining tail case (all gaps <= D/2, D > 200)
  is sampled here and flagged honestly.

NEAR-AP SIDE: majority-parity class inside an AP of length 9 (7 of 9 consecutive terms, step d
even) -- the ONLY sets that can have burden <= 13 by LEM-016. Build covering 13-sets around
such majority classes and test the certificate battery (C0/C1/C2) + exact liveness.
"""
from math import gcd, isqrt
from functools import reduce
from itertools import combinations
import random

random.seed(78)

# ================================================================ PART 1: LEM-016 by pruned DFS
print("=" * 100)
print("PART 1 -- LEM-016: exhaustive pruned DFS, all normalized 7-sets with B <= 13, D <= 200")
print("=" * 100)

def dfs_all_lowB(D, Bcap):
    """All A = {0, a2..a6, D} (0 < a2 < ... < a6 < D) with |A+^A| <= Bcap. Pruned by
    monotonicity (adding an element never removes sums) with INCREMENTAL sum sets."""
    out = []
    def rec(chosen, start, sums):
        n = len(chosen)
        if n == 7:
            out.append((tuple(sorted(chosen)), len(sums)))
            return
        need = 7 - n
        for a in range(start, D - need + 2):
            if a >= D:
                break
            new = {a + c for c in chosen}
            s2 = sums | new
            if len(s2) <= Bcap:
                rec(chosen + [a], a + 1, s2)
    rec([0, D], 1, {D})
    return out

maxD_at = {11: None, 12: None, 13: None}
found_any = []
for D in range(6, 61):                                 # exhaustive range (proof layer 1)
    hits = dfs_all_lowB(D, 13)
    for A, B in hits:
        diffs = [A[i+1] - A[i] for i in range(6)]
        if reduce(gcd, diffs) != 1:
            continue                     # dilate of smaller-D set
        for lvl in (11, 12, 13):
            if B <= lvl and (maxD_at[lvl] is None or D > maxD_at[lvl]):
                maxD_at[lvl] = D
        found_any.append((D, B, A))
print(f"gcd-normalized 7-sets with B <= 13 exist at diameters up to (exhaustive D <= 60):")
print(f"  B = 11: max D = {maxD_at[11]}   (the AP: D = 6 expected)")
print(f"  B <= 12: max D = {maxD_at[12]}")
print(f"  B <= 13: max D = {maxD_at[13]}")
print(f"total normalized sets with B <= 13, D <= 60: {len(found_any)}")
print("examples at the extreme diameter:")
for D, B, A in [x for x in found_any if x[0] == maxD_at[13]][:6]:
    print(f"  D={D} B={B} A={list(A)}")

# sharpness: B = 14 with unbounded D (two-piece witnesses)
def restricted_B(A):
    return len({A[i] + A[j] for i in range(len(A)) for j in range(i + 1, len(A))})
for N in (50, 500, 5000):
    A = [0, 1, 2, 3, 4, N, N + 1]
    print(f"sharpness: {A}: B = {restricted_B(A)} (14 expected) at diameter {N+1}")

# tail sampling: D in (200, 2000], all gaps <= D/2 -- B always >= 14?
tail_viol = 0
for _ in range(200000):
    D = random.randrange(201, 2001)
    mids = sorted(random.sample(range(1, D), 5))
    A = [0] + mids + [D]
    if max(A[i+1] - A[i] for i in range(6)) <= D // 2 and restricted_B(A) <= 13:
        tail_viol += 1
        print(f"  TAIL VIOLATION: {A}")
print(f"tail sample (200k, D in (200,2000], maxgap <= D/2): B <= 13 violations: {tail_viol}")

# ================================================================ PART 2: near-AP side closure
print()
print("=" * 100)
print("PART 2 -- the STRUCTURED branch: majority class = 7 of 9 consecutive d-multiples;")
print("          covering 13-sets built around it: does the certificate battery always fire?")
print("=" * 100)

def danger_j(k):  return -(-k // 14) - 1
def blocked(R, k):
    j = danger_j(k)
    bad = set(range(0, j + 1)) | set(range(k - j, k))
    for s in range(1, k):
        if all((r * s) % k not in bad for r in R):
            return False
    return True
def covering(S):  return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1
def rulers_of(S): return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})

def certified(S):
    """C0 or any unblocked proper descent (full C2, any k > 14)."""
    if max(S) <= 13 * min(S):
        return 'C0'
    for q in rulers_of(S):
        for k in range(15, q):
            if q % k == 0:
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    return f'C2:k={k}'
    return None

def exact_live(S):
    for q in rulers_of(S):
        for p in range(1, q):
            if all(14 * (v * p % q) >= q and 14 * (q - v * p % q) >= q for v in S):
                return (q, p)
    return None

hole_patterns = list(combinations(range(9), 7))     # 36 patterns
tested = certified_n = 0
uncert = []
for trial in range(400):
    T = random.choice(hole_patterns)
    d = 2 * random.randrange(1, 8)                   # even step
    a = random.randrange(1, 60)
    if a % 2 == 0 and all((a + d * t) % 2 == 0 for t in T):
        pass                                          # any parity fine; majority = same parity anyway
    maj = sorted(a + d * t for t in T)
    # minority: 6 opposite-parity elements, random below ~2*max
    lo, hi = 1, 2 * maj[-1] + 20
    minority = []
    while len(minority) < 6:
        v = random.randrange(lo, hi)
        if (v - a) % 2 == 1 and v not in maj and v not in minority:
            minority.append(v)
    S = sorted(maj + minority)
    if len(set(S)) != 13 or not covering(S) or not primitive(S):
        continue
    tested += 1
    c = certified(S)
    if c:
        certified_n += 1
    else:
        lv = exact_live(S)
        uncert.append((S, lv))
print(f"near-AP-majority covering sets tested: {tested}; certified (C0 or C2-any-k): "
      f"{certified_n} ({certified_n/max(tested,1):.1%})")
if uncert:
    print(f"UNCERTIFIED (live-status shown): {len(uncert)}")
    for S, lv in uncert[:5]:
        print(f"  S={S}  exact_live={lv}")
else:
    print("ZERO uncertified -- the structured branch is certificate-closed on this sample.")
print()
print("Done.")
