#!/usr/bin/env python3
"""
c3_spectrum_kpc1.py -- THREAD C (HYP-2368, kind-pasteur-2026-06-10-S1 sub-agent kpc1)

PART 1 (SANITY): verify the score-sequence identity for the number of 3-cycles
    c3(T) = C(n,3) - sum_v C(s_v, 2)
(Kendall & Babington Smith 1940; also in Moon, Topics on Tournaments; often quoted
via Goodman) against DIRECT 3-cycle counting:
  - ALL 2^C(n,2) labeled tournaments for n = 4, 5  (64 and 1024 tournaments)
  - 500 random labeled tournaments each for n = 6, 7, 8  (fixed seed, reproducible)

PART 2 (BRUTE SPECTRUM): enumerate ALL Landau score sequences (nondecreasing,
prefix sums >= C(k,2), total = C(n,2)) by backtracking for n = 3..13 and collect
the achievable set of f(s) = sum C(s_i,2).  By Landau's theorem every such
sequence is realized by a tournament, and by the identity above
    achievable c3 values on n vertices  =  { C(n,3) - f(s) : s Landau }.
Report: max c3 (vs the classical formula (n^3-n)/24 odd / (n^3-4n)/24 even),
spectrum size, and ALL gaps (missing values) in [0, max].

EXACT integer arithmetic throughout. No floats.
"""

import sys
import random
from itertools import combinations

def comb2(x):
    return x * (x - 1) // 2

def comb3(x):
    return x * (x - 1) * (x - 2) // 6

# ----------------------------------------------------------------------------
# PART 1: sanity of the identity
# ----------------------------------------------------------------------------

def c3_direct(n, beats):
    """beats[i][j] = 1 if i -> j. Count cyclic triples directly:
    a triple is cyclic iff no vertex of it beats the other two."""
    cnt = 0
    for i, j, k in combinations(range(n), 3):
        # out-degrees within the triple
        di = beats[i][j] + beats[i][k]
        dj = beats[j][i] + beats[j][k]
        dk = beats[k][i] + beats[k][j]
        if di != 2 and dj != 2 and dk != 2:
            cnt += 1
    return cnt

def c3_score(n, beats):
    scores = [sum(beats[i]) for i in range(n)]
    return comb3(n) - sum(comb2(s) for s in scores)

def all_tournaments_check(n):
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    bad = 0
    seen_values = set()
    for mask in range(1 << m):
        beats = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (mask >> b) & 1:
                beats[i][j] = 1
            else:
                beats[j][i] = 1
        d = c3_direct(n, beats)
        s = c3_score(n, beats)
        if d != s:
            bad += 1
            if bad <= 5:
                print(f"  MISMATCH n={n} mask={mask}: direct={d} score-formula={s}")
        seen_values.add(d)
    print(f"n={n}: ALL {1 << m} labeled tournaments checked, mismatches={bad}, "
          f"c3 values seen={sorted(seen_values)}")
    return bad == 0, seen_values

def random_tournaments_check(n, trials, rng):
    pairs = list(combinations(range(n), 2))
    bad = 0
    for _ in range(trials):
        beats = [[0] * n for _ in range(n)]
        for (i, j) in pairs:
            if rng.getrandbits(1):
                beats[i][j] = 1
            else:
                beats[j][i] = 1
        d = c3_direct(n, beats)
        s = c3_score(n, beats)
        if d != s:
            bad += 1
            if bad <= 5:
                print(f"  MISMATCH n={n}: direct={d} score-formula={s}")
    print(f"n={n}: {trials} random tournaments checked, mismatches={bad}")
    return bad == 0

# ----------------------------------------------------------------------------
# PART 2: brute-force Landau enumeration
# ----------------------------------------------------------------------------

def landau_fset(n):
    """Set of f(s) = sum C(s_i,2) over all Landau score sequences on n vertices.
    Backtracking with prefix-sum pruning. Also counts sequences."""
    N2 = comb2(n)
    fset = set()
    count = [0]
    sys.setrecursionlimit(10000)

    def rec(k, last, psum, f):
        if k == n:
            if psum == N2:
                fset.add(f)
                count[0] += 1
            return
        r = n - k  # scores still to choose (including this one)
        for s in range(last, n):
            ps = psum + s
            if ps + (r - 1) * s > N2:      # remaining scores >= s, monotone in s
                break
            if ps < comb2(k + 1):           # Landau prefix condition at k+1
                continue
            if ps + (r - 1) * (n - 1) < N2:  # cannot reach the total
                continue
            rec(k + 1, s, ps, f + comb2(s))

    rec(0, 0, 0, 0)
    return fset, count[0]

def max_c3_formula(n):
    if n % 2 == 1:
        return (n**3 - n) // 24
    return (n**3 - 4 * n) // 24

def spectrum_report(n):
    fset, nseq = landau_fset(n)
    C3 = comb3(n)
    c3set = {C3 - f for f in fset}
    mx = max(c3set)
    mn = min(c3set)
    full = set(range(0, mx + 1))
    gaps = sorted(full - c3set)
    formula = max_c3_formula(n)
    ok_max = (mx == formula)
    print(f"n={n:2d}: #LandauSeqs={nseq:7d}  min c3={mn}  max c3={mx} "
          f"(formula {formula} {'OK' if ok_max else 'MISMATCH!'})  "
          f"spectrum size={len(c3set)}  gaps={gaps if gaps else 'NONE'}")
    return c3set, gaps, nseq

def main():
    print("=" * 78)
    print("PART 1: SANITY -- c3(T) = C(n,3) - sum C(s_v,2) vs direct count")
    print("=" * 78)
    ok = True
    for n in (4, 5):
        good, _ = all_tournaments_check(n)
        ok = ok and good
    rng = random.Random(20260610)
    for n in (6, 7, 8):
        ok = ok and random_tournaments_check(n, 500, rng)
    print(f"PART 1 RESULT: {'ALL IDENTITY CHECKS PASS' if ok else '*** FAILURES ***'}")
    print()
    print("=" * 78)
    print("PART 2: BRUTE SPECTRUM via full Landau-sequence enumeration, n=3..13")
    print("=" * 78)
    seq_counts = []
    sizes = []
    maxes = []
    for n in range(3, 14):
        c3set, gaps, nseq = spectrum_report(n)
        seq_counts.append(nseq)
        sizes.append(len(c3set))
        maxes.append(max(c3set))
    print()
    print("OEIS lookup material:")
    print(f"  #Landau sequences (n=3..13)   : {seq_counts}")
    print(f"    (expect A000571: 2,4,9,22,59,167,490,1486,4639,14805,48107)")
    print(f"  max c3 sequence  (n=3..13)    : {maxes}")
    print(f"  spectrum size    (n=3..13)    : {sizes}")

if __name__ == "__main__":
    main()
