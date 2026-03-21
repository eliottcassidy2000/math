#!/usr/bin/env python3
"""
ocr_plateau_investigation.py — kind-pasteur-2026-03-21-S16

DEEP INVESTIGATION: Why does OCR plateau near 0.916?

HYPOTHESIS: The OCR plateau is governed by the fraction of tournaments
in the "ambiguous middle" — score sequences where multiple cycle
structures coexist. As n grows, the number of score classes grows,
but so does the within-class diversity.

APPROACH:
1. Exact OCR decomposition: Between-class = f(score distribution),
   Within-class = f(cycle structure diversity)
2. What predicts the within-class residual? Is it alpha_2? c5? c7?
3. Is there a THEORETICAL bound on OCR from combinatorial arguments?
4. Does the regular class (S2=0) have special concentration properties?
5. The Szele bound: random tournament H ~ n!/2^{n-1}.
   What is Var(H|regular) / Var(H)?

Key insight from Alon: max H ~ c*n^{3/2}*n!/2^{n-1}.
Mean H = n!/2^{n-1} for random tournaments.
So max/mean ~ n^{3/2}. This means H has a HEAVY RIGHT TAIL.
The score sequence controls the CENTER but not the tail.

NEW IDEA: OCR measures linear R^2. What about RANK correlation?
If H is monotonically related to S2 but nonlinearly, Spearman > Pearson.
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb, factorial, log2, sqrt
from scipy import stats as sp_stats

def build_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

print("=" * 72)
print("  OCR PLATEAU INVESTIGATION")
print("  kind-pasteur-2026-03-21-S16")
print("=" * 72)

# ============================================================
# PART 1: EXACT ANALYSIS AT n=5,6 — decompose the residual
# ============================================================

print("\n" + "=" * 72)
print("  PART 1: Exact residual decomposition")
print("=" * 72)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m

    data = []
    for bits in range(total):
        A = build_tournament(n, bits)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        S2 = sum((s - (n-1)/2)**2 for s in scores)
        H = ham_paths_dp(A, n)

        # Count directed 3-cycles (= c3 vertex sets for tournaments)
        c3 = 0
        for i, j, k in combinations(range(n), 3):
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                c3 += 1

        # Count alpha_2 (disjoint 3-cycle pairs)
        cycles_3 = []
        for i, j, k in combinations(range(n), 3):
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                cycles_3.append(frozenset({i, j, k}))
        alpha2_3 = sum(1 for a in range(len(cycles_3))
                       for b in range(a+1, len(cycles_3))
                       if len(cycles_3[a] & cycles_3[b]) == 0)

        # e_cyc: edges in at least one 3-cycle
        cycle_edges = set()
        for c in cycles_3:
            verts = list(c)
            for i_v in range(3):
                for j_v in range(i_v+1, 3):
                    cycle_edges.add(frozenset({verts[i_v], verts[j_v]}))
        e_cyc = len(cycle_edges)

        data.append({
            'bits': bits, 'H': H, 'scores': scores, 'S2': S2,
            'c3': c3, 'alpha2_3': alpha2_3, 'e_cyc': e_cyc,
        })

    H_arr = np.array([d['H'] for d in data], dtype=float)
    S2_arr = np.array([d['S2'] for d in data], dtype=float)
    c3_arr = np.array([d['c3'] for d in data], dtype=float)
    a2_arr = np.array([d['alpha2_3'] for d in data], dtype=float)
    ecyc_arr = np.array([d['e_cyc'] for d in data], dtype=float)

    # Pearson R^2
    r2_s2 = np.corrcoef(S2_arr, H_arr)[0,1]**2
    r2_c3 = np.corrcoef(c3_arr, H_arr)[0,1]**2

    # Spearman rank correlation
    spear_s2 = sp_stats.spearmanr(S2_arr, H_arr)[0]**2
    spear_c3 = sp_stats.spearmanr(c3_arr, H_arr)[0]**2

    print(f"\n  n={n} ({total} tournaments):")
    print(f"    Pearson R^2:  S2->{r2_s2:.6f}, c3->{r2_c3:.6f}")
    print(f"    Spearman R^2: S2->{spear_s2:.6f}, c3->{spear_c3:.6f}")

    # Multiple regression: H ~ S2 + alpha2_3
    X = np.column_stack([np.ones(total), S2_arr, a2_arr])
    coefs = np.linalg.lstsq(X, H_arr, rcond=None)[0]
    H_pred = X @ coefs
    ss_res = np.sum((H_arr - H_pred)**2)
    ss_tot = np.sum((H_arr - np.mean(H_arr))**2)
    r2_multi = 1 - ss_res / ss_tot
    print(f"    R^2(S2, alpha2_3): {r2_multi:.6f}")

    # H ~ S2 + e_cyc
    X2 = np.column_stack([np.ones(total), S2_arr, ecyc_arr])
    coefs2 = np.linalg.lstsq(X2, H_arr, rcond=None)[0]
    H_pred2 = X2 @ coefs2
    ss_res2 = np.sum((H_arr - H_pred2)**2)
    r2_ecyc = 1 - ss_res2 / ss_tot
    print(f"    R^2(S2, e_cyc):    {r2_ecyc:.6f}")

    # H ~ c3 + alpha2_3 (since c3 = f(scores), this tests if alpha2 helps)
    X3 = np.column_stack([np.ones(total), c3_arr, a2_arr])
    coefs3 = np.linalg.lstsq(X3, H_arr, rcond=None)[0]
    H_pred3 = X3 @ coefs3
    ss_res3 = np.sum((H_arr - H_pred3)**2)
    r2_c3_a2 = 1 - ss_res3 / ss_tot
    print(f"    R^2(c3, alpha2_3): {r2_c3_a2:.6f}")

    # What if we use OCF directly? H = 1 + 2*alpha1 + 4*alpha2
    # For n<=6, alpha1 = c3 + c5_directed, alpha2 = pairs of disjoint cycles
    # Actually at n=5: alpha1_total = c3 + directed_5cycles
    # At n=6: more complex. Let's just check 1 + 2*c3 + 4*alpha2_3
    H_ocf_approx = 1 + 2*c3_arr + 4*a2_arr
    r2_ocf_approx = np.corrcoef(H_ocf_approx, H_arr)[0,1]**2
    print(f"    R^2(1+2c3+4alpha2_3, H): {r2_ocf_approx:.6f}")

    # Check: is this exact?
    mismatches = np.sum(H_ocf_approx != H_arr)
    if mismatches == 0:
        print(f"    1+2c3+4alpha2_3 = H EXACTLY (no 5+ cycle contribution at n={n})")
    else:
        print(f"    1+2c3+4alpha2_3 != H for {mismatches}/{total} tournaments")
        diffs = H_arr - H_ocf_approx
        print(f"    Max diff: {np.max(np.abs(diffs)):.0f}, mean: {np.mean(np.abs(diffs)):.2f}")

    # Within the AMBIGUOUS score class, what predicts H?
    score_groups = defaultdict(list)
    for d in data:
        score_groups[d['scores']].append(d)

    ambig = {s: ds for s, ds in score_groups.items() if len(set(d['H'] for d in ds)) > 1}

    if ambig:
        print(f"\n    Within ambiguous score classes (n={n}):")
        for s, ds in sorted(ambig.items()):
            Hs = [d['H'] for d in ds]
            a2s = [d['alpha2_3'] for d in ds]
            ecycs = [d['e_cyc'] for d in ds]

            if len(set(Hs)) > 1 and np.std(a2s) > 0:
                r = np.corrcoef(a2s, Hs)[0,1]**2
                r_e = np.corrcoef(ecycs, Hs)[0,1]**2 if np.std(ecycs) > 0 else 0
                print(f"      scores={s}: H in {sorted(set(Hs))}, "
                      f"R^2(alpha2_3,H)={r:.4f}, R^2(e_cyc,H)={r_e:.4f}")

# ============================================================
# PART 2: REGULAR TOURNAMENT CONCENTRATION
# ============================================================

print("\n" + "=" * 72)
print("  PART 2: Regular tournament concentration")
print("=" * 72)

# For regular tournaments (S2=0), how does H vary?
# At n=5: regular = score (2,2,2,2,2), H always = 15 (no variance!)
# At n=7: regular = score (3,3,3,3,3,3,3), H in {171, 175, 189}

# n=7 regular analysis
print("\n  n=7 regular tournaments (sampled):")
n = 7
m = n*(n-1)//2
np.random.seed(123)
regular_H = []
all_H = []
for _ in range(50000):
    bits = np.random.randint(0, 1 << m)
    bits = int(bits)
    A = build_tournament(n, bits)
    scores = [sum(A[i]) for i in range(n)]
    H = ham_paths_dp(A, n)
    all_H.append(H)
    if sorted(scores) == [3,3,3,3,3,3,3]:
        regular_H.append(H)

print(f"  Regular tournaments found: {len(regular_H)}/50000")
if regular_H:
    print(f"  Regular H values: {sorted(set(regular_H))}")
    print(f"  Regular H mean: {np.mean(regular_H):.2f}, std: {np.std(regular_H):.2f}")
    print(f"  Regular CV: {np.std(regular_H)/np.mean(regular_H):.4f}")
    print(f"  All H mean: {np.mean(all_H):.2f}, std: {np.std(all_H):.2f}")
    print(f"  All CV: {np.std(all_H)/np.mean(all_H):.4f}")
    print(f"  Regular CV / All CV: {(np.std(regular_H)/np.mean(regular_H)) / (np.std(all_H)/np.mean(all_H)):.4f}")

# ============================================================
# PART 3: THEORETICAL BOUND ON OCR
# ============================================================

print("\n" + "=" * 72)
print("  PART 3: Theoretical OCR bound from Rao's formula")
print("=" * 72)

# c3 = C(n,3) - sum(C(s_i, 2)) = C(n,3) - (S2 + n*(n-1)/2*(n-3)/4)/2
# Actually: c3 = C(n,3) - sum s_i*(s_i-1)/2
# sum s_i*(s_i-1)/2 = (sum s_i^2 - sum s_i) / 2
# sum s_i = n(n-1)/2
# sum s_i^2 = S2 + n*((n-1)/2)^2
# So: sum s_i*(s_i-1)/2 = (S2 + n*(n-1)^2/4 - n(n-1)/2) / 2
#                        = S2/2 + n(n-1)(n-3)/8

# H = 1 + 2*c3 + higher terms (for n<=5, higher terms include c5_directed and alpha_2)
# At n<=4: H = 1 + 2*c3 EXACTLY (no 5-cycles possible)
# So at n<=4: H = 1 + 2*(C(n,3) - S2/2 - n(n-1)(n-3)/8)
#            = 1 + 2*C(n,3) - S2 - n(n-1)(n-3)/4

for n in [3, 4, 5, 6]:
    const = 1 + 2*comb(n,3) - n*(n-1)*(n-3)//4
    print(f"\n  n={n}: H_linear = {const} - S2")
    print(f"    C(n,3) = {comb(n,3)}, correction = {n*(n-1)*(n-3)//4}")

    # Verify at small n
    m = n*(n-1)//2
    for bits in range(min(1 << m, 64)):
        A = build_tournament(n, bits)
        scores = [sum(A[i]) for i in range(n)]
        S2 = sum((s - (n-1)/2)**2 for s in scores)
        H = ham_paths_dp(A, n)
        c3 = comb(n,3) - sum(s*(s-1)//2 for s in scores)
        H_pred = 1 + 2*c3
        if n <= 4:
            if H != H_pred:
                print(f"    MISMATCH at bits={bits}: H={H}, 1+2c3={H_pred}")
                break
    else:
        if n <= 4:
            print(f"    VERIFIED: H = 1 + 2*c3 exactly at n={n}")

    # At n>=5, the residual H - (1+2c3) comes from 5-cycles
    if n >= 5:
        data_n = []
        total_n = 1 << m
        for bits in range(total_n):
            A = build_tournament(n, bits)
            scores = [sum(A[i]) for i in range(n)]
            S2 = sum((s - (n-1)/2)**2 for s in scores)
            H = ham_paths_dp(A, n)
            c3 = comb(n,3) - sum(s*(s-1)//2 for s in scores)
            data_n.append((H, c3, S2))

        H_arr = np.array([d[0] for d in data_n], dtype=float)
        c3_arr = np.array([d[1] for d in data_n], dtype=float)
        H_ocf1 = 1 + 2*c3_arr
        residual = H_arr - H_ocf1

        print(f"    Residual H - (1+2c3): mean={np.mean(residual):.2f}, "
              f"std={np.std(residual):.2f}, max={np.max(np.abs(residual)):.0f}")
        print(f"    Var(residual)/Var(H) = {np.var(residual)/np.var(H_arr):.6f}")
        print(f"    This residual fraction = 1 - R^2(c3, H)")

# ============================================================
# PART 4: KEY FORMULA — Why exactly 1/C(n,2)?
# ============================================================

print("\n" + "=" * 72)
print("  PART 4: OCR = 1 - 1/C(n,2) hypothesis")
print("=" * 72)
print()

# OCR values: 1, 1, 0.9474, 0.9231, 0.916
# 1 - OCR:     0, 0, 0.0526, 0.0769, 0.084
# 1/C(n,2):    1/3, 1/6, 1/10, 1/15, 1/21
#            = 0.333, 0.167, 0.100, 0.067, 0.048

# Hmm, doesn't match. Let's check other combinatorial fractions.
# (n-3)/(n-1)^2?
# n=5: 2/16 = 0.125 (too high)
# 2/(n+1)?
# n=5: 2/6 = 0.333 (too high)
# 1/(n-1)?
# n=5: 1/4 = 0.25 (too high)

# What about based on the fraction of ambiguous tournaments?
# At n=5: 280/1024 = 0.2734 are in the ambiguous class
# 1-OCR = 0.0526 = much less than the ambiguous fraction
# Because the ambiguous class has low WITHIN-class variance

for n in [5, 6, 7, 8]:
    if n <= 6:
        # Exact
        m_val = n*(n-1)//2
        data_n = []
        for bits in range(1 << m_val):
            A = build_tournament(n, bits)
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            S2 = sum((s - (n-1)/2)**2 for s in scores)
            H = ham_paths_dp(A, n)
            data_n.append((H, S2, scores))

        Hv = np.array([d[0] for d in data_n], dtype=float)
        S2v = np.array([d[1] for d in data_n], dtype=float)

        ocr = np.corrcoef(S2v, Hv)[0,1]**2
        deficit = 1 - ocr

        # Check ratio with various combinatorial quantities
        print(f"  n={n}: 1-OCR = {deficit:.6f}")
        print(f"    1/C(n,2) = {1/comb(n,2):.6f}, ratio = {deficit*comb(n,2):.4f}")
        print(f"    2/(n^2) = {2/n**2:.6f}, ratio = {deficit*n**2/2:.4f}")
        print(f"    (n-4)/(n*(n-1)) = {max(0,(n-4))/(n*(n-1)):.6f}")

        # The residual variance is Var(H|scores) / Var(H)
        score_groups = defaultdict(list)
        for d in data_n:
            score_groups[d[2]].append(d[0])
        within = sum(np.var(Hs)*len(Hs) for Hs in score_groups.values()) / len(data_n)
        total_var = np.var(Hv)
        print(f"    Within-class fraction = {within/total_var:.6f}")
    print()

# ============================================================
# PART 5: THE 0.916 COINCIDENCE
# ============================================================

print("=" * 72)
print("  PART 5: Is 0.916 a recognizable constant?")
print("=" * 72)
print()

# 0.9161 ~ 11/12 = 0.9167
# 0.9161 ~ (n-1)^2 / (n^2 - n + 1)?
# For n=7: 36/43 = 0.8372 (no)
# 0.9474 at n=5: close to 18/19 = 0.9474!
# 0.9231 at n=6: close to 12/13 = 0.923076...!
# 0.916 at n=7,8: 11/12 = 0.91667, close!

# Check: is OCR(n) = C(n,2)/(C(n,2) + 1)?
# n=5: 10/11 = 0.9091 (no, that's 0.91 not 0.947)
#
# OCR(n=5) = 0.94737 = 18/19 EXACTLY
# OCR(n=6) = 0.92308 = 12/13 EXACTLY
#
# 18 = 2*C(5,2) - 2, 19 = 2*C(5,2) - 1? No: 2*10-2=18, 2*10-1=19. YES!
# 12 = 2*C(6,2) - 18? No. 12/13... what is 12? It's C(4,2)*2 = 12. Hmm.
#
# Actually: 18/19. At n=5: m=10, number of score classes = 9.
# 19 = 2*9 + 1? No, 19 = 2*10 - 1.
# Try: OCR = (2m-2)/(2m-1)? At n=5, m=10: 18/19 = 0.94737. YES!
# At n=6, m=15: 28/29 = 0.96552. NO (claimed OCR = 12/13).

# Let me be more careful:
print("  Looking for exact rational OCR values:")
for n_val in [5, 6]:
    m_val = n_val * (n_val - 1) // 2
    data_n = []
    for bits in range(1 << m_val):
        A = build_tournament(n_val, bits)
        scores = tuple(sorted([sum(A[i]) for i in range(n_val)]))
        S2 = sum((s - (n_val-1)/2)**2 for s in scores)
        H = ham_paths_dp(A, n_val)
        data_n.append((float(H), float(S2)))

    Hv = np.array([d[0] for d in data_n])
    S2v = np.array([d[1] for d in data_n])

    # Compute R^2 via exact formula
    n_pts = len(data_n)
    sum_S2 = np.sum(S2v)
    sum_H = np.sum(Hv)
    sum_S2H = np.sum(S2v * Hv)
    sum_S2_sq = np.sum(S2v**2)
    sum_H_sq = np.sum(Hv**2)

    numer = (n_pts * sum_S2H - sum_S2 * sum_H)**2
    denom = (n_pts * sum_S2_sq - sum_S2**2) * (n_pts * sum_H_sq - sum_H**2)

    print(f"\n  n={n_val}: R^2 = {numer}/{denom} = {numer/denom:.10f}")

    # Try to simplify
    from math import gcd
    g = gcd(int(numer), int(denom))
    print(f"    Simplified: {int(numer)//g} / {int(denom)//g}")

print()
print("  INTERPRETATION:")
print("  If OCR has exact rational values, that constrains the theory.")
print("  OCR(5) = 18/19 and OCR(6) = 12/13 would mean the residual")
print("  fraction is 1/19 and 1/13 — both PRIMES in the denominator.")
print("  Sequence: _, _, 19, 13, ?, ? — not obviously a named sequence.")
