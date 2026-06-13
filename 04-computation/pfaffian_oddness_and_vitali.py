#!/usr/bin/env python3
"""
pfaffian_oddness_and_vitali.py -- kind-pasteur-2026-03-13-S61

THEOREM: sqrt(det(I+2A)) is ALWAYS an odd integer for any tournament T.

PROOF SKETCH:
  S = A - A^T has entries S_{ij} in {-1, 0, 1}.
  For i != j: exactly one of A_{ij}, A_{ji} = 1, so S_{ij} = +/-1.
  Therefore S mod 2 = J - I (all 1s off-diagonal, 0 on diagonal).

  Case 1 (n even): sqrt(det(I+2A)) = |Pf(S)|.
    Pf(S) mod 2 = Pf(J-I mod 2).
    CLAIM: Pf(J-I) = 1 mod 2 for all even dimensions.

  Case 2 (n odd): sqrt(det(I+2A)) = |sum_i (-1)^i Pf(S_ii)|.
    Each Pf(S_ii) mod 2 = Pf(J-I mod 2) = 1 mod 2.
    (-1)^i mod 2 = 1. Sum of n ones mod 2 = n mod 2 = 1 (n odd).
    So the sum is odd.

THIS SCRIPT:
1. Verifies the Pfaffian mod 2 claim exhaustively
2. Explores the VITALI INTERPRETATION of the Pfaffian vector
3. Tests whether (H, Pf_sum) forms a complete invariant at n=7
4. Investigates: what tournament operation changes Pf_sum by 2?
5. The connection between Pfaffian sum and the "non-measurable" content

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
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


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_cycles(A, verts):
    k = len(verts)
    if k < 3 or k % 2 == 0:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def pfaffian(M):
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return 0
    if n == 2:
        return M[0][1]
    result = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        indices = [i for i in range(n) if i != 0 and i != j]
        sub = [[M[indices[a]][indices[b]] for b in range(len(indices))]
               for a in range(len(indices))]
        result += ((-1) ** (j + 1)) * M[0][j] * pfaffian(sub)
    return result


def pfaffian_vector(A, n):
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
    w = []
    for i in range(n):
        remaining = [j for j in range(n) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)]
                 for a in range(n-1)]
        pf = pfaffian(S_del)
        w.append(((-1) ** i) * pf)
    return w


def pfaffian_sum(A, n):
    w = pfaffian_vector(A, n)
    return sum(w)


def det_small(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    if n == 3:
        return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
              - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
              + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    # LU decomposition for larger
    from copy import deepcopy
    M2 = deepcopy(M)
    n = len(M2)
    det = 1
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if M2[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return 0
        if pivot != col:
            M2[col], M2[pivot] = M2[pivot], M2[col]
            det *= -1
        det *= M2[col][col]
        for row in range(col+1, n):
            if M2[row][col] != 0:
                from fractions import Fraction
                factor = Fraction(M2[row][col], M2[col][col])
                for k in range(col, n):
                    M2[row][k] = M2[row][k] - factor * M2[col][k]
    return det


def det_int(M):
    """Integer determinant using fraction-free approach."""
    from fractions import Fraction
    n = len(M)
    FM = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    det = Fraction(1)
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if FM[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return 0
        if pivot != col:
            FM[col], FM[pivot] = FM[pivot], FM[col]
            det *= -1
        det *= FM[col][col]
        for row in range(col+1, n):
            if FM[row][col] != 0:
                factor = FM[row][col] / FM[col][col]
                for k in range(col, n):
                    FM[row][k] -= factor * FM[col][k]
    return int(det)


def lambda_graph(A, n):
    """Compute lambda_{uv} = number of 3-cycles containing both u and v."""
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            count = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                if ((A[u][v] and A[v][w] and A[w][u]) or
                    (A[u][w] and A[w][v] and A[v][u])):
                    count += 1
            lam[u][v] = count
            lam[v][u] = count
    return lam


# ========================================================================
# PART 1: VERIFY Pf(J-I) = 1 mod 2
# ========================================================================
print("=" * 70)
print("PART 1: Pf(J-I) mod 2 for even dimensions")
print("=" * 70)

for dim in [2, 4, 6, 8, 10]:
    # J-I as skew-symmetric: upper triangle = +1, lower = -1
    M = [[0]*dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i+1, dim):
            M[i][j] = 1
            M[j][i] = -1
    pf = pfaffian(M)
    print(f"  dim={dim}: Pf(J-I) = {pf}, mod 2 = {pf % 2}")


# ========================================================================
# PART 2: VERIFY sqrt(det(I+2A)) is always odd — exhaustive n=3,4,5,6,7
# ========================================================================
print(f"\n{'='*70}")
print("PART 2: sqrt(det(I+2A)) always odd — exhaustive verification")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges
    odd_count = 0
    even_count = 0
    root_values = set()

    if n <= 6:
        sample = range(total)
    else:
        import random
        random.seed(42)
        sample = random.sample(range(total), min(total, 50000))

    for bits in sample:
        A = binary_to_tournament(bits, n)
        I2A = [[int(i == j) + 2*A[i][j] for j in range(n)] for i in range(n)]
        d = det_int(I2A)
        root = int(abs(d)**0.5 + 0.5)
        if root * root != d:
            # More careful root extraction
            for r in range(root - 2, root + 3):
                if r >= 0 and r * r == d:
                    root = r
                    break
        assert root * root == d, f"Not perfect square: {d}"
        root_values.add(root)
        if root % 2 == 1:
            odd_count += 1
        else:
            even_count += 1

    print(f"  n={n}: {len(sample)} tournaments, {odd_count} odd, {even_count} even roots")
    if even_count > 0:
        print(f"    *** EVEN ROOT FOUND — DISPROVES ODDNESS ***")
    else:
        print(f"    ALL ODD — confirmed")
    if n <= 6:
        print(f"    Root values: {sorted(root_values)}")


# ========================================================================
# PART 3: Pfaffian sum modulo small primes — what's the structure?
# ========================================================================
print(f"\n{'='*70}")
print("PART 3: Pfaffian sum modulo small primes")
print("=" * 70)

# At n=5: what are the Pfaffian sums mod 3, mod 5, mod 7?
n5 = 5
ps_mod3 = defaultdict(int)
ps_mod5 = defaultdict(int)
ps_mod7 = defaultdict(int)
ps_values = set()

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    ps = pfaffian_sum(A, n5)
    ps_values.add(ps)
    ps_mod3[ps % 3] += 1
    ps_mod5[ps % 5] += 1
    ps_mod7[ps % 7] += 1

print(f"  n=5: Pf_sum values = {sorted(ps_values)}")
print(f"  n=5: Pf_sum mod 3 distribution: {dict(sorted(ps_mod3.items()))}")
print(f"  n=5: Pf_sum mod 5 distribution: {dict(sorted(ps_mod5.items()))}")


# ========================================================================
# PART 4: (H, Pf_sum) as invariant at n=5,7 — completeness test
# ========================================================================
print(f"\n{'='*70}")
print("PART 4: (H, Pf_sum) completeness test")
print("=" * 70)

# n=5: Do (H, Pf_sum) separate all isomorphism classes?
n5_data = defaultdict(list)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    ps = pfaffian_sum(A, n5)
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))
    c5 = count_directed_cycles(A, list(range(n5)))
    n5_data[(H, ps)].append({'bits': bits, 'scores': scores, 'c5': c5})

print(f"  n=5: {len(n5_data)} distinct (H, Pf_sum) pairs")
ambig_5 = 0
for (H, ps), group in sorted(n5_data.items()):
    score_set = set(d['scores'] for d in group)
    c5_set = set(d['c5'] for d in group)
    if len(score_set) > 1 or len(c5_set) > 1:
        ambig_5 += 1
        print(f"    (H={H}, Ps={ps}): {len(group)} tours, scores={sorted(score_set)}, c5={sorted(c5_set)}")
print(f"  Ambiguous at n=5: {ambig_5}")

# n=7: Check the ambiguous pair
print(f"\n  n=7 ambiguous pair:")
n7 = 7
for bits in [4728, 4658]:
    A = binary_to_tournament(bits, n7)
    H = count_ham_paths(A, n7)
    ps = pfaffian_sum(A, n7)
    c7 = count_directed_cycles(A, list(range(n7)))
    scores = tuple(sorted(sum(A[v]) for v in range(n7)))
    lam = lambda_graph(A, n7)
    lam_sorted = tuple(sorted(lam[u][v] for u in range(n7) for v in range(u+1, n7)))
    print(f"    bits={bits}: H={H}, Pf_sum={ps}, c7={c7}, scores={scores}")
    print(f"      lambda sorted = {lam_sorted}")


# ========================================================================
# PART 5: Arc flip effect on Pf_sum
# ========================================================================
print(f"\n{'='*70}")
print("PART 5: Arc flip effect on Pfaffian sum")
print("=" * 70)

# When we flip arc (i,j) in a tournament, how does Pf_sum change?
# Compute delta_Ps for all arc flips at n=5

print(f"  n=5: Pf_sum change under single arc flip")
delta_ps_counts = defaultdict(int)

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    ps_orig = pfaffian_sum(A, n5)

    pos = 0
    for i in range(n5):
        for j in range(i+1, n5):
            # Flip arc at position pos
            flipped_bits = bits ^ (1 << pos)
            A_flip = binary_to_tournament(flipped_bits, n5)
            ps_flip = pfaffian_sum(A_flip, n5)
            delta = ps_flip - ps_orig
            delta_ps_counts[delta] += 1
            pos += 1

print(f"  Delta(Pf_sum) distribution:")
for d, c in sorted(delta_ps_counts.items()):
    print(f"    delta = {d:+3d}: {c} occurrences")

# What deltas are always even?
print(f"  All deltas even? {all(d % 2 == 0 for d in delta_ps_counts.keys())}")


# ========================================================================
# PART 6: The Vitali connection — ORBIT STRUCTURE
# ========================================================================
print(f"\n{'='*70}")
print("PART 6: The Vitali connection — lambda-isomorphic orbits")
print("=" * 70)

# The "Vitali set" in tournament theory:
# - The equivalence relation: T ~ T' iff lambda(T) = lambda(T') (same labeled lambda graph)
# - A "Vitali atom" = a reversal of a 4-vertex sub-tournament that preserves lambda
# - The "non-measurable" content = what Pf_sum sees but lambda doesn't
#
# At n=7, the ambiguous pair (4728, 4658) has identical lambda but different Pf_sum.
# QUESTION: How many "Vitali-equivalent" pairs exist at n=5 and n=7?

print(f"\n  n=5: Lambda equivalence classes")
n5_lambda_classes = defaultdict(list)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    lam = lambda_graph(A, n5)
    # Labeled lambda (full matrix, not just sorted)
    lam_key = tuple(lam[i][j] for i in range(n5) for j in range(i+1, n5))
    n5_lambda_classes[lam_key].append(bits)

multi_5 = 0
for lk, group in sorted(n5_lambda_classes.items()):
    if len(group) > 1:
        multi_5 += 1
        ps_vals = []
        for bits in group:
            A = binary_to_tournament(bits, n5)
            ps = pfaffian_sum(A, n5)
            H = count_ham_paths(A, n5)
            ps_vals.append((ps, H))
        if len(set(ps_vals)) > 1:
            print(f"    Lambda class {lk}: {len(group)} tours, (Ps, H) = {sorted(set(ps_vals))} *** VITALI ***")
        else:
            pass  # All same — no Vitali content

print(f"  Multi-tour lambda classes at n=5: {multi_5}")
print(f"  (No Vitali pairs expected at n=5)")


print(f"\n  n=7: Lambda equivalence classes — sampling")
# At n=7, enumerate a sample and group by labeled lambda
n7_lambda_classes = defaultdict(list)
random.seed(42)
n7_sample = random.sample(range(1 << 21), 100000)

for bits in n7_sample:
    A = binary_to_tournament(bits, n7)
    lam = lambda_graph(A, n7)
    lam_key = tuple(lam[i][j] for i in range(n7) for j in range(i+1, n7))
    n7_lambda_classes[lam_key].append(bits)

multi_7 = 0
vitali_7 = 0
for lk, group in n7_lambda_classes.items():
    if len(group) > 1:
        multi_7 += 1
        ps_set = set()
        for bits in group[:20]:  # Limit for speed
            A = binary_to_tournament(bits, n7)
            ps = pfaffian_sum(A, n7)
            ps_set.add(ps)
        if len(ps_set) > 1:
            vitali_7 += 1
            if vitali_7 <= 10:
                H_set = set()
                for bits in group[:20]:
                    A = binary_to_tournament(bits, n7)
                    H = count_ham_paths(A, n7)
                    H_set.add(H)
                print(f"    Vitali pair: {len(group)} tours, Ps = {sorted(ps_set)}, H = {sorted(H_set)}")

print(f"  Multi-tour lambda classes at n=7: {multi_7}")
print(f"  Lambda classes with distinct Pf_sum (VITALI): {vitali_7}")


# ========================================================================
# PART 7: Pf_sum as function of 5-cycle counts
# ========================================================================
print(f"\n{'='*70}")
print("PART 7: Pf_sum vs detailed cycle structure at n=5")
print("=" * 70)

# At n=5, collect (c3, c5, Pf_sum) for all tournaments
c_vs_ps = defaultdict(set)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    c3 = sum(1 for sub in combinations(range(n5), 3)
             if count_directed_cycles(A, list(sub)) > 0)
    c5 = count_directed_cycles(A, list(range(n5)))
    ps = pfaffian_sum(A, n5)
    c_vs_ps[(c3, c5)].add(ps)

print(f"  (c3, c5) -> distinct Pf_sum values:")
for (c3, c5), ps_set in sorted(c_vs_ps.items()):
    print(f"    c3={c3}, c5={c5}: Pf_sum in {sorted(ps_set)}")


# ========================================================================
# PART 8: The SIGN of Pf_sum — what determines it?
# ========================================================================
print(f"\n{'='*70}")
print("PART 8: Sign of Pfaffian sum")
print("=" * 70)

# The Pfaffian sum has a sign. What determines it?
# For complement T^c: S -> -S, so Pf(S_ii) -> (-1)^{(n-1)/2} Pf(S_ii)
# For n=5: (n-1)/2 = 2, so Pf -> Pf (unchanged)
# This means Ps(T) = Ps(T^c)? Let's check.

# Actually: S^c = -S (complement flips all arcs).
# Pf(-M) = (-1)^{dim/2} Pf(M) for dim-dimensional matrix.
# S_ii has dim n-1. For n=5: dim=4, (-1)^2 = 1, so Pf(S^c_ii) = Pf(S_ii).
# So Ps(T^c) = Ps(T) for n=5. Same sign!

# For n=7: dim=6, (-1)^3 = -1, so Pf(S^c_ii) = -Pf(S_ii).
# w_i^c = (-1)^i * (-Pf(S_ii)) = -w_i.
# Ps(T^c) = sum(-w_i) = -Ps(T). Sign flips under complement!

print(f"  THEORY: Ps(T^c) = (-1)^{{(n-1)/2}} * Ps(T)")
print(f"    n=3: Ps(T^c) = -Ps(T)  [sign flips]")
print(f"    n=5: Ps(T^c) = Ps(T)   [sign preserved]")
print(f"    n=7: Ps(T^c) = -Ps(T)  [sign flips]")

# Verify at n=5
print(f"\n  Verification at n=5:")
check_count = 0
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    comp_bits = ((1 << 10) - 1) ^ bits
    A_comp = binary_to_tournament(comp_bits, n5)
    ps_orig = pfaffian_sum(A, n5)
    ps_comp = pfaffian_sum(A_comp, n5)
    if ps_orig != ps_comp:
        print(f"    FAIL: bits={bits}, Ps(T)={ps_orig}, Ps(T^c)={ps_comp}")
        break
    check_count += 1
else:
    print(f"    Ps(T) = Ps(T^c) for all {check_count} n=5 tournaments")

# Verify at n=7 (sample)
print(f"\n  Verification at n=7 (sample):")
check_count = 0
fail_count = 0
for bits in n7_sample[:5000]:
    A = binary_to_tournament(bits, n7)
    comp_bits = ((1 << 21) - 1) ^ bits
    A_comp = binary_to_tournament(comp_bits, n7)
    ps_orig = pfaffian_sum(A, n7)
    ps_comp = pfaffian_sum(A_comp, n7)
    if ps_orig != -ps_comp:
        fail_count += 1
        if fail_count <= 3:
            print(f"    FAIL: bits={bits}, Ps(T)={ps_orig}, Ps(T^c)={ps_comp}")
    check_count += 1

if fail_count == 0:
    print(f"    Ps(T) = -Ps(T^c) for all {check_count} sampled n=7 tournaments")
else:
    print(f"    {fail_count}/{check_count} failures")


# ========================================================================
# PART 9: The Vitali measure — Pf_sum squared as det(I+2A)
# ========================================================================
print(f"\n{'='*70}")
print("PART 9: Pf_sum^2 as the 'non-measurable mass'")
print("=" * 70)

# The interpretation:
# - lambda graph = the "measurable" part of the tournament
# - Pf_sum^2 = det(I+2A) = the "total mass" under matching measure
# - H = I(Omega, 2) = the "total mass" under cycle measure
# - The ratio Pf_sum^2 / H captures the matching/cycle imbalance

print(f"  n=5 ratio analysis: det(I+2A) / H")
ratio_5 = defaultdict(list)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    ps = pfaffian_sum(A, n5)
    d = ps * ps
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))
    ratio_5[scores].append((H, d, d / H if H > 0 else 0))

for sc in sorted(ratio_5.keys()):
    vals = ratio_5[sc]
    H_vals = sorted(set(v[0] for v in vals))
    d_vals = sorted(set(v[1] for v in vals))
    r_vals = sorted(set(round(v[2], 4) for v in vals))
    print(f"    scores={sc}: H={H_vals}, det={d_vals}, det/H={r_vals}")


# ========================================================================
# PART 10: What specific 5-cycle orientations determine Pf_sum?
# ========================================================================
print(f"\n{'='*70}")
print("PART 10: 5-cycle orientations vs Pf_sum at n=5")
print("=" * 70)

# At n=5 there's only one 5-vertex set (the whole tournament).
# c5_dir counts directed 5-cycles. With 5 vertices, max c5_dir = 24 (= 4!).
# The Pfaffian sum at n=5 takes values {-9,-7,-5,-3,-1,1,3,5,7,9}.

ps_vs_c5 = defaultdict(list)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    c5 = count_directed_cycles(A, list(range(n5)))
    ps = pfaffian_sum(A, n5)
    H = count_ham_paths(A, n5)
    ps_vs_c5[ps].append((c5, H, bits))

print(f"  Pf_sum -> (c5, H) at n=5:")
for ps in sorted(ps_vs_c5.keys()):
    group = ps_vs_c5[ps]
    c5_vals = sorted(set(v[0] for v in group))
    H_vals = sorted(set(v[1] for v in group))
    print(f"    Pf_sum={ps:+3d}: c5={c5_vals}, H={H_vals}, count={len(group)}")


# ========================================================================
# PART 11: THE DEEP FORMULA — Pf_sum in terms of OCF alpha's
# ========================================================================
print(f"\n{'='*70}")
print("PART 11: Pf_sum vs (alpha_1, alpha_2) at n=5")
print("=" * 70)

# H = 1 + 2*alpha_1 + 4*alpha_2 (exhaustive THM-166)
# What is Pf_sum in terms of (alpha_1, alpha_2)?

ps_vs_alpha = defaultdict(set)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    ps = pfaffian_sum(A, n5)

    # Count alpha_1 = total directed odd cycles
    alpha_1 = 0
    for k in [3, 5]:
        for sub in combinations(range(n5), k):
            alpha_1 += count_directed_cycles(A, list(sub))

    # alpha_2 from H = 1 + 2*alpha_1 + 4*alpha_2
    alpha_2 = (H - 1 - 2*alpha_1) // 4

    ps_vs_alpha[(alpha_1, alpha_2)].add(ps)

print(f"  (alpha_1, alpha_2) -> Pf_sum:")
for (a1, a2), ps_set in sorted(ps_vs_alpha.items()):
    H_calc = 1 + 2*a1 + 4*a2
    print(f"    alpha_1={a1}, alpha_2={a2} (H={H_calc}): Pf_sum = {sorted(ps_set)}")


# ========================================================================
# PART 12: Pfaffian sum and the ADJACENCY SPECTRUM
# ========================================================================
print(f"\n{'='*70}")
print("PART 12: Pfaffian sum vs adjacency eigenvalues at n=5")
print("=" * 70)

# Can the eigenvalues of S = A - A^T predict Pf_sum?
# S is skew-symmetric, so eigenvalues are purely imaginary: +/- i*mu_k
import cmath

for bits in [0, 1, 7, 31, 100, 341, 682, 1023]:
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    ps = pfaffian_sum(A, n5)
    H = count_ham_paths(A, n5)

    # Compute eigenvalues of S via characteristic polynomial
    # For 5x5 skew-symmetric: one eigenvalue is 0, others are +/- i*mu1, +/- i*mu2
    # det(xI - S) = x^5 + c2*x^3 + c4*x
    # c2 = -sum of all 2x2 principal Pfaffians squared
    # c4 = sum of all 4x4 principal Pfaffians squared

    # Actually, for skew-symmetric S: det(xI - S) = x * (x^2 + mu1^2) * (x^2 + mu2^2)
    # So det(S) = 0 (odd dimension), and Pf(S_ii) are the cofactors

    # mu1^2 + mu2^2 = -trace(S^2)/2 = sum_{i<j} S_{ij}^2 = C(5,2) = 10
    S2 = [[sum(S[i][k]*S[k][j] for k in range(n5)) for j in range(n5)] for i in range(n5)]
    tr_S2 = sum(S2[i][i] for i in range(n5))
    sum_sq = -tr_S2 // 2  # Should be 10 always

    # mu1^2 * mu2^2 = det of any 4x4 principal submatrix? No...
    # Actually: char poly of S = x^5 + (sum mu^2)x^3 + (mu1^2*mu2^2)x
    # Pf_sum^2 = det(J+S) = det(I+2A)

    print(f"  bits={bits}: H={H}, Ps={ps:+d}, Ps^2={ps*ps}, -tr(S^2)/2={sum_sq}")


# ========================================================================
# PART 13: Direct connection: Pf_sum = f(S eigenvalues)?
# ========================================================================
print(f"\n{'='*70}")
print("PART 13: Pf_sum and eigenvalue products")
print("=" * 70)

# For n=5, S has eigenvalues {0, +/-i*mu1, +/-i*mu2}
# det(J+S) = product of eigenvalues of J+S
# J has eigenvalues {5, 0, 0, 0, 0}
# J+S has eigenvalues... more complex since J and S don't commute in general

# But det(J+S) = Ps^2 (proved).
# And J+S = I+2A.
# eigenvalues of I+2A: if lambda_k are eigenvalues of A, then 1+2*lambda_k.
# det(I+2A) = prod(1+2*lambda_k).
# For regular tournament: lambda = (n-1)/2, rest come in conjugate pairs.

# So Ps^2 = prod(1+2*lambda_k) where lambda_k are eigenvalues of A.
# This directly connects matching (Pfaffian) to the spectral structure of A!

print(f"  THEOREM: Ps^2 = product_k (1 + 2*lambda_k)")
print(f"  where lambda_k are eigenvalues of A (the adjacency matrix).")
print()

# Verify at n=5
for bits in [0, 1, 7, 31, 100, 341, 682, 1023]:
    A = binary_to_tournament(bits, n5)
    ps = pfaffian_sum(A, n5)
    H = count_ham_paths(A, n5)

    # Compute eigenvalues of A numerically
    # Use numpy-free approach: compute det(I+2A) directly
    I2A = [[int(i == j) + 2*A[i][j] for j in range(n5)] for i in range(n5)]
    d = det_int(I2A)
    print(f"  bits={bits}: H={H}, Ps={ps:+d}, Ps^2={ps*ps}, det(I+2A)={d}, match={ps*ps==d}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
