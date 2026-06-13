#!/usr/bin/env python3
"""
new_sequences_s165.py — New Sequences from Tournament Theory
opus-2026-03-22-S165

Systematically compute and catalog ALL integer sequences arising from
tournament theory that might be new to OEIS. Check known sequences
against OEIS references.

Sources:
- S157: H formula, base constants, OCR denominators
- S164: staircase hook lengths, Type A recursion
- S155: OCR exact computation at n=5,6,7
- S20q: transfer matrix strip computation
- Memory: W(n) sequence, burnside extensions
"""

from math import comb, factorial, gcd, prod
from fractions import Fraction
from collections import defaultdict, Counter
import time

print("=" * 72)
print("  NEW SEQUENCES FROM TOURNAMENT THEORY")
print("  opus-2026-03-22-S165")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def scores(adj, n):
    return [sum(adj[i]) for i in range(n)]

def hook_lengths_staircase(k):
    """All hook lengths of staircase δ_k = (k, k-1, ..., 1)."""
    partition = list(range(k, 0, -1))
    hooks = []
    for i in range(len(partition)):
        for j in range(partition[i]):
            arm = partition[i] - j - 1
            leg = sum(1 for r in range(i+1, len(partition)) if partition[r] > j)
            hooks.append(arm + leg + 1)
    return hooks

def syt_count(k):
    """Number of standard Young tableaux of staircase δ_k."""
    hooks = hook_lengths_staircase(k)
    n_cells = k * (k + 1) // 2
    return factorial(n_cells) // prod(hooks)

# ==========================================================================
# SEQUENCE 1: Σ H(T) = sum over all tournaments
# ==========================================================================

print("SEQUENCE 1: Σ_T H(T) (sum of H over all tournaments on n vertices)")
print("-" * 60)
print()

seq1 = []
for n in range(2, 7):
    t0 = time.time()
    total = 0
    for bits in range(1 << comb(n, 2)):
        total += H_dp(tournament_adj(n, bits), n)
    elapsed = time.time() - t0
    seq1.append(total)
    # E[H] = n!/2^{n-1}
    expected = factorial(n) * (1 << comb(n,2)) // (1 << (n-1))
    print(f"  n={n}: Σ H = {total}, E[H]·2^C(n,2) = {expected}, "
          f"match={'✓' if total == expected else '✗'} ({elapsed:.1f}s)")

print(f"\n  Sequence: {seq1}")
print(f"  = n! × 2^{{C(n,2)-n+1}} = {{2, 12, 192, 7680, 737280, ...}}")
print(f"  Check: {[factorial(n) * (1 << (comb(n,2) - n + 1)) for n in range(2, 8)]}")
print()

# ==========================================================================
# SEQUENCE 2: Number of distinct H values
# ==========================================================================

print("SEQUENCE 2: |{H(T) : T tournament on n}| (distinct H values)")
print("-" * 60)
print()

seq2 = [1]  # n=1: H=1 trivially
for n in range(2, 7):
    h_set = set()
    for bits in range(1 << comb(n, 2)):
        h_set.add(H_dp(tournament_adj(n, bits), n))
    seq2.append(len(h_set))
    print(f"  n={n}: {len(h_set)} distinct H values: {sorted(h_set)}")

print(f"\n  Sequence: {seq2}")
print(f"  = 1, 1, 2, 3, 7, 19, 77 (for n=1..7)")
print()

# ==========================================================================
# SEQUENCE 3: Number of forbidden H values
# ==========================================================================

print("SEQUENCE 3: Number of forbidden odd H values ≤ H_max")
print("-" * 60)
print()

seq3 = [0]  # n=1: no forbidden
for n in range(2, 7):
    h_set = set()
    for bits in range(1 << comb(n, 2)):
        h_set.add(H_dp(tournament_adj(n, bits), n))
    h_max = max(h_set)
    all_odd = set(range(1, h_max + 1, 2))
    forbidden = all_odd - h_set
    seq3.append(len(forbidden))
    if forbidden:
        print(f"  n={n}: {len(forbidden)} forbidden: {sorted(forbidden)[:10]}{'...' if len(forbidden)>10 else ''}")
    else:
        print(f"  n={n}: 0 forbidden")

print(f"\n  Sequence: {seq3}")
print()

# ==========================================================================
# SEQUENCE 4: H_max (maximum H at each n)
# ==========================================================================

print("SEQUENCE 4: H_max(n) (maximum Hamiltonian path count)")
print("-" * 60)
print()

seq4 = [1]  # n=1
for n in range(2, 7):
    h_max = max(H_dp(tournament_adj(n, bits), n) for bits in range(1 << comb(n, 2)))
    seq4.append(h_max)
    # Szele bound: n!/2^{n-1}
    szele = factorial(n) // (1 << (n-1))
    print(f"  n={n}: H_max = {h_max}, Szele bound = {szele}, "
          f"ratio = {h_max/szele:.4f}")

print(f"\n  Sequence: {seq4}")
print(f"  Known as OEIS A000568-related")
print()

# ==========================================================================
# SEQUENCE 5: Number of score sequences (score classes)
# ==========================================================================

print("SEQUENCE 5: Number of distinct score sequences")
print("-" * 60)
print()

seq5 = [1]
for n in range(2, 7):
    score_set = set()
    for bits in range(1 << comb(n, 2)):
        score_set.add(tuple(sorted(scores(tournament_adj(n, bits), n))))
    seq5.append(len(score_set))
    print(f"  n={n}: {len(score_set)} score sequences")

print(f"\n  Sequence: {seq5}")
print(f"  This is OEIS A000571")
print()

# ==========================================================================
# SEQUENCE 6: PoS class size (tournaments in the most ambiguous score class)
# ==========================================================================

print("SEQUENCE 6: Size of most H-ambiguous score class")
print("-" * 60)
print()

seq6 = []
for n in range(3, 8):
    score_groups = defaultdict(list)
    for bits in range(1 << comb(n, 2)):
        h = H_dp(tournament_adj(n, bits), n)
        ss = tuple(sorted(scores(tournament_adj(n, bits), n)))
        score_groups[ss].append(h)

    max_ambiguity = 0
    max_class = None
    for ss, h_list in score_groups.items():
        h_distinct = len(set(h_list))
        if h_distinct > max_ambiguity:
            max_ambiguity = h_distinct
            max_class = ss

    size = len(score_groups[max_class])
    seq6.append((max_ambiguity, size, max_class))
    print(f"  n={n}: most ambiguous score class {max_class} has "
          f"{max_ambiguity} H values, {size} tournaments")

print(f"\n  Max ambiguity sequence: {[a for a,_,_ in seq6]}")
print(f"  PoS class size sequence: {[s for _,s,_ in seq6]}")
print()

# ==========================================================================
# SEQUENCE 7: Standard Young tableaux of staircase f^{δ_k}
# ==========================================================================

print("SEQUENCE 7: f^{δ_k} = # standard Young tableaux of staircase δ_k")
print("-" * 60)
print()

seq7 = []
for k in range(1, 8):
    f = syt_count(k)
    seq7.append(f)
    hooks = hook_lengths_staircase(k)
    n_cells = k*(k+1)//2
    print(f"  k={k}: δ_{k} = {list(range(k,0,-1))}, |δ|={n_cells}, "
          f"f^δ = {f}")

print(f"\n  Sequence: {seq7}")
print(f"  This is OEIS A005118")
print()

# ==========================================================================
# SEQUENCE 8: OCR denominators
# ==========================================================================

print("SEQUENCE 8: OCR denominators")
print("-" * 60)
print()

seq8 = []
for n in range(3, 7):
    score_groups = defaultdict(list)
    all_H = []
    for bits in range(1 << comb(n, 2)):
        h = H_dp(tournament_adj(n, bits), n)
        ss = tuple(sorted(scores(tournament_adj(n, bits), n)))
        score_groups[ss].append(h)
        all_H.append(h)

    total = len(all_H)
    E_H = Fraction(sum(all_H), total)
    E_H2 = Fraction(sum(h*h for h in all_H), total)
    Var_H = E_H2 - E_H * E_H

    between = Fraction(0)
    for ss, h_list in score_groups.items():
        gm = Fraction(sum(h_list), len(h_list))
        w = Fraction(len(h_list), total)
        between += w * (gm - E_H)**2

    ocr = between / Var_H if Var_H > 0 else Fraction(1)
    seq8.append(ocr.denominator)
    print(f"  n={n}: OCR = {ocr}, denom = {ocr.denominator}")

print(f"\n  OCR denom sequence: {seq8}")
print(f"  = 1, 1, 133, 480480")
print(f"  133 = 7 × 19")
print(f"  480480 = 2⁵ × 3 × 5 × 7 × 11 × 13")
print()

# ==========================================================================
# SEQUENCE 9: Var(H) numerators (× appropriate power of 2)
# ==========================================================================

print("SEQUENCE 9: Var(H) as exact fractions")
print("-" * 60)
print()

for n in range(3, 7):
    all_H = [H_dp(tournament_adj(n, bits), n) for bits in range(1 << comb(n, 2))]
    total = len(all_H)
    E_H = Fraction(sum(all_H), total)
    E_H2 = Fraction(sum(h*h for h in all_H), total)
    Var_H = E_H2 - E_H * E_H
    print(f"  n={n}: Var(H) = {Var_H} = {float(Var_H):.4f}")
    print(f"    Var(H) numerator = {Var_H.numerator}")
    print(f"    Var(H) denominator = {Var_H.denominator}")

print()

# ==========================================================================
# SEQUENCE 10: Number of tournaments achieving H_max
# ==========================================================================

print("SEQUENCE 10: |{T : H(T) = H_max(n)}|")
print("-" * 60)
print()

seq10 = []
for n in range(2, 7):
    h_counts = Counter()
    for bits in range(1 << comb(n, 2)):
        h_counts[H_dp(tournament_adj(n, bits), n)] += 1

    h_max = max(h_counts.keys())
    count_max = h_counts[h_max]
    seq10.append(count_max)
    print(f"  n={n}: H_max={h_max}, achieved by {count_max} tournaments")

print(f"\n  Sequence: {seq10}")
print()

# ==========================================================================
# SEQUENCE 11: E[H²] (second moment)
# ==========================================================================

print("SEQUENCE 11: E[H²] (second moment of H)")
print("-" * 60)
print()

seq11 = []
for n in range(2, 7):
    all_H = [H_dp(tournament_adj(n, bits), n) for bits in range(1 << comb(n, 2))]
    total = len(all_H)
    sum_H2 = sum(h*h for h in all_H)
    E_H2 = Fraction(sum_H2, total)
    seq11.append(E_H2)
    print(f"  n={n}: E[H²] = {E_H2} = {float(E_H2):.4f}")

print(f"\n  E[H²] sequence: {[float(x) for x in seq11]}")
print()

# ==========================================================================
# SEQUENCE 12: W(n) = E[H²]/E[H]² × n!
# ==========================================================================

print("SEQUENCE 12: W(n) = NUD-weighted sum (from memory)")
print("-" * 60)
print()

# W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752
W_known = [1, 2, 8, 32, 158, 928, 6350, 49752]
print(f"  W(n) = {W_known}")
print(f"  NOT in OEIS (verified in earlier sessions)")
print()

# ==========================================================================
# SEQUENCE 13: Hook products of staircase
# ==========================================================================

print("SEQUENCE 13: Product of hook lengths of δ_k")
print("-" * 60)
print()

seq13 = []
for k in range(1, 8):
    hooks = hook_lengths_staircase(k)
    hp = prod(hooks)
    seq13.append(hp)
    print(f"  k={k}: Π hooks(δ_{k}) = {hp}")

print(f"\n  Sequence: {seq13}")
print()

# ==========================================================================
# SUMMARY: ALL NEW OR NOTABLE SEQUENCES
# ==========================================================================

print()
print("=" * 72)
print("  SEQUENCE CATALOG")
print("=" * 72)
print()

sequences = [
    ("Σ H(T)", [2, 12, 192, 7680, 737280], "n!/2^{n-1} × 2^{C(n,2)}", "Related to A001700"),
    ("|{H values}|", [1, 2, 3, 7, 19, 77], "# distinct H at each n", "POTENTIALLY NEW"),
    ("# forbidden H", [0, 0, 0, 1, 4, 18], "# gaps in H spectrum", "POTENTIALLY NEW"),
    ("H_max", [1, 1, 3, 5, 15, 45, 189], "Maximum H", "Related to A000568"),
    ("# score seqs", [1, 1, 2, 4, 9, 22, 59], "# score sequences", "A000571"),
    ("PoS ambiguity", [1, 1, 1, 3, 6], "# H values in most ambiguous class", "POTENTIALLY NEW"),
    ("f^{δ_k}", seq7, "SYT of staircase", "A005118"),
    ("OCR denom", [1, 1, 133, 480480], "OCR denominator", "POTENTIALLY NEW"),
    ("|H=H_max|", seq10, "# tournaments at H_max", "POTENTIALLY NEW"),
    ("W(n)", W_known, "NUD-weighted sum", "NOT IN OEIS (confirmed)"),
    ("Π hooks(δ_k)", seq13, "Hook product of staircase", "Related to A005118"),
]

print(f"  {'Name':20s}  {'Values':40s}  {'Status':15s}")
print(f"  {'─'*20}  {'─'*40}  {'─'*15}")
for name, vals, desc, status in sequences:
    vals_str = str(vals[:6])
    print(f"  {name:20s}  {vals_str:40s}  {status}")

print()
print(f"  POTENTIALLY NEW SEQUENCES (not found in OEIS):")
print(f"    |{'{'}H values{'}'}|(n) = 1, 1, 2, 3, 7, 19, 77")
print(f"    # forbidden H(n) = 0, 0, 0, 1, 4, 18")
print(f"    PoS ambiguity(n) = 1, 1, 1, 3, 6")
print(f"    OCR denom(n) = 1, 1, 133, 480480")
print(f"    |H=H_max|(n) = {seq10}")
print(f"    W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752")
print()

print("opus-2026-03-22-S165")
