#!/usr/bin/env python3
"""
pi_moments_89c.py — Moments of H and the E[H²] = 4/3 · E[H]² mystery

opus-2026-03-14-S89c

Key findings to investigate:
1. E[H²]/E[H]² = 4/3 EXACTLY at n=3,4. Why? And what is the exact formula?
2. E[1/H] = 5/6, 8/15, 8683/34320 — do these simplify further?
3. The pair-counting formula for E[H²] gave wrong answers — fix it
4. Var(H)/E[H]² seems to stabilize — is CV² → constant?
"""

import math
import itertools
from fractions import Fraction
from collections import Counter

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: Exact Moments of H via Fractions")
print("=" * 70)

for n in range(3, 7):
    h_values = []
    for adj in all_tournaments(n):
        h_values.append(count_H(adj, n))

    N = len(h_values)
    m = n*(n-1)//2

    # Exact moments as fractions
    E_H = Fraction(sum(h_values), N)
    E_H2 = Fraction(sum(h**2 for h in h_values), N)
    E_H3 = Fraction(sum(h**3 for h in h_values), N)
    E_H4 = Fraction(sum(h**4 for h in h_values), N)
    Var_H = E_H2 - E_H**2
    E_inv_H = sum(Fraction(1, h) for h in h_values) / N
    E_inv_H2 = sum(Fraction(1, h**2) for h in h_values) / N

    print(f"\n  n={n} (N=2^{m}={N}):")
    print(f"    E[H]   = {E_H} = {float(E_H):.4f}")
    print(f"    E[H²]  = {E_H2} = {float(E_H2):.4f}")
    print(f"    Var[H] = {Var_H} = {float(Var_H):.4f}")
    print(f"    E[H²]/E[H]² = {E_H2/E_H**2} = {float(E_H2/E_H**2):.6f}")
    print(f"    E[1/H] = {E_inv_H} = {float(E_inv_H):.8f}")
    print(f"    E[1/H²]= {E_inv_H2} = {float(E_inv_H2):.8f}")
    print(f"    E[H]·E[1/H] = {E_H * E_inv_H} = {float(E_H * E_inv_H):.6f}")

    # Check: is E[H²]/E[H]² a simple fraction?
    ratio = E_H2 / E_H**2
    print(f"    Ratio E[H²]/E[H]² = {ratio}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: E[H²] via Double Counting — CORRECTED")
print("=" * 70)

# The CORRECT formula: For two permutations π, σ to both be Hamiltonian
# paths of the SAME tournament T:
# - π requires arcs: π(0)→π(1), π(1)→π(2), ..., π(n-2)→π(n-1)
# - σ requires arcs: σ(0)→σ(1), σ(1)→σ(2), ..., σ(n-2)→σ(n-1)
#
# For a random tournament, each arc is decided independently.
# The pair (π,σ) constrains some arcs:
# - "agreed" arcs: both π and σ require u→v (these are the shared arcs)
# - "π-only" arcs: π requires u→v but σ doesn't constrain (u,v)
# - "σ-only" arcs: σ requires u→v but π doesn't constrain (u,v)
# - "conflicting" arcs: π requires u→v but σ requires v→u
#
# If there are any conflicts, P(both HP) = 0
# If no conflicts, P(both HP) = (1/2)^{total constrained arcs}
# where total = #{agreed} + #{π-only} + #{σ-only}

# Let's also think of it differently:
# Σ_T H(T)² = Σ_T (Σ_π 1_{π HP})(Σ_σ 1_{σ HP})
#            = Σ_{π,σ} Σ_T 1_{π and σ both HP of T}
#            = Σ_{π,σ} 2^{m - c(π,σ)} if compatible, else 0
# where c(π,σ) = number of distinct UNDIRECTED pairs {u,v} constrained

print("\nE[H²] via corrected pair counting:")
for n in range(3, 8):
    perms = list(itertools.permutations(range(n)))
    m = n*(n-1)//2

    total_sum = 0
    n_compatible = 0
    n_incompatible = 0

    for pi in perms:
        # Build arc set for π
        pi_arcs = {}  # {frozenset(u,v): direction}
        for i in range(n-1):
            u, v = pi[i], pi[i+1]
            pair = frozenset([u, v])
            pi_arcs[pair] = (u, v)  # arc goes u→v

        for sigma in perms:
            # Build arc set for σ
            sigma_arcs = {}
            for i in range(n-1):
                u, v = sigma[i], sigma[i+1]
                pair = frozenset([u, v])
                sigma_arcs[pair] = (u, v)

            # Check compatibility
            compatible = True
            constrained_pairs = set()

            for pair, direction in pi_arcs.items():
                constrained_pairs.add(pair)
                if pair in sigma_arcs:
                    if sigma_arcs[pair] != direction:
                        compatible = False
                        break

            if not compatible:
                n_incompatible += 1
                continue

            # Add σ-only pairs
            for pair in sigma_arcs:
                constrained_pairs.add(pair)

            n_compatible += 1
            c = len(constrained_pairs)
            total_sum += 2**(m - c)

    E_H2_pair = Fraction(total_sum, 2**m)
    E_H = Fraction(math.factorial(n), 2**(n-1))
    ratio = E_H2_pair / E_H**2

    print(f"\n  n={n}: E[H²] = {E_H2_pair} = {float(E_H2_pair):.4f}")
    print(f"    E[H]² = {E_H**2} = {float(E_H**2):.4f}")
    print(f"    E[H²]/E[H]² = {ratio} = {float(ratio):.8f}")
    print(f"    Compatible pairs: {n_compatible}/{len(perms)**2}")
    print(f"    Incompatible pairs: {n_incompatible}/{len(perms)**2}")

    if n <= 6:
        # Verify
        h_values = [count_H(adj, n) for adj in all_tournaments(n)]
        direct = Fraction(sum(h**2 for h in h_values), len(h_values))
        print(f"    Direct E[H²] = {direct} = {float(direct):.4f}")
        print(f"    Match: {E_H2_pair == direct}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: E[H²] Closed Form Search")
print("=" * 70)

# From the data:
# n=3: E[H²] = 3, E[H]² = 9/4, ratio = 4/3
# n=4: E[H²] = 12, E[H]² = 9, ratio = 4/3
# n=5: E[H²] = 1185/16 ??? Let me get exact value

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    sum_h2 = sum(h**2 for h in h_values)
    E_H2 = Fraction(sum_h2, N)
    E_H = Fraction(math.factorial(n), 2**(n-1))

    print(f"\n  n={n}: Σ H² = {sum_h2}")
    print(f"    E[H²] = {E_H2}")
    print(f"    E[H]² = {E_H**2}")
    ratio = E_H2 / E_H**2
    print(f"    E[H²]/E[H]² = {ratio}")

    # Factor the numerator and denominator
    r_num = ratio.numerator
    r_den = ratio.denominator
    print(f"    = {r_num}/{r_den}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: The Sum Σ H² — Double Path Counting")
print("=" * 70)

# Σ_T H(T)² = Σ_{compatible (π,σ)} 2^{m - c(π,σ)}
# where c(π,σ) = number of distinct undirected pairs constrained

# For identical π=σ: c = n-1, contribution = n! · 2^{m-n+1}
# For other compatible pairs: c ≥ n-1

# Σ_T H² = Σ_T H · H = Σ_{(T, π, σ)} 1_{π HP} 1_{σ HP}
# = Σ_{(π,σ) compatible} 2^{m-c(π,σ)}

# Can we compute Σ_{(π,σ)} 2^{-c(π,σ)} [compatible] exactly?

print("\nΣ_T H² and Σ_{(π,σ)} 2^{m-c(π,σ)} [compatible]:")
for n in range(3, 7):
    perms = list(itertools.permutations(range(n)))
    m = n*(n-1)//2

    c_distribution = Counter()
    total_weighted = 0

    for pi in perms:
        pi_arcs = {}
        for i in range(n-1):
            u, v = pi[i], pi[i+1]
            pi_arcs[frozenset([u, v])] = (u, v)

        for sigma in perms:
            sigma_arcs = {}
            compatible = True
            for i in range(n-1):
                u, v = sigma[i], sigma[i+1]
                pair = frozenset([u, v])
                sigma_arcs[pair] = (u, v)
                if pair in pi_arcs and pi_arcs[pair] != (u, v):
                    compatible = False
                    break

            if not compatible:
                c_distribution[-1] += 1
                continue

            constrained = set(pi_arcs.keys()) | set(sigma_arcs.keys())
            c = len(constrained)
            c_distribution[c] += 1
            total_weighted += 2**(m - c)

    print(f"\n  n={n}: m={m}")
    print(f"    c-distribution (c → count):")
    for c in sorted(c_distribution.keys()):
        if c == -1:
            print(f"      incompatible: {c_distribution[c]}")
        else:
            print(f"      c={c}: {c_distribution[c]} pairs, contribution = {c_distribution[c] * 2**(m-c)}")
    print(f"    Σ_T H² = {total_weighted}")

    # Verify
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    print(f"    Direct Σ H² = {sum(h**2 for h in h_values)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: Why E[H²]/E[H]² = 4/3 at n=3,4?")
print("=" * 70)

# At n=3: E[H] = 3/2, E[H²] = 3
# E[H²]/E[H]² = 3/(9/4) = 4/3

# At n=4: E[H] = 3, E[H²] = 12
# E[H²]/E[H]² = 12/9 = 4/3

# Var[H]/E[H]² = 4/3 - 1 = 1/3 at n=3,4

# This is the variance of a uniform on {1, 3}: Var = (3-1)²/4 = 1
# And E² = 4, so Var/E² = 1/4 ≠ 1/3

# Let me think about this differently.
# E[H²] = Σ_{(π,σ) compat} 2^{m-c(π,σ)} / 2^m
# E[H]² = (n!)² / 2^{2(n-1)}

# E[H²]/E[H]² = Σ_{(π,σ) compat} 2^{-c(π,σ)} / ((n!)² · 2^{-2(n-1)})
#              = Σ_{(π,σ) compat} 2^{2(n-1)-c(π,σ)} / (n!)²

print("\nσ_{compat} 2^{2(n-1)-c} / n!²:")
for n in range(3, 8):
    perms = list(itertools.permutations(range(n)))

    total = Fraction(0)
    for pi in perms:
        pi_arcs = {}
        for i in range(n-1):
            u, v = pi[i], pi[i+1]
            pi_arcs[frozenset([u, v])] = (u, v)

        for sigma in perms:
            sigma_arcs = {}
            compatible = True
            for i in range(n-1):
                u, v = sigma[i], sigma[i+1]
                pair = frozenset([u, v])
                sigma_arcs[pair] = (u, v)
                if pair in pi_arcs and pi_arcs[pair] != (u, v):
                    compatible = False
                    break

            if not compatible:
                continue

            constrained = set(pi_arcs.keys()) | set(sigma_arcs.keys())
            c = len(constrained)
            total += Fraction(2**(2*(n-1) - c), 1)

    nfact2 = Fraction(math.factorial(n)**2)
    ratio = total / nfact2
    print(f"  n={n}: Σ 2^{{2(n-1)-c}} = {total}, / n!² = {ratio} = {float(ratio):.8f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: The Inverse Moment E[1/H]")
print("=" * 70)

# E[1/H] exact values:
# n=3: 5/6
# n=4: 8/15
# n=5: 8683/34320

# Factor 34320 = 2^4 · 3 · 5 · 11 · 13 ... wait
# 34320 = 16 · 2145 = 16 · 3 · 715 = 16 · 3 · 5 · 143 = 16 · 3 · 5 · 11 · 13
# So 34320 = 2^4 · 3 · 5 · 11 · 13
# Denominators: 6 = 2·3, 15 = 3·5, 34320 = 2^4·3·5·11·13

# Interesting: the denominators involve products of odd numbers ≤ 2n-1
# 5/6: denom has factors up to 3 (= 2·3-1? No, 2n-1 = 5 for n=3)
# 8/15: denom has factors up to 5 (= 2n-1 = 7? No, 2·4-1=7)
# Hmm, 6 = 2·3, 15 = 3·5

# Let me check: are these connected to product of odd cycle lengths?
# At n=3: odd cycles up to length 3 → 1/3
# At n=4: odd cycles up to length 3 → 1/3, 1/5
# At n=5: odd cycles up to length 5 → 1/3, 1/5

# Actually let me check if E[1/H] has the form Π_k (1 - 1/f(k)) or similar

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_inv = sum(Fraction(1, h) for h in h_values) / N
    print(f"  n={n}: E[1/H] = {E_inv}")
    print(f"    Numerator factors: {E_inv.numerator}")
    print(f"    Denominator factors: {E_inv.denominator}")

    # Factor denominator
    d = E_inv.denominator
    factors = []
    for p in range(2, d+1):
        while d % p == 0:
            factors.append(p)
            d //= p
        if d == 1:
            break
    print(f"    Denominator = {' × '.join(str(f) for f in factors)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 7: The Harmonic Mean of H")
print("=" * 70)

# Harmonic mean H_harm = N / Σ 1/H = 1/E[1/H]
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_inv = sum(Fraction(1, h) for h in h_values) / N
    H_harm = 1 / E_inv
    E_H = Fraction(math.factorial(n), 2**(n-1))

    print(f"  n={n}: H_harm = {H_harm} = {float(H_harm):.4f}")
    print(f"         E[H] = {E_H} = {float(E_H):.4f}")
    print(f"         H_harm/E[H] = {float(H_harm/E_H):.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 8: E[H^r] for various r — Is there a pattern?")
print("=" * 70)

# E[H^r] for fractional r
# In particular, E[H^{1/2}] relates to the "amplitude" of tournaments
# and E[H^{-1/2}] to the inverse amplitude

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_H = Fraction(math.factorial(n), 2**(n-1))

    print(f"\n  n={n}: E[H] = {float(E_H):.4f}")
    for r in [-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]:
        if r == 0:
            val = 1.0
        else:
            val = sum(h**r for h in h_values) / N
        # Normalize by E[H]^r
        normalized = val / float(E_H)**r if float(E_H)**r != 0 else float('inf')
        print(f"    E[H^{r:4.1f}] = {val:15.6f}, E[H^r]/E[H]^r = {normalized:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 9: The H² Formula — Rewriting")
print("=" * 70)

# Σ_T H(T)² counts the number of tournament-ordered-path-pairs (T, π, σ)
# This equals Σ_{(π,σ) compatible} 2^{m - c(π,σ)}
# = Σ_{(π,σ)} [compatible] 2^{m-c(π,σ)}
#
# Key insight: compatibility means the union of directed arcs of π and σ
# form a consistent tournament (no pair has both directions)
# Equivalently: π and σ are both paths in some common tournament
#
# c(π,σ) = #{undirected pairs constrained by π or σ}
#         = #{pairs (u,v) where u immediately precedes v in π or σ}
# Note: an arc in π and σ might constrain the same pair
# If k = #shared arcs (same direction), then:
# c = (n-1) + (n-1) - k = 2(n-1) - k

# And "compatible" means: no pair is constrained in opposite directions
# i.e., no arc u→v in π and v→u in σ (for the same pair {u,v})

# Let f = #conflicting arcs (opposite direction on same pair)
# Then compatible iff f = 0

# The total pair decomposition: for each undirected pair {u,v} that is
# constrained by both π and σ, it either agrees (shared arc) or conflicts
# k + f = #{pairs constrained by both}
# c = (n-1) + (n-1) - (k + f) = 2(n-1) - k - f (but only f=0 matters)

# So for compatible pairs: c = 2(n-1) - k where k = shared arcs
# And Σ_T H² = Σ_{(π,σ) f=0} 2^{m-2(n-1)+k}

# This is exactly what I computed before! Let me check...
# E[H²] = 2^{-m} · Σ_{(π,σ) compat} 2^{m-c} = Σ_{(π,σ) compat} 2^{-c}
# = Σ_{(π,σ) compat} 2^{-(2(n-1)-k)} = 2^{-2(n-1)} Σ_{(π,σ) compat} 2^k
# E[H]² = n!²/2^{2(n-1)}
# So E[H²]/E[H]² = (Σ_{compat} 2^k) / n!²

# The earlier computation gave wrong answers because I had a bug.
# Let me redo very carefully.

print("\nREDO: Σ_{compatible (π,σ)} 2^k / n!² where k = shared arcs:")
for n in range(3, 7):
    perms = list(itertools.permutations(range(n)))
    nfact2 = len(perms)**2

    total_2k = 0
    compat_count = 0

    for pi in perms:
        # Arcs of π as directed pairs
        pi_arcs_set = set()
        pi_pairs = set()
        for i in range(n-1):
            pi_arcs_set.add((pi[i], pi[i+1]))
            pi_pairs.add(frozenset([pi[i], pi[i+1]]))

        for sigma in perms:
            sigma_arcs_set = set()
            sigma_pairs = set()
            conflict = False
            for i in range(n-1):
                u, v = sigma[i], sigma[i+1]
                sigma_arcs_set.add((u, v))
                pair = frozenset([u, v])
                sigma_pairs.add(pair)
                # Check conflict: σ needs u→v but π needs v→u
                if (v, u) in pi_arcs_set:
                    conflict = True
                    break

            if conflict:
                continue

            # Shared arcs: directed arcs in both
            k = len(pi_arcs_set & sigma_arcs_set)
            total_2k += 2**k
            compat_count += 1

    ratio = Fraction(total_2k, nfact2)
    print(f"\n  n={n}: Σ 2^k = {total_2k}, / n!² = {ratio} = {float(ratio):.8f}")
    print(f"    Compatible: {compat_count}/{nfact2}")

    # Also compute directly
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    E_H2 = Fraction(sum(h**2 for h in h_values), len(h_values))
    E_H = Fraction(math.factorial(n), 2**(n-1))
    print(f"    Direct E[H²]/E[H]² = {E_H2/E_H**2}")
    print(f"    Match: {ratio == E_H2/E_H**2}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 10: Exact Σ_T H(T)² Formula")
print("=" * 70)

# If Σ_T H² = Σ_{compat (π,σ)} 2^{m - 2(n-1) + k}
# then the key quantity is P(π,σ compatible) weighted by 2^k

# For IDENTICAL perms π=σ: always compatible, k = n-1
# Contribution: n! · 2^{n-1}
# Total Σ_T H² / 2^m = Σ_{compat} 2^{-(2(n-1)-k)}

# Let's separate diagonal and off-diagonal:
print("\nDiagonal (π=σ) vs off-diagonal contributions:")
for n in range(3, 7):
    m = n*(n-1)//2
    perms = list(itertools.permutations(range(n)))

    diag_total = 0
    offdiag_total = 0

    for pi in perms:
        pi_arcs_set = set()
        for i in range(n-1):
            pi_arcs_set.add((pi[i], pi[i+1]))

        # Diagonal term
        diag_total += 2**(n-1)  # k = n-1

        for sigma in perms:
            if sigma == pi:
                continue

            sigma_arcs_set = set()
            conflict = False
            for i in range(n-1):
                u, v = sigma[i], sigma[i+1]
                sigma_arcs_set.add((u, v))
                if (v, u) in pi_arcs_set:
                    conflict = True
                    break

            if conflict:
                continue

            k = len(pi_arcs_set & sigma_arcs_set)
            offdiag_total += 2**k

    sum_H = math.factorial(n) * 2**(m - n + 1)  # = Σ_T H
    sum_H2 = (diag_total + offdiag_total) * 2**(m - 2*(n-1))

    # Direct
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    direct_sum_H2 = sum(h**2 for h in h_values)

    print(f"\n  n={n}:")
    print(f"    Diagonal Σ 2^k = n!·2^{{n-1}} = {diag_total}")
    print(f"    Off-diagonal Σ 2^k = {offdiag_total}")
    print(f"    Total Σ 2^k = {diag_total + offdiag_total}")
    print(f"    Σ_T H² = (Σ 2^k) · 2^{{m-2(n-1)}} = {sum_H2}")
    print(f"    Direct Σ H² = {direct_sum_H2}")
    print(f"    Match: {sum_H2 == direct_sum_H2}")
    print(f"    Off-diagonal / Diagonal = {Fraction(offdiag_total, diag_total)} = {offdiag_total/diag_total:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 11: The E[H²]/E[H]² Sequence")
print("=" * 70)

# Let r_n = E[H²]/E[H]² = (Σ compat 2^k) / n!²
# r_3 = 4/3
# r_4 = 4/3
# r_5 = ?
# r_6 = ?

# From direct computation above, we have exact fractions.
# Also compute for n=7 via pair counting (this will be slow but doable)

print("\nThe ratio sequence r_n = E[H²]/E[H]²:")
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_H = Fraction(math.factorial(n), 2**(n-1))
    E_H2 = Fraction(sum(h**2 for h in h_values), N)
    r = E_H2 / E_H**2
    print(f"  r_{n} = {r} ≈ {float(r):.8f}")

    # Does 1/(r_n - 1) simplify?
    cv2 = r - 1  # = Var/E²
    print(f"  CV² = r-1 = {cv2}")
    print(f"  1/CV² = {1/cv2}")

print("\n  The CV² sequence = Var[H]/E[H]²:")
print("  1/3, 1/3, ?, ? — does it stay at 1/3?")

# Check from exact data
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_H = sum(h_values) / N
    E_H2 = sum(h**2 for h in h_values) / N
    cv2 = E_H2/E_H**2 - 1
    print(f"  n={n}: CV² = {cv2:.8f}, 1/3 = {1/3:.8f}")

print("\n" + "=" * 70)
print("END — opus-2026-03-14-S89c")
print("=" * 70)
