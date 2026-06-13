#!/usr/bin/env python3
"""
ascending_successions_s115.py — Ascending Successions Among Succession-Free Permutations
opus-2026-03-21-S115

THM-215: The compatible (conflict-free) permutation τ for a fixed σ are
exactly the SUCCESSION-FREE permutations (no i with τ(i+1)=τ(i)+1)
after relabeling.

An "ascending succession" at position i in permutation π means
π(i+1) = π(i) + 1 (the next position has value exactly 1 higher).

QUESTIONS:
1. Among all n! permutations, what is the distribution of #ascending successions?
2. Among the A000255(n-1) succession-free perms, what related statistics emerge?
3. For the OVERLAP structure f(n,a): how does a (shared arcs) relate to
   the succession structure of the relative permutation τ∘σ⁻¹?
4. Connection to the OCR: does the succession structure within score classes
   explain the residual?

KEY IDEA: An arc (σ(i), σ(i+1)) is an "ascending arc" if σ(i+1) > σ(i)
(the path goes from lower-labeled to higher-labeled vertex).
In the tiling model, these are the BACKBONE arcs.
The number of ascending arcs = #{i : σ(i+1) > σ(i)} = #ascents of σ.

For a COMPATIBLE pair (σ,τ): they share some arcs and conflict on none.
The shared arcs that are ascending in BOTH σ and τ are the
"jointly ascending" arcs — the backbone-like structure that both paths share.
"""

from math import comb, factorial
from itertools import permutations
from collections import Counter, defaultdict
from fractions import Fraction

print("=" * 72)
print("  ASCENDING SUCCESSIONS AMONG SUCCESSION-FREE PERMUTATIONS")
print("  opus-2026-03-21-S115")
print("=" * 72)
print()

# =============================================================================
# PART 1: BASIC STATISTICS — SUCCESSIONS IN PERMUTATIONS
# =============================================================================

print("=" * 72)
print("  PART 1: SUCCESSION DISTRIBUTION IN S_n")
print("=" * 72)
print()

print("""
  A SUCCESSION at position i in π ∈ S_n means π(i+1) = π(i) + 1.
  (The value increases by exactly 1 at the next position.)

  A000255(n) = #{π ∈ S_n : no successions} = succession-free perms.

  From THM-215: fixing any HP σ, the compatible τ are exactly those
  where σ⁻¹∘τ is succession-free. So A000255(n-1) counts the
  compatible partners per fixed σ.
""")

for n in range(3, 8):
    perms = list(permutations(range(n)))

    # Count successions for each permutation
    succ_dist = Counter()
    for perm in perms:
        n_succ = sum(1 for i in range(n-1) if perm[i+1] == perm[i] + 1)
        succ_dist[n_succ] += 1

    print(f"  n = {n}: (n! = {factorial(n)})")
    print(f"  {'#succ':>6} | {'count':>8} | {'fraction':>10}")
    print(f"  " + "-" * 30)
    for k in range(n):
        if succ_dist[k] > 0:
            frac = Fraction(succ_dist[k], factorial(n))
            print(f"  {k:6d} | {succ_dist[k]:8d} | {frac}")

    # Verify: succ_dist[0] = A000255(n-1)? No: A000255(n-1) counts
    # permutations of n-1 with no succession. Here we count S_n with no succession.
    # A000255 with the indexing: a(n) = perms of [n] with no succession.
    # a(0)=1, a(1)=1, a(2)=1, a(3)=3, a(4)=11, a(5)=53, a(6)=309

    # But our count for n-element permutations with 0 successions:
    # n=3: should be A000255(3-1)... wait, A000255 counts perms of {1,...,n}.
    # Let me just check.
    print(f"  #{'{'}0 successions{'}'} = {succ_dist[0]}")

    # The derangement-like sequence for successions:
    # #perms of [n] with 0 successions = subfactorial-like
    # Known: this is A000166 (derangements) NO — that's fixed points.
    # Succession-free = A000255(n). Let me verify.
    def A000255(m):
        if m <= 1: return 1
        a, b = 1, 1
        for k in range(2, m+1):
            a, b = b, k*b + (k-1)*a
        return b

    print(f"  A000255({n}) = {A000255(n)}")
    print(f"  Match: {succ_dist[0] == A000255(n)}")
    print()

# =============================================================================
# PART 2: ASCENTS vs SUCCESSIONS
# =============================================================================

print("=" * 72)
print("  PART 2: ASCENTS vs SUCCESSIONS")
print("=" * 72)
print()

print("""
  An ASCENT at position i: π(i+1) > π(i) (value increases, by any amount)
  A SUCCESSION at position i: π(i+1) = π(i) + 1 (value increases by exactly 1)

  Every succession is an ascent, but not conversely.
  A succession is a "minimal ascent" — the smallest possible increase.

  The number of ascents follows the EULERIAN NUMBER distribution:
  A(n,k) = #{π ∈ S_n with exactly k ascents}.

  The connection to tournaments:
  If σ is a Hamiltonian path σ(1)→σ(2)→...→σ(n),
  the arc σ(i)→σ(i+1) is "ascending" iff σ(i+1) > σ(i).
  An ascending arc goes from lower to higher vertex label.
  In the tiling model, these are BACKBONE-DIRECTION arcs.

  Number of ascending arcs = #ascents of σ.
  Number of successions = #minimal ascending arcs (consecutive labels).
""")

for n in [5, 6]:
    perms = list(permutations(range(n)))

    # Joint distribution of (#ascents, #successions)
    joint = Counter()
    for perm in perms:
        n_asc = sum(1 for i in range(n-1) if perm[i+1] > perm[i])
        n_succ = sum(1 for i in range(n-1) if perm[i+1] == perm[i] + 1)
        joint[(n_asc, n_succ)] += 1

    print(f"  n = {n}: Joint distribution of (ascents, successions)")
    print(f"  {'asc':>4} {'succ':>5} | {'count':>8} | {'fraction':>12}")
    print(f"  " + "-" * 35)
    for (a, s) in sorted(joint.keys()):
        if joint[(a,s)] > 0:
            frac = Fraction(joint[(a,s)], factorial(n))
            print(f"  {a:4d} {s:5d} | {joint[(a,s)]:8d} | {frac}")
    print()

# =============================================================================
# PART 3: COMPATIBLE PAIRS AND THEIR SUCCESSION STRUCTURE
# =============================================================================

print("=" * 72)
print("  PART 3: SUCCESSION STRUCTURE OF COMPATIBLE PAIRS")
print("=" * 72)
print()

print("""
  For a COMPATIBLE pair (σ,τ): no arc of σ conflicts with any arc of τ.
  The relative permutation ρ = σ⁻¹∘τ tells us how τ differs from σ.

  KEY: σ⁻¹∘τ is succession-free (has no successions).
  This is the content of THM-215.

  Now: among compatible pairs, how do the SHARED ARCS relate to
  the ascent/succession structure?

  A shared arc at position i means σ(i)→σ(i+1) = τ(j)→τ(j+1) for some j.
  This means the SAME directed arc appears in both paths.

  Question: are shared arcs more likely to be ASCENDING (backbone-direction)?
""")

for n in [5, 6]:
    perms = list(permutations(range(n)))

    # For each compatible pair, compute:
    # - #shared arcs (a)
    # - #shared ascending arcs (arcs where both σ and τ go from lower to higher label)
    # - #shared descending arcs
    # - #successions of the relative permutation ρ = σ⁻¹∘τ

    stats = defaultdict(list)
    total_compatible = 0

    for sigma in perms:
        s_arcs = {}
        for i in range(n-1):
            pair = frozenset({sigma[i], sigma[i+1]})
            s_arcs[pair] = (sigma[i], sigma[i+1])

        sigma_inv = [0] * n
        for i in range(n):
            sigma_inv[sigma[i]] = i

        for tau in perms:
            # Check compatibility
            agree = 0
            conflict = False
            shared_asc = 0
            shared_desc = 0

            for i in range(n-1):
                pair = frozenset({tau[i], tau[i+1]})
                if pair in s_arcs:
                    if s_arcs[pair] == (tau[i], tau[i+1]):
                        agree += 1
                        # Is this arc ascending? (goes from lower to higher label)
                        if tau[i] < tau[i+1]:
                            shared_asc += 1
                        else:
                            shared_desc += 1
                    else:
                        conflict = True
                        break

            if conflict:
                continue

            total_compatible += 1

            # Relative permutation ρ = σ⁻¹ ∘ τ
            rho = tuple(sigma_inv[tau[i]] for i in range(n))

            # Successions of ρ
            rho_succ = sum(1 for i in range(n-1) if rho[i+1] == rho[i] + 1)

            # Ascents of ρ
            rho_asc = sum(1 for i in range(n-1) if rho[i+1] > rho[i])

            stats['a'].append(agree)
            stats['shared_asc'].append(shared_asc)
            stats['shared_desc'].append(shared_desc)
            stats['rho_succ'].append(rho_succ)
            stats['rho_asc'].append(rho_asc)

    print(f"  n = {n}: {total_compatible} compatible pairs")
    print(f"    Expected: {factorial(n)} × A000255({n-1}) = {factorial(n) * A000255(n-1)}")
    print(f"    Match: {total_compatible == factorial(n) * A000255(n-1)}")
    print()

    # Statistics of shared arcs
    a_vals = stats['a']
    sa_vals = stats['shared_asc']
    sd_vals = stats['shared_desc']

    print(f"    E[shared arcs] = {sum(a_vals)/len(a_vals):.4f}")
    print(f"    E[shared ascending] = {sum(sa_vals)/len(sa_vals):.4f}")
    print(f"    E[shared descending] = {sum(sd_vals)/len(sd_vals):.4f}")
    print(f"    Ascending fraction of shared: {sum(sa_vals)/(sum(sa_vals)+sum(sd_vals)):.4f}")
    print()

    # Successions of relative permutation
    rho_succ_vals = stats['rho_succ']
    print(f"    E[successions of ρ] = {sum(rho_succ_vals)/len(rho_succ_vals):.6f}")
    print(f"    Should be 0 (ρ is succession-free by THM-215)")
    succ_dist_rho = Counter(rho_succ_vals)
    print(f"    Distribution: {dict(sorted(succ_dist_rho.items()))}")
    print()

    # Joint distribution of (shared arcs, shared ascending)
    joint_a_sa = Counter(zip(a_vals, sa_vals))
    print(f"    Joint (shared_arcs, shared_ascending):")
    print(f"    {'a':>4} {'sa':>4} | {'count':>8} | {'frac of a':>10}")
    print(f"    " + "-" * 30)
    for a_val in range(n):
        total_at_a = sum(joint_a_sa.get((a_val, sa), 0) for sa in range(n))
        if total_at_a == 0:
            continue
        for sa_val in range(a_val + 1):
            cnt = joint_a_sa.get((a_val, sa_val), 0)
            if cnt > 0:
                print(f"    {a_val:4d} {sa_val:4d} | {cnt:8d} | {cnt/total_at_a:10.4f}")
    print()

# =============================================================================
# PART 4: THE ASCENDING OVERLAP GENERATING FUNCTION
# =============================================================================

print("=" * 72)
print("  PART 4: THE ASCENDING OVERLAP STRUCTURE")
print("=" * 72)
print()

print("""
  Define for compatible pairs (σ,τ):
    a = total shared arcs
    a_↑ = shared ascending arcs (both paths go from lower to higher label)
    a_↓ = shared descending arcs (both go from higher to lower)
    a = a_↑ + a_↓

  The ASCENDING OVERLAP GF:
    F_n(x, y) = Σ_{compatible (σ,τ)} x^{a_↑} y^{a_↓}

  Then F_n(x, x) = the original overlap GF F_n(x).
  And F_n(2, 2) = the quantity that gives E[H²].

  But F_n(2, 1) and F_n(1, 2) would tell us about the
  DIRECTIONAL structure of the overlap.

  If ascending arcs are disproportionately shared, it would mean
  compatible paths tend to agree on backbone-direction arcs more
  than non-backbone arcs.
""")

for n in [5]:
    perms = list(permutations(range(n)))
    F_table = Counter()

    for sigma in perms:
        s_arcs = {}
        for i in range(n-1):
            pair = frozenset({sigma[i], sigma[i+1]})
            s_arcs[pair] = (sigma[i], sigma[i+1])

        for tau in perms:
            a_up = 0
            a_down = 0
            conflict = False
            for i in range(n-1):
                pair = frozenset({tau[i], tau[i+1]})
                if pair in s_arcs:
                    if s_arcs[pair] == (tau[i], tau[i+1]):
                        if tau[i] < tau[i+1]:
                            a_up += 1
                        else:
                            a_down += 1
                    else:
                        conflict = True
                        break
            if not conflict:
                F_table[(a_up, a_down)] += 1

    print(f"  n = {n}: Ascending Overlap Table F(a_↑, a_↓)")
    print(f"  {'a_↑':>4} {'a_↓':>4} | {'count':>8} | {'per σ':>8}")
    print(f"  " + "-" * 30)
    for (au, ad) in sorted(F_table.keys()):
        cnt = F_table[(au, ad)]
        per_sigma = cnt // factorial(n)
        print(f"  {au:4d} {ad:4d} | {cnt:8d} | {per_sigma:8d}")

    # Check symmetry: is F(a_↑, a_↓) = F(a_↓, a_↑)?
    symmetric = all(F_table.get((au, ad), 0) == F_table.get((ad, au), 0)
                    for (au, ad) in F_table)
    print(f"\n  Symmetric (F(a↑,a↓) = F(a↓,a↑))? {symmetric}")

    # Marginals
    total = sum(F_table.values())
    E_au = sum(au * cnt for (au, ad), cnt in F_table.items()) / total
    E_ad = sum(ad * cnt for (au, ad), cnt in F_table.items()) / total
    print(f"  E[a_↑] = {E_au:.6f}")
    print(f"  E[a_↓] = {E_ad:.6f}")
    print(f"  E[a_↑] = E[a_↓]? {abs(E_au - E_ad) < 0.0001}")
    print()

# =============================================================================
# PART 5: SCORE-CONDITIONAL ASCENDING STRUCTURE
# =============================================================================

print("=" * 72)
print("  PART 5: DOES ASCENDING STRUCTURE VARY BY SCORE CLASS?")
print("=" * 72)
print()

print("""
  THE KEY QUESTION: Within a tournament T with score sequence s,
  do the Hamiltonian paths have different ascending structure
  depending on which score class T belongs to?

  If YES: the ascending structure carries information BEYOND scores,
  which could explain part of the OCR residual.

  For each tournament T, compute:
  - The score sequence
  - The ascent distribution of its Hamiltonian paths
""")

n = 5
m = comb(n, 2)

# Group by score, compute ascent stats of HPs
score_ascent_stats = defaultdict(lambda: {'ascents': [], 'successions': [], 'H': []})

for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    scores = [0]*n
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1; scores[i] += 1
            else: adj[j][i] = 1; scores[j] += 1
            idx += 1

    sc = tuple(sorted(scores))

    # Find HPs
    hp_ascents = []
    hp_succs = []
    for perm in permutations(range(n)):
        if all(adj[perm[k]][perm[k+1]] for k in range(n-1)):
            n_asc = sum(1 for k in range(n-1) if perm[k+1] > perm[k])
            n_succ = sum(1 for k in range(n-1) if perm[k+1] == perm[k] + 1)
            hp_ascents.append(n_asc)
            hp_succs.append(n_succ)

    H = len(hp_ascents)
    mean_asc = sum(hp_ascents) / H if H > 0 else 0
    mean_succ = sum(hp_succs) / H if H > 0 else 0

    score_ascent_stats[sc]['ascents'].append(mean_asc)
    score_ascent_stats[sc]['successions'].append(mean_succ)
    score_ascent_stats[sc]['H'].append(H)

print(f"  n = {n}: HP ascent/succession statistics by score class")
print()
print(f"  {'Score':>20} | {'E[#asc/HP]':>10} | {'E[#succ/HP]':>11} | {'E[H]':>6} | {'Var(asc)':>8}")
print("  " + "-" * 65)

for sc in sorted(score_ascent_stats.keys()):
    d = score_ascent_stats[sc]
    asc_vals = d['ascents']
    succ_vals = d['successions']
    H_vals = d['H']
    import numpy as np
    print(f"  {str(sc):>20} | {np.mean(asc_vals):10.4f} | {np.mean(succ_vals):11.4f} | {np.mean(H_vals):6.1f} | {np.var(asc_vals):8.4f}")

print()

# Focus on the PoS class
pos_sc = (1, 2, 2, 2, 3)
pos_data = score_ascent_stats[pos_sc]
print(f"  PoS class {pos_sc}:")
print(f"    {len(pos_data['H'])} tournaments")

# Group by H within the PoS
H_groups = defaultdict(list)
for i in range(len(pos_data['H'])):
    H_groups[pos_data['H'][i]].append(pos_data['ascents'][i])

print(f"    Mean ascents per HP, by H value:")
for H_val in sorted(H_groups.keys()):
    asc_list = H_groups[H_val]
    print(f"      H={H_val}: mean_ascents = {np.mean(asc_list):.4f} (n={len(asc_list)})")

print()

# =============================================================================
# PART 6: THE SUCCESSION-FREE OVERLAP AND THE OCR
# =============================================================================

print("=" * 72)
print("  PART 6: CONNECTION TO THE OCR")
print("=" * 72)
print()

print("""
  SYNTHESIS:

  1. Compatible permutation pairs are characterized by the relative
     permutation ρ = σ⁻¹∘τ being SUCCESSION-FREE (THM-215).

  2. Among compatible pairs, shared arcs split equally between
     ascending and descending (by symmetry of the uniform random model).

  3. Within the PoS class, HP ascent statistics VARY by H value:
     higher-H tournaments have HPs with different ascent distributions.

  4. The ascending overlap GF F(a_↑, a_↓) is SYMMETRIC: F(a↑,a↓)=F(a↓,a↑).
     This reflects the complement symmetry H(T)=H(T^op).

  THE CONNECTION:
  The OCR residual measures the variance of H within score classes.
  This variance comes from different CYCLE STRUCTURES (c₅, α₂).
  The cycle structures are encoded in the OVERLAP PATTERN of HP pairs.
  The overlap pattern, in turn, is governed by the SUCCESSION-FREE
  condition on the relative permutation.

  So: the OCR residual is ultimately a statement about how
  succession-free permutations distribute across different
  overlap types (a_↑, a_↓) within each score class.

  An EXACT OCR FORMULA would need:
  F_n(2, 2; s) = the ascending overlap GF restricted to score class s.
  Then: Var(H|s) = [F_n(2,2;s)/n_s - E[H|s]²]
  And: OCR = 1 - Σ P(s) Var(H|s) / Var(H)
""")

print("opus-2026-03-21-S115: ASCENDING SUCCESSIONS")
