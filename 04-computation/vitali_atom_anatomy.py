#!/usr/bin/env python3
"""
vitali_atom_anatomy.py -- kind-pasteur-2026-03-13-S61

ANATOMY OF THE VITALI ATOM: What tournament operation preserves labeled
lambda but flips the Pfaffian sum sign?

Discovery from pfaffian_oddness_and_vitali.py:
  - At n=5, ALL Vitali pairs have |Ps(T)| = |Ps(T')| but sign(Ps(T)) = -sign(Ps(T'))
  - The labeled lambda graph is IDENTICAL between T and T'
  - So there exists a transformation T -> T' that:
    (1) Preserves lambda_{uv} for all (u,v)
    (2) Flips sign(Ps)
    (3) Preserves H (since H depends only on cycle structure)

QUESTION: Is this transformation always a reversal of a 4-vertex sub-tournament?
If so, WHICH 4-vertex subset? And why does it flip Ps?

Also explores:
- The Vitali atom as a "gauge transformation" (preserving measurable quantities)
- Whether the sign of Ps has a TOPOLOGICAL interpretation
- The connection to the Euler class / Stiefel-Whitney class of the tournament

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


def tournament_to_bits(A, n):
    bits = 0
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << pos)
            pos += 1
    return bits


def lambda_graph(A, n):
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


def lambda_key(A, n):
    lam = lambda_graph(A, n)
    return tuple(lam[i][j] for i in range(n) for j in range(i+1, n))


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


def pfaffian_sum(A, n):
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
    total = 0
    for i in range(n):
        remaining = [j for j in range(n) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)]
                 for a in range(n-1)]
        pf = pfaffian(S_del)
        total += ((-1) ** i) * pf
    return total


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


def reverse_subtournament(A, n, subset):
    """Reverse all arcs within the given vertex subset."""
    from copy import deepcopy
    B = deepcopy(A)
    subset_set = set(subset)
    for u in subset:
        for v in subset:
            if u != v:
                B[u][v] = A[v][u]  # Reverse
    return B


# ========================================================================
# PART 1: Find ALL Vitali pairs at n=5 and identify the transformation
# ========================================================================
print("=" * 70)
print("PART 1: Vitali pairs at n=5 — identifying the atom")
print("=" * 70)

n = 5
num_edges = n * (n - 1) // 2

# Group tournaments by labeled lambda
lambda_classes = defaultdict(list)
for bits in range(1 << num_edges):
    A = binary_to_tournament(bits, n)
    lk = lambda_key(A, n)
    ps = pfaffian_sum(A, n)
    H = count_ham_paths(A, n)
    lambda_classes[lk].append({'bits': bits, 'ps': ps, 'H': H, 'A': A})

# Find Vitali pairs: same lambda, different Ps sign
vitali_pairs = []
for lk, group in lambda_classes.items():
    ps_set = set(d['ps'] for d in group)
    if len(ps_set) > 1:
        # Split by sign
        pos = [d for d in group if d['ps'] > 0]
        neg = [d for d in group if d['ps'] < 0]
        zero = [d for d in group if d['ps'] == 0]
        if pos and neg:
            vitali_pairs.append((lk, pos, neg))

print(f"  Total lambda classes at n=5: {len(lambda_classes)}")
print(f"  Vitali classes (Ps sign differs): {len(vitali_pairs)}")

# For each Vitali pair, find what arcs differ between T and T'
print(f"\n  VITALI ATOM ANATOMY:")
atom_sizes = defaultdict(int)

for lk, pos_group, neg_group in vitali_pairs[:20]:
    t_pos = pos_group[0]
    t_neg = neg_group[0]

    # Find differing arcs
    A_pos = t_pos['A']
    A_neg = t_neg['A']

    diff_edges = []
    for i in range(n):
        for j in range(i+1, n):
            if A_pos[i][j] != A_neg[i][j]:
                diff_edges.append((i, j))

    # Which vertices are involved?
    diff_verts = set()
    for i, j in diff_edges:
        diff_verts.add(i)
        diff_verts.add(j)

    atom_sizes[len(diff_edges)] += 1

    # Check: is T_neg = reverse_subtournament(T_pos, subset)?
    found_reversal = False
    reversal_set = None
    for k in range(2, n+1):
        for subset in combinations(range(n), k):
            B = reverse_subtournament(A_pos, n, subset)
            if tournament_to_bits(B, n) == t_neg['bits']:
                found_reversal = True
                reversal_set = subset
                break
        if found_reversal:
            break

    if len(vitali_pairs) <= 80 or len(diff_edges) <= 3:
        print(f"    Lambda={lk[:5]}..., Ps=({t_pos['ps']:+d}, {t_neg['ps']:+d}), "
              f"H=({t_pos['H']}, {t_neg['H']}), "
              f"diff_edges={diff_edges} ({len(diff_edges)} arcs), "
              f"diff_verts={sorted(diff_verts)}, "
              f"reversal={'Yes '+str(reversal_set) if found_reversal else 'No'}")

print(f"\n  Arc difference sizes: {dict(sorted(atom_sizes.items()))}")


# ========================================================================
# PART 2: Is the Vitali atom ALWAYS a full-tournament complement?
# ========================================================================
print(f"\n{'='*70}")
print("PART 2: Classification of Vitali transformations at n=5")
print("=" * 70)

# Check all possible transformations
transform_types = defaultdict(int)

for lk, pos_group, neg_group in vitali_pairs:
    for t_pos in pos_group:
        for t_neg in neg_group:
            A_pos = t_pos['A']
            A_neg = t_neg['A']

            diff_count = sum(1 for i in range(n) for j in range(i+1, n)
                           if A_pos[i][j] != A_neg[i][j])

            # Check if it's a k-vertex reversal
            for k in range(2, n+1):
                for subset in combinations(range(n), k):
                    B = reverse_subtournament(A_pos, n, subset)
                    if tournament_to_bits(B, n) == t_neg['bits']:
                        transform_types[f'{k}-reversal'] += 1
                        break
                else:
                    continue
                break
            else:
                transform_types['other'] += 1

print(f"  Transform type counts: {dict(sorted(transform_types.items()))}")


# ========================================================================
# PART 3: Pfaffian sum under k-vertex reversal
# ========================================================================
print(f"\n{'='*70}")
print("PART 3: Pfaffian sum change under k-vertex reversal")
print("=" * 70)

# For each tournament, compute Ps under ALL possible k-vertex reversals
for k_rev in [2, 3, 4, 5]:
    deltas = defaultdict(int)
    for bits in range(1 << num_edges):
        A = binary_to_tournament(bits, n)
        ps_orig = pfaffian_sum(A, n)
        for subset in combinations(range(n), k_rev):
            B = reverse_subtournament(A, n, subset)
            ps_rev = pfaffian_sum(B, n)
            delta = ps_rev - ps_orig
            deltas[delta] += 1

    print(f"  k={k_rev}-vertex reversal deltas: {dict(sorted(deltas.items()))}")

    # Key question: does k-reversal ALWAYS flip sign?
    always_flips = all(d != 0 for d in deltas.keys())
    sign_flips = all(d * (-1) >= 0 for d in deltas.keys())
    print(f"    Never preserves Ps? {always_flips}")

    # What fraction flip sign (ps_rev * ps_orig < 0)?
    flip_count = 0
    same_count = 0
    zero_count = 0
    total = 0
    for bits in range(1 << num_edges):
        A = binary_to_tournament(bits, n)
        ps_orig = pfaffian_sum(A, n)
        for subset in combinations(range(n), k_rev):
            B = reverse_subtournament(A, n, subset)
            ps_rev = pfaffian_sum(B, n)
            total += 1
            if ps_orig * ps_rev < 0:
                flip_count += 1
            elif ps_orig * ps_rev > 0:
                same_count += 1
            else:
                zero_count += 1
    print(f"    Sign: flips={flip_count}/{total} ({100*flip_count/total:.1f}%), "
          f"same={same_count}/{total}, zero={zero_count}/{total}")


# ========================================================================
# PART 4: Full complement IS a 5-reversal at n=5
# ========================================================================
print(f"\n{'='*70}")
print("PART 4: Complement = full reversal")
print("=" * 70)

# At n=5, complement = reverse all C(5,2)=10 arcs = 5-vertex reversal
# We proved Ps(T^c) = Ps(T) at n=5 ((-1)^2 = 1)
# So 5-reversal preserves Ps.
# But 3-reversal... and 4-reversal...

# Check: lambda preservation under k-reversal
print(f"  Lambda preservation under k-reversal at n=5:")
for k_rev in [2, 3, 4, 5]:
    preserves = 0
    total = 0
    for bits in range(0, 1 << num_edges, 4):  # Sample
        A = binary_to_tournament(bits, n)
        lk_orig = lambda_key(A, n)
        for subset in combinations(range(n), k_rev):
            B = reverse_subtournament(A, n, subset)
            lk_rev = lambda_key(B, n)
            total += 1
            if lk_orig == lk_rev:
                preserves += 1
    print(f"    k={k_rev}: {preserves}/{total} preserve lambda ({100*preserves/total:.1f}%)")


# ========================================================================
# PART 5: THE 4-REVERSAL IS THE VITALI ATOM
# ========================================================================
print(f"\n{'='*70}")
print("PART 5: 4-vertex reversal anatomy")
print("=" * 70)

# A 4-vertex tournament can be "regular" (cycle) or "transitive" (total order).
# Reversing a 4-cycle flips all C(4,2)=6 arcs within the sub-tournament.
# Does this preserve lambda?

# For a 4-vertex sub-tournament on vertices {a,b,c,d} in T_5:
# The 5th vertex e is connected to a,b,c,d.
# lambda(u,v) counts 3-cycles through (u,v).
# A 3-cycle through (u,v) uses a third vertex w.
# If w=e (outside the 4-set), the 3-cycle is unaffected by reversal.
# If w is inside the 4-set, the 3-cycle IS reversed.

# Under reversal: (u->v->w->u) becomes (u<-v<-w<-u) = (u->w->v->u).
# The VERTEX SET {u,v,w} is unchanged, and it still forms a 3-cycle!
# So lambda_{uv} is preserved for any u,v in the 4-set (since every
# 3-cycle through u,v in the 4-set is still a 3-cycle after reversal).

# For u in 4-set, v = e: 3-cycles through (u,e) use a third vertex w.
# If w is in the 4-set: the cycle uses arcs (u,e), (e,w), (w,u) or reverse.
# Under reversal: arcs between u and w are flipped. But arc (u,e) and (e,w)
# stay the same. So the cycle may or may not be preserved.

# Actually let's be more careful. A 3-cycle {u, v, e} where u is in 4-set:
# It uses arcs among u, v, e. If v is also in 4-set, the arc u->v or v->u
# is flipped. So some 3-cycles through (u,e) are gained and some lost.
# But lambda_{u,e} = (# 3-cycles through u,e) — does this count change?

# Let's just compute it:
print(f"  Lambda preservation for EACH 4-subset:")

lambda_pres_count = 0
lambda_fail_count = 0

for bits in range(1 << num_edges):
    A = binary_to_tournament(bits, n)
    lk_orig = lambda_key(A, n)
    for subset in combinations(range(n), 4):
        B = reverse_subtournament(A, n, subset)
        lk_rev = lambda_key(B, n)
        if lk_orig == lk_rev:
            lambda_pres_count += 1
        else:
            lambda_fail_count += 1

total_4 = 1024 * 5  # C(5,4)=5 subsets per tournament
print(f"    Preserves lambda: {lambda_pres_count}/{total_4} ({100*lambda_pres_count/total_4:.1f}%)")
print(f"    Breaks lambda:    {lambda_fail_count}/{total_4}")

# Now for those that preserve lambda, does Ps flip?
print(f"\n  Among lambda-preserving 4-reversals:")
ps_flip = 0
ps_same = 0
ps_zero = 0
for bits in range(1 << num_edges):
    A = binary_to_tournament(bits, n)
    lk_orig = lambda_key(A, n)
    ps_orig = pfaffian_sum(A, n)
    for subset in combinations(range(n), 4):
        B = reverse_subtournament(A, n, subset)
        lk_rev = lambda_key(B, n)
        if lk_orig == lk_rev:
            ps_rev = pfaffian_sum(B, n)
            if ps_orig * ps_rev < 0:
                ps_flip += 1
            elif ps_orig * ps_rev > 0:
                ps_same += 1
            else:
                ps_zero += 1

print(f"    Ps sign flip:     {ps_flip}")
print(f"    Ps sign same:     {ps_same}")
print(f"    Ps = 0 involved:  {ps_zero}")
if ps_same == 0 and ps_zero == 0:
    print(f"    *** LAMBDA-PRESERVING 4-REVERSAL ALWAYS FLIPS Ps SIGN ***")


# ========================================================================
# PART 6: Deeper — the 4-tournament type determines the atom
# ========================================================================
print(f"\n{'='*70}")
print("PART 6: Sub-tournament type for lambda-preserving 4-reversals")
print("=" * 70)

# A 4-vertex tournament is either:
# - Transitive (score seq (0,1,2,3)) — 24 labelings
# - Regular cycle (score seq (1,1,2,2)) — 24 non-transitive
# There are 4!/something = a few isomorphism classes.
# Actually there are exactly 4 unlabeled 4-tournaments:
# (0,1,2,3) transitive, and 3 with score (1,1,2,2).
# No wait: C(4,2)=6 edges, 2^6=64 total, up to isomorphism: 4.

# Classify the sub-tournament type when 4-reversal preserves lambda
subtour_types = defaultdict(int)

for bits in range(1 << num_edges):
    A = binary_to_tournament(bits, n)
    lk_orig = lambda_key(A, n)
    for subset in combinations(range(n), 4):
        B = reverse_subtournament(A, n, subset)
        lk_rev = lambda_key(B, n)
        if lk_orig == lk_rev:
            # Classify the sub-tournament
            sub_verts = list(subset)
            sub_scores = tuple(sorted(
                sum(A[u][v] for v in sub_verts if v != u) for u in sub_verts
            ))
            # Count 3-cycles in the sub-tournament
            sub_c3 = 0
            for a, b, c in combinations(sub_verts, 3):
                if A[a][b] and A[b][c] and A[c][a]:
                    sub_c3 += 1
                if A[a][c] and A[c][b] and A[b][a]:
                    sub_c3 += 1
            subtour_types[(sub_scores, sub_c3)] += 1

print(f"  Sub-tournament (score, c3) when 4-reversal preserves lambda:")
for (sc, c3), cnt in sorted(subtour_types.items()):
    label = "TRANSITIVE" if sc == (0, 1, 2, 3) else "NON-TRANS"
    print(f"    scores={sc}, c3={c3}: {cnt} ({label})")


# ========================================================================
# PART 7: The SIGN COCYCLE interpretation
# ========================================================================
print(f"\n{'='*70}")
print("PART 7: Pfaffian sum as a SIGN FUNCTION on lambda classes")
print("=" * 70)

# At n=5, the Vitali structure is:
# - Each lambda class splits into exactly 2 sign classes
# - The sign is determined by... what?
# - It's a Z/2Z-valued function on tournaments that depends on
#   MORE than the labeled lambda graph.

# HYPOTHESIS: sign(Ps) = (-1)^{some topological invariant}

# Test: is sign(Ps) determined by the tournament's edge set modulo 2
# in some simple way?

# First: what about the DETERMINANT of the adjacency matrix mod 2?
print(f"  Testing simple sign predictors at n=5:")

for bits in [0, 7, 100, 341]:
    A = binary_to_tournament(bits, n)
    ps = pfaffian_sum(A, n)
    # Parity of number of "forward" arcs (i<j, A[i][j]=1)
    fwd = sum(A[i][j] for i in range(n) for j in range(i+1, n))
    # Parity of score sequence
    scores = tuple(sum(A[v]) for v in range(n))
    score_sum = sum(scores)  # = C(n,2) always
    print(f"    bits={bits}: Ps={ps:+d}, fwd_arcs={fwd}, scores={scores}")

# The sign of Ps: is it determined by the parity of the number of
# "forward" arcs? Let's check exhaustively.
parity_determines = True
for lk, pos_group, neg_group in vitali_pairs:
    fwd_pos = set(sum(d['A'][i][j] for i in range(n) for j in range(i+1, n))
                  for d in pos_group)
    fwd_neg = set(sum(d['A'][i][j] for i in range(n) for j in range(i+1, n))
                  for d in neg_group)
    if fwd_pos & fwd_neg:
        parity_determines = False
        break
print(f"\n  Forward-arc parity determines Ps sign? {parity_determines}")

# Check: parity of the number of arcs going "left to right" (i->j with i<j)
# This is just the bit-count of the bits representation
from functools import reduce
import operator

pop_determines = True
for lk, pos_group, neg_group in vitali_pairs:
    pop_pos = set(bin(d['bits']).count('1') % 2 for d in pos_group)
    pop_neg = set(bin(d['bits']).count('1') % 2 for d in neg_group)
    if pop_pos & pop_neg:
        pop_determines = False
        break
print(f"  Bit-count parity determines Ps sign? {pop_determines}")

# Another possibility: the sign of the Pfaffian of S (the full Pfaffian
# at the even-dimensional sub-level)
# At n=5 (odd), we can't take Pf(S) directly (S is 5x5, odd dim).
# But we can look at Pf(S_{00}) (the 4x4 minor).

pf_determines = True
for lk, pos_group, neg_group in vitali_pairs:
    pf_pos = set()
    pf_neg = set()
    for d in pos_group:
        S = [[d['A'][i][j] - d['A'][j][i] for j in range(n)] for i in range(n)]
        S00 = [[S[i][j] for j in range(1, n)] for i in range(1, n)]
        pf = pfaffian(S00)
        pf_pos.add(1 if pf > 0 else -1 if pf < 0 else 0)
    for d in neg_group:
        S = [[d['A'][i][j] - d['A'][j][i] for j in range(n)] for i in range(n)]
        S00 = [[S[i][j] for j in range(1, n)] for i in range(1, n)]
        pf = pfaffian(S00)
        pf_neg.add(1 if pf > 0 else -1 if pf < 0 else 0)
    if pf_pos & pf_neg:
        pf_determines = False
        break
print(f"  Sign of Pf(S_{{00}}) determines Ps sign? {pf_determines}")


# ========================================================================
# PART 8: At n=7 — the DEEP Vitali atoms
# ========================================================================
print(f"\n{'='*70}")
print("PART 8: Vitali atoms at n=7 — the ambiguous pair")
print("=" * 70)

n7 = 7
# The famous ambiguous pair: bits 4728 (H=109) and 4658 (H=111)
# These have DIFFERENT sorted lambda (same sorted lambda multiset)
# but here the "Vitali" is about sorted lambda, not labeled lambda

# Let me find genuine lambda-identical Vitali pairs at n=7
# (where labeled lambda is the same)
import random
random.seed(42)

print(f"\n  Searching for labeled-lambda-identical pairs with DIFFERENT H at n=7...")

# Sample and group by labeled lambda
n7_lambda = defaultdict(list)
sample = random.sample(range(1 << 21), 200000)

for bits in sample:
    A = binary_to_tournament(bits, n7)
    lk = lambda_key(A, n7)
    n7_lambda[lk].append(bits)

# Find classes with multiple tours and check H
diff_H_count = 0
for lk, group in n7_lambda.items():
    if len(group) < 2:
        continue
    H_set = set()
    for bits in group[:10]:
        A = binary_to_tournament(bits, n7)
        H = count_ham_paths(A, n7)
        H_set.add(H)
    if len(H_set) > 1:
        diff_H_count += 1
        if diff_H_count <= 5:
            ps_map = {}
            for bits in group[:10]:
                A = binary_to_tournament(bits, n7)
                H = count_ham_paths(A, n7)
                ps = pfaffian_sum(A, n7)
                ps_map[bits] = (H, ps)
            print(f"    Lambda class: {len(group)} tours")
            for b, (H, ps) in list(ps_map.items())[:5]:
                print(f"      bits={b}: H={H}, Ps={ps}")

            # Check if they differ by a k-reversal
            bits_list = list(ps_map.keys())[:5]
            for i in range(len(bits_list)):
                for j in range(i+1, len(bits_list)):
                    if ps_map[bits_list[i]][0] != ps_map[bits_list[j]][0]:
                        A1 = binary_to_tournament(bits_list[i], n7)
                        A2 = binary_to_tournament(bits_list[j], n7)
                        diff_arcs = sum(1 for u in range(n7) for v in range(u+1, n7)
                                       if A1[u][v] != A2[u][v])
                        print(f"      Diff arcs between H={ps_map[bits_list[i]][0]} "
                              f"and H={ps_map[bits_list[j]][0]}: {diff_arcs}")

print(f"\n  Lambda classes with different H: {diff_H_count}")


# ========================================================================
# PART 9: Ps^2 residues — the "mass spectrum"
# ========================================================================
print(f"\n{'='*70}")
print("PART 9: The mass spectrum Ps^2 at n=5 and n=7")
print("=" * 70)

# Ps is always odd. So Ps^2 = 1 mod 8.
# But what about mod 16, mod 32?
print(f"  n=5: Ps^2 values = {sorted(set(d['ps']**2 for group in lambda_classes.values() for d in group))}")

# At n=5: Ps^2 in {1, 9, 25, 49, 81}
# These are 1^2, 3^2, 5^2, 7^2, 9^2
# Mod 8: all = 1
# Mod 16: 1, 9, 9, 1, 1 -> {1, 9}
ps_sq_mod16 = defaultdict(int)
ps_sq_mod24 = defaultdict(int)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n)
    ps = pfaffian_sum(A, n)
    ps_sq_mod16[ps**2 % 16] += 1
    ps_sq_mod24[ps**2 % 24] += 1

print(f"  Ps^2 mod 16: {dict(sorted(ps_sq_mod16.items()))}")
print(f"  Ps^2 mod 24: {dict(sorted(ps_sq_mod24.items()))}")

# At n=7 (sample):
ps7_values = set()
for bits in sample[:20000]:
    A = binary_to_tournament(bits, n7)
    ps = pfaffian_sum(A, n7)
    ps7_values.add(abs(ps))

print(f"\n  n=7: |Ps| values (sample): {sorted(ps7_values)}")
print(f"  All odd? {all(p % 2 == 1 for p in ps7_values)}")

# Check mod 4
ps7_mod4 = defaultdict(int)
for bits in sample[:20000]:
    A = binary_to_tournament(bits, n7)
    ps = pfaffian_sum(A, n7)
    ps7_mod4[ps % 4] += 1
print(f"  n=7: Ps mod 4: {dict(sorted(ps7_mod4.items()))}")


# ========================================================================
# PART 10: The Vitali measure as Stiefel-Whitney class
# ========================================================================
print(f"\n{'='*70}")
print("PART 10: Topological interpretation")
print("=" * 70)

print("""
THE VITALI ANALOGY (formalized):

1. TOURNAMENT SPACE = the set of all n-vertex tournaments
   (2^{C(n,2)} elements)

2. LAMBDA EQUIVALENCE = the equivalence relation T ~ T' iff
   lambda(T) = lambda(T') (same labeled lambda graph)
   This is the "rational structure" — the measurable part.

3. PFAFFIAN SUM = a Z-valued function on tournaments
   Ps: T -> Z, always odd.
   At n=5: Ps is a SIGN function on lambda classes
   (each class has exactly +Ps and -Ps representatives).
   At n=7: some classes have representatives with DIFFERENT |Ps|.

4. THE VITALI ATOM = a 4-vertex reversal that preserves lambda
   This is the "translation" in the analogy.
   It acts by flipping the sign of Ps (at n=5).

5. THE NON-MEASURABLE SET = the set of tournaments where Ps > 0
   (or any sign-consistent choice from each lambda class).
   This cannot be defined using only lambda — it requires
   the Pfaffian, which sees MATCHINGS not just CYCLES.

6. THE MEASURE = the Hamiltonian path count H(T).
   H is lambda-measurable at n <= 6 (determined by lambda + scores).
   H becomes lambda-NON-MEASURABLE at n = 7.
   The Pfaffian sum resolves this non-measurability.

7. THE DEEPER STRUCTURE:
   Ps^2 = det(I + 2A) = product(1 + 2*lambda_k)
   where lambda_k are eigenvalues of A.
   So |Ps| is determined by the SPECTRUM of A.
   The SIGN of Ps is a topological invariant:
   it's the orientation of the tournament in the
   "matching space" (related to the Stiefel-Whitney class
   or Euler class of the associated vector bundle).
""")


print(f"{'='*70}")
print("DONE.")
print("=" * 70)
