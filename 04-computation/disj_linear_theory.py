#!/usr/bin/env python3
"""
disj_linear_theory.py -- kind-pasteur-2026-03-13-S60

WHY is disj(k1,k2) LINEAR in c_k for circulant tournaments?

THEORY: For circulant tournament on Z_p with connection set S:
- c_k = number of directed k-cycles
- Every k-cycle on vertex set V has a "type" determined by the gap sequence
  (differences between consecutive vertices mod p in the cyclic ordering)
- The number of directed cycles on a specific k-element subset V depends only
  on the "difference set" of V (multiset of pairwise differences)

KEY OBSERVATION: In a circulant tournament, two cycles on vertex sets V1, V2
are disjoint iff V1 ∩ V2 = ∅. This depends on the VERTEX SETS, not on
which directed cycles they support.

So: disj(k1,k2) = sum over (V1, V2) with |V1|=k1, |V2|=k2, V1∩V2=∅
                   of n_cycles(V1) * n_cycles(V2)

where n_cycles(V) = number of directed Hamiltonian cycles on the subtournament
induced by V.

This is a SUM OF PRODUCTS. For it to be linear in total c_k, we need:

  sum_{V1∩V2=∅} n(V1)*n(V2) = linear function of {sum_V n(V) : |V|=k}

This is NOT obvious. Let's test: define n(V) for each k-vertex set, and see
if the "disjoint product sum" is determined by the "marginal sums" (total c_k).

At p=7: all C(7,3)=35 three-element subsets have SOME directed 3-cycles.
Since c3=14 is constant, each 3-element set has exactly 14/35 = 0.4 cycles on avg.
But actually each 3-element set has 0 or 2 directed 3-cycles (a 3-cycle and its reverse).
So 7 sets have 2 cycles and 28 sets have 0 cycles.

Wait, c3=14 and there are C(7,3)=35 possible 3-subsets.
14 directed 3-cycles / 2 (each set supports 0 or 2) = 7 subsets with cycles.
35 - 7 = 28 subsets without.

disj(3,3) counts pairs (C1, C2) of DIRECTED 3-cycles with disjoint vertex sets.
Since each active 3-subset has exactly 2 directed cycles, and we're counting
ordered-within-set pairs:
  For Paley p=7: disj(3,3) = 7 (from H_ck_theory output)
  7 active 3-subsets, C(7,2)=21 pairs of 3-subsets total
  Of these, how many are disjoint? 7 vertices, choose 3+3=6, so C(7,6)*...
  Actually: 7 active 3-subsets. Each pair of disjoint active subsets contributes
  2*2=4 to directed disj count (C1,C2 unordered but directed cycles ordered).
  Wait, for SAME length, we count unordered pairs of DIRECTED cycles.
  So each pair of disjoint active 3-subsets contributes C(2,2) + 2*1 = ...
  No: 2 directed cycles on V1, 2 on V2. Unordered pairs from different sets:
  2*2 = 4 ordered pairs, but we want unordered: it's just 2*2 = 4 pairs where
  we pick one from each set. But since k1=k2=3, we need unordered:
  Number of unordered pairs of DIRECTED 3-cycles from two disjoint sets = 2*2/1 = 4?
  No, if V1 has cycles {A, A'} and V2 has {B, B'}, unordered pairs are:
  (A,B), (A,B'), (A',B), (A',B') = 4 pairs. But we're counting unordered
  pairs of DIRECTED cycles: these are 4 distinct unordered pairs.

  So disj(3,3) = (number of disjoint active 3-subset pairs) * 4.
  Paley: disj(3,3)=7, so disjoint active pairs = 7/4 = 1.75. Not integer!

  Something's wrong. Let me recount from the actual script output.
  Actually, let me just compute this directly.
"""

from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def deep_analysis(p, S, name):
    """Detailed analysis of n(V) distribution and disjoint structure."""
    A = build_adj(p, S)

    print(f"\n{'='*70}")
    print(f"  {name}, p={p}, S={S}")
    print(f"{'='*70}")

    for k in range(3, p + 1, 2):
        # Count n(V) for each k-subset
        nV = {}
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            nV[fs] = nc

        total_ck = sum(nV.values())
        active = {V: n for V, n in nV.items() if n > 0}

        # Distribution of n(V)
        dist = defaultdict(int)
        for n in nV.values():
            dist[n] += 1

        print(f"\n  k={k}: c_{k}={total_ck}, {len(active)}/{len(nV)} active subsets")
        print(f"    n(V) distribution: {dict(sorted(dist.items()))}")

        # By circulant symmetry, n(V) depends only on the difference multiset
        # Check: do all translates V+t have the same n(V)?
        translate_check = True
        for subset in combinations(range(p), k):
            fs0 = frozenset(subset)
            nc0 = nV[fs0]
            for t in range(1, p):
                fs_t = frozenset((v + t) % p for v in subset)
                if nV[fs_t] != nc0:
                    translate_check = False
                    break
            if not translate_check:
                break

        print(f"    Translation-invariant n(V)? {translate_check}")

        if translate_check:
            # n(V) is constant on orbits of Z_p.
            # Each orbit has p members (since p is prime and k < p).
            # Number of orbits = C(p,k)/p
            n_orbits = len(nV) // p
            print(f"    {n_orbits} orbits of size {p}")

            # Show n(V) per orbit
            seen_orbits = set()
            orbit_nV = []
            for subset in combinations(range(p), k):
                fs = frozenset(subset)
                # Canonical form: minimal translate
                canon = min(frozenset((v + t) % p for v in subset) for t in range(p))
                if canon not in seen_orbits:
                    seen_orbits.add(canon)
                    orbit_nV.append((canon, nV[fs]))

            orbit_nV.sort(key=lambda x: (-x[1], x[0]))
            if len(orbit_nV) <= 20:
                for orb, nc in orbit_nV:
                    print(f"      orbit {sorted(orb)}: n(V)={nc}")

    # Now the KEY test: for disjoint pairs, does the structure factor?
    # disj(k1,k2) = sum_{V1 cap V2 = 0} n(V1) * n(V2)
    # If n(V) were constant (= c_k / C(p,k)), then
    #   disj(k1,k2) = (c_k1/C(p,k1)) * (c_k2/C(p,k2)) * D(p,k1,k2)
    # where D(p,k1,k2) = number of disjoint (k1,k2)-subset pairs.
    # This would make disj quadratic in c_k.
    #
    # But n(V) is NOT constant — different orbits have different n(V).
    # The linearity must come from the CORRELATION between n(V1) and
    # the disjoint structure being special for circulant tournaments.

    print(f"\n  --- DISJOINT PAIR DECOMPOSITION (p={p}) ---")

    # For (3,3) pairs specifically
    if 3 + 3 <= p:
        k = 3
        nV = {}
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            nV[fs] = nc

        active_sets = [V for V, n in nV.items() if n > 0]

        # Count disjoint active pairs and their contribution
        disj_count = 0  # count of disjoint (directed cycle) pairs
        for i in range(len(active_sets)):
            for j in range(i + 1, len(active_sets)):
                if not (active_sets[i] & active_sets[j]):
                    n1 = nV[active_sets[i]]
                    n2 = nV[active_sets[j]]
                    disj_count += n1 * n2
        # Also self-pairs for same vertex set: C(n(V), 2)
        # No: different vertex sets only for disjoint
        # For same vertex set, V1=V2 means V1 ∩ V2 = V1 ≠ ∅
        # So we only consider V1 ≠ V2

        # For directed cycles FROM the same vertex set:
        # These always conflict (same vertices), so contribute 0 to disj

        # Wait: the original code stores each directed cycle as a frozenset
        # and counts unordered pairs of directed cycles.
        # If V1 has n1 directed cycles and V2 has n2, and V1∩V2=∅,
        # then the number of disjoint directed cycle pairs = n1*n2.
        # For V1=V2 (same set), disjoint count = 0.
        # For k1=k2, unordered pairs: n1*n2/1 since V1≠V2 and |V1|=|V2|=k.
        # But we need to be careful about the factor of 2.

        # The original code loops i < j over the LIST of directed cycles.
        # If V1 has 2 cycles (call them a,a') and V2 has 2 cycles (b,b'):
        # The list is [a(V1), a'(V1), ..., b(V2), b'(V2), ...]
        # Pairs: (a,b), (a,b'), (a',b), (a',b') = 4 unordered pairs
        # All with disjoint vertex sets.
        # So contribution = n1*n2 = 2*2 = 4 per disjoint pair of active sets.

        # For k=3, n(V) is always 0 or 2 (cycle + reverse).
        # So contribution per disjoint active pair = 2*2 = 4.

        n_disj_active_pairs = sum(
            1 for i in range(len(active_sets))
            for j in range(i + 1, len(active_sets))
            if not (active_sets[i] & active_sets[j])
        )

        print(f"\n  (3,3) pairs:")
        print(f"    Active 3-subsets: {len(active_sets)}/{len(nV)}")
        print(f"    Disjoint active pairs: {n_disj_active_pairs}")
        print(f"    disj(3,3) = sum n1*n2 = {disj_count}")
        print(f"    Contribution per pair: {disj_count/n_disj_active_pairs if n_disj_active_pairs else 'N/A'}")

        # Now the key: how many disjoint 3-subset pairs are there in Z_p?
        all_disj_pairs = sum(
            1 for V1 in nV for V2 in nV
            if V1 < V2 and not (V1 & V2)
        )

        # By orbit: group 3-subsets by their orbit
        seen_orbits = {}
        orbit_of = {}
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            canon = min(frozenset((v + t) % p for v in subset) for t in range(p))
            orbit_of[fs] = canon
            if canon not in seen_orbits:
                seen_orbits[canon] = nV[fs]

        # Cross-orbit disjoint counts
        orbit_ids = sorted(seen_orbits.keys())
        print(f"\n    Orbits: {len(orbit_ids)}")
        for oid in orbit_ids:
            print(f"      orbit {sorted(oid)}: n(V)={seen_orbits[oid]}")

        # For each pair of orbits, count disjoint pairs
        print(f"\n    Cross-orbit disjoint pair counts:")
        for i, o1 in enumerate(orbit_ids):
            for j, o2 in enumerate(orbit_ids):
                if i > j:
                    continue
                members1 = [V for V in nV if orbit_of[V] == o1]
                members2 = [V for V in nV if orbit_of[V] == o2]
                count = 0
                if o1 == o2:
                    for a in range(len(members1)):
                        for b in range(a + 1, len(members1)):
                            if not (members1[a] & members1[b]):
                                count += 1
                else:
                    for V1 in members1:
                        for V2 in members2:
                            if not (V1 & V2):
                                count += 1
                n1, n2 = seen_orbits[o1], seen_orbits[o2]
                contrib = count * n1 * n2
                if count > 0:
                    print(f"      ({sorted(o1)}, {sorted(o2)}): "
                          f"n1={n1}, n2={n2}, disj_pairs={count}, contrib={contrib}")


# Run analysis
for p_val in [7]:
    m = (p_val - 1) // 2
    S_qr = sorted(j for j in range(1, p_val) if pow(j, (p_val - 1) // 2, p_val) == 1)
    S_int = list(range(1, m + 1))

    deep_analysis(p_val, S_qr, "Paley")
    deep_analysis(p_val, S_int, "Interval")

for p_val in [11]:
    m = (p_val - 1) // 2
    S_qr = sorted(j for j in range(1, p_val) if pow(j, (p_val - 1) // 2, p_val) == 1)
    deep_analysis(p_val, S_qr, "Paley")

print("\nDONE.")
