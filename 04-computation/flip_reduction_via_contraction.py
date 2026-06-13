"""
flip_reduction_via_contraction.py
kind-pasteur-2026-03-07-S35

From THM-082 deletion-contraction: H(T) = H(T\e) + H(T/e).
Key corollary: For an arc flip e=(u->v) to e'=(v->u) giving T -> T':
  H(T) - H(T') = H(T/e) - H(T'/e')

since T\e = T'\e' (same digraph: T without the u-v edge).

This script verifies this identity and studies what H(T/e) - H(T'/e') looks like.

The two contractions T/e and T'/e' differ ONLY in the edges involving w:
  T/e:   (x,w) iff (x,u),  (w,x) iff (v,x)
  T'/e': (x,w) iff (x,v),  (w,x) iff (u,x)

So T/e and T'/e' are related by swapping u's in-role and v's out-role.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
from itertools import permutations
from collections import defaultdict


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def ham_paths_dp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def contract_edge(adj, n, u, v):
    """Contract e=(u->v): w gets IN from u, OUT from v."""
    others = [x for x in range(n) if x != u and x != v]
    new_n = n - 1
    # w = index 0, others = indices 1,2,...
    omap = {x: i+1 for i, x in enumerate(others)}

    new_adj = [[0]*new_n for _ in range(new_n)]
    for x in others:
        for y in others:
            new_adj[omap[x]][omap[y]] = adj[x][y]
    for x in others:
        if adj[x][u]:  # (x,w) iff (x,u)
            new_adj[omap[x]][0] = 1
        if adj[v][x]:  # (w,x) iff (v,x)
            new_adj[0][omap[x]] = 1
    return new_adj, new_n


def flip_arc(adj, n, u, v):
    """Return tournament with arc u->v flipped to v->u."""
    new_adj = [row[:] for row in adj]
    new_adj[u][v] = 0
    new_adj[v][u] = 1
    return new_adj


def test_flip_reduction(n):
    """Verify H(T) - H(T') = H(T/e) - H(T'/e') for all tournaments and all arcs."""
    num_bits = n * (n - 1) // 2
    num_T = 1 << num_bits

    total_tests = 0
    identity_holds = 0
    max_delta = 0
    delta_dist = defaultdict(int)

    # Also track: what does T/e look like? Is it always a tournament?
    contraction_is_tournament = 0
    contraction_not_tournament = 0

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                # Flip arc u->v to v->u
                adj_flip = flip_arc(adj, n, u, v)
                H_T_flip = ham_paths_dp(adj_flip, n)

                # Contract e=(u->v) in T
                con_adj, cn = contract_edge(adj, n, u, v)
                H_con = ham_paths_dp(con_adj, cn)

                # Contract e'=(v->u) in T'
                con_adj_flip, cn2 = contract_edge(adj_flip, n, v, u)
                H_con_flip = ham_paths_dp(con_adj_flip, cn2)

                # Check identity
                lhs = H_T - H_T_flip
                rhs = H_con - H_con_flip

                total_tests += 1
                if lhs == rhs:
                    identity_holds += 1

                delta_dist[lhs] += 1
                if abs(lhs) > max_delta:
                    max_delta = abs(lhs)

                # Check if contraction is a tournament
                is_tourn = True
                for i in range(cn):
                    for j in range(i+1, cn):
                        if con_adj[i][j] + con_adj[j][i] != 1:
                            is_tourn = False
                            break
                    if not is_tourn:
                        break
                if is_tourn:
                    contraction_is_tournament += 1
                else:
                    contraction_not_tournament += 1

    print(f"\n=== n={n}: Flip Reduction via Contraction ===")
    print(f"Total tests: {total_tests}")
    print(f"H(T)-H(T') = H(T/e)-H(T'/e'): {identity_holds}/{total_tests} ({100*identity_holds/total_tests:.1f}%)")
    print(f"Max |H(T)-H(T')|: {max_delta}")
    print(f"Contraction is tournament: {contraction_is_tournament}/{total_tests}")
    print(f"Contraction NOT tournament: {contraction_not_tournament}/{total_tests}")

    # Show delta distribution
    print(f"\nDelta = H(T)-H(T') distribution:")
    for d in sorted(delta_dist.keys()):
        print(f"  delta={d:+4d}: {delta_dist[d]} cases")

    return identity_holds == total_tests


def study_contraction_pair(n):
    """
    Study the two contractions T/e and T'/e' more closely.
    They differ only in how w connects to other vertices.
    """
    num_bits = n * (n - 1) // 2

    print(f"\n=== n={n}: Contraction pair structure ===")

    # How do T/e and T'/e' relate? They have the same vertex set and
    # same edges among non-w vertices. They differ only in w's edges.

    # For x != u,v:
    #   T/e:   (x,w) iff (x,u),  (w,x) iff (v,x)
    #   T'/e': (x,w) iff (x,v),  (w,x) iff (u,x)
    #
    # Recall T' = T with u<->v flipped. So:
    #   T'/e': (x,w) iff T'(x,v)=T(x,v),  (w,x) iff T'(u,x)=T(u,x) for x!=u,v
    #   WAIT: T' flips only the u-v edge. For x != u,v:
    #     T'(x,v) = T(x,v), T'(u,x) = T(u,x)
    #   So T'/e': (x,w) iff T(x,v), (w,x) iff T(u,x)
    #
    # Compare with T/e: (x,w) iff T(x,u), (w,x) iff T(v,x)
    #
    # So the difference between T/e and T'/e' is:
    #   IN to w: T/e uses u's in-profile, T'/e' uses v's in-profile
    #   OUT from w: T/e uses v's out-profile, T'/e' uses u's out-profile
    #
    # In other words: T/e merges (IN=u, OUT=v), T'/e' merges (IN=v, OUT=u).

    # Interesting: if u and v have the same in/out profiles to other vertices,
    # then T/e = T'/e' and H(T) = H(T').

    # Count how often in-profile equality implies H-equality
    profile_eq = 0
    profile_eq_H = 0
    profile_neq = 0

    for bits in range(1 << num_bits):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                # Check if u and v have same profiles to others
                others = [x for x in range(n) if x != u and x != v]
                in_same = all(adj[x][u] == adj[x][v] for x in others)
                out_same = all(adj[u][x] == adj[v][x] for x in others)

                if in_same and out_same:
                    profile_eq += 1
                    adj_flip = flip_arc(adj, n, u, v)
                    H_flip = ham_paths_dp(adj_flip, n)
                    if H_T == H_flip:
                        profile_eq_H += 1
                else:
                    profile_neq += 1

    print(f"Profile u=v (in AND out to others): {profile_eq}")
    if profile_eq > 0:
        print(f"  Of these, H(T)=H(T'): {profile_eq_H}/{profile_eq}")
    print(f"Profile u!=v: {profile_neq}")


def study_contraction_structure(n):
    """
    At what point does T/e stop being a tournament?
    A tournament contraction T/e is a tournament iff for every x != u,v:
      exactly one of (x,w) or (w,x) holds.
    This means: for every x, exactly one of T(x,u) or T(v,x) holds.
    But T(x,u) + T(u,x) = 1 and T(x,v) + T(v,x) = 1.
    So (x,w) = T(x,u), (w,x) = T(v,x).
    For T/e to be a tournament: T(x,u) + T(v,x) = 1 for all x != u,v.
    This means T(x,u) = T(x,v) for all x (since T(v,x) = 1 - T(x,v)).
    I.e., u and v have the SAME in-profile from all other vertices.
    Equivalently: u and v have the same SCORE among others, with identical beaters.
    """
    num_bits = n * (n - 1) // 2

    tour_count = 0
    non_tour_count = 0

    for bits in range(1 << num_bits):
        adj = tournament_from_bits(n, bits)
        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                others = [x for x in range(n) if x != u and x != v]
                is_tourn = all(adj[x][u] + adj[v][x] == 1 for x in others)
                if is_tourn:
                    tour_count += 1
                else:
                    non_tour_count += 1

    total = tour_count + non_tour_count
    print(f"\n=== n={n}: Contraction tournament check ===")
    print(f"T/e is tournament: {tour_count}/{total} ({100*tour_count/total:.1f}%)")
    print(f"T/e is NOT tournament: {non_tour_count}/{total} ({100*non_tour_count/total:.1f}%)")
    print(f"(Tournament iff T(x,u)=T(x,v) for all other x)")


# Run tests
for n in [4, 5]:
    test_flip_reduction(n)

study_contraction_pair(4)
study_contraction_structure(4)
study_contraction_structure(5)

print("\n=== SUMMARY ===")
print("THM-082: H(T) = H(T\\e) + H(T/e) — proved by bijection")
print("Corollary: H(T) - H(T') = H(T/e) - H(T'/e') — arc flip reduces to contraction difference")
print("T/e and T'/e' differ only in w's edges: T/e = (IN=u, OUT=v), T'/e' = (IN=v, OUT=u)")
