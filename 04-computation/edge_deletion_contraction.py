#!/usr/bin/env python3
"""
Edge deletion-contraction for the Redei-Berge function of tournaments.

From Mitrovic (arXiv:2504.20968), the noncommutative Redei-Berge function satisfies:
  W_X = W_{X\e} - W_{X/e}^up

For the commutative version (Mitrovic-Stojadinovic arXiv:2407.18608):
  U_X = sum_{S subset E, S nonempty} (-1)^{|S|-1} U_{X\S}

For a tournament T with edge e = (u,v):
  T\e is T minus the edge u->v (no longer a tournament!)
  T/e is the contraction: merge u,v into a single vertex

Key question: what does edge deletion-contraction look like in terms of
the OCF / independence polynomial?

When we delete edge e=(u,v) from tournament T:
  - T\e still has the edge v->u (since T has both directions between u,v)
    Wait, no: T is a tournament, so either u->v OR v->u, not both.
    If e=(u,v) is an edge, deleting it makes u,v have NO edge between them.
    So T\e is NOT a tournament - it's missing one edge.

When we contract edge e=(u,v):
  - Merge u,v into vertex w
  - For any other vertex x: edge w->x iff v->x in T; edge x->w iff x->u in T
  - The result is a tournament on n-1 vertices! (assuming u,v agree on most neighbors)
    Actually: for each x != u,v, we get w->x iff v->x, and x->w iff x->u.
    Since u->v, x->u, v->x gives a consistent path x->u->v->... through w.
    The contraction T/e IS a tournament on n-1 vertices.

This is interesting: T/e is a tournament, so H(T/e) is well-defined.
And H(T\e) is the number of Hamiltonian paths in the non-tournament T\e.

From the formula: ps_1(U_T)(1) = H(T)
  ps_1(U_{T\e})(1) - ps_1(U_{T/e})(1) = ps_1(U_T)(1) = H(T)
Wait, the formula is U_X = sum_S (-1)^{|S|-1} U_{X\S}, which is different from
the deletion-contraction W = W_{X\e} - W_{X/e}^up.

Let me just compute ps_1 values and check relationships.

opus-2026-03-07-S38
"""
from itertools import permutations, combinations
from collections import defaultdict


def ham_paths(n, adj):
    """Count Hamiltonian paths in digraph given as adjacency dict/set."""
    count = 0
    for p in permutations(range(n)):
        if all((p[i], p[i+1]) in adj for i in range(n-1)):
            count += 1
    return count


def make_tournament_edges(n, A):
    """Get edge set from adjacency matrix."""
    return {(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]}


def contract_edge(n, edges, u, v):
    """Contract edge (u,v): merge u into v.
    For other vertex x:
      w->x iff v->x in original
      x->w iff x->u in original
    Return (n-1, new_edges) with vertices relabeled to 0..n-2.
    """
    # New vertex set: all except u, with v acting as merged vertex
    new_verts = [x for x in range(n) if x != u]
    relabel = {x: i for i, x in enumerate(new_verts)}

    new_edges = set()
    for x in new_verts:
        for y in new_verts:
            if x == y:
                continue
            if x == v and y != v:
                # merged vertex -> y: use v->y
                if (v, y) in edges:
                    new_edges.add((relabel[x], relabel[y]))
            elif y == v and x != v:
                # x -> merged vertex: use x->u
                if (x, u) in edges:
                    new_edges.add((relabel[x], relabel[y]))
            else:
                # neither endpoint is merged
                if (x, y) in edges:
                    new_edges.add((relabel[x], relabel[y]))

    return n - 1, new_edges


def delete_edge(n, edges, u, v):
    """Delete edge (u,v) from digraph."""
    return n, edges - {(u, v)}


def count_odd_cycles(n, edges):
    """Count directed odd cycles by length."""
    counts = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            seen = set()
            for p in permutations(verts):
                if all((p[i], p[(i+1) % length]) in edges for i in range(length)):
                    min_idx = list(p).index(min(p))
                    canonical = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if canonical not in seen:
                        seen.add(canonical)
                        counts[length] += 1
    return dict(counts)


def count_alpha(n, edges):
    """Count alpha_1 (total odd cycles) and alpha_2 (disjoint pairs)."""
    all_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for p in permutations(verts):
                if all((p[i], p[(i+1) % length]) in edges for i in range(length)):
                    min_idx = list(p).index(min(p))
                    canonical = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if frozenset(verts) not in {frozenset(c) for c in all_cycles}:
                        all_cycles.append(canonical)
                    break

    # Actually need to track vertex sets
    cycle_vsets = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for p in permutations(verts):
                if all((p[i], p[(i+1) % length]) in edges for i in range(length)):
                    cycle_vsets.append(frozenset(verts))
                    break

    alpha1 = len(cycle_vsets)
    alpha2 = 0
    for i in range(len(cycle_vsets)):
        for j in range(i+1, len(cycle_vsets)):
            if len(cycle_vsets[i] & cycle_vsets[j]) == 0:
                alpha2 += 1

    return alpha1, alpha2


def main():
    # Cyclic tournament C_5
    n = 5
    edges = set()
    for i in range(n):
        edges.add((i, (i+1) % n))
        edges.add((i, (i+2) % n))

    print(f"=== Cyclic Tournament C_5 ===")
    H = ham_paths(n, edges)
    print(f"H(C_5) = {H}")
    a1, a2 = count_alpha(n, edges)
    print(f"alpha_1 = {a1}, alpha_2 = {a2}")
    print(f"OCF: 1 + 2*{a1} + 4*{a2} = {1 + 2*a1 + 4*a2}")

    # Edge deletion-contraction for each edge
    print(f"\n--- Edge deletion-contraction ---")
    for u, v in sorted(edges):
        # Delete edge (u,v)
        n_del, edges_del = delete_edge(n, edges, u, v)
        H_del = ham_paths(n_del, edges_del)

        # Contract edge (u,v)
        n_con, edges_con = contract_edge(n, edges, u, v)
        H_con = ham_paths(n_con, edges_con)

        # Check: is T/e a tournament?
        is_tourn = True
        for x in range(n_con):
            for y in range(x+1, n_con):
                has_xy = (x, y) in edges_con
                has_yx = (y, x) in edges_con
                if has_xy == has_yx:  # both or neither
                    is_tourn = False

        print(f"  e=({u},{v}): H(T\\e)={H_del}, H(T/e)={H_con}, "
              f"H(T\\e)-H(T/e)={H_del-H_con}, tournament={is_tourn}")

    print(f"\n  H(T) = {H}, should match H(T\\e) - H(T/e) for all e? "
          f"(Not necessarily - ps_1 doesn't commute with deletion-contraction simply)")

    # Paley T_7
    print(f"\n=== Paley Tournament T_7 ===")
    n = 7
    QR = {1, 2, 4}
    edges = {(i, j) for i in range(n) for j in range(n)
             if i != j and (j - i) % n in QR}
    H = ham_paths(n, edges)
    print(f"H(T_7) = {H}")

    # Sample a few edges
    edge_list = sorted(edges)[:7]  # first 7 edges
    print(f"\n--- Edge deletion-contraction (first 7 edges) ---")
    for u, v in edge_list:
        n_del, edges_del = delete_edge(n, edges, u, v)
        H_del = ham_paths(n_del, edges_del)

        n_con, edges_con = contract_edge(n, edges, u, v)
        H_con = ham_paths(n_con, edges_con)

        is_tourn = True
        for x in range(n_con):
            for y in range(x+1, n_con):
                has_xy = (x, y) in edges_con
                has_yx = (y, x) in edges_con
                if has_xy == has_yx:
                    is_tourn = False

        a1_con, a2_con = count_alpha(n_con, edges_con) if is_tourn else (None, None)

        print(f"  e=({u},{v}): H(T\\e)={H_del}, H(T/e)={H_con}, "
              f"diff={H_del-H_con}, tourn={is_tourn}", end="")
        if is_tourn:
            print(f", OCF(T/e)=1+2*{a1_con}+4*{a2_con}={1+2*a1_con+4*a2_con}")
        else:
            print()

    # Key question: does H(T\e) - H(T/e) = H(T) or something else?
    # The formula from Mitrovic is about U_X (symmetric function), not just ps_1.
    # At the ps_1(1) level: ps_1(U_T)(1) = H(T) for tournaments.
    # But ps_1(U_{T\e})(1) ≠ H(T\e) in general (T\e is not a tournament).
    # Actually, ps_1(U_X)(1) counts listings of V with no X-descent for any digraph X.
    # For tournament T, this equals H(T^op) = H(T).
    # For T\e: a listing with no (T\e)-descent means no consecutive pair (a,b)
    # with a->b in T\e. So u,v can appear consecutively with u before v
    # (since u->v was deleted). This is NOT the same as H(T\e).
    # ps_1(U_{T\e})(1) = #{listings with no T\e-descent} = #{HPs in complement of T\e}.

    print(f"\n=== Complement interpretation ===")
    n = 5
    edges = set()
    for i in range(n):
        edges.add((i, (i+1) % n))
        edges.add((i, (i+2) % n))

    comp_edges = {(i, j) for i in range(n) for j in range(n) if i != j} - edges
    print(f"C_5 complement (= C_5^op) has H = {ham_paths(n, comp_edges)}")

    for u, v in [(0, 1)]:
        n_del, edges_del = delete_edge(n, edges, u, v)
        comp_del = {(i, j) for i in range(n) for j in range(n) if i != j} - edges_del
        H_comp_del = ham_paths(n, comp_del)
        print(f"  After deleting ({u},{v}): comp has H = {H_comp_del}")
        print(f"  ps_1(U_(T\\e))(1) should be = {H_comp_del}")


if __name__ == "__main__":
    main()
