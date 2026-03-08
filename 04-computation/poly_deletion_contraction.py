"""
poly_deletion_contraction.py
kind-pasteur-2026-03-07-S35

Polynomial-level deletion-contraction for the forward-edge polynomial F(T,x).

For tournament T with edge e=(u->v), define:
  G_{u,v}(x) = sum_{P: u immed. before v in P} x^{asc(P)-1}
    (asc counts all forward T-edges including u->v; we subtract 1 for the u->v step)

Then: F_T(x) = F_{T\\e}(x) + (x-1) * G_{u,v}(x)

where T\\e is T with edge u->v removed (a non-tournament digraph).

At x=1: F_T(1) = n! = F_{T\\e}(1)  (consistent: (x-1) factor vanishes)
At x^{n-1} coeff: [x^{n-1}]F_T = H(T), [x^{n-1}]F_{T\\e} = H(T\\e),
  [x^{n-1}](x-1)*G = G_{n-2} = H(T/e) = # Ham paths using u->v

So the polynomial identity generalizes THM-082 to the polynomial level.

This script verifies the polynomial identity and studies G_{u,v}(x).
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


def compute_F(adj, n):
    """F(T,x) = sum_P x^{asc_T(P)}. Returns coefficient list [F_0, F_1, ..., F_{n-1}]."""
    F = [0] * n
    for P in permutations(range(n)):
        asc = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[asc] += 1
    return F


def compute_F_digraph(adj, n):
    """For general digraph: asc counts edges present in direction of traversal."""
    F = [0] * n
    for P in permutations(range(n)):
        asc = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[asc] += 1
    return F


def compute_G_uv(adj, n, u, v):
    """
    G_{u,v}(x) = sum over permutations P with u immediately before v,
    of x^{asc(P)-1}, where asc counts T-forward edges INCLUDING u->v.
    Returns coefficient list [G_0, ..., G_{n-2}].
    """
    G = [0] * (n - 1)
    for P in permutations(range(n)):
        # Check if u is immediately before v in P
        for i in range(n - 1):
            if P[i] == u and P[i + 1] == v:
                asc = sum(1 for j in range(n - 1) if adj[P[j]][P[j + 1]])
                # asc includes the u->v step (which is always forward since adj[u][v]=1)
                # We want x^{asc-1}
                G[asc - 1] += 1
                break
    return G


def poly_sub(A, B):
    """A - B as polynomial coefficient lists."""
    n = max(len(A), len(B))
    result = [0] * n
    for i in range(len(A)):
        result[i] += A[i]
    for i in range(len(B)):
        result[i] -= B[i]
    return result


def poly_mul_xm1(G):
    """Multiply polynomial G by (x-1). Returns list of length len(G)+1."""
    n = len(G)
    result = [0] * (n + 1)
    for i in range(n):
        result[i + 1] += G[i]  # x * G
        result[i] -= G[i]      # -G
    return result


def test_polynomial_dc(n):
    """Verify F_T(x) = F_{T\\e}(x) + (x-1) * G_{u,v}(x) for all tournaments and edges."""
    num_bits = n * (n - 1) // 2
    num_T = 1 << num_bits

    total_tests = 0
    pass_count = 0

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F_T = compute_F(adj, n)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                # Compute F(T\e)
                adj_del = [row[:] for row in adj]
                adj_del[u][v] = 0
                F_del = compute_F_digraph(adj_del, n)

                # Compute G_{u,v}
                G = compute_G_uv(adj, n, u, v)

                # Check F_T = F_del + (x-1)*G
                xm1_G = poly_mul_xm1(G)
                # xm1_G has length n, F_del has length n, F_T has length n
                # (x-1)*G has degree n-1 same as F_T

                rhs = [0] * n
                for i in range(n):
                    rhs[i] = F_del[i]
                for i in range(min(n, len(xm1_G))):
                    rhs[i] += xm1_G[i]

                total_tests += 1
                if F_T == rhs:
                    pass_count += 1
                else:
                    print(f"FAIL: bits={bits}, e=({u},{v})")
                    print(f"  F_T = {F_T}")
                    print(f"  F_del = {F_del}")
                    print(f"  G = {G}")
                    print(f"  (x-1)*G = {xm1_G}")
                    print(f"  rhs = {rhs}")

    print(f"\nn={n}: Polynomial DC: {pass_count}/{total_tests} pass")
    return pass_count == total_tests


def study_G_structure(n):
    """Study the structure of G_{u,v}(x) across tournaments."""
    num_bits = n * (n - 1) // 2
    num_T = 1 << num_bits

    # Is G palindromic? Anti-palindromic?
    palindromic_count = 0
    anti_palindromic_count = 0
    total_tests = 0

    # Sample G values
    G_values = set()

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                G = compute_G_uv(adj, n, u, v)
                total_tests += 1
                G_values.add(tuple(G))

                # Check palindrome: G_k = G_{n-2-k}
                deg = n - 2
                is_pal = all(G[k] == G[deg - k] for k in range(deg + 1))
                if is_pal:
                    palindromic_count += 1

                # Check anti-palindrome: G_k = -G_{n-2-k}
                is_anti = all(G[k] == -G[deg - k] for k in range(deg + 1))
                if is_anti:
                    anti_palindromic_count += 1

    print(f"\nn={n}: G_{'{u,v}'} structure:")
    print(f"  Total (T,u,v) triples: {total_tests}")
    print(f"  Distinct G polynomials: {len(G_values)}")
    print(f"  Palindromic: {palindromic_count}/{total_tests}")
    print(f"  Anti-palindromic: {anti_palindromic_count}/{total_tests}")


def study_flip_polynomial(n):
    """
    For arc flip e=(u->v) -> e'=(v->u):
    F_T(x) - F_{T'}(x) = (x-1) * [G_{u,v}^T(x) - G_{v,u}^{T'}(x)]

    where G_{u,v}^T is computed in T, G_{v,u}^{T'} is computed in T'.
    But G_{v,u}^{T'} sums over perms with v immediately before u in T'.

    Can we relate G_{u,v}^T and G_{v,u}^{T'} more directly?
    """
    num_bits = n * (n - 1) // 2
    num_T = 1 << num_bits

    total = 0
    g_equal = 0
    diff_palindrome = 0

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                # T has u->v
                G_uv_T = compute_G_uv(adj, n, u, v)

                # T' = flip u->v to v->u
                adj_flip = [row[:] for row in adj]
                adj_flip[u][v] = 0
                adj_flip[v][u] = 1
                G_vu_Tp = compute_G_uv(adj_flip, n, v, u)

                total += 1
                if G_uv_T == G_vu_Tp:
                    g_equal += 1

                # Check if D = G_uv_T - G_vu_Tp is anti-palindromic
                D = poly_sub(G_uv_T, G_vu_Tp)
                deg = n - 2
                is_anti = all(D[k] == -D[deg - k] for k in range(deg + 1))
                if is_anti:
                    diff_palindrome += 1

    print(f"\nn={n}: Flip polynomial analysis:")
    print(f"  Total arc flips: {total}")
    print(f"  G_uv^T = G_vu^T': {g_equal}/{total}")
    print(f"  D = G_uv - G_vu anti-palindromic: {diff_palindrome}/{total}")


def connect_G_to_contraction_F(n):
    """
    Is G_{u,v}(x) related to F(T/e, x) in any clean way?

    T/e has n-1 vertices with w = merged (u,v).
    G_{u,v}(x) = sum_{P: ...u,v...} x^{asc(P)-1}

    The permutations with u immediately before v biject with
    permutations of T/e: replace ...u,v... with ...w...

    But the asc count in T vs T/e differs because:
    - Edges to/from w in T/e depend on u's in-profile and v's out-profile
    - The edge BEFORE u (say z->u) maps to z->w, which has T/e[z][w] = T[z][u]
    - The edge AFTER v (say v->y) maps to w->y, which has T/e[w][y] = T[v][y]

    So asc(P_T) - 1 = asc(P_{T/e}) iff the before/after edges have the same
    forward/backward status in both T and T/e. Since T/e inherits those edges
    directly, this IS the case!

    Therefore: G_{u,v}(x) = F(T/e, x)  (exactly!)
    """
    num_bits = n * (n - 1) // 2
    num_T = 1 << num_bits

    total = 0
    match = 0

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue

                G = compute_G_uv(adj, n, u, v)

                # Contract T/e
                others = [x for x in range(n) if x != u and x != v]
                cn = n - 1
                con_adj = [[0]*cn for _ in range(cn)]
                omap = {x: i+1 for i, x in enumerate(others)}
                for x in others:
                    for y in others:
                        con_adj[omap[x]][omap[y]] = adj[x][y]
                for x in others:
                    if adj[x][u]:
                        con_adj[omap[x]][0] = 1
                    if adj[v][x]:
                        con_adj[0][omap[x]] = 1

                F_con = compute_F_digraph(con_adj, cn)

                total += 1
                if G == F_con:
                    match += 1

    print(f"\nn={n}: G_{{u,v}}(x) = F(T/e, x)?  {match}/{total}")
    if match != total:
        print("  NOT equal — the forward-edge count depends on vertex labeling")


# Run tests
print("=" * 60)
print("Polynomial Deletion-Contraction for F(T,x)")
print("=" * 60)

for n in [4, 5]:
    test_polynomial_dc(n)

print("\n" + "=" * 60)
print("G structure analysis")
print("=" * 60)
for n in [4, 5]:
    study_G_structure(n)

print("\n" + "=" * 60)
print("Flip polynomial analysis")
print("=" * 60)
for n in [4, 5]:
    study_flip_polynomial(n)

print("\n" + "=" * 60)
print("G vs F(T/e) comparison")
print("=" * 60)
for n in [4, 5]:
    connect_G_to_contraction_F(n)

print("\nDONE")
