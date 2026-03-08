#!/usr/bin/env python3
"""
f_poly_flip_formula.py — How F(T,x) changes under a single arc flip.

For tournament T with arc u->v, let T' = flip(T, u->v) = T with arc reversed.
Then T' is still a tournament.

F(T,x) - F(T',x) depends only on paths using u,v consecutively.
Path ...u,v...: in T, this step is forward (x). In T', this step is backward (1).
  So contribution to F(T) is x^{fwd}, contribution to F(T') is x^{fwd-1}.
  Difference: x^{fwd} - x^{fwd-1} = x^{fwd-1}(x-1).

Path ...v,u...: in T, this step is backward (1). In T', this step is forward (x).
  So contribution to F(T) is x^{fwd}, contribution to F(T') is x^{fwd+1}.
  Difference: x^{fwd} - x^{fwd+1} = -x^{fwd}(x-1).

So: F(T,x) - F(T',x) = (x-1) * [sum_{P: ...u,v...} x^{fwd(P)-1} - sum_{P: ...v,u...} x^{fwd(P)}]

Let G_uv(x) = sum_{P: ...u,v...} x^{fwd(P)-1}  (shift by removing the u->v forward step)
Let G_vu(x) = sum_{P: ...v,u...} x^{fwd(P)}

Then F(T,x) - F(T',x) = (x-1) * [G_uv(x) - G_vu(x)]

By palindrome: F(T,x) = x^{n-1} F(T, 1/x), and similarly for T'.
So the difference is also palindromic in some sense.

QUESTION: What is G_uv(x) - G_vu(x)?

Author: opus-2026-03-07-S45
"""
from itertools import permutations, combinations
import math

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
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def flip(adj, n, u, v):
    B = [row[:] for row in adj]
    B[u][v] = 0
    B[v][u] = 1
    return B

def poly_mul(A, B):
    """Multiply two polynomials represented as coefficient lists."""
    n = len(A) + len(B) - 1
    C = [0]*n
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C[i+j] += a*b
    return C

# ============================================================
# TEST: F(T) - F(T') = (x-1) * (G_uv - G_vu)
# ============================================================
print("=" * 60)
print("F(T,x) - F(T',x) = (x-1) * [G_uv(x) - G_vu(x)]")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    pass_count = 0
    total = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F_T = compute_F(adj, n)

        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                total += 1

                adj_flip = flip(adj, n, u, v)
                F_flip = compute_F(adj_flip, n)

                diff = [F_T[k] - F_flip[k] for k in range(n)]

                # Compute G_uv and G_vu
                G_uv = [0]*(n-1)  # degree n-2 (shifted)
                G_vu = [0]*(n-1)  # degree n-2
                for P in permutations(range(n)):
                    fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
                    uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
                    uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
                    if uses_uv:
                        G_uv[fwd - 1] += 1  # shift down by 1
                    elif uses_vu:
                        G_vu[fwd] += 1

                # G_uv - G_vu
                max_deg = max(len(G_uv), len(G_vu))
                G_diff = [(G_uv[k] if k < len(G_uv) else 0) - (G_vu[k] if k < len(G_vu) else 0)
                          for k in range(max_deg)]

                # (x-1) * G_diff = [-G_diff[0], G_diff[0]-G_diff[1], ..., G_diff[-1]]
                prod = [0]*n
                for k in range(len(G_diff)):
                    prod[k] -= G_diff[k]      # -1 * G_diff
                    prod[k+1] += G_diff[k]    # x * G_diff

                if diff == prod:
                    pass_count += 1

    print(f"  {pass_count}/{total} pass")

# ============================================================
# EXPLORE: What is G_uv(x) - G_vu(x)?
# ============================================================
print("\n" + "=" * 60)
print("G_uv(x) - G_vu(x) structure")
print("=" * 60)

n = 4
m = n*(n-1)//2

seen = set()
for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)
    key = tuple(F_T)
    if key in seen:
        continue
    seen.add(key)

    for u in range(n):
        for v in range(n):
            if adj[u][v]:
                break
        else:
            continue
        break

    G_uv = [0]*(n-1)
    G_vu = [0]*(n-1)
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
        uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
        if uses_uv:
            G_uv[fwd - 1] += 1
        elif uses_vu:
            G_vu[fwd] += 1

    G_diff = [G_uv[k] - G_vu[k] for k in range(n-1)]

    # Is G_diff palindromic?
    d = n - 2
    palindrome = all(G_diff[k] == G_diff[d-k] for k in range((d+1)//2 + 1))
    # Is G_diff anti-palindromic?
    anti_palindrome = all(G_diff[k] == -G_diff[d-k] for k in range((d+1)//2 + 1))

    print(f"  bits={bits} arc {u}->{v}: G_uv={G_uv}, G_vu={G_vu}, diff={G_diff}, "
          f"palindrome={palindrome}, anti_palindrome={anti_palindrome}")

# ============================================================
# EXPLORE: G_uv(x) as contraction polynomial?
# ============================================================
print("\n" + "=" * 60)
print("G_uv(x) vs F(T/e, x) — is G_uv a contraction?")
print("=" * 60)

# G_uv[k] = # perms with ...u,v... and fwd = k+1 (after shift)
# This is the "contracted" polynomial if fwd counts match.
# From v2 analysis: they DON'T match in general due to edge direction changes.
# But G_uv IS well-defined as the shifted F_uv.

# The TOTAL sum(G_uv) + sum(G_vu) = # paths using u,v consecutively
# = H(T/e) by the H deletion-contraction identity.

n = 4
m = n*(n-1)//2

def ham_paths_dp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def contract_arc(adj, n, u, v):
    others = sorted([x for x in range(n) if x != u and x != v])
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    for i, x in enumerate(others):
        if adj[x][u]:
            B[i+1][0] = 1
        if adj[v][x]:
            B[0][i+1] = 1
    for i, x in enumerate(others):
        for j, y in enumerate(others):
            B[i+1][j+1] = adj[x][y]
    return B, new_n

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue

            G_uv = [0]*(n-1)
            G_vu = [0]*(n-1)
            for P in permutations(range(n)):
                fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
                uses_uv = any(P[i]==u and P[i+1]==v for i in range(n-1))
                uses_vu = any(P[i]==v and P[i+1]==u for i in range(n-1))
                if uses_uv:
                    G_uv[fwd - 1] += 1
                elif uses_vu:
                    G_vu[fwd] += 1

            con_adj, con_n = contract_arc(adj, n, u, v)
            is_t = all(con_adj[i][j] + con_adj[j][i] == 1
                       for i in range(con_n) for j in range(i+1, con_n))

            if is_t:
                F_con = compute_F(con_adj, con_n)
            else:
                F_con = None

            total_through = sum(G_uv) + sum(G_vu)
            H_con = ham_paths_dp(con_adj, con_n)

            if total_through != H_con:
                print(f"  MISMATCH: bits={bits} arc {u}->{v}: through={total_through}, H_con={H_con}")
                break
    else:
        continue
    break
else:
    print(f"  sum(G_uv) + sum(G_vu) = H(T/e) for ALL arcs at n=4: PASS")

    # Now check: is G_uv + G_vu = F(T/e) as polynomials?
    ok = 0
    total = 0
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        for u in range(n):
            for v in range(n):
                if u == v or not adj[u][v]:
                    continue
                con_adj, con_n = contract_arc(adj, n, u, v)
                is_t = all(con_adj[i][j] + con_adj[j][i] == 1
                           for i in range(con_n) for j in range(i+1, con_n))
                if not is_t:
                    continue

                G_uv = [0]*(n-1)
                G_vu = [0]*(n-1)
                for P in permutations(range(n)):
                    fwd = sum(1 for i_p in range(n-1) if adj[P[i_p]][P[i_p+1]])
                    uses_uv = any(P[i_p]==u and P[i_p+1]==v for i_p in range(n-1))
                    uses_vu = any(P[i_p]==v and P[i_p+1]==u for i_p in range(n-1))
                    if uses_uv:
                        G_uv[fwd - 1] += 1
                    elif uses_vu:
                        G_vu[fwd] += 1

                F_con = compute_F(con_adj, con_n)
                G_sum = [G_uv[k] + G_vu[k] for k in range(con_n)]

                total += 1
                if G_sum == F_con:
                    ok += 1

    print(f"  G_uv + G_vu = F(T/e) as polynomials: {ok}/{total}")

# ============================================================
# KEY DISCOVERY: If G_uv + G_vu != F(T/e), what IS the correction?
# ============================================================
print("\n" + "=" * 60)
print("CORRECTION TERM: F(T/e) vs G_uv + G_vu")
print("=" * 60)

n = 4
m = n*(n-1)//2

for bits in range(min(16, 1 << m)):
    adj = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            con_adj, con_n = contract_arc(adj, n, u, v)
            is_t = all(con_adj[i][j] + con_adj[j][i] == 1
                       for i in range(con_n) for j in range(i+1, con_n))
            if not is_t:
                continue

            G_uv = [0]*(n-1)
            G_vu = [0]*(n-1)
            for P in permutations(range(n)):
                fwd = sum(1 for i_p in range(n-1) if adj[P[i_p]][P[i_p+1]])
                uses_uv = any(P[i_p]==u and P[i_p+1]==v for i_p in range(n-1))
                uses_vu = any(P[i_p]==v and P[i_p+1]==u for i_p in range(n-1))
                if uses_uv:
                    G_uv[fwd - 1] += 1
                elif uses_vu:
                    G_vu[fwd] += 1

            F_con = compute_F(con_adj, con_n)
            G_sum = [G_uv[k] + G_vu[k] for k in range(con_n)]
            correction = [F_con[k] - G_sum[k] for k in range(con_n)]

            if any(c != 0 for c in correction):
                print(f"  bits={bits} arc {u}->{v}: correction={correction}")
                print(f"    G_uv={G_uv[:con_n]}, G_vu={G_vu[:con_n]}, F_con={F_con}")

                # Count edge disagreements
                others = sorted([x for x in range(n) if x != u and x != v])
                disagree_in = [(x, adj[x][u], adj[x][v]) for x in others if adj[x][u] != adj[x][v]]
                disagree_out = [(x, adj[u][x], adj[v][x]) for x in others if adj[u][x] != adj[v][x]]
                print(f"    IN-disagree: {disagree_in}")
                print(f"    OUT-disagree: {disagree_out}")
