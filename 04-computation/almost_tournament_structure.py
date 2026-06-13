"""
almost_tournament_structure.py
kind-pasteur-2026-03-07-S38

Deep analysis of the "almost-tournament claim" needed for THM-086 proof:
  c_j(T\e) = 0 mod 3 for j < val(n) - 1

T\e is a tournament minus one arc. Can we express c_j(T\e) in terms of
tournament invariants to prove the claim?

Key observations:
1. T\e = T minus arc (u->v). Now neither u->v nor v->u exists.
2. F(T\e, x) counts permutations weighted by forward edges in T\e.
3. c_j(T\e) = sum_k C(k,j) F_k(T\e)

Strategy: Express F(T\e, x) in terms of F(T, x) and some correction.
- A permutation P in T\e has the same forward edges as in T, EXCEPT
  if P has u immediately before v (or v immediately before u).
- If u is at position i and v at position i+1 in P:
  In T: this contributes to F_{fwd+1} if u->v, or F_fwd if v->u.
  In T\e: u->v arc is deleted, so this contributes to F_fwd (no forward edge).
  Difference: each such P loses 1 from its forward count.

Let's formalize and test this.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from math import comb


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


def compute_F_dp_digraph(adj, n):
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if adj[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]
    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


def compute_N_uv(adj, n, u, v):
    """Count permutations where u is immediately before v, by forward count.
    N_uv[k] = #{P : P has u at position i, v at position i+1, and fwd(P in T) = k}
    """
    # DP tracking: last vertex, forward count, and mask
    # But we also need to know if the (u,v) adjacency contributed a forward edge
    # Simpler: compute F(T,x) restricted to perms with u immediately before v

    # Build all perms with u immediately before v
    # Treat (u,v) as a single block in positions (i, i+1)
    # The remaining n-2 elements go in the other n-2 positions
    # The block can be at positions (0,1), (1,2), ..., (n-2, n-1)

    from itertools import permutations
    N = [0] * n
    others = [x for x in range(n) if x != u and x != v]

    for pos in range(n - 1):  # block at (pos, pos+1)
        for perm in permutations(others):
            # Build full permutation
            P = list(perm[:pos]) + [u, v] + list(perm[pos:])
            # Count forward edges in T
            fwd = 0
            for i in range(n - 1):
                if adj[P[i]][P[i+1]]:
                    fwd += 1
            N[fwd] += 1

    return N


# ============================================================
# Test: F(T\e, x) = F(T, x) - shift of N_uv
# ============================================================
print("=" * 70)
print("Relationship: F(T\\e, x) vs F(T, x)")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    random.seed(42)

    for trial in range(5):
        bits = random.randint(0, (1 << m) - 1)
        adj = tournament_from_bits(n, bits)

        # Find arc u->v
        u, v = (0, 1) if adj[0][1] else (1, 0)

        # F(T, x)
        F_T = compute_F_dp_digraph(adj, n)

        # F(T\e, x) — delete arc u->v
        adj_del = [row[:] for row in adj]
        adj_del[u][v] = 0
        F_del = compute_F_dp_digraph(adj_del, n)

        # N_uv: perms with u immediately before v
        N_uv = compute_N_uv(adj, n, u, v)

        # Hypothesis: F_del[k] = F_T[k] + N_uv[k+1] - N_uv[k+1]?
        # Actually: when u->v is deleted, perms with u imm. before v lose 1 fwd edge
        # So F_del[k] = F_T[k] - N_uv[k] (lost from k) + N_uv[k+1] (gained from k+1 -> k)
        # Wait: the perms that had u before v with fwd=k in T now have fwd=k-1 in T\e
        # (since the u->v edge no longer counts). But only if adj[u][v]=1.
        # Since we chose u->v, adj[u][v]=1, so yes: each such perm loses 1.
        # F_del[k] = F_T[k] - N_uv[k] + N_uv[k+1]  (shift N down by 1)

        # But also: perms with v immediately before u are unaffected (adj[v][u]=0 in both).
        # Wait, in T\e, adj[u][v]=0 AND adj[v][u]=0 (neither arc exists).
        # In T, adj[v][u]=0 (since u->v). So v->u was already not a forward edge.
        # In T\e, v->u is STILL not a forward edge (adj_del[v][u] = 0).
        # Actually wait: in the original tournament adj[v][u] = 0 since u->v.
        # In T\e: adj_del[u][v] = 0, adj_del[v][u] = 0.
        # So perms with v before u: in T, adj[v][u]=0 so no fwd edge. In T\e, same. No change.
        # Perms with u before v: in T, adj[u][v]=1 so fwd edge. In T\e, adj_del[u][v]=0. Lost 1.

        # Therefore: F_del[k] = F_T[k] - N_uv[k] + N_uv[k+1]  for k >= 0
        # (where N_uv[n] = 0)

        F_pred = [0] * n
        for k in range(n):
            F_pred[k] = F_T[k] - N_uv[k] + (N_uv[k+1] if k+1 < n else 0)

        ok = (F_pred == F_del)
        if trial < 3 or not ok:
            print(f"\n  n={n}, trial {trial}: {'OK' if ok else 'FAIL'}")
            if not ok:
                print(f"    F_T   = {F_T}")
                print(f"    F_del = {F_del}")
                print(f"    F_pred= {F_pred}")
                print(f"    N_uv  = {N_uv}")


# ============================================================
# Express c_j(T\e) in terms of c_j(T) and correction from N_uv
# ============================================================
print("\n" + "=" * 70)
print("c_j(T\\e) = c_j(T) - correction from N_uv")
print("=" * 70)

# c_j(T\e) = sum_k C(k,j) F_del[k]
#           = sum_k C(k,j) [F_T[k] - N_uv[k] + N_uv[k+1]]
#           = c_j(T) - sum_k C(k,j) N_uv[k] + sum_k C(k,j) N_uv[k+1]
#           = c_j(T) - n_j + sum_k C(k,j) N_uv[k+1]
# where n_j = sum_k C(k,j) N_uv[k]
#
# sum_k C(k,j) N_uv[k+1] = sum_{k'=1}^{n-1} C(k'-1, j) N_uv[k']
#                         = sum_k C(k-1, j) N_uv[k]
#
# So c_j(T\e) = c_j(T) - sum_k [C(k,j) - C(k-1,j)] N_uv[k]
#             = c_j(T) - sum_k C(k-1, j-1) N_uv[k]   [by Pascal's identity]
#             = c_j(T) - n_{j-1}^{shifted}
# where n_{j-1}^{shifted} = sum_k C(k-1, j-1) N_uv[k] = sum_{k>=1} C(k-1,j-1) N_uv[k]
#                          = sum_{m>=0} C(m, j-1) N_uv[m+1]
#
# This is the (j-1)-th Taylor coefficient of N_uv SHIFTED!
# In other words, if we define N_uv^*(x) = sum_k N_uv[k+1] x^k (shift down by 1),
# then its j-1 Taylor coeff is sum_m C(m,j-1) N_uv[m+1].
# So: c_j(T\e) = c_j(T) - [Taylor_{j-1} of N_uv^*]

print("\nFormula: c_j(T\\e) = c_j(T) - sum_{k>=1} C(k-1, j-1) * N_uv[k]")
print("         = c_j(T) - [c_{j-1} of the shifted N_uv polynomial]")

# Verify this
for n in [5, 6]:
    m = n * (n - 1) // 2
    random.seed(42)

    for trial in range(5):
        bits = random.randint(0, (1 << m) - 1)
        adj = tournament_from_bits(n, bits)
        u, v = (0, 1) if adj[0][1] else (1, 0)

        F_T = compute_F_dp_digraph(adj, n)
        c_T = [sum(comb(k, j) * F_T[k] for k in range(n)) for j in range(n)]

        adj_del = [row[:] for row in adj]
        adj_del[u][v] = 0
        F_del = compute_F_dp_digraph(adj_del, n)
        c_del = [sum(comb(k, j) * F_del[k] for k in range(n)) for j in range(n)]

        N_uv = compute_N_uv(adj, n, u, v)

        ok = True
        for j in range(n):
            if j == 0:
                correction = N_uv[0]  # C(-1,-1) is ill-defined; but c_0(T\e) = c_0(T) - N_uv[0]... let's handle manually
                # Actually for j=0: c_0(T\e) = sum_k F_del[k] = sum_k [F_T[k] - N_uv[k] + N_uv[k+1]]
                # = c_0(T) - sum N_uv[k] + sum N_uv[k+1] = c_0(T) - N_uv[0]
                correction = N_uv[0]
            else:
                correction = sum(comb(k-1, j-1) * N_uv[k] for k in range(1, n))
            predicted = c_T[j] - correction
            if predicted != c_del[j]:
                ok = False
                print(f"  n={n}, trial {trial}, j={j}: FAIL")

        if trial < 2:
            print(f"  n={n}, trial {trial}: formula verified for all j")


# ============================================================
# What is N_uv mod 3? Is c_{j-1}(N_uv^*) = 0 mod 3?
# ============================================================
print("\n" + "=" * 70)
print("N_uv Taylor coefficients mod 3")
print("=" * 70)

# If c_j(T) = 0 mod 3 for j < val(n), and c_j(T\e) = 0 mod 3 for j < val(n)-1,
# then the correction sum_{k>=1} C(k-1,j-1) N_uv[k] must be 0 mod 3 for j < val(n)-1.
#
# But we KNOW c_j(T) = 0 mod 3 for j < val(n) (THM-086).
# So the claim c_j(T\e) = 0 for j < val(n)-1 is equivalent to:
#   sum_{k>=1} C(k-1, j-1) N_uv[k] = 0 mod 3 for j < val(n)-1
# i.e., c_{j-1}(N_uv^*) = 0 mod 3 for j-1 < val(n)-2
# i.e., the shifted N_uv polynomial has Taylor zeros for j < val(n)-2.

for n in [5, 6, 7]:
    m = n * (n - 1) // 2
    val_n = 2 * ((n - 1) // 2)
    random.seed(42)

    if n <= 6:
        num_samples = min(500, 1 << m)
        exhaustive = (num_samples == 1 << m)
    else:
        num_samples = 500
        exhaustive = False

    failures = [0] * n

    for trial in range(num_samples):
        if exhaustive:
            bits = trial
        else:
            bits = random.randint(0, (1 << m) - 1)

        adj = tournament_from_bits(n, bits)
        u, v = (0, 1) if adj[0][1] else (1, 0)
        N_uv = compute_N_uv(adj, n, u, v)

        # Taylor coefficients of shifted N_uv: c_j^* = sum_{k>=1} C(k-1, j) N_uv[k]
        for j in range(n):
            cj_star = sum(comb(k-1, j) * N_uv[k] for k in range(1, n))
            if cj_star % 3 != 0:
                failures[j] += 1

    method = "exhaustive" if exhaustive else f"{num_samples} samples"
    print(f"\n  n={n}, val(n)={val_n} ({method}):")
    print(f"  Need c_j^*(N_uv) = 0 mod 3 for j < val(n)-2 = {val_n - 2}")
    for j in range(min(n, val_n + 1)):
        status = "ALWAYS 0" if failures[j] == 0 else f"{failures[j]} FAIL"
        marker = ""
        if j == val_n - 2:
            marker = " <-- need up to here"
        print(f"    c_{j}^*(N_uv): {status}{marker}")


# ============================================================
# What IS N_uv? Can we express it in terms of simpler objects?
# ============================================================
print("\n" + "=" * 70)
print("Structure of N_uv")
print("=" * 70)

# N_uv[k] = #{perms with u imm. before v and fwd(P,T)=k}
# This is related to F(T, x) with the constraint that u is immediately before v.
# If we condition on u at position i and v at position i+1, the remaining
# n-2 elements form two groups: those before position i and those after position i+1.
# This looks like a "descent at position i" problem in the restricted permutation.
#
# Actually: N_uv is related to the transfer matrix!
# N_uv[k] = sum over positions i of M[*, u, i] * t(u,v) * M[v, *, n-1-i]
# where M counts path weights. But t(u,v) = 1 (arc exists in T).
# More precisely: N_uv = sum_i (paths ending at u at position i) * (paths starting at v at position i+1)
# This is exactly the "factored" path count.

# For now, let's check: is N_uv tournament-independent mod 3?
print("\nIs N_uv (mod 3) tournament-independent?")

for n in [4, 5]:
    m = n * (n - 1) // 2
    N_mod3_set = set()

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        u, v = (0, 1) if adj[0][1] else (1, 0)
        N_uv = compute_N_uv(adj, n, u, v)
        N_mod3 = tuple(x % 3 for x in N_uv)
        N_mod3_set.add(N_mod3)

    print(f"  n={n}: {len(N_mod3_set)} distinct N_uv mod 3 patterns")
    for pat in sorted(N_mod3_set):
        print(f"    {pat}")


print("\nDONE")
