"""
c4_induction_test.py
kind-pasteur-2026-03-07-S37

Test inductive approach: c_j(T) = c_j(T\e) + c_{j-1}(T/e)
(from THM-083 polynomial DC).

If c_{j-1}(T/e) = 0 mod 3 for j-1 < val(n-1), and c_j(T\e) = 0 mod 3
for j < val(n), then c_j(T) = 0 mod 3 for j < min(val(n-1)+1, val(n)).

For n=7 (odd): val(7)=6, val(6)=4.
DC gives: c_j(T) = c_j(T\e) + c_{j-1}(T/e).
From T/e (n-1=6 verts): c_{j-1}=0 mod 3 for j-1 < 4, i.e. j < 5.
So c_j(T) = c_j(T\e) mod 3 for j=4,5.
The question becomes: is c_j(T\e) = 0 mod 3 for j=4,5?

T\e is a digraph missing one arc. Does the THM-086 pattern extend to
"almost-tournaments" (complete minus one arc)?
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
    """F(D,x) for a general digraph (not necessarily tournament).
    asc counts positions where adj[P[i]][P[i+1]] = 1."""
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


def compute_contraction(adj, n, u, v):
    """Contract arc u->v: merge u,v into single vertex w.
    w inherits in-arcs from v, out-arcs from u."""
    # New tournament on n-1 vertices: relabel to [0, n-2]
    # w takes the position of min(u,v), others shifted
    verts = [i for i in range(n) if i != v]  # keep u, remove v
    new_n = n - 1
    new_adj = [[0]*new_n for _ in range(new_n)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            if i == u and j == u:
                continue
            # Determine arc direction
            src = i if i != u else u  # u is the merged vertex
            dst = j if j != u else u
            # For merged vertex w=u:
            # out-arcs of w: arcs from u to others (not v)
            # in-arcs of w: arcs from others to v
            if src == u and dst != u:
                # w -> dst: use v's out-arc (w inherits OUT from head v)
                new_adj[i_idx][j_idx] = adj[v][dst]
            elif src != u and dst == u:
                # src -> w: use u's in-arc (w inherits IN from tail u)
                new_adj[i_idx][j_idx] = adj[src][u]
            else:
                new_adj[i_idx][j_idx] = adj[src][dst]
    return new_adj, new_n


def delete_arc(adj, n, u, v):
    """Delete arc u->v from tournament. Returns digraph with missing arc."""
    new_adj = [row[:] for row in adj]
    new_adj[u][v] = 0  # arc u->v removed
    # Note: adj[v][u] = 0 too (was already 0 since u->v was in T)
    # So now neither u->v nor v->u exists
    return new_adj


# Test the DC identity: c_j(T) = c_j(T\e) + c_{j-1}(T/e)
print("=" * 70)
print("Verify DC for Taylor coefficients")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
random.seed(42)

for trial in range(20):
    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_from_bits(n, bits)
    F_T = compute_F_dp_digraph(adj, n)
    c_T = [sum(comb(k, j) * F_T[k] for k in range(n)) for j in range(n)]

    # Pick an arc u->v
    u, v = 0, 1
    if adj[0][1] == 0:
        u, v = 1, 0

    # T\e
    adj_del = delete_arc(adj, n, u, v)
    F_del = compute_F_dp_digraph(adj_del, n)
    c_del = [sum(comb(k, j) * F_del[k] for k in range(n)) for j in range(n)]

    # T/e
    adj_con, new_n = compute_contraction(adj, n, u, v)
    F_con = compute_F_dp_digraph(adj_con, new_n)
    c_con = [sum(comb(k, j) * F_con[k] for k in range(new_n)) for j in range(new_n)]

    # Check: c_j(T) = c_j(T\e) + c_{j-1}(T/e)
    ok = True
    for j in range(1, n):
        lhs = c_T[j]
        rhs = c_del[j] + c_con[j - 1]
        if lhs != rhs:
            ok = False
            print(f"  Trial {trial}, j={j}: c_T={lhs}, c_del+c_con={rhs} MISMATCH")

    if ok and trial < 3:
        print(f"  Trial {trial}: DC identity verified for all j=1..{n-1}")


# Now test: for T\e (tournament minus one arc), is c_j(T\e) = 0 mod 3 for j < val(n)?
print("\n" + "=" * 70)
print("c_j(T\\e) mod 3 for tournament minus one arc")
print("=" * 70)

for n in [6, 7]:
    m_bits = n * (n - 1) // 2
    val = 2 * ((n - 1) // 2)
    random.seed(42)
    num_samples = 3000 if n == 7 else 5000

    failures = [0] * n
    for _ in range(num_samples):
        bits = random.randint(0, (1 << m_bits) - 1)
        adj = tournament_from_bits(n, bits)
        # Pick arc 0->1 or 1->0 (whichever exists)
        u, v = (0, 1) if adj[0][1] else (1, 0)
        adj_del = delete_arc(adj, n, u, v)
        F_del = compute_F_dp_digraph(adj_del, n)
        c_del = [sum(comb(k, j) * F_del[k] for k in range(n)) for j in range(n)]

        for j in range(n):
            if c_del[j] % 3 != 0:
                failures[j] += 1

    print(f"\n  n={n}, val(n)={val}, {num_samples} samples:")
    for j in range(n):
        status = "ALWAYS 0" if failures[j] == 0 else f"{failures[j]} failures"
        print(f"    c_{j}(T\\e): {status}")


print("\nDONE")
