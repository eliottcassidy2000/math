"""
fk_mod2_proof.py
kind-pasteur-2026-03-07-S38

THEOREM (THM-094): F_k(T) = A(n,k) mod 2 for ALL tournaments T on n vertices.

Equivalently: F_k(T) = C(n-1, k) mod 2 for all T.
(Since A(n,k) = C(n-1,k) mod 2, a classical fact.)

PROOF by induction on n using deletion-contraction (THM-083):

Base case: n=1,2. Only one tournament (up to labeling direction), F = A_n(x).

Inductive step: Assume for all (n-1)-vertex tournaments T', F(T',x) = A_{n-1}(x) mod 2.
Let T_1, T_2 be n-vertex tournaments differing by one arc flip (e -> e').
Then T_1\\e = T_2\\e' (same digraph D), and T_1/e, T_2/e' are (n-1)-vertex tournaments.

DC identity (THM-083):
  F_{T_1}(x) = F_D(x) + (x-1)*F_{T_1/e}(x)
  F_{T_2}(x) = F_D(x) + (x-1)*F_{T_2/e'}(x)

Difference: F_{T_1}(x) - F_{T_2}(x) = (x-1)*[F_{T_1/e}(x) - F_{T_2/e'}(x)]

By inductive hypothesis: F_{T_1/e}(x) = A_{n-1}(x) = F_{T_2/e'}(x) mod 2.
Therefore: F_{T_1}(x) = F_{T_2}(x) mod 2.

Since any tournament connects to the transitive via arc flips:
  F(T,x) = F(T_trans, x) = A_n(x) mod 2. QED.

This script verifies the theorem exhaustively for small n and checks the
proof structure computationally.
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


def compute_F_dp(adj, n):
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


def eulerian_numbers(n):
    if n == 1:
        return [1]
    prev = eulerian_numbers(n - 1)
    A = [0] * n
    for k in range(n):
        A[k] = (k + 1) * prev[k] if k < len(prev) else 0
        if k > 0:
            A[k] += (n - k) * prev[k - 1]
    return A


# ============================================================
# Verification: F_k(T) = A(n,k) mod 2 for all T
# ============================================================
print("=" * 70)
print("THM-094: F_k(T) = A(n,k) mod 2 for all tournaments T")
print("=" * 70)

all_pass = True

for n in range(2, 8):
    A = eulerian_numbers(n)
    A_mod2 = [a % 2 for a in A]
    C_mod2 = [comb(n-1, k) % 2 for k in range(n)]

    # Check that A(n,k) = C(n-1,k) mod 2
    assert A_mod2 == C_mod2, f"Eulerian != binomial mod 2 at n={n}"

    m_bits = n * (n - 1) // 2
    num_T = 1 << m_bits

    random.seed(42)
    if n <= 7:
        samples = range(num_T) if n <= 6 else [random.randint(0, num_T-1) for _ in range(5000)]
    else:
        samples = [random.randint(0, num_T-1) for _ in range(3000)]

    failures = 0
    for bits in samples:
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        F_mod2 = [f % 2 for f in F]
        if F_mod2 != A_mod2:
            failures += 1
            if failures <= 3:
                print(f"  FAIL at n={n}: F mod 2 = {F_mod2}, expected {A_mod2}")

    method = "exhaustive" if n <= 6 else f"{len(list(samples))} samples"
    status = "PASS" if failures == 0 else f"FAIL ({failures})"
    print(f"  n={n} ({method}): {status}")
    print(f"    A(n,k) mod 2 = {A_mod2}")
    print(f"    C(n-1,k) mod 2 = {C_mod2}")
    if failures > 0:
        all_pass = False

if all_pass:
    print("\n  ALL TESTS PASSED")


# ============================================================
# Verify the DC induction step explicitly
# ============================================================
print("\n" + "=" * 70)
print("DC induction step verification")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    random.seed(42)
    A_n = eulerian_numbers(n)
    A_nm1 = eulerian_numbers(n - 1)

    verified = 0
    for trial in range(min(100, 1 << m)):
        bits = random.randint(0, (1 << m) - 1) if trial > 0 else 0
        adj = tournament_from_bits(n, bits)

        # Find an arc u->v
        u, v = (0, 1) if adj[0][1] else (1, 0)

        # F(T)
        F_T = compute_F_dp(adj, n)

        # T\e (delete arc u->v)
        adj_del = [row[:] for row in adj]
        adj_del[u][v] = 0
        F_del = compute_F_dp(adj_del, n)

        # T/e (contract arc u->v)
        # Merged vertex w keeps position of u, removes v
        verts = [i for i in range(n) if i != v]
        new_n = n - 1
        new_adj = [[0]*new_n for _ in range(new_n)]
        for i_idx, i in enumerate(verts):
            for j_idx, j in enumerate(verts):
                if i_idx == j_idx:
                    continue
                if i == u and j != u:
                    # w -> j: inherit OUT from head v
                    new_adj[i_idx][j_idx] = adj[v][j]
                elif i != u and j == u:
                    # i -> w: inherit IN from tail u
                    new_adj[i_idx][j_idx] = adj[i][u]
                else:
                    new_adj[i_idx][j_idx] = adj[i][j]
        F_con = compute_F_dp(new_adj, new_n)

        # Check DC at polynomial level: F_T(x) = F_del(x) + (x-1)*F_con(x)
        # Coefficient of x^k in (x-1)*F_con(x): F_con[k-1] - F_con[k]
        # (where F_con[-1] = 0 and F_con[new_n] = 0)
        for k in range(n):
            lhs = F_T[k]
            con_k = F_con[k] if k < new_n else 0
            con_km1 = F_con[k-1] if 0 < k <= new_n - 1 else (F_con[k-1] if k-1 >= 0 and k-1 < new_n else 0)
            rhs = F_del[k] + (con_km1 if k >= 1 else 0) - con_k
            if lhs != rhs and trial == 0:
                # Try the other convention
                pass
        # Actually let's just verify at the polynomial level
        # F_T(x) = F_del(x) + (x-1)*F_con(x)
        # Evaluate at a few x values
        for x in [0, 1, 2, 3, -1]:
            ft = sum(F_T[k] * x**k for k in range(n))
            fd = sum(F_del[k] * x**k for k in range(n))
            fc = sum(F_con[k] * x**k for k in range(new_n))
            assert ft == fd + (x-1)*fc, f"DC poly fails at n={n}, trial {trial}, x={x}"

        # Check that F_con mod 2 = A_{n-1} mod 2
        F_con_mod2 = [f % 2 for f in F_con]
        A_nm1_mod2 = [a % 2 for a in A_nm1]
        assert F_con_mod2 == A_nm1_mod2, f"Inductive hypothesis fails at n={n}"

        verified += 1

    print(f"  n={n}: DC identity + inductive hypothesis verified for {verified} tournaments")


# ============================================================
# Consequence: F_k(T) is odd iff C(n-1,k) is odd (Lucas' theorem)
# ============================================================
print("\n" + "=" * 70)
print("Consequence: F_k(T) parity from Lucas' theorem")
print("=" * 70)

for n in range(3, 10):
    odd_positions = [k for k in range(n) if comb(n-1, k) % 2 == 1]
    even_positions = [k for k in range(n) if comb(n-1, k) % 2 == 0]
    print(f"  n={n}: F_k odd at k={odd_positions}, F_k even at k={even_positions}")


# ============================================================
# Connection to mod 3: does the DC proof also give mod 3?
# ============================================================
print("\n" + "=" * 70)
print("Why DC proof works for mod 2 but not mod 3")
print("=" * 70)

# For mod 2: F(T/e) = A_{n-1}(x) mod 2 for ALL T/e (by induction).
# So the correction (x-1)*[F(T_1/e) - F(T_2/e')] = 0 mod 2.

# For mod 3: F(T/e) is NOT constant mod 3. Different (n-1)-vertex tournaments
# have different F mod 3 (e.g., one free parameter for odd n-1).
# So (x-1)*[F(T_1/e) - F(T_2/e')] != 0 mod 3 in general.
# We can only get TAYLOR COEFFICIENT vanishing, not full F equality.

print("  Mod 2: F(T/e) = A_{n-1}(x) mod 2 for ALL T/e => full F equality mod 2")
print("  Mod 3: F(T/e) mod 3 varies (1 free param for odd n) => only Taylor zeros")

# Verify: how many distinct F mod 3 patterns at each n?
for n in range(3, 7):
    m = n * (n - 1) // 2
    patterns = set()
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        patterns.add(tuple(f % 3 for f in F))
    print(f"  n={n}: {len(patterns)} distinct F mod 3 (vs 1 mod 2)")


print("\nDONE")
