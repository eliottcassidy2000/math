#!/usr/bin/env python3
"""
VERIFY: W(T,r) = (H(T) / 2^{n-1}) * (2r+1)^{n-1}

This follows from the type B Worpitzky identity if our universal constants
are exactly the type B Eulerian numbers.

Type B Worpitzky identity:
  sum_{k=1}^{n} T_B(n,k) * C(x+k-1, n-1) = (2x+1)^{n-1}

Our finding:
  W(T,r) = (H/2^{n-1}) * sum_{k=0}^{n-1} c_k * C(r+k, n-1)

where c_k = T_B(n, k+1) (type B Eulerian numbers with shifted index).

If the identity holds, then:
  W(T,r) = (H/2^{n-1}) * (2r+1)^{n-1}

Let's verify this directly!
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def W_eval(A, n, r):
    """Evaluate W(T,r) directly"""
    total = 0.0
    for perm in permutations(range(n)):
        if not all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            continue
        prod = 1.0
        for i in range(n-1):
            s_e = A[perm[i]][perm[i+1]] - 0.5
            prod *= (r + s_e)
        total += prod
    return total

def W_eval_dp(A, n, r):
    """Evaluate W(T,r) using DP, faster"""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    s_e = 0.5  # forward edge: A[v][u] = 1
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)] * (r + s_e)
                elif A[u][v]:
                    s_e = -0.5  # backward edge: A[v][u] = 0
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)] * (r + s_e)
    # Wait, only edges from v to u where A[v][u]=1 are traversed
    # In a tournament, either A[v][u]=1 or A[u][v]=1
    # The HP goes v -> u, so the edge IS v->u if A[v][u]=1 (forward)
    # or the HP goes against the edge if A[v][u]=0 (backward)
    # Actually... in our HP enumeration, we require A[P_i][P_{i+1}]=1
    # So ALL edges in the path are forward. That means s_e = +0.5 always!
    # No wait, "forward" means P_i < P_{i+1} (in the labeling order)
    # not in the tournament direction.
    #
    # Let me re-derive: W(T,r) = sum_P prod_{i} (r + s_e)
    # where s_e = A[P_i][P_{i+1}] - 1/2
    # Since P is a HP of T, we have A[P_i][P_{i+1}] = 1 for all i
    # So s_e = 1 - 1/2 = 1/2 always!!
    # Therefore W(T,r) = sum_P (r + 1/2)^{n-1} = H * (r + 1/2)^{n-1}
    pass

# WAIT. I think I have a fundamental confusion about the definition.
# Let me re-read the definition.
#
# The "forward edge polynomial" F(T,x) counts HPs by number of
# "forward edges" where forward means P_i < P_{i+1} in the natural labeling.
#
# So F_k = #{HPs P where exactly k edges go from lower to higher label}
#
# In a HP P = (P_0, P_1, ..., P_{n-1}), the edge P_i -> P_{i+1}
# EXISTS in the tournament (A[P_i][P_{i+1}] = 1).
# It is a "forward edge" if P_i < P_{i+1} (ascending in the labeling).
#
# W(T,r) = sum_P prod_i (r + s_e(P_i, P_{i+1}))
# where s_e = +1/2 if P_i < P_{i+1} (forward), -1/2 if P_i > P_{i+1} (backward).
#
# NOT s_e = A[P_i][P_{i+1}] - 1/2 (which would always be +1/2).
#
# So s_e depends on the LABEL ORDER, not the arc direction.

print("=" * 70)
print("CORRECTED VERIFICATION: W(T,r) = (H/2^{n-1}) * (2r+1)^{n-1} ?")
print("=" * 70)

random.seed(42)

def W_eval_correct(A, n, r):
    """W(T,r) = sum_P prod_i (r + s_e) where s_e = +1/2 if P_i < P_{i+1}, else -1/2"""
    total = 0.0
    for perm in permutations(range(n)):
        if not all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            continue
        prod = 1.0
        for i in range(n-1):
            if perm[i] < perm[i+1]:
                s_e = 0.5
            else:
                s_e = -0.5
            prod *= (r + s_e)
        total += prod
    return total

for n in [3, 4, 5, 6, 7]:
    print(f"\nn={n}:")
    for trial in range(5 if n <= 5 else 3):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)

        for r in [0, 0.5, 1, 2, -1]:
            W_actual = W_eval_correct(A, n, r)
            W_predicted = (H / 2**(n-1)) * (2*r + 1)**(n-1)
            diff = abs(W_actual - W_predicted)
            match = diff < 1e-8

            if not match and r == 0.5:
                print(f"  Trial {trial}, r={r}: W={W_actual:.6f}, H/2^{n-1}*(2r+1)^{n-1}={W_predicted:.6f}, MISMATCH!")
            elif r == 0.5:
                pass  # Expected to match at r=1/2 since both give H

        # Check at non-trivial r
        r = 1.7
        W_actual = W_eval_correct(A, n, r)
        W_predicted = (H / 2**(n-1)) * (2*r + 1)**(n-1)
        diff = abs(W_actual - W_predicted)
        match = diff < 1e-6 * max(abs(W_actual), 1)

        if trial < 3:
            print(f"  Trial {trial}: H={H}, r=1.7: W={W_actual:.4f}, predicted={W_predicted:.4f}, match={match}")

# This probably DOESN'T work — W(T,r) depends on the tournament structure,
# not just H. Let me verify more carefully.

print("\n\n" + "=" * 70)
print("DETAILED CHECK: Does W(T,r) = (H/2^{n-1})*(2r+1)^{n-1}?")
print("=" * 70)

for n in [4, 5]:
    print(f"\nn={n}:")
    # Compare two tournaments with same H but different structure
    from itertools import permutations as perms

    tournaments = []
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A_t = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A_t[i][j] = 1
            else: A_t[j][i] = 1
        H = ham_path_count_dp(A_t, n)
        tournaments.append((H, A_t))

    # Group by H
    from collections import defaultdict
    by_H = defaultdict(list)
    for H, A_t in tournaments:
        by_H[H].append(A_t)

    for H in sorted(by_H.keys()):
        if len(by_H[H]) < 2:
            continue
        # Compare W(T,r) for two tournaments with same H
        A1 = by_H[H][0]
        A2 = by_H[H][1]

        r = 2.3
        W1 = W_eval_correct(A1, n, r)
        W2 = W_eval_correct(A2, n, r)
        pred = (H / 2**(n-1)) * (2*r + 1)**(n-1)

        if abs(W1 - W2) > 1e-8:
            print(f"  H={H}: W1(r=2.3) = {W1:.6f}, W2(r=2.3) = {W2:.6f}, predicted = {pred:.6f}")
            print(f"    -> W differs between tournaments with same H!")
            break
        else:
            if H <= 5:
                print(f"  H={H}: W1=W2={W1:.6f}, predicted={pred:.6f}, match={abs(W1-pred)<1e-6}")

# ===== So what's really going on? =====
# The Worpitzky expansion gave us W = (H/2^{n-1}) * sum c_k C(r+k, n-1)
# and the c_k are type B Eulerian numbers T_B(n-1, k+1).
# But sum T_B(n-1, k+1) C(r+k, n-1) = (2r+1)^{n-2} ... let me check the indexing.

print("\n\n" + "=" * 70)
print("TYPE B WORPITZKY IDENTITY CHECK")
print("=" * 70)

def typeB_eulerian(n, k):
    """T_B(n,k) for k=1,...,n. Using the formula from OEIS A060187."""
    # T(n,k) = sum_{i=1}^{k} (-1)^{k-i} C(n, k-i) (2i-1)^{n-1}
    # Note: k is 1-indexed here
    return sum((-1)**(k-i) * comb(n, k-i) * (2*i-1)**(n-1) for i in range(1, k+1))

for n in range(2, 8):
    row = [typeB_eulerian(n, k) for k in range(1, n+1)]
    print(f"  T_B({n},k) = {row}")

    # Check Worpitzky identity: sum_k T_B(n,k) C(x+k-1, n-1) = (2x+1)^{n-1}
    for x in [0, 1, 2, 3]:
        lhs = sum(typeB_eulerian(n, k) * comb(x + k - 1, n - 1) for k in range(1, n+1))
        rhs = (2*x + 1)**(n-1)
        print(f"    x={x}: LHS={lhs}, (2x+1)^{n-1}={rhs}, match={lhs==rhs}")

# So sum_k T_B(n,k) C(x+k-1, n-1) = (2x+1)^{n-1}
# In our notation (0-indexed): sum_{k=0}^{n-1} c_k C(r+k, n-1) = (2r+1)^{n-1}
# where c_k = T_B(n, k+1)
#
# But this is an IDENTITY in r — it holds for ALL r and ALL n.
# So our theorem W(T,r) = (H/2^{n-1}) * sum_k c_k C(r+k, n-1)
# automatically simplifies to W(T,r) = (H/2^{n-1}) * (2r+1)^{n-1}
#
# But this should be easy to verify directly!

print("\n\n" + "=" * 70)
print("FINAL DIRECT VERIFICATION")
print("=" * 70)

random.seed(42)

for n in [3, 4, 5, 6]:
    print(f"\nn={n}:")

    for trial in range(10 if n <= 5 else 3):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)

        # Check at several r values
        all_match = True
        for r in [0, 0.5, 1.0, -0.5, 2.0, 3.5, -1.7]:
            W_actual = W_eval_correct(A, n, r)
            W_predicted = (H / 2**(n-1)) * (2*r + 1)**(n-1)

            if abs(W_actual - W_predicted) > 1e-6 * max(abs(W_actual), 1):
                all_match = False
                if trial < 3:
                    print(f"  Trial {trial}, r={r}: W={W_actual:.6f}, predicted={W_predicted:.6f}, MISMATCH!")

        if all_match and trial < 3:
            print(f"  Trial {trial}: H={H}, W(T,r) = (H/2^{n-1})*(2r+1)^{n-1} for all tested r. EXACT!")
