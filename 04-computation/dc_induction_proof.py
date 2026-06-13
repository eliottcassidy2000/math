"""
dc_induction_proof.py
kind-pasteur-2026-03-07-S37

Verify the inductive proof of THM-086 via deletion-contraction.

PROOF SKETCH:
  c_j(T) = c_j(T\\e) + c_{j-1}(T/e)   [THM-083 in Taylor form]

  If c_{j-1}(T/e) = 0 mod 3 for all (n-1)-vertex tournaments when j-1 < val(n-1),
  and c_j(T\\e) = 0 mod 3 for all n-vertex "almost-tournaments" when j < val(n)-1,
  then c_j(T) = 0 mod 3 for j < val(n-1)+1 = val(n)-1.
  Then palindrome gives one more: c_{val(n)-1} = 0.

KEY QUESTIONS:
  1. Does c_j(T\\e) = 0 mod 3 for j < val(n)-1?
  2. Does the palindrome give the last step?
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


def delete_arc(adj, n, u, v):
    new_adj = [row[:] for row in adj]
    new_adj[u][v] = 0
    return new_adj


# Test: for T\\e at various n, find the threshold where c_j stops being universal 0 mod 3
print("=" * 70)
print("c_j(T\\e) universal zeros: precise threshold")
print("=" * 70)

for n in [5, 6, 7, 8]:
    val_n = 2 * ((n - 1) // 2)
    m_bits = n * (n - 1) // 2
    random.seed(42)
    num_samples = min(5000, 1 << m_bits)
    exhaustive = (num_samples == 1 << m_bits)

    failures = [0] * n

    if exhaustive:
        for bits in range(1 << m_bits):
            adj = tournament_from_bits(n, bits)
            u, v = (0, 1) if adj[0][1] else (1, 0)
            adj_del = delete_arc(adj, n, u, v)
            F = compute_F_dp_digraph(adj_del, n)
            c = [sum(comb(k, j) * F[k] for k in range(n)) for j in range(n)]
            for j in range(n):
                if c[j] % 3 != 0:
                    failures[j] += 1
    else:
        for _ in range(num_samples):
            bits = random.randint(0, (1 << m_bits) - 1)
            adj = tournament_from_bits(n, bits)
            u, v = (0, 1) if adj[0][1] else (1, 0)
            adj_del = delete_arc(adj, n, u, v)
            F = compute_F_dp_digraph(adj_del, n)
            c = [sum(comb(k, j) * F[k] for k in range(n)) for j in range(n)]
            for j in range(n):
                if c[j] % 3 != 0:
                    failures[j] += 1

    method = "exhaustive" if exhaustive else f"{num_samples} samples"
    print(f"\n  n={n}, val(n)={val_n} ({method}):")

    threshold = n
    for j in range(n):
        status = "ALWAYS 0" if failures[j] == 0 else f"{failures[j]} FAIL"
        marker = ""
        if failures[j] == 0 and (j + 1 < n and failures[j + 1] > 0):
            threshold = j + 1
            marker = " <-- threshold"
        elif failures[j] > 0 and (j == 0 or failures[j - 1] == 0):
            marker = " <-- first failure"
        print(f"    c_{j}(T\\e): {status}{marker}")

    val_del = threshold
    print(f"    => val(T\\e) = {val_del} vs val(T) = {val_n}, val(T)-1 = {val_n - 1}")
    if val_del >= val_n - 1:
        print(f"    DC induction step: c_j(T\\e)=0 for j < val(n)-1: HOLDS")
    else:
        print(f"    DC induction step: FAILS at this n")


# Summary of the induction
print("\n" + "=" * 70)
print("INDUCTION SUMMARY")
print("=" * 70)
print("""
For tournament T on n vertices, arc e = (u->v):
  c_j(T) = c_j(T\\e) + c_{j-1}(T/e)

Step 1: c_{j-1}(T/e) = 0 mod 3 for j-1 < val(n-1) [by induction on n-1 verts]
Step 2: c_j(T\\e) = 0 mod 3 for j < val(n)-1 [separate claim for almost-tournaments]
Step 3: Combined gives c_j(T) = 0 for j < min(val(n-1)+1, val(n)-1) = val(n)-1
Step 4: Palindrome: c_{val(n)-1} = 0 mod 3 [from palindrome + c_j=0 for j < val(n)-1]

For n odd: val(n) = n-1, val(n-1) = n-3.
  Step 1 gives j-1 < n-3, i.e. j < n-2.
  Step 2 gives j < n-2 (= val(n)-1).
  Step 3: j < n-2. Need j < n-1.
  Step 4: palindrome gives c_{n-2} = 0 (since n-1 is even, palindrome forces it).

For n even: val(n) = n-2, val(n-1) = n-2.
  Step 1 gives j-1 < n-2, i.e. j < n-1.
  Step 2 gives j < n-3 (= val(n)-1).
  Step 3: j < n-3. Need j < n-2.
  Step 4: palindrome gives c_{n-3} = 0 (one more from palindrome).
""")


print("DONE")
