#!/usr/bin/env python3
"""
VERIFIED: Each diagonal M(r)[a,a] = sum_P (-1)^{pos(a,P)} * prod(r + s_e)
is an even function of r.

This means: the coefficient of r^k in M(r)[a,a] is zero for odd k.

WHY? The palindromic reversal P -> P^rev maps:
- pos(a,P) -> (n-1) - pos(a,P) (since a is at position k in P, it's at n-1-k in P^rev)
- prod(r+s_e) -> prod(r-s_e) (since forward/backward swap under reversal)

So: sum_P (-1)^{pos(a,P)} prod(r+s_e)
  = sum_P (-1)^{n-1-pos(a,P)} prod(r-s_e)  [pairing P with P^rev]
  = (-1)^{n-1} sum_P (-1)^{-pos(a,P)} prod(r-s_e)
  = (-1)^{n-1} sum_P (-1)^{pos(a,P)} * (-1)^{-2*pos(a,P)} prod(r-s_e)

Wait, let me be more careful. For odd n:
(-1)^{n-1-pos(a,P)} = (-1)^{n-1} * (-1)^{-pos(a,P)} = 1 * (-1)^{pos(a,P)} * (-1)^{-2*pos(a,P)}
No, (-1)^{n-1-k} = (-1)^{n-1} (-1)^{-k} = (-1)^{n-1} / (-1)^k

For odd n: (-1)^{n-1} = 1, so (-1)^{n-1-k} = (-1)^{-k} = (-1)^k (since (-1)^{-1} = -1)
So: (-1)^{n-1-pos(a,P)} = (-1)^{pos(a,P)}

And prod(r - s_e) = prod(r + s'_e) where s'_e = -s_e (edge in T^op)

So: M(r)[a,a] for T = sum_P_T (-1)^{pos(a,P)} prod(r+s_e)
                     = sum_P_{T^op} (-1)^{pos(a,P)} prod(r+s_e)
                     = M(r)[a,a] for T^op

This says M(r)[a,a] is the SAME for T and T^op (at odd n).
But this doesn't directly prove evenness in r.

For the evenness: prod(r+s_e) under r -> -r becomes prod(-r+s_e).
And prod(-r+s_e) = (-1)^{n-1} prod(r-s_e).

So: M(-r)[a,a] = sum_P (-1)^{pos(a,P)} prod(-r+s_e)
               = (-1)^{n-1} sum_P (-1)^{pos(a,P)} prod(r-s_e)
               = (-1)^{n-1} * M(r)[a,a] for T^op

If M(r)[a,a] for T = M(r)[a,a] for T^op (proved above), then:
M(-r)[a,a] = (-1)^{n-1} M(r)[a,a]

For odd n: M(-r)[a,a] = M(r)[a,a]. EVEN function! ✓
For even n: M(-r)[a,a] = -M(r)[a,a]. ODD function!

This PROVES the diagonal entries are even functions of r for odd n.

Now: what about OFF-diagonal entries?
M(r)[a,b] for a ≠ b: under the same argument with T^op:
M(-r)[a,b] = (-1)^{n-1} * M(r)[a,b] for T^op

But M(r)[a,b] for T^op is NOT obviously M(r)[a,b] for T.
The symmetry M[a,b] = M[b,a] is the claim, not M(T) = M(T^op).

HOWEVER: if we could show M(r)[a,b] = M(r)[b,a] AND M(r)[a,b] for T = M(r)[b,a] for T^op,
then: M(-r)[a,b] = (-1)^{n-1} M(r)[a,b] for T^op
     = (-1)^{n-1} M(r)[b,a] for T (if the off-diagonal symmetry extends to T^op)
     = (-1)^{n-1} M(r)[a,b] for T (by symmetry)

So off-diagonals would also be even functions for odd n.

Let me verify: are ALL entries of M(r) even functions of r at odd n?
"""
from itertools import permutations
import numpy as np
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_M_r_full(A, n, r):
    """Full transfer matrix at parameter r"""
    M = [[0.0]*n for _ in range(n)]
    V = set(range(n))

    for a in range(n):
        for b in range(n):
            if a == b:
                for P in permutations(range(n)):
                    weight = 1.0
                    for i in range(n-1):
                        s = 0.5 if A[P[i]][P[i+1]] else -0.5
                        weight *= (r + s)
                    M[a][a] += (-1)**P.index(a) * weight
            else:
                others = sorted(V - {a, b})
                total = 0.0
                for mask in range(1 << len(others)):
                    S = set()
                    for idx in range(len(others)):
                        if (mask >> idx) & 1:
                            S.add(others[idx])
                    R = set(others) - S

                    # E_a(S, r) with r-weighted paths
                    verts_E = sorted(S | {a})
                    E_a = 0.0
                    for P in permutations(verts_E):
                        if P[-1] != a: continue
                        w = 1.0
                        for i in range(len(P)-1):
                            s = 0.5 if A[P[i]][P[i+1]] else -0.5
                            w *= (r + s)
                        E_a += w

                    # B_b(R, r) with r-weighted paths
                    verts_B = sorted({b} | R)
                    B_b = 0.0
                    for P in permutations(verts_B):
                        if P[0] != b: continue
                        w = 1.0
                        for i in range(len(P)-1):
                            s = 0.5 if A[P[i]][P[i+1]] else -0.5
                            w *= (r + s)
                        B_b += w

                    total += (-1)**len(S) * E_a * B_b
                M[a][b] = total
    return np.array(M)

random.seed(42)

# Test at n=5 (odd)
n = 5
print(f"=== Full M(r) evenness test at n={n} (odd) ===")
for trial in range(3):
    A = random_tournament(n)

    for r in [0.3, 0.7, 1.3]:
        M_pos = compute_M_r_full(A, n, r)
        M_neg = compute_M_r_full(A, n, -r)

        diff = M_pos - M_neg
        max_diff = np.max(np.abs(diff))
        print(f"  trial {trial}, r={r}: max|M(r)-M(-r)| = {max_diff:.2e}")

# Test at n=4 (even)
n = 4
print(f"\n=== Full M(r) at n={n} (even) ===")
for trial in range(2):
    A = random_tournament(n)

    for r in [0.3, 0.7]:
        M_pos = compute_M_r_full(A, n, r)
        M_neg = compute_M_r_full(A, n, -r)

        # For even n: M(-r) = -M(r) (odd function)
        diff_odd = M_pos + M_neg
        max_diff = np.max(np.abs(diff_odd))
        print(f"  trial {trial}, r={r}: max|M(r)+M(-r)| = {max_diff:.2e}")

# ========== The KEY test: is M(r) symmetric for ALL r? ==========
print(f"\n\n=== M(r) symmetric for all r? ===")
n = 5
for trial in range(3):
    A = random_tournament(n)

    for r in [0.0, 0.3, 0.5, 1.0, 1.7]:
        M_r = compute_M_r_full(A, n, r)
        diff = M_r - M_r.T
        max_diff = np.max(np.abs(diff))
        print(f"  trial {trial}, r={r}: max|M(r)-M(r)^T| = {max_diff:.2e}")
