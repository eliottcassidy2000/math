#!/usr/bin/env python3
"""
PROOF ATTEMPT: M(r)[a,b] = M(r)[b,a] via path reversal

M(r)[a,b] = sum_S (-1)^|S| * E_a(S,r) * B_b(R,r)

where S ⊆ V\{a,b}, R = V\{a,b}\S

E_a(S,r) = sum_{paths P through S∪{a} ending at a} prod(r + s_e(P))
B_b(R,r) = sum_{paths Q through R∪{b} starting at b} prod(r + s_e(Q))

Under reversal of BOTH P and Q:
- P^rev through S∪{a} starts at a (instead of ending there)
- Q^rev through R∪{b} ends at b (instead of starting there)

And the weights: prod(r+s_e) maps to prod(r-s_e) = prod(r+s_e') for T^op.

So: E_a(S,r) in T = B_a(S,r') in T^op where r' maps via some transformation.
Actually: sum_P_ending_at_a prod(r+s_e) in T
        = sum_P^rev_starting_at_a prod(r-s_e) in T
        = sum_Q_starting_at_a prod(r-s_e) in T

This is NOT B_a(S,r) (which uses r+s, not r-s). It's B_a(S,-r) in... no.

Wait: prod(r - s_e) = prod(r + (-s_e)).
For each edge, -s_e reverses the forward/backward assignment.
So prod(r - s_e) evaluated on path P in T = prod(r + s_e) evaluated on P in T^op.

Therefore:
E_a(S, r) in T = B_a(S, r) in T^op  [reversing path direction]

And: M(r)[a,b] in T = sum_S (-1)^|S| E_a(S,r) B_b(R,r)
                     = sum_S (-1)^|S| B_a(S,r) E_b(R,r) in T^op [by above]
                     = M(r)[b,a] in T^op [since B_b(R) role becomes E_b(R)]

Wait, that's exactly M(r)[b,a] in T^op!

So: M(r)[a,b] for T = M(r)[b,a] for T^op.

This is an ANTI-TRANSPOSITION relation: M(T) = M(T^op)^T.

But we also know (from the diagonal analysis) that M(T) and M(T^op) have
the same diagonal entries for odd n.

And if M(T) is symmetric (M = M^T), then M(T) = M(T^op)^T = M(T^op).

So symmetry of M(T) would IMPLY M(T) = M(T^op), which we should verify.

Let me check: does M(T) = M(T^op)?
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

def reverse_tournament(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]

def compute_M(A, n):
    M = [[0]*n for _ in range(n)]
    V = set(range(n))
    for a in range(n):
        for b in range(n):
            if a == b:
                total = 0
                for P in permutations(range(n)):
                    if all(A[P[i]][P[i+1]] for i in range(n-1)):
                        total += (-1)**P.index(a)
                M[a][a] = total
            else:
                others = sorted(V - {a, b})
                total = 0
                for mask in range(1 << len(others)):
                    S = set()
                    for idx in range(len(others)):
                        if (mask >> idx) & 1:
                            S.add(others[idx])
                    R = set(others) - S
                    verts_E = sorted(S | {a})
                    E_a = sum(1 for P in permutations(verts_E)
                              if P[-1] == a and all(A[P[i]][P[i+1]] for i in range(len(P)-1)))
                    verts_B = sorted({b} | R)
                    B_b = sum(1 for P in permutations(verts_B)
                              if P[0] == b and all(A[P[i]][P[i+1]] for i in range(len(P)-1)))
                    total += (-1)**len(S) * E_a * B_b
                M[a][b] = total
    return np.array(M, dtype=float)

random.seed(42)

print("=== M(T) vs M(T^op) ===")
for n in [3, 4, 5]:
    print(f"\nn={n}:")
    for trial in range(5 if n < 5 else 3):
        A = random_tournament(n)
        A_op = reverse_tournament(A)

        M_T = compute_M(A, n)
        M_Top = compute_M(A_op, n)

        # Check M(T) = M(T^op)
        diff1 = np.max(np.abs(M_T - M_Top))
        # Check M(T) = M(T^op)^T
        diff2 = np.max(np.abs(M_T - M_Top.T))

        H = int(np.trace(M_T)) if n % 2 == 1 else sum(M_T[i][j] for i in range(n) for j in range(n))

        print(f"  trial {trial}: |M(T)-M(T^op)|={diff1:.2e}, |M(T)-M(T^op)^T|={diff2:.2e}")

        if diff1 < 0.01:
            print(f"    M(T) = M(T^op) !")
        elif diff2 < 0.01:
            print(f"    M(T) = M(T^op)^T !")
        else:
            print(f"    Neither M=M^op nor M=M^op^T")
            if n <= 4:
                print(f"    M(T) = \n{M_T.astype(int)}")
                print(f"    M(T^op) = \n{M_Top.astype(int)}")

# ========== THE PROOF PATH ==========
print("\n\n=== PROOF PATH ===")
print("""
If M(T) = M(T^op)^T (anti-transposition), AND M(T) is symmetric (M = M^T),
then M(T) = M(T^op).

Conversely, if M(T) = M(T^op) and M(T) = M(T^op)^T,
then M(T) = M(T)^T (symmetric).

So: M(T) symmetric ⟺ M(T) = M(T^op)
given that M(T) = M(T^op)^T always holds.

If we can prove M(T) = M(T^op)^T algebraically (path reversal),
then proving M(T) = M(T^op) (complement invariance) would complete
the symmetry proof!

And we already know W(T) = W(T^op) for odd n (trace level).
The question is whether this extends to the FULL MATRIX level.
""")
