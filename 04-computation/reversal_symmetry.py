#!/usr/bin/env python3
"""
HYPOTHESIS T: Transfer matrix symmetry from path reversal

M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(R)
where S ⊆ V\{a,b}, R = V\{a,b}\S

E_a(S) = #{HPs in T[S∪{a}] ending at a}
B_b(R) = #{HPs in T[R∪{b}] starting from b}

Under reversal of paths in T^op:
E_a^{op}(S) = #{HPs in T^op[S∪{a}] ending at a}
            = #{HPs in T[S∪{a}] starting from a} = B_a(S)

So: M^{op}[a,b] = sum_S (-1)^|S| B_a(S) * E_b(R)

If T has a "reversal symmetry" (T ≅ T^op via some permutation),
then M[a,b] = M^op[perm(a),perm(b)] = M[perm(a),perm(b)] with swapped roles.

But we want M[a,b] = M[b,a] for ALL tournaments, not just self-converse ones.

The palindromic structure gives: for each HP P = (v_0,...,v_{n-1}) in T,
P^rev = (v_{n-1},...,v_0) is a HP in T^op.

So H(T) = H(T^op) always.

Now: M[b,a] = sum_S (-1)^|S| E_b(S) * B_a(R)

For symmetry M[a,b] = M[b,a], we need:
sum_S (-1)^|S| [E_a(S)*B_b(R) - E_b(S)*B_a(R)] = 0

This is a SUM over 2^{n-2} subsets of a SIGNED PRODUCT.

INSIGHT: Under path reversal, E_a(S) in T maps to B_a(S) in T^op.
If we could map between T and T^op while preserving the inclusion-exclusion
structure, we'd get the symmetry.

HYPOTHESIS: The "r-polynomial" version M(r) might make the symmetry transparent.

Let M(r)[a,b] = sum_P r^{pos(a,P)} * (-1)^{pos(a,P)} * (stuff)

Actually, the W-polynomial version:

W(r) = tr(M(r)) where M(r) has entries
M(r)[a,b] = sum_{P through a...b} prod_{edges} (r + s_e)

Wait, I need to think about what the r-parameterized version actually is.

Let me just compute M[a,b] directly and verify symmetry, then look for
patterns in the actual matrix structure.
"""
from itertools import permutations, combinations
import numpy as np
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_M(A, n):
    """Compute transfer matrix M[a][b] using definition"""
    M = [[0]*n for _ in range(n)]
    V = set(range(n))

    for a in range(n):
        for b in range(n):
            if a == b:
                # Diagonal: M[a,a] = sum_P (-1)^{pos(a,P)}
                total = 0
                for P in permutations(range(n)):
                    if all(A[P[i]][P[i+1]] for i in range(n-1)):
                        pos_a = P.index(a)
                        total += (-1)**pos_a
                M[a][a] = total
            else:
                # Off-diagonal: M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R)
                others = sorted(V - {a, b})
                total = 0
                for mask in range(1 << len(others)):
                    S = set()
                    for idx in range(len(others)):
                        if (mask >> idx) & 1:
                            S.add(others[idx])
                    R = set(others) - S

                    # E_a(S) = #{HPs in T[S∪{a}] ending at a}
                    verts_E = sorted(S | {a})
                    E_a = 0
                    for P in permutations(verts_E):
                        if P[-1] != a: continue
                        valid = True
                        for i in range(len(P)-1):
                            if not A[P[i]][P[i+1]]:
                                valid = False
                                break
                        if valid: E_a += 1

                    # B_b(R) = #{HPs in T[{b}∪R] starting from b}
                    verts_B = sorted({b} | R)
                    B_b = 0
                    for P in permutations(verts_B):
                        if P[0] != b: continue
                        valid = True
                        for i in range(len(P)-1):
                            if not A[P[i]][P[i+1]]:
                                valid = False
                                break
                        if valid: B_b += 1

                    total += (-1)**len(S) * E_a * B_b
                M[a][b] = total
    return M

random.seed(42)

print("=== Transfer matrix M[a,b] verification ===")
for n in [3, 4, 5]:
    print(f"\nn={n}:")
    for trial in range(3 if n < 5 else 2):
        A = random_tournament(n)
        M = compute_M(A, n)

        # Check symmetry
        sym = True
        for a in range(n):
            for b in range(a+1, n):
                if M[a][b] != M[b][a]:
                    sym = False

        H = sum(1 for P in permutations(range(n))
                if all(A[P[i]][P[i+1]] for i in range(n-1)))
        tr = sum(M[a][a] for a in range(n))

        print(f"  trial {trial}: H={H}, tr(M)={tr}, symmetric={'YES' if sym else 'NO'}")
        if n <= 4:
            print(f"  M = ")
            for row in M:
                print(f"    {row}")

# ========== Look for structure in M ==========
print("\n\n=== Structure of M at n=5 ===")
n = 5
for trial in range(5):
    A = random_tournament(n)
    M = compute_M(A, n)
    H = sum(1 for P in permutations(range(n))
            if all(A[P[i]][P[i+1]] for i in range(n-1)))

    M_np = np.array(M, dtype=float)
    eigs = np.round(np.sort(np.linalg.eigvals(M_np))[::-1], 4)
    scores = [sum(A[i]) for i in range(n)]

    print(f"\n  trial {trial}: H={H}, scores={scores}")
    print(f"  M = {M_np.astype(int).tolist()}")
    print(f"  eigenvalues of M: {eigs}")
    print(f"  tr(M) = {int(np.trace(M_np))}")
    print(f"  M symmetric: {np.allclose(M_np, M_np.T)}")

    # Check: is M related to A in a simple way?
    # M[a,b] as function of scores, adjacency...
    for a in range(n):
        for b in range(a+1, n):
            edge = "→" if A[a][b] else "←"
            print(f"    ({a}{edge}{b}): M[{a},{b}]={M[a][b]:>3}, M[{b},{a}]={M[b][a]:>3}, "
                  f"diff={M[a][b]-M[b][a]}")
