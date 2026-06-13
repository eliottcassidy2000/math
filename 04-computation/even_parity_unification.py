#!/usr/bin/env python3
"""
UNIFYING EVEN PARITY: M(r) = M(-r) and W(z)·W(-z) = 1

Two "even parity" symmetries in our framework:
1. THM-030 Cor 2: M(r) = M(-r) — transfer matrix even in r
2. IO reciprocity: W(z)·W(-z) = 1 — walk GF reciprocity

Are these related? Can one be derived from the other?

CONNECTION HYPOTHESIS:
W(z,r) is the walk GF with arc weights t(i,j) = r + s_{ij}.
At the multilinear level, extracting the Ham-path coefficient:
  H(r) = multilinear coeff of W(z,r,x_1,...,x_n)
  M(r)[a,b] = some extraction from W

If W(z,r)·W(-z,-r) = W(z,r) (commutative IO reciprocity),
and this holds at the multilinear level, it might imply M(r) = M(-r).

Let's test whether W(z,r)·W(-z,-r) = 1 at the (z,r)-parameter level.

Also: explore the relationship between H(r) and W at specific z values.
"""

from fractions import Fraction
from itertools import permutations, combinations
import numpy as np

def adjacency_weighted(A, r):
    """Weighted adjacency: t(i,j) = r + s_{ij} where s_{ij} = A[i][j] - 1/2."""
    n = len(A)
    T = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                T[i][j] = r + (A[i][j] - 0.5)
    return T

def W_gf(A, z, r):
    """IO walk GF: W(z,r) = det(I + z*T^T) / det(I - z*T) at parameters z, r."""
    n = len(A)
    T = adjacency_weighted(A, r)
    num = np.linalg.det(np.eye(n) + z * T.T)
    den = np.linalg.det(np.eye(n) - z * T)
    if abs(den) < 1e-15:
        return float('inf')
    return num / den

def H_weighted(A, r):
    """Total Ham path weight at parameter r."""
    n = len(A)
    total = 0.0
    for perm in permutations(range(n)):
        w = 1.0
        for k in range(n-1):
            w *= r + (A[perm[k]][perm[k+1]] - 0.5)
        total += w
    return total

# =====================================================================
# Test W(z,r) · W(-z,-r) = 1?
# =====================================================================
print("=" * 70)
print("W(z,r) · W(-z,-r) = 1? (Joint reciprocity)")
print("=" * 70)

n = 3
A3 = [[0,1,0],[0,0,1],[1,0,0]]  # 3-cycle
A3_trans = [[0,1,1],[0,0,1],[0,0,0]]  # transitive

for A, name in [(A3, "3-cycle"), (A3_trans, "transitive")]:
    print(f"\n  {name} (n=3):")
    for z, r in [(0.3, 0.4), (0.2, 0.3), (0.1, 0.6), (-0.2, 0.5)]:
        Wp = W_gf(A, z, r)
        Wm = W_gf(A, -z, -r)
        product = Wp * Wm
        print(f"    z={z:+.1f}, r={r:+.1f}: W(z,r)={Wp:.6f}, W(-z,-r)={Wm:.6f}, "
              f"product={product:.6f}")

# Also test W(z,r) · W(-z,r) = ? (negate z only, keep r)
print()
print("=" * 70)
print("W(z,r) · W(-z,r) = 1? (z-only reciprocity)")
print("=" * 70)

for A, name in [(A3, "3-cycle"), (A3_trans, "transitive")]:
    print(f"\n  {name} (n=3):")
    for z, r in [(0.3, 0.4), (0.2, 0.3), (0.1, 0.6)]:
        Wp = W_gf(A, z, r)
        Wm = W_gf(A, -z, r)  # negate z only
        product = Wp * Wm
        print(f"    z={z:+.1f}, r={r:+.1f}: product={product:.6f}")

# =====================================================================
# H(r) and W at specific z
# =====================================================================
print()
print("=" * 70)
print("H(r) vs W(z,r) AT SPECIFIC z VALUES")
print("=" * 70)

for A, name in [(A3, "3-cycle"), (A3_trans, "transitive")]:
    print(f"\n  {name}:")
    for r in [0.2, 0.3, 0.5, 0.7]:
        Hr = H_weighted(A, r)
        # Try various z values
        for z in [0.1, 0.3, 0.5, 0.8]:
            Wz = W_gf(A, z, r)
            if Wz != float('inf'):
                print(f"    r={r:.1f}, z={z:.1f}: H(r)={Hr:.4f}, W(z,r)={Wz:.4f}")

# =====================================================================
# n=5: Same tests
# =====================================================================
print()
print("=" * 70)
print("n=5: JOINT RECIPROCITY TESTS")
print("=" * 70)

n = 5
A5_paley = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j-i)%n in [1, 2]:
            A5_paley[i][j] = 1

A5_trans = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        A5_trans[i][j] = 1

for A, name in [(A5_paley, "Paley"), (A5_trans, "transitive")]:
    print(f"\n  {name} (n=5):")

    # W(z,r) · W(-z,-r) = ?
    for z, r in [(0.1, 0.2), (0.2, 0.3), (0.3, 0.1)]:
        Wp = W_gf(A, z, r)
        Wm = W_gf(A, -z, -r)
        print(f"    z={z}, r={r}: W*W(-z,-r) = {Wp*Wm:.6f}")

    # W(z,r) · W(-z,r) = ?
    for z, r in [(0.1, 0.2), (0.2, 0.3), (0.3, 0.1)]:
        Wp = W_gf(A, z, r)
        Wm = W_gf(A, -z, r)
        print(f"    z={z}, r={r}: W*W(-z,r) = {Wp*Wm:.6f}")

# =====================================================================
# H(r) = H(-r)? (even in r)
# =====================================================================
print()
print("=" * 70)
print("H(r) = H(-r)? (even in r)")
print("=" * 70)

for n, A, name in [(3, A3, "3-cycle"), (3, A3_trans, "trans n=3"),
                    (5, A5_paley, "Paley n=5"), (5, A5_trans, "trans n=5")]:
    for r in [0.1, 0.2, 0.3, 0.4]:
        Hr = H_weighted(A, r)
        Hnr = H_weighted(A, -r)
        print(f"  {name}: H({r})={Hr:.6f}, H({-r})={Hnr:.6f}, "
              f"H(r)=H(-r)? {abs(Hr-Hnr)<1e-10}")

# H(r) should contain only ODD powers of r (n-1 is even at odd n)
# Actually: H(r) = sum_P prod_{edges in P} (r + s_e)
# = sum_P (r + s_1)(r + s_2)...(r + s_{n-1})
# Each path has n-1 edges, so H(r) is degree n-1 polynomial in r.
# THM-030 says H(r) = tr(M(r)) at odd n. Since M(r) is even in r,
# tr(M(r)) is also even in r. So H(r) should be even at ODD n!

# At EVEN n: tr(M(r)) = 0 always, so this doesn't constrain H(r).
# H(r) at even n could be odd, even, or neither.

print()
print("=" * 70)
print("CONCLUSION: PARITY STRUCTURE")
print("=" * 70)
print("""
RESULTS:

1. W(z,r) · W(-z,-r) = 1 (JOINT reciprocity) — FAILS!
   Negating BOTH z and r does NOT give reciprocal in general.

2. W(z,r) · W(-z,r) = 1 (z-only reciprocity) — TRUE
   This follows from det(M^T) = det(M) regardless of r.

3. H(r) is EVEN in r at odd n:
   H(r) = H(-r) confirmed for 3-cycle, trans n=3, Paley n=5, trans n=5.

4. The IO reciprocity W(z)·W(-z)=1 and the M parity M(r)=M(-r)
   are INDEPENDENT symmetries:
   - W reciprocity: from det(T^T) = det(T) (matrix algebra, z-variable only)
   - M parity: from THM-030 inclusion-exclusion (r-variable)
   Both are "even parity" but in different variables.
""")
