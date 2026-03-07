#!/usr/bin/env python3
"""
PROOF: |det(M)| = F(n+1) for the transitive tournament T_n.

The transfer matrix for T_n is tridiagonal:
  M[i,i] = (-1)^i
  M[i,i+1] = M[i+1,i] = (-1)^i

Conjugating by D = diag((-1)^i):
  D*M*D = I + tridiag(S)

where S is the skew-adjacency matrix (S[i,j] = 1 if i<j, -1 if i>j).
But tridiag(S) = U - L (upper/lower bidiagonals of 1s and -1s).

So D*M*D = I + U - L, where:
  (I+U-L)[i,i] = 1
  (I+U-L)[i,i+1] = 1  (for i < n-1)
  (I+U-L)[i+1,i] = -1  (for i < n-1)

The determinant of this n x n matrix should be F(n+1).

PROOF by cofactor expansion along the last row:
Let A_n = I + U - L (n x n matrix).
det(A_n) = 1 * det(A_{n-1}) + 1 * det(B_{n-1})

where B_{n-1} is obtained by deleting row n and column n-1.

Actually, let's use the standard tridiagonal determinant recurrence.
For A_n with diagonal a, super-diagonal b, sub-diagonal c:

det(A_n) = a_n * det(A_{n-1}) - b_{n-1} * c_{n-1} * det(A_{n-2})

Here a = 1 (constant), b = 1 (constant), c = -1 (constant).
So: det(A_n) = 1 * det(A_{n-1}) - 1*(-1) * det(A_{n-2})
            = det(A_{n-1}) + det(A_{n-2})

This is the FIBONACCI recurrence!

With det(A_1) = 1 = F(2) and det(A_2) = 1*1 - 1*(-1) = 2 = F(3):
  det(A_3) = F(3) + F(2) = 3 = F(4)
  det(A_4) = F(4) + F(3) = 5 = F(5)
  det(A_5) = F(5) + F(4) = 8 = F(6)

So det(D*M*D) = det(I+U-L) = F(n+1).

Since det(D) = prod(-1)^i = (-1)^{0+1+...+(n-1)} = (-1)^{n(n-1)/2}:
det(M) = det(D)^{-2} * det(D*M*D) = det(D*M*D) = F(n+1)

Wait, det(D) = (-1)^{n(n-1)/2}, so det(D)^2 = (-1)^{n(n-1)}.
And D*M*D has det = det(D)^2 * det(M) = (-1)^{n(n-1)} * det(M).

So det(M) = (-1)^{n(n-1)} * F(n+1).

Check: n=3: (-1)^6 * 3 = 3, but det(M) = -3. Hmm.

Wait, D*M is NOT D*M*D. Let me recheck.

D * M has entries (D*M)[i,j] = (-1)^i * M[i,j].
And D*M != D*M*D unless D is self-inverse (which it is since D^2 = I).

But D*M is NOT similar to M. D*M*D IS similar to M:
det(D*M*D) = det(D)^2 * det(M).

Let me verify numerically.

kind-pasteur-2026-03-06-S25b (continuation)
"""

import numpy as np

def fibonacci(n):
    """Return F(n) where F(1)=F(2)=1."""
    a, b = 1, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return a


# ============================================================
# Verify: det(D*M*D) = F(n+1) for transitive tournament
# ============================================================
print("=" * 70)
print("det(M) for transitive tournament T_n")
print("=" * 70)

for n in range(2, 12):
    # Construct M directly
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    D = np.diag([(-1) ** i for i in range(n)])

    # D*M*D
    DMD = D @ M @ D

    # D*M (what we checked was I+U-L)
    DM = D @ M

    det_M = np.linalg.det(M)
    det_DM = np.linalg.det(DM)
    det_DMD = np.linalg.det(DMD)
    det_D = np.linalg.det(D)

    fn1 = fibonacci(n + 1)

    print(f"  n={n}: det(M) = {det_M:10.1f}  det(D*M) = {det_DM:10.1f}  det(D*M*D) = {det_DMD:10.1f}  F(n+1) = {fn1}")

    # Check what equals F(n+1)
    # det(D*M) = det(D)*det(M)
    # det(D*M*D) = det(D)^2 * det(M)


# ============================================================
# Analyze the relationship
# ============================================================
print("\n" + "=" * 70)
print("Relationships")
print("=" * 70)

for n in range(2, 10):
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    D = np.diag([(-1) ** i for i in range(n)])

    det_M = round(np.linalg.det(M))
    det_D = round(np.linalg.det(D))
    det_DM = round(np.linalg.det(D @ M))
    det_DMD = round(np.linalg.det(D @ M @ D))

    fn1 = fibonacci(n + 1)
    sign_n = (-1) ** (n * (n - 1) // 2)

    print(f"  n={n}: det(M)={det_M:6d}  det(D)={det_D:3d}  det(D*M)={det_DM:6d}  det(D*M*D)={det_DMD:6d}  F(n+1)={fn1}  (-1)^(n(n-1)/2)={sign_n}")

    # Check: det(D*M) should be det(D)*det(M)
    assert det_DM == det_D * det_M, f"Failed at n={n}"

    # det(D*M*D) should be det(D)^2 * det(M) = det(M) (since det(D)^2 = 1)
    assert det_DMD == det_M, f"det(DMD) != det(M) at n={n}"


# ============================================================
# So det(M) = det(D*M*D). Now what is D*M*D?
# ============================================================
print("\n" + "=" * 70)
print("D*M*D matrix structure")
print("=" * 70)

for n in [3, 5]:
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    D = np.diag([(-1) ** i for i in range(n)])
    DMD = D @ M @ D

    print(f"\n  n={n}: D*M*D =")
    for row in DMD:
        print(f"    {[int(x) for x in row]}")

    # D*M*D[i,j] = (-1)^i * M[i,j] * (-1)^j = (-1)^{i+j} * M[i,j]
    # M[i,i] = (-1)^i, so DMD[i,i] = (-1)^{2i} * (-1)^i = (-1)^i
    # M[i,i+1] = (-1)^i, so DMD[i,i+1] = (-1)^{2i+1} * (-1)^i = (-1)^{3i+1}

    # Hmm that doesn't look right. Let me compute directly.
    for i in range(n):
        for j in range(max(0,i-1), min(n,i+2)):
            val = int(DMD[i,j])
            expected = ((-1)**(i+j)) * int(M[i,j])
            print(f"    DMD[{i},{j}] = (-1)^{i+j} * M[{i},{j}] = (-1)^{i+j} * {int(M[i,j])} = {expected} (actual {val})")


# ============================================================
# D*M (not D*M*D) is the matrix I + U - L
# ============================================================
print("\n" + "=" * 70)
print("D*M = I + U - L (confirmed)")
print("=" * 70)

for n in range(2, 8):
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    D = np.diag([(-1) ** i for i in range(n)])
    DM = D @ M

    # Build I + U - L
    IUL = np.eye(n, dtype=int)
    for i in range(n - 1):
        IUL[i, i + 1] = 1
        IUL[i + 1, i] = -1

    match = np.array_equal(DM.astype(int), IUL)
    det_DM = round(np.linalg.det(DM))

    print(f"  n={n}: D*M = I+U-L: {match}, det(D*M) = {det_DM}")

    # det(D*M) = det(D)*det(M)
    det_D = round(np.linalg.det(D))
    det_M = round(np.linalg.det(M))
    print(f"         det(D) = {det_D}, det(M) = {det_M}, det(D)*det(M) = {det_D*det_M}")

    # det(I+U-L) satisfies Fibonacci recurrence
    # det(A_n) = 1*det(A_{n-1}) - 1*(-1)*det(A_{n-2}) = det(A_{n-1}) + det(A_{n-2})
    fn1 = fibonacci(n + 1)
    print(f"         det(I+U-L) = F({n+1}) = {fn1}: {abs(det_DM) == fn1 or det_DM == fn1}")


# ============================================================
# FORMULA: det(M) = det(D) * F(n+1) = (-1)^{n(n-1)/2} * F(n+1)
# ============================================================
print("\n" + "=" * 70)
print("FORMULA: det(M) = (-1)^{n(n-1)/2} * F(n+1)")
print("=" * 70)

for n in range(2, 12):
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    det_M = round(np.linalg.det(M))
    fn1 = fibonacci(n + 1)
    sign = (-1) ** (n * (n - 1) // 2)
    expected = sign * fn1

    # Actually det(D*M) = det(D)*det(M)
    # det(I+U-L) = F(n+1) (by Fibonacci recurrence)
    # So det(D)*det(M) = F(n+1)
    # det(M) = F(n+1) / det(D) = F(n+1) * det(D) (since det(D)^2 = 1)
    # det(D) = prod_{i=0}^{n-1} (-1)^i = (-1)^{sum_{i=0}^{n-1} i} = (-1)^{n(n-1)/2}

    det_D = (-1) ** (n * (n - 1) // 2)
    formula = det_D * fn1

    print(f"  n={n:2d}: det(M) = {det_M:8d}, (-1)^{{n(n-1)/2}} * F(n+1) = {formula:8d}, match: {det_M == formula}")


print("\n" + "=" * 70)
print("PROVED!")
print("=" * 70)
print("""
For the transitive tournament T_n:
  M is tridiagonal with M[i,i] = M[i,i+1] = (-1)^i
  D*M = I + U - L where D = diag((-1)^i)
  det(I+U-L) = F(n+1) by the Fibonacci recurrence
    (tridiagonal det recurrence: f_n = f_{n-1} + f_{n-2})
  det(M) = det(D)^{-1} * F(n+1) = (-1)^{n(n-1)/2} * F(n+1)

COROLLARY: For the transitive tournament T_n:
  Eigenvalues of M are 2*cos(k*pi/(n+1)) for k = 1,...,n
  (since I+U-L has eigenvalues 1 + 2i*sin(k*pi/(n+1)),
   and |det(D)|=1 doesn't change eigenvalue magnitudes,
   but the transformation M = D^{-1}*(I+U-L) is NOT a similarity transform)

Actually, M and I+U-L are NOT similar (D*M != D*M*D^{-1} since D^{-1}=D but D*M != M*D).
The eigenvalues of M are the eigenvalues of M (which we computed directly).
The eigenvalues of I+U-L = D*M are 1+2i*sin(kpi/(n+1)), which have modulus
sqrt(1 + 4*sin^2(kpi/(n+1))).

The eigenvalues of M are +-sqrt(1 + 4*sin^2(kpi/(n+1))) = +-sqrt(3 - 2cos(2kpi/(n+1))).
For n=5: k=1,2,3,4,5 give sqrt(3-2cos(2pi/6))=sqrt(2), sqrt(3-2cos(4pi/6))=2,
sqrt(3-2cos(pi))=sqrt(5)... that doesn't match.

Let me just verify the eigenvalues numerically.
""")

for n in [5, 7]:
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = (-1) ** i
        if i < n - 1:
            M[i, i + 1] = (-1) ** i
            M[i + 1, i] = (-1) ** i

    evals_M = sorted(np.linalg.eigvalsh(M))[::-1]
    print(f"  n={n}: M eigenvalues = {[round(e, 6) for e in evals_M]}")

    # Chebyshev: 2*cos(k*pi/(n+1)) for k=1,...,n
    cheb = sorted([2 * np.cos(k * np.pi / (n + 1)) for k in range(1, n + 1)])[::-1]
    print(f"         2*cos(k*pi/{n+1}) = {[round(c, 6) for c in cheb]}")
    print(f"         Match: {np.allclose(evals_M, cheb)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
