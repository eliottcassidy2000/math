"""Audit A (THM-4762 (1)): independent check of the SL2(Z) / Pisano facts.  Auditor for opus-2026-10-22-S55.

Checks (exact integer arithmetic, written independently of catmap_mtc.py):
  * A = t s t^-1 s^-1 = Q^2 (sympy, exact);
  * pi(N) via the Fibonacci pair recursion; ord(A mod N) via matrix powers; ord = pi/2 for 3 <= N <= NMAX;
  * pi(N) even for N >= 3; no N <= NMAX with pi(N) = 15; least N with pi(N) = 120;
  * rank of apparition alpha(N), pi/alpha in {1,2,4}; |<A> cap Z(SL2(Z/N))| (central powers of A) for N <= NMAX;
  * the sharper fact |<A> cap Z| in {1,2} for N >= 3 (proof in the report: pi/alpha = 4 iff alpha odd).
Usage: python a_pisano.py [NMAX]"""
import sys
from math import gcd
import sympy as sp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000

s = sp.Matrix([[0, -1], [1, 0]]); t = sp.Matrix([[1, 1], [0, 1]]); Q = sp.Matrix([[1, 1], [1, 0]])
A = t * s * t.inv() * s.inv()
print('t s t^-1 s^-1 =', A.tolist(), '; Q^2 =', (Q ** 2).tolist(), '; equal:', A == Q ** 2, '; s t =', (s * t).tolist(),
      'trace', (s * t).trace(), '; s^4 = I:', s ** 4 == sp.eye(2), '; (s t)^6 = I:', (s * t) ** 6 == sp.eye(2))
print('eigenvalues of A:', [sp.nsimplify(e) for e in A.eigenvals()], '; phi^2 =', sp.nsimplify(((1 + sp.sqrt(5)) / 2) ** 2))


def pisano(N):
    if N == 1:
        return 1
    a, b = 0, 1
    for k in range(1, 6 * N + 1):          # pi(N) <= 6N (classical)
        a, b = b, (a + b) % N
        if (a, b) == (0, 1):
            return k
    raise RuntimeError(N)


def matmul(X, Y, N):
    return ((X[0] * Y[0] + X[1] * Y[2]) % N, (X[0] * Y[1] + X[1] * Y[3]) % N,
            (X[2] * Y[0] + X[3] * Y[2]) % N, (X[2] * Y[1] + X[3] * Y[3]) % N)


def order_and_central(M, N):
    """order of M in SL2(Z/N) and the list of exponents j (1 <= j <= order) with M^j central (scalar mod N)."""
    I = (1 % N, 0, 0, 1 % N)
    X = I; j = 0; cen = []
    while True:
        X = matmul(X, M, N); j += 1
        if X[1] == 0 and X[2] == 0 and X[0] == X[3]:
            cen.append((j, X[0]))
        if X == I:
            return j, cen


def alpha(N):
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % N; k += 1
        if a == 0:
            return k


Am = (2, 1, 1, 1)
bad_ord, bad_even, ratio_set, cent_sizes, fifteen = [], [], {}, {}, []
for N in range(1, NMAX + 1):
    p = pisano(N)
    if p == 15:
        fifteen.append(N)
    if N >= 2:
        o, cen = order_and_central(Am, N)
        if N >= 3 and o * 2 != p:
            bad_ord.append(N)
        if N >= 3 and p % 2:
            bad_even.append(N)
        if N >= 3:
            r = p // alpha(N)
            ratio_set[r] = ratio_set.get(r, 0) + 1
            cent_sizes[len(cen)] = cent_sizes.get(len(cen), 0) + 1
            assert p % alpha(N) == 0
    if N == 2:
        print('N = 2: pi(2) =', p, ', ord(A mod 2) =', order_and_central(Am, 2)[0])
print(f'N <= {NMAX}: ord(A mod N) = pi(N)/2 fails for N >= 3 at: {bad_ord}; pi(N) odd for N >= 3 at: {bad_even}; '
      f'N with pi(N) = 15: {fifteen}')
print('pi(1..30):', [pisano(N) for N in range(1, 31)])
print('pi(4), pi(6), pi(30) =', pisano(4), pisano(6), pisano(30), '; N <= NMAX with pi(N) = 120:',
      [N for N in range(1, NMAX + 1) if pisano(N) == 120])
print('pi/alpha value counts (N >= 3):', ratio_set)
print('number of central powers of A in <A> (incl. A^ord = I), N >= 3, value counts:', cent_sizes)
# odd alpha <=> pi/alpha = 4 check
viol = [N for N in range(3, NMAX + 1) if (pisano(N) // alpha(N) == 4) != (alpha(N) % 2 == 1)]
print('pi/alpha = 4 <=> alpha odd, violations for 3 <= N <= NMAX:', viol)
for N in (16, 80, 176, 5, 3, 8):
    o, cen = order_and_central(Am, N)
    print(f'  N = {N}: alpha = {alpha(N)}, pi = {pisano(N)}, ord A = {o}, central powers (j, lambda) = {cen}')
