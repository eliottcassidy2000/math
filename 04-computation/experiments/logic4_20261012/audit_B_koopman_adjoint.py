"""Audit B (2026-10-12) of THM-4623(3): (a) the 'adjoint remark'.  The transpose of the Koopman operator
theta f = f o T restricted to level L (functions on Z/2^L, averaged back to level L) is the Markov matrix
P_L(a -> b) = (1/2) #{x in {a, a + 2^L} : T x = b mod 2^L}.  By the parity-vector bijection this is the
de Bruijn shift chain, so its spectrum is {1, 0, ..., 0} (P_L^L = J/2^L), NOT 'all nontrivial eigenvalues
of modulus 1/2' (that is THM-4520's 3-adic frequency operator on Z/(2 3^{n-1}), a different object).
(b) the cylinder claim used in the eigenvalue argument: every residue class mod 2^L contains a point p
with T^k p = p for every k >= L (periodically repeated parity words); checked exactly in Q for L = 4,
k in {4, 5, 7, 9, 11} and for L = 6, k in {6, 7, 11}.
"""
import numpy as np
from fractions import Fraction

def T(x): return x // 2 if x % 2 == 0 else (3 * x + 1) // 2   # 2-adic Collatz map on integers (x odd: (3x+1)/2)

print("(a) spectrum of the level-L transfer matrix P_L (transpose of theta restricted to level L):")
for L in range(1, 9):
    N = 2 ** L
    P = np.zeros((N, N))
    for a in range(N):
        for x in (a, a + N):
            P[a, T(x) % N] += 0.5
    ev = np.linalg.eigvals(P)
    mods = sorted(abs(ev), reverse=True)
    PL = np.linalg.matrix_power(P, L)
    print("  L=%d: eigenvalue moduli (top 3) %s; max|lambda| after the leading 1: %.2e; ||P^L - J/2^L||_max = %.1e"
          % (L, ["%.3f" % m for m in mods[:3]], mods[1], abs(PL - 1.0 / N).max()))

print("(b) periodic points in every residue class mod 2^L (exact rational 2-adic points):")
def periodic_point(word):
    """the unique x in Q (odd denominator) with parity vector word^infinity under T."""
    # T_w(x) = (3^s x + c)/2^k affine; fixed point x = c/(2^k - 3^s)
    a, c = Fraction(1), Fraction(0)          # x_j = a x + c after j steps
    for w in word:
        if w == 1: a, c = 3 * a / 2, (3 * c + 1) / 2
        else: a, c = a / 2, c / 2
    x = c / (1 - a)
    return x
def Tq(x):
    # parity of a rational with odd denominator = parity of numerator * inverse(denominator) mod 2
    num, den = x.numerator, x.denominator
    assert den % 2 == 1
    if num % 2 == 0: return x / 2
    return (3 * x + 1) / 2
import itertools
for L, ks in ((4, (4, 5, 7, 9, 11)), (6, (6, 7, 11))):
    N = 2 ** L
    for k in ks:
        classes = set()
        for word in itertools.product((0, 1), repeat=k):
            x = periodic_point(word)
            # check T^k x = x and read the residue of x mod 2^L
            y = x
            for _ in range(k): y = Tq(y)
            assert y == x
            res = (x.numerator * pow(x.denominator, -1, N)) % N
            # the residue must equal the class cut out by the first L parities (Terras bijection); record it
            classes.add(res)
        print("  L=%d, k=%d: residue classes mod 2^L hit by period-k points: %d of %d" % (L, k, len(classes), N))
