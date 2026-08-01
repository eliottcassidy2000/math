"""Large census of the Hadamard multiplier question (= SFC(2) at p=0, all windows):
does A = sum_j C(m,j)(sj)! lam^j share a root with B = sum_j C(m,j)(sj+s)! lam^j ?
Equal-degree pair; modular resultant (nonzero mod one prime is a one-way certificate).
"""
from math import comb, factorial

PRIMES = [(1 << 61) - 1, 2305843009213693951 % ((1 << 61) - 1) or 2147483647, 1000000007]
PRIMES = [(1 << 61) - 1, 2147483647, 1000000007]


def res_mod(a, b, mod):
    """resultant of polys given ascending-coefficient lists, over F_mod (Euclidean)."""
    A = [x % mod for x in a][::-1]
    B = [x % mod for x in b][::-1]
    def strip(P):
        i = 0
        while i < len(P) and P[i] == 0:
            i += 1
        return P[i:]
    A, B = strip(A), strip(B)
    if not A or not B:
        return 0
    res = 1
    while True:
        dA, dB = len(A) - 1, len(B) - 1
        if dB == 0:
            res = res * pow(B[0], dA, mod) % mod
            return res % mod
        if dA < dB:
            A, B = B, A
            if (dA % 2) and (dB % 2):
                res = (-res) % mod
            continue
        inv = pow(B[0], mod - 2, mod)
        # A := A - (A[0]/B[0]) x^{dA-dB} B, repeatedly
        while len(A) - 1 >= dB and any(A):
            if A[0] == 0:
                A = A[1:]
                continue
            f = A[0] * inv % mod
            sh = (len(A) - 1) - dB
            for i in range(len(B)):
                A[i + 0] = (A[i] - f * B[i]) % mod if i + 0 < len(A) else 0
            A = A[1:]
        A = strip(A)
        if not A:
            return 0
        dAn = len(A) - 1
        res = res * pow(B[0], dA - dAn, mod) % mod
        if ((dA % 2) and (dB % 2)):
            pass
        A, B = B, A
        if (dB % 2) and (dAn % 2):
            res = (-res) % mod


import sympy as sp
lam = sp.Symbol('lam')

def nonzero_res(m, s):
    a = [comb(m, j) * factorial(s * j) for j in range(m + 1)]
    b = [comb(m, j) * factorial(s * j + s) for j in range(m + 1)]
    for mod in PRIMES:
        P1 = sp.Poly([x % mod for x in a[::-1]], lam, modulus=mod)
        P2 = sp.Poly([x % mod for x in b[::-1]], lam, modulus=mod)
        R = sp.resultant(P1, P2)
        if R % mod != 0:
            return True
    return False


print("Hadamard census: does A share a root with its w-Hadamard transform B?")
print("(equivalently SFC(2) at p=0, window k = m-1)\n")
fails, tested = [], 0
for s in range(1, 11):
    bad = []
    for m in range(1, 26):
        if not nonzero_res(m, s):
            bad.append(m)
        tested += 1
    print(f"   s={s:2d}: m=1..25 -> {'no shared root anywhere' if not bad else 'SHARED at m='+str(bad)}",
          flush=True)
    fails += [(s, m) for m in bad]
print(f"\n   TOTAL {tested} (s,m) cells, s<=10, m<=25 (windows k<=24)")
print(f"   cells with a shared root: {fails if fails else 'NONE'}")
