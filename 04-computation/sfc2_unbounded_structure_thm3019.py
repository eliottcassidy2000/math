"""SFC(2) unbounded: structure + a large modular census.

Facts established here:
 (S1) L(f^m) = int_0^inf f(t)^m e^{-t} dt.
 (S2) For REAL coefficients and m EVEN, L(f^m) > 0.  Since any two (or three)
      consecutive integers contain an even one, a common zero of a consecutive
      moment window is necessarily NON-REAL.  (Applies to SFC(N), any N, any k.)
 (S3) Consequently Res(I_{k+1}, I_{k+2}) >= 0 always: the even-index member has
      no real root, so the resultant taken from that side is a product over
      conjugate pairs of |.|^2.  Hence SFC(2) at window k  <=>  Res > 0.
 (S4) For (p,s) = (0,1) there is an exact recurrence  I_m = 1 + m lam I_{m-1},
      so I_m = I_{m+1} = 0 gives I_{m+1} = 1 + (m+1) lam I_m = 1, a
      contradiction: SFC(2) holds for that family at EVERY window, for free.
Census: Res(I_{k+1},I_{k+2}) mod several primes (nonzero mod one prime is a
one-way certificate of nonvanishing).
"""
from math import comb, factorial
import sympy as sp

lam = sp.Symbol('lam')
PRIMES = [(1 << 61) - 1, 2147483647, 1000000007]


def I_coeffs(m, p, s):
    return [comb(m, j) * factorial(p * m + s * j) for j in range(m + 1)]


def res_mod(c1, c2, mod):
    """resultant of two polynomials given by coefficient lists (ascending) mod prime"""
    A = [x % mod for x in c1][::-1]          # descending
    B = [x % mod for x in c2][::-1]
    while A and A[0] == 0: A = A[1:]
    while B and B[0] == 0: B = B[1:]
    if not A or not B: return 0
    res = 1
    while len(B) > 1:
        d = len(A) - len(B)
        if d < 0:
            A, B = B, A
            if (len(A) - 1) % 2 and (len(B) - 1) % 2:
                res = (-res) % mod
            continue
        inv = pow(B[0], mod - 2, mod)
        f = A[0] * inv % mod
        A = [(A[i] - f * (B[i - d] if 0 <= i - d < len(B) else 0)) % mod
             for i in range(len(A))]
        while A and A[0] == 0: A = A[1:]
        if not A: return 0
        res = res * pow(B[0], (len(A) - 1) - (len(A) - 1), mod) % mod
        # standard Euclidean resultant bookkeeping
        degA_old = None
    return None


# use sympy over ZZ modulo a prime instead (robust)
def resultant_nonzero(p, s, k, primes=PRIMES):
    c1 = I_coeffs(k + 1, p, s)
    c2 = I_coeffs(k + 2, p, s)
    for mod in primes:
        P1 = sp.Poly([x % mod for x in c1[::-1]], lam, modulus=mod)
        P2 = sp.Poly([x % mod for x in c2[::-1]], lam, modulus=mod)
        try:
            R = sp.resultant(P1, P2)
        except Exception:
            continue
        if R % mod != 0:
            return True, mod
    return False, None


print("(S4) exact recurrence check for (p,s)=(0,1): I_m = 1 + m*lam*I_{m-1}")
ok = True
for m in range(1, 9):
    lhs = sum(I_coeffs(m, 0, 1)[j] * lam ** j for j in range(m + 1))
    rhs = 1 + m * lam * sum(I_coeffs(m - 1, 0, 1)[j] * lam ** j for j in range(m))
    ok &= (sp.expand(lhs - rhs) == 0)
print(f"   verified m=1..8: {ok}")

print("\nSFC(2) census: Res(I_{k+1},I_{k+2}) != 0 (modular one-way certificates)")
fails = []
tested = 0
for p in range(0, 9):
    for s in range(1, 9):
        for k in range(0, 9):
            good, mod = resultant_nonzero(p, s, k)
            tested += 1
            if not good:
                fails.append((p, s, k))
    print(f"   p={p}: done (cumulative tested {tested}, failures {len(fails)})", flush=True)
print(f"\n   TOTAL tested {tested} cells (p<=8, s<=8, k<=8, i.e. exponents up to "
      f"q = p+s <= 16 and windows to 8)")
print(f"   cells with vanishing resultant: {fails if fails else 'NONE'}")
