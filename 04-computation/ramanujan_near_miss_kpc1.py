#!/usr/bin/env python3
# ramanujan_near_miss_kpc1.py -- session kind-pasteur-2026-06-10-S1, Thread D (HYP-2370)
# Ramanujan's near-miss family (lost notebook; cf. Hirschhorn, Math. Mag. 68 (1995);
# Ono & Trebat-Leder, "The 1729 K3 surface", Res. Number Theory (2015/16), arXiv:1510.00735).
#
# Generating functions (Ramanujan):
#   sum a_n x^n = (1 + 53x +  9x^2) / D(x)
#   sum b_n x^n = (2 - 26x - 12x^2) / D(x)
#   sum c_n x^n = (2 +  8x - 10x^2) / D(x),   D(x) = 1 - 82x - 82x^2 + x^3.
# CLAIM (Ramanujan):  a_n^3 + b_n^3 = c_n^3 + (-1)^n   for all n >= 0.
#
# This script:
#   (1) expands the series exactly (pure-int convolution recurrence
#       u_n = N_n + 82 u_{n-1} + 82 u_{n-2} - u_{n-3}), prints triples n = 0..10,
#       and checks the cubic identity with bigints;
#   (2) extends the HOMOGENEOUS recurrence backwards (it is reversible: leading and
#       trailing coefficients of D are units) and shows n = -1 gives (-9, 12, 10),
#       i.e. 9^3 + 10^3 = 12^3 + 1 = 1729: the taxicab number IS in the family;
#   (3) PROVES the identity for ALL n (both directions) by pure integer linear
#       algebra: a_n, b_n, c_n are linear functionals of C^n for the 3x3 companion
#       matrix C (det C = -1); hence a_n^3, b_n^3, c_n^3 are linear functionals of
#       M^n, M = C(x)C(x)C (Kronecker cube, 27x27, integer, det = -1); (-1)^n = r^n
#       with r = -1 a root of chi_M (product of the three eigenvalues of C); so
#       d_n = a_n^3 + b_n^3 - c_n^3 - (-1)^n satisfies the order-27 integer linear
#       recurrence with characteristic polynomial chi_M (Cayley-Hamilton), monic
#       with constant term != 0 (reversible).  27 consecutive zeros then force
#       d == 0 for every integer n.  We verify d_n = 0 for n = -10..40 (51 values,
#       covering 27 consecutive) and chi_M(-1) = 0 and chi_M(0) != 0, all exactly.
# Pure python ints + Fraction (for the char poly; result is integer and checked so).

import sys
from fractions import Fraction

failures = 0
def check(label, ok):
    global failures
    if not ok:
        failures += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

# ---------------------------------------------------------------- series expansion
D = [1, -82, -82, 1]                       # D(x) = 1 - 82x - 82x^2 + x^3
NUM_A = [1, 53, 9]
NUM_B = [2, -26, -12]
NUM_C = [2, 8, -10]

def series(num, n_max):
    """Coefficients of num(x)/D(x) for n = 0..n_max, exact ints."""
    u = []
    for n in range(n_max + 1):
        Nn = num[n] if n < len(num) else 0
        un = Nn
        if n >= 1: un += 82 * u[n - 1]
        if n >= 2: un += 82 * u[n - 2]
        if n >= 3: un -= u[n - 3]
        u.append(un)
    return u

N_FWD = 40
a = series(NUM_A, N_FWD)
b = series(NUM_B, N_FWD)
c = series(NUM_C, N_FWD)

# Backwards extension of the homogeneous recurrence u_{n+3} = 82u_{n+2}+82u_{n+1}-u_n,
# solved for u_n:  u_n = 82 u_{n+2} + 82 u_{n+1} - u_{n+3}.
N_BACK = 10
def extend_back(u, k):
    v = list(u)            # v[i] holds u_{i - k}
    for _ in range(k):
        v.insert(0, 82 * v[1] + 82 * v[0] - v[2])
    return v

A = extend_back(a, N_BACK)   # A[i] = a_{i - N_BACK}
B = extend_back(b, N_BACK)
C_ = extend_back(c, N_BACK)
def at(u, n): return u[n + N_BACK]

print("=" * 78)
print("PART 4a -- THE FAMILY, FORWARD (n = 0..10), exact bigints")
print("=" * 78)
print(f"{'n':>3} {'a_n':>22} {'b_n':>22} {'c_n':>22}  a^3+b^3-c^3")
all_ok = True
for n in range(0, 11):
    an, bn, cn = at(A, n), at(B, n), at(C_, n)
    d = an**3 + bn**3 - cn**3
    ok = (d == (-1) ** n)
    all_ok &= ok
    print(f"{n:>3} {an:>22} {bn:>22} {cn:>22}  {d:+d}  {'PASS' if ok else 'FAIL'}")
check("a_n^3 + b_n^3 = c_n^3 + (-1)^n for n = 0..10", all_ok)
print()
print("Sanity vs. literature: (a_0,b_0,c_0) = (1,2,2): 1+8 = 9 = 8+1;")
print("(a_1,b_1,c_1) = (135,138,172): 135^3+138^3 = 172^3 - 1 = 5088447;")
print("(a_2,b_2,c_2) = (11161,11468,14258).  All match Hirschhorn (1995).")

print()
print("=" * 78)
print("PART 4b -- THE FAMILY, BACKWARD: 1729 LIVES AT n = -1")
print("=" * 78)
print("The homogeneous recurrence u_{n+3} = 82u_{n+2} + 82u_{n+1} - u_n has unit")
print("leading AND trailing coefficients, so it runs backwards over Z:")
print("u_n = 82u_{n+2} + 82u_{n+1} - u_{n+3}.  (For n < 0 these are NOT power-series")
print("coefficients -- the series ones vanish -- but the canonical bi-infinite")
print("extension of the homogeneous recurrence, as in Ono--Trebat-Leder.)")
back_ok = True
for n in range(-3, 0):
    an, bn, cn = at(A, n), at(B, n), at(C_, n)
    d = an**3 + bn**3 - cn**3
    back_ok &= (d == (1 if n % 2 == 0 else -1))
    print(f"  n = {n}: (a,b,c) = ({an}, {bn}, {cn}),  a^3+b^3-c^3 = {d:+d}")
check("identity a_n^3+b_n^3-c_n^3 = (-1)^n holds at n = -3,-2,-1", back_ok)
am1, bm1, cm1 = at(A, -1), at(B, -1), at(C_, -1)
check("(a_{-1}, b_{-1}, c_{-1}) = (-9, 12, 10)", (am1, bm1, cm1) == (-9, 12, 10))
check("(-9)^3 + 12^3 = 10^3 - 1  (n = -1, (-1)^n = -1)", am1**3 + bm1**3 == cm1**3 - 1)
check("rearranged: 9^3 + 10^3 = 12^3 + 1", 9**3 + 10**3 == 12**3 + 1)
check("9^3 + 10^3 = 1729 = 12^3 + 1^3 (taxicab)", 9**3 + 10**3 == 1729 and 12**3 + 1 == 1729)
print("So (9,10,12) belongs to the SAME Ramanujan family, at index n = -1 of the")
print("backward extension -- 1729 is not a companion parametrization but the")
print("(-1)-st term of the lost-notebook family.  (Ono & Trebat-Leder 2015 make the")
print("same observation en route to the K3 surface; see citation block below.)")

print()
print("=" * 78)
print("PART 4c -- PROOF FOR ALL n (Cayley-Hamilton on the Kronecker cube)")
print("=" * 78)

# Companion matrix of u_{n+3} = 82u_{n+2} + 82u_{n+1} - u_n acting on (u_n,u_{n+1},u_{n+2}):
Cmat = [[0, 1, 0],
        [0, 0, 1],
        [-1, 82, 82]]

def mat_mul(X, Y):
    n, m, p = len(X), len(Y[0]), len(Y)
    return [[sum(X[i][k] * Y[k][j] for k in range(p)) for j in range(m)] for i in range(n)]

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

check("det C = -1 (recurrence reversible over Z)", det3(Cmat) == -1)

def kron(X, Y):
    rx, cx, ry, cy = len(X), len(X[0]), len(Y), len(Y[0])
    return [[X[i][j] * Y[k][l] for j in range(cx) for l in range(cy)]
            for i in range(rx) for k in range(ry)]

M = kron(kron(Cmat, Cmat), Cmat)          # 27 x 27 integer matrix
print(f"M = C(x)C(x)C built: {len(M)}x{len(M[0])} integer matrix.")

# Characteristic polynomial chi_M via Faddeev-LeVerrier (exact rationals).
def charpoly(Mm):
    n = len(Mm)
    F = [[Fraction(x) for x in row] for row in Mm]
    Ak = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        Ak[i][i] = Fraction(1)
    coeffs = [Fraction(1)]                 # chi(x) = x^n + c1 x^{n-1} + ... + cn
    for k in range(1, n + 1):
        Ak = mat_mul(F, Ak)
        tr = sum(Ak[i][i] for i in range(n))
        ck = -tr / k
        coeffs.append(ck)
        for i in range(n):
            Ak[i][i] += ck
    return coeffs

chi = charpoly(M)
check("chi_M has integer coefficients", all(f.denominator == 1 for f in chi))
chi = [int(f) for f in chi]
print(f"chi_M(x) = x^27 + ... ; constant term chi_M(0) = {chi[-1]}")
check("chi_M(0) != 0  (recurrence for d_n reversible)", chi[-1] != 0)

def chi_eval(t):
    v = 0
    for cf in chi:
        v = v * t + cf
    return v

check("chi_M(-1) = 0  (so (-1)^n satisfies the order-27 recurrence)", chi_eval(-1) == 0)

# Proof skeleton, now made exact:
#   s_n = (u_n, u_{n+1}, u_{n+2})^T satisfies s_{n+1} = C s_n; the homogeneous
#   recurrence holds for a,b,c at all n>=3 and for the backward extension at all
#   n in Z, and C is invertible over Z, so u_n = e_1^T C^n s_0 for ALL n in Z.
#   Hence u_n^3 = (e_1 (x) e_1 (x) e_1)^T M^n (s_0 (x) s_0 (x) s_0), a linear
#   functional of M^n.  By Cayley-Hamilton, chi_M(M) = 0, so every linear
#   functional of M^n -- in particular a_n^3, b_n^3, c_n^3 -- satisfies
#       sum_{i=0}^{27} chi_i * t_{n+27-i} = 0   for all n in Z          (*)
#   (valid for negative n too since M is invertible: det M = det(C)^9 = -1...
#   checked below).  And (-1)^n satisfies (*) because chi_M(-1) = 0.
#   So d_n = a_n^3 + b_n^3 - c_n^3 - (-1)^n satisfies (*).  chi is monic with
#   chi(0) != 0, so 27 consecutive zeros of d force d == 0 on all of Z.

# det M check (Kronecker: det(A(x)B) = det(A)^m det(B)^n):
# det M = det(C)^9 * det(C(x)C)^3 = ... easier: eigenvalue product = (abc)^9 = (-1)^9.
# We verify via chi: det M = (-1)^27 * chi_M(0).
detM = (-1) ** 27 * chi[-1]
check("det M = -1 (M invertible over Z; (*) holds for all n in Z)", detM == -1)

# Verify the recurrence (*) really annihilates d_n on our computed window,
# and that d_n = 0 on >= 27 consecutive n.
def d_at(n):
    return at(A, n) ** 3 + at(B, n) ** 3 - at(C_, n) ** 3 - (1 if n % 2 == 0 else -1)

window = list(range(-N_BACK, N_FWD + 1))
dvals = {n: d_at(n) for n in window}
zeros_ok = all(v == 0 for v in dvals.values())
check(f"d_n = 0 for ALL n in [{-N_BACK}, {N_FWD}]  ({len(window)} values, exact bigints)", zeros_ok)

# (*) check on the window (redundant given zeros, but confirms the machinery):
rec_ok = True
for n0 in range(-N_BACK, N_FWD + 1 - 27):
    acc = 0
    for i, cf in enumerate(chi):           # chi_0 x^27 + ... + chi_27
        acc += cf * dvals[n0 + 27 - i]
    if acc != 0:
        rec_ok = False
check("order-27 recurrence (*) annihilates d_n on the window", rec_ok)

print()
if zeros_ok and rec_ok and chi_eval(-1) == 0 and chi[-1] != 0:
    print("PROVED (all n in Z, hence in particular all n >= 0):")
    print("    a_n^3 + b_n^3 = c_n^3 + (-1)^n")
    print("for Ramanujan's sequences, INCLUDING the backward extension; the n = -1")
    print("term is the taxicab identity 9^3 + 10^3 = 12^3 + 1 = 1729.")
    print("Proof inputs, all exact integer arithmetic: det C = -1; chi_M integral,")
    print("monic, chi_M(0) = " + str(chi[-1]) + ", chi_M(-1) = 0; d_n = 0 at 51 >= 27")
    print("consecutive indices.  No floats anywhere.")
else:
    print("*** PROOF INCOMPLETE -- see FAIL lines above ***")

print()
print("=" * 78)
print("CITATIONS / READING")
print("=" * 78)
print("- Ramanujan, lost notebook (the three generating functions; see also")
print("  M. D. Hirschhorn, 'An amazing identity of Ramanujan', Math. Mag. 68 (1995)")
print("  199-201, and Hirschhorn, Math. Mag. 69 (1996) 267-269 for proofs).")
print("- K. Ono, S. Trebat-Leder, 'The 1729 K3 surface', Res. Number Theory (2015),")
print("  2:26 (2016), arXiv:1510.00735: Ramanujan's near-miss family comes from his")
print("  study of Euler's quartic surface X^3+Y^3 = Z^3+W^3; he had found (decades")
print("  before the term existed) an elliptically-fibered K3 surface of Picard")
print("  number 18 whose cubic twists give infinitely many elliptic curves over Q")
print("  of rank >= 2; 1729 sits in the family (our n = -1).  [Full text not")
print("  re-fetched this session -- web fetch failed; abstract + press coverage")
print("  verified via search results, accessed 2026-06-10.]")
print("- Eisenstein coda: x^3 + y^3 - z^3 = +-1 says the norm-form factorization")
print("  (x+y)(x+wy)(x+w^2y) = z^3 +- 1 in Z[w], w = zeta_6 - 1... precisely:")
print("  x^3+y^3 = N(x+wy)*(x+y) -- misses the Fermat cube z^3 by a UNIT of Z[w].")
print("  Near-misses are unit-distance-from-Fermat statements; the unit group of")
print("  Z[w] is finite (6 roots of unity), which is why +-1 is the only possible")
print("  'miss' and why the family is infinite yet never closes the gap (FLT n=3,")
print("  Euler/Gauss: x^3+y^3=z^3 has no nontrivial solution).")

print()
print("=" * 78)
if failures == 0:
    print("ALL CHECKS PASS (0 failures)")
else:
    print(f"*** {failures} CHECK(S) FAILED ***")
sys.exit(0 if failures == 0 else 1)
