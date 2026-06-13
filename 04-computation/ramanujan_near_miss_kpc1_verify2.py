# =====================================================================
# ADVERSARIAL VERIFIER (independent pass 2) -- session kind-pasteur-2026-06-10-S1
# Independent check of thread D claims D7 / D8 (Ramanujan near-miss family).
# FRESH code; methods deliberately DIFFERENT from the worker:
#   * coefficients via direct rational-function long division (convolution
#     against D(x), numerator taps), NOT recurrence-first;
#   * all-n proof via an ORDER-7 annihilator q(x) built from the exact
#     factorization chi_C(x) = (x+1)(x^2-83x+1) and integer symmetric
#     functions of the quadratic's roots -- NOT the worker's 27x27
#     Cayley-Hamilton/Faddeev-LeVerrier route;
#   * the worker's chi_M(0)=1, chi_M(-1)=0, det M=-1 cross-checked by
#     fraction-free Bareiss integer determinants of the explicit 27x27
#     Kronecker cube (two determinants, no char poly).
# Pure python ints. No floats. Window n = -15..45 (61 values > worker's 51).
# =====================================================================
import sys
from math import isqrt

fails = 0
def check(name, cond, detail=""):
    global fails
    line = ("PASS " if cond else "FAIL ") + name
    if detail:
        line += "  |  " + detail
    print(line)
    if not cond:
        fails += 1

def sgn(n):                       # (-1)^n for any integer n, no floats
    return 1 if n % 2 == 0 else -1

# ----------------------------------------------------------------------
# PART 1: series coefficients by direct long division of N(x)/D(x)
# D(x) = 1 - 82x - 82x^2 + x^3  =>  u_n = N_n + 82u_{n-1} + 82u_{n-2} - u_{n-3}
# (u_j = 0 for j < 0; N_n = 0 for n >= 3)
# ----------------------------------------------------------------------
NA = [1, 53, 9]
NB = [2, -26, -12]
NC = [2, 8, -10]
NHI, NLO = 45, -15

def expand(num):
    u = {}
    for n in range(0, NHI + 1):
        nn = num[n] if n < 3 else 0
        u[n] = (nn + 82 * u.get(n - 1, 0) + 82 * u.get(n - 2, 0)
                - u.get(n - 3, 0))
    # canonical backward extension of the reversible homogeneous recurrence:
    # u_n = -u_{n+3} + 82u_{n+2} + 82u_{n+1}
    for n in range(-1, NLO - 1, -1):
        u[n] = -u[n + 3] + 82 * u[n + 2] + 82 * u[n + 1]
    return u

A = expand(NA)
B = expand(NB)
C = expand(NC)

print("=" * 70)
print("PART 1: forward triples from direct series division")
print("=" * 70)
for n in range(0, 6):
    print("  n=%2d  (a,b,c) = (%d, %d, %d)" % (n, A[n], B[n], C[n]))
check("D7: (a,b,c)_0 == (1,2,2)", (A[0], B[0], C[0]) == (1, 2, 2))
check("D7: (a,b,c)_1 == (135,138,172)", (A[1], B[1], C[1]) == (135, 138, 172))
check("D7: (a,b,c)_2 == (11161,11468,14258)",
      (A[2], B[2], C[2]) == (11161, 11468, 14258))
check("D7: n=10 entries have 20 digits",
      (len(str(A[10])), len(str(B[10])), len(str(C[10]))) == (20, 20, 20),
      "a_10=%d" % A[10])

print()
print("=" * 70)
print("PART 2: cubic identity a^3+b^3-c^3 = (-1)^n on n = %d..%d" % (NLO, NHI))
print("=" * 70)
bad = []
for n in range(NLO, NHI + 1):
    d = A[n]**3 + B[n]**3 - C[n]**3 - sgn(n)
    if d != 0:
        bad.append((n, d))
check("D7: d_n == 0 exactly for all %d indices n in [%d,%d]"
      % (NHI - NLO + 1, NLO, NHI), not bad, "violations: %s" % bad[:3])

print()
print("=" * 70)
print("PART 3: backward extension and 1729 (D8)")
print("=" * 70)
for n in (-1, -2, -3):
    print("  n=%2d  (a,b,c) = (%d, %d, %d)" % (n, A[n], B[n], C[n]))
check("D8: (a,b,c)_{-1} == (-9, 12, 10)", (A[-1], B[-1], C[-1]) == (-9, 12, 10))
check("D8: (a,b,c)_{-2} == (-791, 1010, 812)",
      (A[-2], B[-2], C[-2]) == (-791, 1010, 812))
check("D8: (a,b,c)_{-3} == (-65601, 83802, 67402)",
      (A[-3], B[-3], C[-3]) == (-65601, 83802, 67402))
check("D8: (-9)^3 + 12^3 == 10^3 - 1 == 999",
      (-9)**3 + 12**3 == 10**3 - 1 == 999)
check("D8: rearranged: 9^3 + 10^3 == 1729 == 12^3 + 1",
      9**3 + 10**3 == 1729 and 12**3 + 1 == 1729)

# ----------------------------------------------------------------------
# PART 4: ALL-n proof, independent route (order-7 annihilator).
#
# chi_C(x) = x^3 - 82x^2 - 82x + 1 = (x+1)(x^2 - 83x + 1) =: (x+1) g(x).
# g has discriminant 6885 > 0, not a square => phi != psi irrational,
# and g(-1) = 85 != 0 => the three roots {-1, phi, psi} are distinct
# => companion C is diagonalizable => a_n = P(-1)^n + Q phi^n + R psi^n.
# Cubing: a_n^3 is a Z-bar-linear combination of L^n where L runs over
# triple products of roots. With phi*psi = 1, phi+psi = 83 the complete
# multiset of triple-product VALUES is
#   {-1 (from (-1)^3 and (-1)phi*psi), phi, psi, -phi^2, -psi^2, phi^3, psi^3}
# with exact integer symmetric data:
#   phi^2+psi^2 = 83^2 - 2 = 6887,  (phi*psi)^2 = 1
#   phi^3+psi^3 = 83^3 - 3*83 = 571538,  (phi*psi)^3 = 1
# so EVERY triple product is a root of
#   q(x) = (x+1) (x^2-83x+1) (x^2+6887x+1) (x^2-571538x+1),   deg q = 7.
# q also kills (-1)^n (factor x+1). Hence d_n = a^3+b^3-c^3-(-1)^n obeys
# the order-7 integer recurrence with char poly q; q is monic with
# q(0) = 1*1*1*1 = 1, so the recurrence is a unit at BOTH ends and runs
# in both directions over Z. Seven consecutive exact zeros (PART 2 gives
# 61) therefore force d == 0 on all of Z. QED (textbook steps only).
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("PART 4: order-7 annihilator route for the all-n proof (D7)")
print("=" * 70)

def polmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            r[i + j] += pi * qj
    return r

# exact factorization check: (x+1)(x^2-83x+1) == x^3 - 82x^2 - 82x + 1
chiC = polmul([1, 1], [1, -83, 1])            # ascending-degree coeff lists?
# NOTE: store DESCENDING is error-prone; use ascending: p[i] = coeff of x^i
# rebuild ascending: (1 + x) * (1 - 83x + x^2)
chiC = polmul([1, 1], [1, -83, 1])
check("chi_C factorization: (x+1)(x^2-83x+1) == 1 - 82x - 82x^2 + x^3 "
      "(ascending coeffs)", chiC == [1, -82, -82, 1], "got %s" % chiC)
disc = 83 * 83 - 4
check("g = x^2-83x+1: disc = 6885 > 0 and NOT a perfect square (phi != psi, "
      "irrational)", disc == 6885 and isqrt(disc)**2 != disc)
check("g(-1) = 85 != 0 (so -1 is not a root of g; 3 distinct eigenvalues)",
      1 + 83 + 1 == 85)
check("symmetric data: phi^2+psi^2 == 6887", 83**2 - 2 == 6887)
check("symmetric data: phi^3+psi^3 == 571538", 83**3 - 3 * 83 == 571538)

# q ascending: (1+x)(1-83x+x^2)(1+6887x+x^2)(1-571538x+x^2)
q = polmul(polmul([1, 1], [1, -83, 1]),
           polmul([1, 6887, 1], [1, -571538, 1]))
check("q monic of degree 7 with q(0) = 1 (unit ends => reversible over Z)",
      len(q) == 8 and q[7] == 1 and q[0] == 1, "q = %s" % q)
qm1 = sum(qi * sgn(i) for i, qi in enumerate(q))
check("q(-1) == 0 (annihilates (-1)^n)", qm1 == 0)

# the q-recurrence must annihilate a^3, b^3, c^3, (-1)^n, and d on the window
def annihilates(seq_at, lo, hi, name):
    viol = [n for n in range(lo, hi - 7 + 1)
            if sum(q[i] * seq_at(n + i) for i in range(8)) != 0]
    check("q-recurrence annihilates %s on n in [%d,%d]" % (name, lo, hi),
          not viol, "violations at n=%s" % viol[:3])

annihilates(lambda n: A[n]**3, NLO, NHI, "a_n^3")
annihilates(lambda n: B[n]**3, NLO, NHI, "b_n^3")
annihilates(lambda n: C[n]**3, NLO, NHI, "c_n^3")
annihilates(lambda n: sgn(n), NLO, NHI, "(-1)^n")

# ----------------------------------------------------------------------
# PART 5: cross-check the worker's 27x27 evidence numbers by Bareiss.
# State map: (u_n,u_{n+1},u_{n+2}) -> (u_{n+1},u_{n+2},82u_{n+2}+82u_{n+1}-u_n)
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("PART 5: Bareiss cross-check of chi_M(0), chi_M(-1), det M (worker route)")
print("=" * 70)
Cm = [[0, 1, 0], [0, 0, 1], [-1, 82, 82]]

def matmul(X, Y):
    n = len(X)
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]

# Cayley-Hamilton sanity for C itself: C^3 - 82C^2 - 82C + I == 0
C2 = matmul(Cm, Cm); C3 = matmul(C2, Cm)
ch = all(C3[i][j] - 82 * C2[i][j] - 82 * Cm[i][j] + (1 if i == j else 0) == 0
         for i in range(3) for j in range(3))
check("Cayley-Hamilton for C: C^3 - 82C^2 - 82C + I == 0", ch)
detC = (Cm[0][0]*(Cm[1][1]*Cm[2][2]-Cm[1][2]*Cm[2][1])
        - Cm[0][1]*(Cm[1][0]*Cm[2][2]-Cm[1][2]*Cm[2][0])
        + Cm[0][2]*(Cm[1][0]*Cm[2][1]-Cm[1][1]*Cm[2][0]))
check("det C == -1 (reversible state map)", detC == -1)
# companion really reproduces the sequences (independent of PART 1 route)
for nm, U in (("a", A), ("b", B), ("c", C)):
    v = [U[0], U[1], U[2]]
    okrec = True
    for n in range(0, 20):
        if v[0] != U[n]:
            okrec = False
            break
        v = [v[1], v[2], 82 * v[2] + 82 * v[1] - v[0]]
    check("companion state map reproduces %s_n (n=0..19)" % nm, okrec)

def kron(X, Y):
    nx, ny = len(X), len(Y)
    return [[X[i][j] * Y[k][l] for j in range(nx) for l in range(ny)]
            for i in range(nx) for k in range(ny)]

M27 = kron(kron(Cm, Cm), Cm)
check("M = C tensor C tensor C is 27x27", len(M27) == 27)

def bareiss_det(mat):
    n = len(mat)
    Mm = [row[:] for row in mat]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if Mm[k][k] == 0:
            piv = next((i for i in range(k + 1, n) if Mm[i][k] != 0), None)
            if piv is None:
                return 0
            Mm[k], Mm[piv] = Mm[piv], Mm[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                Mm[i][j] = (Mm[i][j] * Mm[k][k] - Mm[i][k] * Mm[k][j]) // prev
            Mm[i][k] = 0
        prev = Mm[k][k]
    return sign * Mm[n - 1][n - 1]

negM = [[-e for e in row] for row in M27]
chiM0 = bareiss_det(negM)                       # chi_M(0) = det(0*I - M)
check("chi_M(0) == 1 (worker claim; => det M = -1, 27 odd)", chiM0 == 1,
      "det(-M) = %d" % chiM0)
negImM = [[-e - (1 if i == j else 0) for j, e in enumerate(row)]
          for i, row in enumerate(M27)]
chiMm1 = bareiss_det(negImM)                    # chi_M(-1) = det(-I - M)
check("chi_M(-1) == 0 (worker claim; (-1)^n inside the order-27 recurrence)",
      chiMm1 == 0, "det(-I-M) = %d" % chiMm1)
check("det M == -1", -chiM0 == -1)

print()
print("=" * 70)
print("VERIFIER SUMMARY: %s (%d failures)" % ("ALL CHECKS PASS" if fails == 0
      else "FAILURES PRESENT", fails))
print("=" * 70)
sys.exit(0 if fails == 0 else 1)
