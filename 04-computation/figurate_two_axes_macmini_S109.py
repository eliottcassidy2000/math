#!/usr/bin/env python3
"""figurate_two_axes_macmini_S109.py -- mac-mini-2026-07-15-S109 (owner directive).
THE TWO ORTHOGONAL LIFTS OF THE TRIANGULAR NUMBERS: polygonal (shape axis) vs
polyhedral/simplex (dimension axis), their two triangles, and the exact structure of
their difference.

Conventions (left-justified triangles, row r >= 0, column k = 0..r; m := r - k + 1 is
the index-within-column, so column k's first entry (m=1) sits in row k):
  PASCAL / polyhedral: P_hedral(r,k) = C(r,k); column k top-down = k-simplex numbers
                       C(m+k-1, k).
  POLYGONAL:           Q(r,0) = 1; Q(r,k) = (k+1)-gonal number at index m
                       = C(m,1) + (k-1)*C(m,2).

Verified here:
 (1) columns agree for k <= 2 and for m <= 2 (shared frame).
 (2) ROW SUMS: Pascal -> 2^r; polygonal -> 1 + C(r+1,2) + C(r+1,4) = A000127(r+1)
     (Moser's circle problem), via the double hockey stick sum_j (r-j) C(j,2) = C(r+1,4).
 (3) DIAGONAL (skipped-row) SUMS: Pascal -> Fibonacci; polygonal -> G(n) =
     1,1,2,3,5,8,13,21,33,51,76,111,157,218,295,... ; GF, order-7 recurrence,
     quasi-polynomial closed form, asymptotic n^4/96.
 (4) THE DIFFERENCE: Delta(m,k) = C(m+k-1,k) - [C(m,1)+(k-1)C(m,2)]
     = sum_{j>=2} C(m, j+1) * C(k-1, j)     (Vandermonde tail),
     so diagonal m has coefficient row (C(m,3), C(m,4), ..., C(m,m)) on the basis
     C(k-1,j) -- Pascal's row m with its first three entries deleted.
     First diagonal (m=3): C(k-1,2) = triangular numbers (owner's data, with the typo
     36 (not 35) fixed); second diagonal (m=4): 4C(k-1,2) + C(k-1,3) = 4,13,28,50,80,...;
     down-column (k=3): C(m,3) = tetrahedral numbers.
 (5) THREE INTERPRETATIONS of the polygonal entry (support-<=2 truncation):
     (a) multisets of size k from m types with at most 2 distinct types;
     (b) k-subsets of [r] with at most 2 maximal runs;
     (c) lattice points on the 1-SKELETON of the dilated simplex k*Delta_{m-1}
     -- all equal Q(r,k); and exact-run/support-s counts = C(m,s)C(k-1,s-1).
 (6) THE RUN-TRUNCATION TOWER T_j (entries = k-subsets with <= j runs):
     row sums = sum_{i<=j} C(r+1, 2i)  [j=1: 1+C(r+1,2); j=2: Moser; -> 2^r],
     diagonal sums G_j(n): G_1 = 1 + floor(n^2/4) (A033638-type), G_2 = G, -> Fibonacci.
 (7) MASTER FORMULA N(s,d,m) = (s-2) C(m+d-2, d) + C(m+d-2, d-1): both triangles are
     orthogonal slices; s-axis affine (gluing triangles: P(s+1,m) = P(s,m) + T_{m-1});
     d-axis integral (coning: N(s,d,m) = sum_i N(s,d-1,i)).
 (8) NO-HOLES COMPLETENESS (Brown's criterion a_{k+1} <= 1 + sum_{i<=k} a_i):
     row/diagonal sum sequences 2^n (tight), Fibonacci, Moser, G are all COMPLETE
     (distinct-part representable, no holes); the COLUMNS (squares, triangulars) FAIL
     Brown and need Fermat's polygonal theorem (k-fold repetition) instead.
 (9) the classic locker parity fact tau(n) odd <=> n square (involution d <-> n/d),
     gaps of the open set = odd numbers.
"""
import sys
from math import comb, isqrt
from fractions import Fraction
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

R = 40                                   # rows to build

def polygonal(s, m):                     # m-th s-gonal number
    return (s - 2) * comb(m, 2) + m

def Q(r, k):                             # polygonal triangle entry
    if k == 0: return 1
    m = r - k + 1
    return polygonal(k + 1, m)

def Phed(r, k):                          # Pascal / polyhedral entry
    return comb(r, k)

print("=" * 78)
print("(1) shared frame: columns k<=2 and entries m<=2 coincide")
ok = all(Q(r, k) == Phed(r, k) for r in range(R) for k in range(min(r, 2) + 1)) and \
     all(Q(k + m - 1, k) == Phed(k + m - 1, k) for k in range(R) for m in (1, 2) if k + m - 1 < R)
print("   agree on first 3 columns + first 2 entries of every column:", ok)

print("\n(2) ROW SUMS")
moser = lambda n: 1 + comb(n, 2) + comb(n, 4)        # A000127(n), n points on circle
ok2 = True
for r in range(R):
    rs_p = sum(Phed(r, k) for k in range(r + 1))
    rs_q = sum(Q(r, k) for k in range(r + 1))
    ok2 &= rs_p == 2 ** r
    ok2 &= rs_q == 1 + comb(r + 1, 2) + comb(r + 1, 4) == moser(r + 1)
print("   Pascal rows = 2^r AND polygonal rows = 1 + C(r+1,2) + C(r+1,4) = A000127(r+1):", ok2)
print("   polygonal row sums r=0..9:", [sum(Q(r, k) for k in range(r + 1)) for r in range(10)])
hs = all(sum((r - j) * comb(j, 2) for j in range(1, r + 1)) == comb(r + 1, 4) for r in range(R))
print("   double hockey stick sum_j (r-j)C(j,2) = C(r+1,4):", hs)

print("\n(3) DIAGONAL (skipped-row) SUMS")
fib = [1, 1]
for _ in range(2 * R): fib.append(fib[-1] + fib[-2])
okF = all(sum(Phed(n - k, k) for k in range(n // 2 + 1)) == fib[n] for n in range(R))
print("   Pascal shallow diagonals = Fibonacci:", okF)
G = [sum(Q(n - k, k) for k in range(n // 2 + 1)) for n in range(R)]
print("   G(n), n=0..24:", G[:25])
user = [1, 1, 2, 3, 5, 8, 13, 21, 33, 51, 76, 111, 157, 218]
print("   matches owner's 14 terms:", G[:14] == user)
# GF check: (1 - 2x + 3x^3 - 2x^4 + x^6) / ((1-x)^5 (1+x)^2)
num = [1, -2, 0, 3, -2, 0, 1]
den = [1]                                # (1-x)^5 (1+x)^2 = (1-x)^3 (1-x^2)^2
for _ in range(5):
    den = [a - b for a, b in zip(den + [0], [0] + den)]
for _ in range(2):
    den = [a + b for a, b in zip(den + [0], [0] + den)]
ser = []
for n in range(R):                       # long division
    c = (num[n] if n < len(num) else 0) - sum(den[i] * ser[n - i] for i in range(1, min(n, 7) + 1))
    ser.append(c)
print("   GF (1-2x+3x^3-2x^4+x^6)/((1-x)^5(1+x)^2) reproduces G:", ser == G)
print("   denominator coeffs (recurrence):", den, " -> G(n) = 3G(n-1)-G(n-2)-5G(n-3)+5G(n-4)+G(n-5)-3G(n-6)+G(n-7)")
okrec = all(sum(den[i] * G[n - i] for i in range(8)) == 0 for n in range(7, R))
print("   order-7 recurrence verified:", okrec)
# quasi-polynomial closed form: G(n) = p(n) + (-1)^n q(n), deg p = 4, deg q = 1
import itertools as _it
A = []; b = []
for n in range(10, 24):
    A.append([Fraction(n) ** e for e in range(5)] + [Fraction((-1) ** n) * n ** e for e in range(2)])
    b.append(Fraction(G[n]))
# solve least-exact via Gaussian elimination on first 7 independent rows
import copy
M = [row[:] + [bb] for row, bb in zip(A[:7], b[:7])]
for c in range(7):
    piv = next(i for i in range(c, 7) if M[i][c] != 0)
    M[c], M[piv] = M[piv], M[c]
    M[c] = [v / M[c][c] for v in M[c]]
    for i in range(7):
        if i != c and M[i][c] != 0:
            M[i] = [u - M[i][c] * v for u, v in zip(M[i], M[c])]
coef = [M[i][7] for i in range(7)]
def Gclosed(n):
    p = sum(coef[e] * n ** e for e in range(5))
    q = coef[5] + coef[6] * n
    return p + (-1) ** n * q
okcf = all(Gclosed(n) == G[n] for n in range(R))
print("   closed form G(n) = p(n) + (-1)^n q(n) with")
print("     p(n) =", " + ".join(f"({coef[e]}) n^{e}" for e in range(5)))
print("     q(n) =", f"({coef[5]}) + ({coef[6]}) n")
print("   exact on all", R, "terms:", okcf, "   leading coeff:", coef[4], "(predicted 1/96)")
print("   Fibonacci defect F(n+1)-G(n), n=0..14:", [fib[n + 1] - G[n] for n in range(15)])

print("\n(4) THE DIFFERENCE TRIANGLE  Delta(m,k) = Simplex - Polygonal")
def simplex(m, k): return comb(m + k - 1, k)
def poly_col(m, k): return polygonal(k + 1, m) if k >= 1 else 1
def delta(m, k): return simplex(m, k) - poly_col(m, k)
okV = all(delta(m, k) == sum(comb(m, j + 1) * comb(k - 1, j) for j in range(2, max(m, k) + 2))
          for m in range(1, 22) for k in range(1, 22))
print("   Vandermonde tail  Delta(m,k) = sum_{j>=2} C(m,j+1) C(k-1,j):", okV)
print("   first diagonal  (m=3):", [delta(3, k) for k in range(3, 10)], "= C(k-1,2) = triangulars",
      all(delta(3, k) == comb(k - 1, 2) for k in range(1, 22)))
print("     [owner's '35-21=15' corrected: the polyhedral entry over octagonal 21 is C(9,2)=36; 36-21=15]")
print("   second diagonal (m=4):", [delta(4, k) for k in range(3, 12)],
      "= 4C(k-1,2) + C(k-1,3) = (k-1)(k-2)(k+9)/6:",
      all(delta(4, k) == 4 * comb(k - 1, 2) + comb(k - 1, 3) == (k - 1) * (k - 2) * (k + 9) // 6
          for k in range(1, 22)))
print("   third diagonal  (m=5):", [delta(5, k) for k in range(3, 11)],
      "= 10C(k-1,2)+5C(k-1,3)+C(k-1,4):",
      all(delta(5, k) == 10 * comb(k - 1, 2) + 5 * comb(k - 1, 3) + comb(k - 1, 4) for k in range(1, 22)))
print("   coefficient rows = Pascal row m minus its first three entries:")
for m in range(3, 9):
    print(f"     m={m}: {[comb(m, j + 1) for j in range(2, m)]}")
print("   down-column k=3 (squares):", [delta(m, 3) for m in range(1, 9)], "= C(m,3) tetrahedrals:",
      all(delta(m, 3) == comb(m, 3) for m in range(1, 22)))

print("\n(5) THREE INTERPRETATIONS of the polygonal entry")
def multisets_support_le2(m, k):
    c = 0
    for supp in range(1, min(m, k) + 1):
        if supp <= 2: c += comb(m, supp) * comb(k - 1, supp - 1)
    return c
def subsets_runs_le2(r, k):
    if k == 0: return 1
    c = 0
    for sub in combinations(range(r), k):
        runs = 1 + sum(1 for i in range(k - 1) if sub[i + 1] > sub[i] + 1)
        if runs <= 2: c += 1
    return c
def skeleton1_points(m, k):
    # lattice points of k*Delta_{m-1} with support <= 2
    return m + comb(m, 2) * (k - 1) if k >= 1 else 1
ok5 = True
for r in range(1, 13):
    for k in range(1, r + 1):
        m = r - k + 1
        q = Q(r, k)
        ok5 &= q == multisets_support_le2(m, k) == skeleton1_points(m, k) == subsets_runs_le2(r, k)
print("   Q(r,k) = support-<=2 multisets = <=2-run subsets = 1-skeleton Ehrhart:", ok5)
def exact_runs(r, k, s):
    c = 0
    for sub in combinations(range(r), k):
        runs = 1 + sum(1 for i in range(k - 1) if sub[i + 1] > sub[i] + 1)
        if runs == s: c += 1
    return c
ok5b = all(exact_runs(r, k, s) == comb(r - k + 1, s) * comb(k - 1, s - 1)
           for r in range(1, 12) for k in range(1, r + 1) for s in range(1, k + 1))
print("   exact-run law  #(k-subsets of [r], s runs) = C(m,s) C(k-1,s-1):", ok5b)

print("\n(6) THE RUN-TRUNCATION TOWER T_j")
def Tj(r, k, j):
    if k == 0: return 1
    m = r - k + 1
    return sum(comb(m, s) * comb(k - 1, s - 1) for s in range(1, j + 1))
for j in (1, 2, 3, 4):
    okrow = all(sum(Tj(r, k, j) for k in range(r + 1)) == sum(comb(r + 1, 2 * i) for i in range(j + 1))
                for r in range(25))
    diag = [sum(Tj(n - k, k, j) for k in range(n // 2 + 1)) for n in range(16)]
    print(f"   j={j}: row sums = sum_i<= j C(r+1,2i): {okrow};  diagonal sums G_{j}(n) n=0..15: {diag}")
q2 = all(sum(Tj(n - k, k, 1) for k in range(n // 2 + 1)) == 1 + (n * n) // 4 for n in range(30))
print("   G_1(n) = 1 + floor(n^2/4)  (quarter-squares + 1):", q2)
okconv = all(Tj(r, k, r + 2) == comb(r, k) for r in range(20) for k in range(r + 1))
print("   j -> infinity recovers Pascal (rows -> 2^r, diagonals -> Fibonacci):", okconv)

print("\n(7) MASTER FORMULA N(s,d,m) = (s-2) C(m+d-2,d) + C(m+d-2,d-1)")
def N(s, d, m): return (s - 2) * comb(m + d - 2, d) + comb(m + d - 2, d - 1)
ok7 = all(N(3, k, m) == simplex(m, k) for k in range(1, 15) for m in range(1, 15))
ok7 &= all(N(s, 2, m) == polygonal(s, m) for s in range(3, 15) for m in range(1, 15))
ok7 &= all(N(s, d, m) == sum(N(s, d - 1, i) for i in range(1, m + 1))
           for s in range(3, 10) for d in range(2, 8) for m in range(1, 10))
ok7 &= all(N(s + 1, 2, m) - N(s, 2, m) == comb(m, 2) for s in range(3, 12) for m in range(1, 12))
print("   simplex slice (s=3), polygonal slice (d=2), coning recurrence (d-axis),")
print("   triangle-gluing increment (s-axis, +T_{m-1}):", ok7)
print("   classic checks: square pyramidal N(4,3,m):", [N(4, 3, m) for m in range(1, 7)],
      "  n^2 = T_{n-1}+T_n:", all(n * n == comb(n, 2) + comb(n + 1, 2) for n in range(1, 20)))

print("\n(8) NO-HOLES COMPLETENESS (Brown criterion) of the four sum-sequences")
def brown(seq, name, terms=60):
    s = 0; okB = True; tight = True
    for i, a in enumerate(seq[:terms]):
        if i == 0:
            okB &= a == 1
        else:
            okB &= a <= 1 + s
            tight &= a == 1 + s
        s += a
    return okB, tight
rows_p = [2 ** r for r in range(30)]
rows_q = [moser(n) for n in range(1, 31)]
diag_p = fib[1:31]
diag_q = G[:30] if G[0] == 1 else [1] + G[:29]
for nm, sq in (("2^r (Pascal rows)", rows_p), ("Moser A000127 (polygonal rows)", rows_q),
               ("Fibonacci (Pascal diagonals)", diag_p), ("G (polygonal diagonals)", diag_q)):
    okB, tight = brown(sq, nm)
    print(f"   {nm}: complete={okB}" + ("  [TIGHT: a_k = 1 + sum -- the binary numeration]" if tight else ""))
sq_seq = [i * i for i in range(1, 20)]
tri_seq = [comb(i + 1, 2) for i in range(1, 20)]
print("   columns FAIL Brown: squares (4 > 1+1):", not brown([1] + sq_seq[1:], "")[0],
      "; triangulars (3 > 1+1):", not brown([1] + tri_seq[1:], "")[0],
      " -> Fermat polygonal theorem territory (repetition k times)")
def cover_with_k(seq, k, N):             # every n <= N a sum of at most k terms (repetition)?
    reach = {0}
    for _ in range(k):
        reach = {a + b for a in reach for b in [0] + seq if a + b <= N}
    return all(n in reach for n in range(1, N + 1))
print("   Gauss 3-triangular covers 1..500:", cover_with_k(tri_seq + [comb(i + 1, 2) for i in range(19, 40)], 3, 500),
      ";  Lagrange 4-square covers 1..500:", cover_with_k([i * i for i in range(1, 25)], 4, 500))

print("\n(9) the integer locker (classic, for the record)")
okL = all((len([d for d in range(1, x + 1) if x % d == 0]) % 2 == 1) == (isqrt(x) ** 2 == x)
          for x in range(1, 10001))
print("   tau(n) odd <=> n square (n <= 10^4):", okL,
      ";  gaps of open lockers = odd numbers:", all((k + 1) ** 2 - k * k == 2 * k + 1 for k in range(1, 100)))

print("\nOEIS SEARCH STRINGS:")
print("   G:", ",".join(map(str, G[:18])))
print("   second-diagonal differences:", ",".join(str(delta(4, k)) for k in range(3, 13)))
print("   third-diagonal differences:", ",".join(str(delta(5, k)) for k in range(3, 13)))
print("\nALL DONE")
