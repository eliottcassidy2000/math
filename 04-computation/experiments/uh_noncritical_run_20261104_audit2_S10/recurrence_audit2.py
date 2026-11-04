"""recurrence_audit2.py -- independent audit (audit2, 2026-11-04) of THM-4868 (H): the run-length
recurrence, its elimination from the coupled (N, N_2) system, the generating function, the closed forms
of HYP-9458, and comparison with (i) this audit's own class census (classes_audit2 output),
(ii) the note's section-4 table, (iii) HYP-9460's pre-registered rows 14-18, (iv) HYP-9458's published
rows (orders 5-13, P3-P6 at 14-16, and the full order-16 row).

Usage: python3 recurrence_audit2.py [classes_audit2 output file ...]
"""
import re
import sys
from math import comb

import sympy as sp

NMAX = 18
ok_all = True


def check(cond, msg):
    global ok_all
    print(("  OK   " if cond else "  FAIL ") + msg)
    if not cond:
        ok_all = False


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


# ---------------------------------------------------------------- recurrence and coupled system
def rows_from_recurrence(nmax):
    N = {4: {2: 1, 3: 0}, 5: {2: 1, 3: 1, 4: 1}}
    for n in range(6, nmax + 1):
        N[n] = {}
        for r in range(2, n):
            N[n][r] = (N[n - 1].get(r, 0) + 2 * N[n - 1].get(r - 1, 0) - N[n - 2].get(r - 1, 0))
        # the recurrence must not put mass outside 2..n-1
        for r in list(range(-3, 2)) + [n, n + 1, n + 2]:
            v = N[n - 1].get(r, 0) + 2 * N[n - 1].get(r - 1, 0) - N[n - 2].get(r - 1, 0)
            assert v == 0, (n, r, v)
    return N


def rows_from_coupled(nmax):
    """The two lines of the proof of (H):  N(m+1,r) = N(m,r) + N(m,r-1) + N2(m,r-2),
    N2(m+1,r) = N(m,r) + N2(m,r-1), from N(4,.) = N2(4,.) = [r=2] (S_4 has two sole entries)."""
    N = {4: {2: 1}}
    N2 = {4: {2: 1}}
    for m in range(4, nmax):
        N[m + 1], N2[m + 1] = {}, {}
        for r in range(0, m + 3):
            a = N[m].get(r, 0) + N[m].get(r - 1, 0) + N2[m].get(r - 2, 0)
            b = N[m].get(r, 0) + N2[m].get(r - 1, 0)
            if a:
                N[m + 1][r] = a
            if b:
                N2[m + 1][r] = b
    return N, N2


Nrec = rows_from_recurrence(NMAX)
Ncpl, N2cpl = rows_from_coupled(NMAX)

print("1. Elimination: the recurrence reproduces the coupled system of the proof (orders 4..%d)" % NMAX)
for n in range(4, NMAX + 1):
    a = [Nrec[n].get(r, 0) for r in range(2, n)]
    b = [Ncpl[n].get(r, 0) for r in range(2, n)]
    extra = [r for r in Ncpl[n] if not (2 <= r <= n - 1) and Ncpl[n][r]]
    check(a == b and not extra, "n=%d recurrence row == coupled-system row %s" % (n, a))
print("   totals: u(n) = sum_r N(n,r) and b(n) = sum_r N2(n,r) against F(2n-6), F(2n-7)")
for n in range(4, NMAX + 1):
    u = sum(Ncpl[n].values())
    b = sum(N2cpl[n].values())
    check(u == fib(2 * n - 6) and b == fib(2 * n - 7), "n=%d u=%d (F=%d) b=%d (F=%d)" % (n, u, fib(2 * n - 6), b, fib(2 * n - 7)))

# ---------------------------------------------------------------- generating function
print("\n2. Generating function t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2)")
x, t = sp.symbols("x t")
G = t**2 * x**4 * (1 + (t**2 - t) * x) / (1 - (1 + 2 * t) * x + t * x**2)
ser = sp.series(G, x, 0, NMAX + 1).removeO()
ser = sp.expand(ser)
for n in range(4, NMAX + 1):
    coeff = sp.Poly(ser.coeff(x, n), t)
    row = {int(m[0]): int(c) for m, c in zip(coeff.monoms(), coeff.coeffs())}
    want = {r: v for r, v in Nrec[n].items() if v}
    check(row == want, "x^%d coefficient == recurrence row" % n)
check(sp.simplify(G.subs(t, 1) - x**4 / (1 - 3 * x + x**2)) == 0, "G(x, 1) = x^4 / (1 - 3x + x^2)")
# symbolic solution of the coupled system
A, B = sp.symbols("A B")
sol = sp.solve([sp.Eq(A, t**2 * x**4 + x * ((1 + t) * A + t**2 * B)), sp.Eq(B, t**2 * x**4 + x * (A + t * B))], [A, B], dict=True)[0]
check(sp.simplify(sol[A] - G) == 0, "coupled system solved symbolically gives exactly G")
G2 = t**2 * x**4 * (1 - t * x) / (1 - (1 + 2 * t) * x + t * x**2)
check(sp.simplify(sol[B] - G2) == 0,
      "N_2 generating function = t^2 x^4 (1 - t x) / (1 - (1 + 2t) x + t x^2)  [audit-derived, not in the note]")
ser2 = sp.expand(sp.series(G2, x, 0, NMAX + 1).removeO())
for n in range(4, NMAX + 1):
    coeff = sp.Poly(ser2.coeff(x, n), t)
    row = {int(m[0]): int(c) for m, c in zip(coeff.monoms(), coeff.coeffs())}
    want = {r: v for r, v in N2cpl[n].items() if v}
    check(row == want, "N_2: x^%d coefficient == coupled-system N_2 row" % n)

# ---------------------------------------------------------------- tables quoted in the repository
print("\n3. Tables quoted in the repository")
note_table = {  # note section 4, N(n, r) r = 2..n-1
    5: [1, 1, 1], 6: [1, 2, 3, 2], 7: [1, 3, 6, 7, 4], 8: [1, 4, 10, 16, 16, 8], 9: [1, 5, 15, 30, 41, 36, 16],
    10: [1, 6, 21, 50, 85, 102, 80, 32], 11: [1, 7, 28, 77, 155, 231, 248, 176, 64],
    12: [1, 8, 36, 112, 259, 456, 608, 592, 384, 128], 13: [1, 9, 45, 156, 406, 819, 1289, 1560, 1392, 832, 256]}
hyp9460 = {  # HYP-9460 P7 predicted rows r = 2..n-1
    14: [1, 10, 55, 210, 606, 1372, 2471, 3530, 3920, 3232, 1792, 512],
    15: [1, 11, 66, 275, 870, 2178, 4396, 7183, 9420, 9680, 7424, 3840, 1024],
    16: [1, 12, 78, 352, 1210, 3312, 7380, 13504, 20256, 24600, 23552, 16896, 8192, 2048],
    17: [1, 13, 91, 442, 1639, 4862, 11826, 23868, 40081, 55692, 63072, 56576, 38144, 17408, 4096],
    18: [1, 14, 105, 546, 2171, 6930, 18238, 40140, 74313, 115598, 149856, 159168, 134400, 85504, 36864, 8192]}
hyp9458 = {  # HYP-9458 / S74 note section 4, T(n, k) k = 2..n-1
    5: [1, 1, 1], 6: [2, 3, 2, 1], 7: [4, 7, 6, 3, 1], 8: [8, 16, 16, 10, 4, 1], 9: [16, 36, 41, 30, 15, 5, 1],
    10: [32, 80, 102, 85, 50, 21, 6, 1], 11: [64, 176, 248, 231, 155, 77, 28, 7, 1],
    12: [128, 384, 592, 608, 456, 259, 112, 36, 8, 1], 13: [256, 832, 1392, 1560, 1289, 819, 406, 156, 45, 9, 1],
    16: [2048, 8192, 16896, 23552, 24600, 20256, 13504, 7380, 3312, 1210, 352, 78, 12, 1]}
hyp9458_p = {14: (512, 1792, 3232, 1, 10, 55), 15: (1024, 3840, 7424, 1, 11, 66), 16: (2048, 8192, 16896, 1, 12, 78)}


def Trow(n):
    return [Nrec[n].get(n + 1 - k, 0) for k in range(2, n)]


for n, row in note_table.items():
    check([Nrec[n].get(r, 0) for r in range(2, n)] == row, "note section 4 row n=%d == recurrence" % n)
for n, row in hyp9460.items():
    check([Nrec[n].get(r, 0) for r in range(2, n)] == row, "HYP-9460 P7 row n=%d == recurrence" % n)
for n, row in hyp9458.items():
    check(Trow(n) == row, "HYP-9458 published T(%d, k) row == recurrence read through T(n,k) = N(n, n+1-k)" % n)
for n, vals in hyp9458_p.items():
    T = Trow(n)
    got = (T[0], T[1], T[2], T[n - 3], T[n - 4], T[n - 5])
    check(got == vals, "HYP-9458 P3-P6 at n=%d: %s" % (n, got))

# ---------------------------------------------------------------- closed forms
print("\n4. HYP-9458 closed forms from the recurrence (n up to 60)")
Nbig = rows_from_recurrence(60)


def T_(n, k):
    return Nbig[n].get(n + 1 - k, 0) if n >= 4 else 0


lo = {}
for name, k, f, nmin in [("T(n,2) = 2^(n-5)", 2, lambda n: sp.Rational(2) ** (n - 5), 5),
                         ("T(n,3) = n 2^(n-7)", 3, lambda n: n * sp.Rational(2) ** (n - 7), 6),
                         ("T(n,4) = (n^2+n-8) 2^(n-10)", 4, lambda n: (n * n + n - 8) * sp.Rational(2) ** (n - 10), 7)]:
    good = all(T_(n, k) == f(n) for n in range(nmin, 61))
    fails_below = [n for n in range(4, nmin) if T_(n, k) != f(n)]
    check(good, "%s for %d <= n <= 60 (fails at n = %s, below the stated range)" % (name, nmin, fails_below))
for name, r, f in [("N(n,2) = 1", 2, lambda n: 1), ("N(n,3) = n - 4", 3, lambda n: n - 4), ("N(n,4) = C(n-3,2)", 4, lambda n: comb(n - 3, 2))]:
    check(all(Nbig[n].get(r, 0) == f(n) for n in range(5, 61)), "top diagonal %s for 5 <= n <= 60" % name)
# the note's intermediate one-step recurrences, with their true ranges
r3 = [n for n in range(6, 61) if T_(n, 3) != 2 ** (n - 7) + 2 * T_(n - 1, 3)]
r4 = [n for n in range(7, 61) if T_(n, 4) != sp.Rational(n) * sp.Rational(2) ** (n - 9) + 2 * T_(n - 1, 4)]
print("  note's T(n,3) = 2^(n-7) + 2 T(n-1,3) fails at n in %s (valid for n >= 7)" % r3)
print("  note's T(n,4) = n 2^(n-9) + 2 T(n-1,4) fails at n in %s (valid for n >= 8)" % r4)
# symbolic: closed forms satisfy the T-recurrence
nn = sp.symbols("n", integer=True)
T2 = 2 ** (nn - 5)
T3 = nn * 2 ** (nn - 7)
T4 = (nn**2 + nn - 8) * 2 ** (nn - 10)
check(sp.simplify(T2 - (0 + 2 * T2.subs(nn, nn - 1) - 0)) == 0, "symbolic: T2 satisfies T(n,2) = 2T(n-1,2) (T(.,1) = 0)")
check(sp.simplify(T3 - (T2.subs(nn, nn - 1) + 2 * T3.subs(nn, nn - 1) - T2.subs(nn, nn - 2))) == 0,
      "symbolic: T3 satisfies T(n,3) = T(n-1,2) + 2T(n-1,3) - T(n-2,2)")
check(sp.simplify(T4 - (T3.subs(nn, nn - 1) + 2 * T4.subs(nn, nn - 1) - T3.subs(nn, nn - 2))) == 0,
      "symbolic: T4 satisfies T(n,4) = T(n-1,3) + 2T(n-1,4) - T(n-2,3)")

# ---------------------------------------------------------------- this audit's census
print("\n5. This audit's own class census (classes_audit2 output) against the recurrence")
for fn in sys.argv[1:]:
    txt = open(fn).read()
    blocks = re.findall(r"n=\s*(\d+) classes=(\d+).*?N\(n,r\) r=2\.\.n-1:([\d ]+)\n\s+N2\(n,r\) r=2\.\.n-1:([\d ]+)\n\s+T\(n,k\) k=2\.\.n-1:([\d ]+)", txt, re.S)
    for n, c, nr, n2r, tk in blocks:
        n = int(n)
        nr = [int(v) for v in nr.split()]
        n2r = [int(v) for v in n2r.split()]
        tk = [int(v) for v in tk.split()]
        check(nr == [Nrec[n].get(r, 0) for r in range(2, n)], "%s: n=%d N row %s" % (fn, n, nr))
        check(n2r == [N2cpl[n].get(r, 0) for r in range(2, n)], "%s: n=%d N2 row %s" % (fn, n, n2r))
        check(tk == list(reversed(nr)), "%s: n=%d cut-arc row is the run row reversed" % (fn, n))
        if n in hyp9458:
            check(tk == hyp9458[n], "%s: n=%d cut-arc row == HYP-9458 published row" % (fn, n))
        check(int(c) == fib(2 * n - 6), "%s: n=%d classes = F(2n-6)" % (fn, n))
    cblocks = re.findall(r"childN\((\d+),r\) r=2\.\.\d+:([\d ]+)\n\s+childN2\(\d+,r\) r=2\.\.\d+:([\d ]+)\n\s+childT\(\d+,k\) k=2\.\.\d+:([\d ]+)", txt)
    for n, nr, n2r, tk in cblocks:
        n = int(n)
        nr = [int(v) for v in nr.split()]
        n2r = [int(v) for v in n2r.split()]
        tk = [int(v) for v in tk.split()]
        check(nr == [Nrec[n].get(r, 0) for r in range(2, n)], "%s: children of order %d, N row %s" % (fn, n, nr))
        check(n2r == [N2cpl[n].get(r, 0) for r in range(2, n)], "%s: children of order %d, N2 row" % (fn, n))
        check(tk == list(reversed(nr)), "%s: children of order %d, cut-arc row is the run row reversed" % (fn, n))
        if n in hyp9458:
            check(tk == hyp9458[n], "%s: children of order %d, cut-arc row == HYP-9458 published row %s" % (fn, n, hyp9458[n]))
        if n in hyp9458_p:
            check((tk[0], tk[1], tk[2], tk[n - 3], tk[n - 4], tk[n - 5]) == hyp9458_p[n], "%s: order %d, HYP-9458 P3-P6 values measured" % (fn, n))
        check(sum(nr) == fib(2 * n - 6), "%s: children of order %d number F(2n-6) = %d" % (fn, n, fib(2 * n - 6)))

print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
