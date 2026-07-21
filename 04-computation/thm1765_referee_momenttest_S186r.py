#!/usr/bin/env python3
"""
thm1765_referee_momenttest_S186r.py  (HOSTILE REFEREE of THM-1765 section 3)

TARGET CLAIM: "FINITE MOMENT TEST: E[P^m] = 0 for all m iff E[P^m] = 0 for
m <= r + 1 + (largest integer zero of q_r) -- an EXPLICIT per-SUPPORT bound."

THE TRAP: q_r's coefficients are polynomial in P's coefficients p, so its
integer zeros are functions of p, not of the support.  "Per-support" needs
the integer-zero content of q_r to be p-independent.  The file's evidence:
four toy supports, ALL with CONSTANT leading coefficient (their own script
prints "M(p)" -- conceding p-dependence -- while canon says "per-support").

PROBE: a fifth support with two monomials per charge side,
    P = a Z + b Z^2 + c W + d W^2   (charges 1, 2, -1, -2),
exact moments from the circular-pair functional E[Z^i W^j] = delta_ij 2^i i!:
    mu_m = sum_{i+2j = k+2l, i+j+k+l = m} m!/(i!j!k!l!) a^i b^j c^k d^l
           * 2^{i+2j} (i+2j)!.
Fit the minimal P-recurrence (exact rational linear algebra) at several
coefficient tuples and inspect the leading coefficient q_r(m): does its
structure (degree, roots, integer zeros) move with the COEFFICIENTS?

Also: reproduce the file's H4 fit (validation of the fitter), and check the
DEGENERATE-STRATUM hazard: whether the generic-orderfit recurrence remains
valid/nonzero when coefficients specialize (d = 0: support collapses).

boxeph-referee-2026-07-20-S186r.  Pure python, exact Fractions.
"""
from fractions import Fraction
import math

fact = math.factorial


def mom_ab_cd(a, b, c, d, m):
    """P = aZ + bZ^2 + cW + dW^2, exact E[P^m]."""
    tot = 0
    for i in range(m + 1):
        for j in range(m - i + 1):
            for k in range(m - i - j + 1):
                l = m - i - j - k
                if i + 2 * j != k + 2 * l:
                    continue
                N = i + 2 * j
                tot += (fact(m) // (fact(i) * fact(j) * fact(k) * fact(l))) \
                    * a ** i * b ** j * c ** k * d ** l * 2 ** N * fact(N)
    return tot


def mom_zsw(a, b, c, m):
    """H4 control: P = aZ + bZW + cW (Fractions ok)."""
    tot = Fraction(0)
    for y in range(m + 1):
        rem = m - y
        if rem % 2:
            continue
        x = rem // 2
        coef = Fraction(fact(m), fact(x) * fact(y) * fact(x))
        tot += coef * a ** x * b ** y * c ** x * (Fraction(2) ** (x + y) * fact(x + y))
    return tot


def fit_recurrence(seq, order, deg, nrows_extra=10):
    nunk = (order + 1) * (deg + 1)
    nrows = min(len(seq) - order, nunk + nrows_extra)
    rows = []
    for m in range(nrows):
        row = []
        for k in range(order + 1):
            for dd in range(deg + 1):
                row.append(Fraction(m) ** dd * Fraction(seq[m + k]))
        rows.append(row)
    R = [r[:] for r in rows]
    ncols = nunk
    piv_of_col = [-1] * ncols
    r = 0
    for c in range(ncols):
        piv = None
        for rr in range(r, len(R)):
            if R[rr][c] != 0:
                piv = rr
                break
        if piv is None:
            continue
        R[r], R[piv] = R[piv], R[r]
        pv = R[r][c]
        R[r] = [x / pv for x in R[r]]
        for rr in range(len(R)):
            if rr != r and R[rr][c] != 0:
                f = R[rr][c]
                R[rr] = [x - f * y for x, y in zip(R[rr], R[r])]
        piv_of_col[c] = r
        r += 1
        if r == len(R):
            break
    free = [c for c in range(ncols) if piv_of_col[c] < 0]
    if not free:
        return None
    fc = free[0]
    sol = [Fraction(0)] * ncols
    sol[fc] = Fraction(1)
    for c in range(ncols):
        pr = piv_of_col[c]
        if pr >= 0:
            sol[c] = -R[pr][fc]
    Q = []
    i = 0
    for k in range(order + 1):
        Q.append(sol[i:i + deg + 1])
        i += deg + 1
    return Q


def check_recurrence(seq, Q, deg):
    order = len(Q) - 1
    for m in range(len(seq) - order):
        tot = Fraction(0)
        for k in range(order + 1):
            qk = sum(Q[k][dd] * Fraction(m) ** dd for dd in range(deg + 1))
            tot += qk * Fraction(seq[m + k])
        if tot != 0:
            return False, m
    return True, None


def poly_str(coeffs):
    parts = []
    for dd, cc in enumerate(coeffs):
        if cc != 0:
            parts.append("(%s) m^%d" % (cc, dd))
    return " + ".join(parts) if parts else "0"


def minimal_fit(seq, orders=(2, 3, 4, 5), degs=(0, 1, 2, 3, 4)):
    for order in orders:
        for deg in degs:
            Q = fit_recurrence(seq, order, deg)
            if Q:
                ok, bad = check_recurrence(seq, Q, deg)
                if ok:
                    return order, deg, Q
    return None


print("=" * 78)
print("REFEREE S186r: does the finite-moment-test bound depend only on the SUPPORT?")
print("=" * 78)

# ---- validation: reproduce the file's H4 fit --------------------------------
a, b, c = Fraction(3, 2), Fraction(-2, 3), Fraction(5, 4)
seq = [mom_zsw(a, b, c, m) for m in range(60)]
res = minimal_fit(seq, orders=(2, 3), degs=(1, 2, 3))
if res:
    order, deg, Q = res
    print("\n(V) H4 control P = aZ+bZW+cW @ (3/2,-2/3,5/4): order %d deg %d;"
          " leading q_%d(m) = %s   [file: order 3, deg 2, leading 1]"
          % (order, deg, order, poly_str(Q[order])))

# ---- the fifth support -------------------------------------------------------
tuples = [(1, 1, 1, 1), (2, 1, 1, 1), (1, 3, 1, 2), (5, 1, 7, 1), (1, 1, 1, 0)]
MMAX = 58
print("\n(P) P = aZ + bZ^2 + cW + dW^2, exact moments to m = %d:" % MMAX)
for tup in tuples:
    a, b, c, d = tup
    seq = [mom_ab_cd(a, b, c, d, m) for m in range(MMAX + 1)]
    res = minimal_fit(seq)
    if res is None:
        print("  (a,b,c,d)=%s: NO recurrence at order<=5, deg<=4 "
              "(support-effectivity itself not visible at these bounds)" % (tup,))
        continue
    order, deg, Q = res
    lead = Q[order]
    zs = [m for m in range(0, 400)
          if sum(lead[dd] * Fraction(m) ** dd for dd in range(deg + 1)) == 0]
    print("  (a,b,c,d)=%s: order %d, coeff-deg %d; verified all m <= %d"
          % (tup, order, deg, MMAX - order))
    print("      leading q_%d(m) = %s" % (order, poly_str(lead)))
    print("      integer zeros in [0,400): %s  => M = %d"
          % (zs, order + 1 + (max(zs) if zs else -1)))

print("\nReading: if the leading coefficient's structure (degree/zeros) moves")
print("with (a,b,c,d), the bound M is per-(support, coefficients), refuting the")
print("'per-support' wording.  If it stays constant here too, the claim remains")
print("UNPROVEN (all evidence from constant-leading examples; no uniformity")
print("argument exists in the file, whose own script writes 'M(p)').")
print("DONE.")
