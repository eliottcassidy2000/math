#!/usr/bin/env python3
"""
holonomic_moments_boxeph_S186.py  (HYP-8580, THM-1765 section 3)

HOLONOMIC MOMENTS: Ghat(s,t) is algebraic, so A(t) = int e^{-s} Ghat ds is
holonomic and the moment sequence mu_m = E[P^m] is P-RECURSIVE: there exist
polynomials q_0..q_r (effective in the support) with
    sum_k q_k(m) mu_{m+k} = 0  for all m.
The Gamma-graded ladder = Birkhoff-Trjitzinsky asymptotics of P-recursive
sequences; nullcone membership per support = an EXPLICIT FINITE MOMENT TEST
(vanishing up to order + last-integer-zero of the leading coefficient
propagates forever). m-side effectivity, complementing opus THM-1740's
coefficient-side Groebner decidability.

CHECKS (all moments computed EXACTLY from the circular-pair functional
E[Z^i W^j] = delta_ij 2^i i!):
(H1) P = aZ + b + cW: EGF e^{bt + 2ac t^2}: mu_m = b mu_{m-1} +
     4ac (m-1) mu_{m-2}. Verify to m = 60. FINITE TEST: nullcone iff
     mu_1 = mu_2 = 0 iff b = 0 and ac = 0 iff ONE-SIDED (2-line closure of
     the linear span by recurrence alone).
(H2) P = aZ + b W^2 + c (charges 1,-2,0): EGF e^{ct + 4a^2 b t^3}:
     mu_m = c mu_{m-1} + 12 a^2 b (m-1)(m-2) mu_{m-3}. Verify. FINITE TEST:
     nullcone iff mu_1 = mu_2 = mu_3 = 0 iff c = 0 and a^2 b = 0: one-sided.
(H3) P = aZ^2 + b + cW^2 (charges 2,0,-2): mu_m = sum_i m!/(i!^2 (m-2i)!)
     (ac)^i b^{m-2i} 4^i (2i)!: FIT a P-recurrence (order 2, coeff degree
     <= 2) by exact rational linear algebra from m <= 20; verify m <= 60.
(H4) P = aZ + b ZW + cW (linear charge-0 in s, no constant): exact moments
     via E[Z^{i+j} W^{i+j}]; FIT order <= 3, degree <= 3; verify m <= 60;
     report the leading coefficient q_r(m) and its integer zeros (the
     finite-test bound M(p) in action).

boxeph-2026-07-20-S186. Pure python, exact Fraction arithmetic.
"""
from fractions import Fraction
import math

def fact(n):
    return math.factorial(n)


# ---------- exact moment generators ----------

def mom_linear(a, b, c, m):
    # P = aZ + b + cW
    tot = Fraction(0)
    for k in range(m // 2 + 1):
        tot += (Fraction(fact(m), fact(k) * fact(m - 2 * k)) *
                (Fraction(2) * a * c) ** k * b ** (m - 2 * k))
    return tot


def mom_zww(a, b, c, m):
    # P = aZ + bW^2 + c: i = #Z, j = #W^2, k = #c: i = 2j; i+j+k = m
    tot = Fraction(0)
    j = 0
    while 3 * j <= m:
        i = 2 * j
        k = m - 3 * j
        coef = Fraction(fact(m), fact(i) * fact(j) * fact(k))
        tot += coef * a ** i * b ** j * c ** k * (Fraction(2) ** (2 * j) * fact(2 * j))
        j += 1
    return tot


def mom_z2w2(a, b, c, m):
    # P = aZ^2 + b + cW^2
    tot = Fraction(0)
    for i in range(m // 2 + 1):
        if 2 * i > m:
            break
        coef = Fraction(fact(m), fact(i) * fact(i) * fact(m - 2 * i))
        tot += coef * (a * c) ** i * b ** (m - 2 * i) * (Fraction(4) ** i * fact(2 * i))
    return tot


def mom_zsw(a, b, c, m):
    # P = aZ + b ZW + cW: #Z = i + j... monomials: Z^x (ZW)^y W^z: charge x - z = 0
    # x + y + z = m: E = 2^{x+y} (x+y)! with x = z
    tot = Fraction(0)
    for y in range(m + 1):
        rem = m - y
        if rem % 2:
            continue
        x = rem // 2
        coef = Fraction(fact(m), fact(x) * fact(y) * fact(x))
        tot += coef * a ** x * b ** y * c ** x * (Fraction(2) ** (x + y) * fact(x + y))
    return tot


# ---------- exact recurrence fitting ----------

def fit_recurrence(seq, order, deg):
    """find q_k(m) = sum_d q_{k,d} m^d, k=0..order, not all zero, with
    sum_k q_k(m) seq[m+k] = 0 for all m in the fitting window; exact."""
    nunk = (order + 1) * (deg + 1)
    rows = []
    M = len(seq) - order - 1
    for m in range(0, min(M, nunk + 8)):
        row = []
        for k in range(order + 1):
            for d in range(deg + 1):
                row.append(Fraction(m) ** d * seq[m + k])
        rows.append(row)
    # find nullspace vector by exact Gaussian elimination
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
    # repackage
    Q = []
    i = 0
    for k in range(order + 1):
        Q.append(sol[i:i + deg + 1])
        i += deg + 1
    return Q


def check_recurrence(seq, Q, deg, start, end):
    order = len(Q) - 1
    for m in range(start, end - order):
        tot = Fraction(0)
        for k in range(order + 1):
            qk = sum(Q[k][d] * Fraction(m) ** d for d in range(deg + 1))
            tot += qk * seq[m + k]
        if tot != 0:
            return False, m
    return True, None


print("=" * 78)
print("HOLONOMIC MOMENTS: P-recursive E[P^m] + the finite moment test")
print("=" * 78)

a, b, c = Fraction(3, 2), Fraction(-2, 3), Fraction(5, 4)

print("\n(H1) P = aZ + b + cW: predicted mu_m = b mu_{m-1} + 4ac(m-1) mu_{m-2}:")
seq = [mom_linear(a, b, c, m) for m in range(61)]
ok = all(seq[m] == b * seq[m - 1] + 4 * a * c * (m - 1) * seq[m - 2]
         for m in range(2, 61))
print("     verified m=2..60: %s" % ok)
print("     FINITE TEST: mu_1 = b, mu_2 = b^2 + 4ac: nullcone <=> b=0 & ac=0")
print("     <=> one-sided. Order-2, leading coeff 1: 2 moments suffice. QED-toy.")

print("\n(H2) P = aZ + bW^2 + c: predicted mu_m = c mu_{m-1} + 12a^2 b (m-1)(m-2) mu_{m-3}:")
seq = [mom_zww(a, b, c, m) for m in range(61)]
ok = all(seq[m] == c * seq[m - 1] + 12 * a * a * b * (m - 1) * (m - 2) * seq[m - 3]
         for m in range(3, 61))
print("     verified m=3..60: %s" % ok)
print("     FINITE TEST: 3 moments suffice: c = 0 & a^2 b = 0 <=> one-sided.")

print("\n(H3) P = aZ^2 + b + cW^2: FIT search order<=4, degree<=4:")
seq = [mom_z2w2(a, b, c, m) for m in range(75)]
found3 = None
for order in (2, 3, 4):
    for deg in (2, 3, 4):
        Q = fit_recurrence(seq[:order + (deg + 1) * (order + 1) + 10], order, deg)
        if Q:
            okk, bad = check_recurrence(seq, Q, deg, 0, 75)
            if okk:
                found3 = (order, deg, Q)
                break
    if found3:
        break
if found3:
    order, deg, Q = found3
    print("     recurrence: order %d, coeff degree %d; verified m<=%d" %
          (order, deg, 75 - order - 1))
    lead = Q[order]
    print("     leading q_%d(m) = %s" % (order, " + ".join(
        "(%s) m^%d" % (lead[d], d) for d in range(deg + 1) if lead[d] != 0)))
    zs = [m for m in range(0, 200)
          if sum(lead[d] * Fraction(m) ** d for d in range(deg + 1)) == 0]
    print("     integer zeros of leading coeff in [0,200): %s -> finite test M = %d" %
          (zs, order + (max(zs) if zs else -1) + 1))
else:
    print("     no recurrence at order<=4, deg<=4 (raise bounds)")

print("\n(H4) P = aZ + b ZW + cW: FIT order <=3, degree <=3 from m<=28:")
seq = [mom_zsw(a, b, c, m) for m in range(75)]
found = None
for order in (2, 3):
    for deg in (1, 2, 3):
        Q = fit_recurrence(seq[:order + deg * 4 + 14], order, deg)
        if Q:
            okk, bad = check_recurrence(seq, Q, deg, 0, 75)
            if okk:
                found = (order, deg, Q)
                break
    if found:
        break
if found:
    order, deg, Q = found
    print("     recurrence found: order %d, coeff degree %d; verified m<=%d" %
          (order, deg, 75 - order - 1))
    lead = Q[order]
    print("     leading q_%d(m) = %s" % (order, " + ".join(
        "%s m^%d" % (lead[d], d) for d in range(deg + 1) if lead[d] != 0)))
    # integer zeros of leading coeff (the finite-test bound in action)
    zs = [m for m in range(0, 200)
          if sum(lead[d] * Fraction(m) ** d for d in range(deg + 1)) == 0]
    print("     integer zeros of leading coeff in [0,200): %s" % zs)
    print("     => finite moment test bound M = order + max(zeros, -1) + 1 = %d" %
          (order + (max(zs) if zs else -1) + 1))
else:
    print("     no recurrence found at order<=3, deg<=3 (raise bounds)")

print("\nDONE.")
