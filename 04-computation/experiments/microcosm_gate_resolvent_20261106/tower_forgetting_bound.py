#!/usr/bin/env python3
# THM-4900(7) checks (mac-mini-2026-11-06-microcosm): exact tower weights U_(e,r) = (L^r q^e)(1) by the first-step
# recursion (R); (b)2 U_(e,r) >= (r+1) U_(e-2,r) for all e, r <= 120; the slack S >= U_(e-2,r); (c) the explicit
# first-order error bound on all ordered pairs of isomorphism classes of orders 4-6 (gentourng) for r <= 60; and the
# audit-H1 example (order-5 classes with P = q^5 + 8q^3 + 2q and q^5 + 8q^3 + 6q: same c3, towers differ at r = 20).
# Exact integer / rational arithmetic. Output: tower_forgetting_bound.out
import itertools, os, subprocess
from fractions import Fraction
from math import comb

GENTOURNG = os.environ.get("GENTOURNG", "/opt/homebrew/bin/gentourng")
RMAX, EMAX = 120, 260
U = [[1] * (EMAX + 1)]
for r in range(1, RMAX + 1):
    prev = U[-1]; row = []
    for e in range(EMAX + 1 - 2 * r):
        v = prev[e + 2]; k = 0
        while 2 * k + 1 <= e:
            v += 2 * comb(e, 2 * k + 1) * prev[e - 2 * k]; k += 1
        row.append(v)
    U.append(row)

# (b)2 and the slack
viol = slack_bad = 0
minratio = []
for r in range(RMAX + 1):
    m = min(Fraction(U[r][e], (r + 1) * U[r][e - 2]) for e in range(2, len(U[r])))
    minratio.append(m)
    viol += sum(1 for e in range(2, len(U[r])) if U[r][e] < (r + 1) * U[r][e - 2])
    if r >= 1:
        for e in range(2, len(U[r])):
            S = sum(2 * (comb(e, 2*k+1) - comb(e-2, 2*k+1)) * U[r-1][e-2*k] for k in range(e // 2 + 1))
            slack_bad += S < U[r][e-2]
print(f"(b)2 U_(e,r) >= (r+1) U_(e-2,r) for r <= {RMAX} and 2 <= e <= {EMAX} - 2r: violations = {viol}")
print(f"     slack S >= U_(e-2,r): violations = {slack_bad}")
print("     min over e of U_(e,r)/((r+1)U_(e-2,r)) at r = 0, 1, 2, 5, 10, 20, 40, 80, 120:",
      ", ".join(f"{float(minratio[r]):.3f}" for r in (0, 1, 2, 5, 10, 20, 40, 80, 120)))

def classes(n):
    out = subprocess.run([GENTOURNG, "-q", str(n)], capture_output=True, text=True, check=True).stdout.split()
    res = []
    for s in out:
        A = [[0] * n for _ in range(n)]; t = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[t] == '1': A[i][j] = 1
                else: A[j][i] = 1
                t += 1
        res.append(A)
    return res

def sw_coeffs(A):
    # P_T(q) = sum_pi C(q + f(pi), n) (Session210 (7)); monomial coefficients by exact Lagrange interpolation
    n = len(A)
    fd = [0] * n
    for p in itertools.permutations(range(n)):
        fd[sum(A[p[i]][p[i + 1]] for i in range(n - 1))] += 1
    xs = list(range(n + 1))
    ys = [sum(c * comb(q + f, n) for f, c in enumerate(fd)) for q in xs]
    coeffs = [Fraction(0)] * (n + 1)
    for i, xi in enumerate(xs):
        basis = [Fraction(1)]; denom = 1
        for j, xj in enumerate(xs):
            if j == i: continue
            basis = [Fraction(0)] + basis
            for t in range(len(basis) - 1):
                basis[t] -= xj * basis[t + 1]
            denom *= (xi - xj)
        for t in range(n + 1):
            coeffs[t] += Fraction(ys[i], denom) * basis[t]
    assert all(c.denominator == 1 for c in coeffs)
    return [int(coeffs[n - 2 * j]) for j in range(n // 2 + 1)]

def c3(A):
    return sum(1 for a, b, c in itertools.combinations(range(len(A)), 3)
               if (A[a][b] and A[b][c] and A[c][a]) or (A[b][a] and A[c][b] and A[a][c]))

checked = bad = 0; worst = Fraction(0)
example = {}
for n in (4, 5, 6):
    data = []
    for A in classes(n):
        sw = sw_coeffs(A); data.append((sw, c3(A)))
        assert sw[0] == 1 and sw[1] == 2 * c3(A)
    for (swT, cT), (swS, cS) in itertools.product(data, data):
        sT = sum(abs(x) for x in swT[1:]); sS = sum(abs(x) for x in swS[1:])
        for r in range(1, 61):
            if r + 1 <= sS: continue
            HT = sum(swT[j] * U[r][n - 2 * j] for j in range(len(swT)))
            HS = sum(swS[j] * U[r][n - 2 * j] for j in range(len(swS)))
            rho = Fraction(U[r][n - 2], U[r][n])
            E = Fraction(HT, HS) - 1 - 2 * (cT - cS) * rho
            B = rho * (sT + sS + 2 * abs(cT - cS) * sS) / (r + 1 - sS)
            checked += 1; bad += abs(E) > B
            if B: worst = max(worst, abs(E) / B)
    print(f"order {n}: {len(data)} classes; all sw_j >= 0: {all(min(sw) >= 0 for sw, _ in data)}")
    if n == 5:
        example = {tuple(sw): cc for sw, cc in data}
print(f"(c) explicit error bound on {checked} (T, T', r) triples (orders 4-6, r <= 60): violations = {bad}; "
      f"max |E|/bound = {float(worst):.4f}")

a, b = (1, 8, 2), (1, 8, 6)
assert a in example and b in example and example[a] == example[b] == 4
r = 20
Ha = sum(a[j] * U[r][5 - 2 * j] for j in range(3)); Hb = sum(b[j] * U[r][5 - 2 * j] for j in range(3))
print(f"audit-H1 example: P = q^5+8q^3+2q vs q^5+8q^3+6q (c3 = 4 both): H ratio - 1 at r = 20 = {float(Fraction(Ha, Hb) - 1):.3e}")
print("all done")
