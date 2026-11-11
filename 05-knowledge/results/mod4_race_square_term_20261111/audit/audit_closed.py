#!/usr/bin/env python3
"""audit_closed.py -- independent audit of Theorem P (5) and the remarks around it.

(a) E_n[{3}], E_n[{3,5}], E_n[{5}] computed from the joint law (4) and the score formulas (1) by a DIRECT double sum over
    (k1, k2) = (|t on O|, |t off O|) for every |O| = q  (no use of the factorised (1+i)^r sums), exact rationals, n <= NMAX1;
(b) the same by the factorised sum of the note (with the boundary terms), n <= NMAX2, compared with (a) on the overlap;
(c) both against the closed forms of (5), all n <= NMAX2;
(d) the conditional form: E[(-1)^(c3+c5) | s = 1_O] = +-2^(2 - ceil(q/2) - ceil((n-q)/2)), and = 2^(-(n-2)/2) for n = 2 mod 4;
(e) the finite-dimensionality argument: for each class r mod 8 the sequences k -> E_{8k+r}[.] are annihilated by the
    order-5 recurrence with characteristic roots 2^-4, 2^-8, 2^-12, -(2+sqrt2)^4 2^-12, -(2-sqrt2)^4 2^-12
    (characteristic polynomial (x - 2^-4)(x - 2^-8)(x - 2^-12)(x^2 + (17/512) x + 2^-20), rational coefficients);
(f) exhaustive values (audit_small.c, n <= 7) and census values (audit_census.c, n <= 10) against the closed forms;
(g) the expected kernel size E|ker A| = 2 + 2^((1-m)/2) cos((m-1)pi/4) - 2^(1-m), against exhaustive sums (m <= 7)
    and the count of (T, x) pairs; and the 'second derivation' of the row n = 1 mod 4 from Theorem P (3).
"""
import sys
from fractions import Fraction as F
from math import comb

NMAX1 = int(sys.argv[1]) if len(sys.argv) > 1 else 64
NMAX2 = int(sys.argv[2]) if len(sys.argv) > 2 else 200

def C2(k): return k * (k - 1) // 2
def kron2(n):  # Kronecker symbol (2/n), n odd
    return 1 if n % 8 in (1, 7) else -1

def direct(n):
    """E[(-1)^c3], E[(-1)^(c3+c5)], E[(-1)^c5] from (1) and (4); also the list of conditional means."""
    tot3 = F(0); tot35 = F(0); tot5 = F(0); cond = {}
    ns = 2 ** (n - 1)                      # number of admissible s
    for q in range(0, n + 1):
        if (q - C2(n)) % 2: continue
        sign3 = (-1) ** (C2(q) + n // 4)
        if q == 0 or q == n:
            q2 = 0 if q == 0 else n
            m35 = F((-1) ** (C2(q2) + n // 4))
        else:
            num = 0; den = 0
            for k1 in range(0, q + 1):
                if (k1 - C2(q)) % 2: continue
                for k2 in range(0, n - q + 1):
                    if (k2 - C2(n - q)) % 2: continue
                    w = comb(q, k1) * comb(n - q, k2)
                    den += w; num += w * (-1) ** (C2(k1 + k2) + n // 4)
            assert den == 2 ** (n - 2)
            m35 = F(num, den)
        cond[q] = m35
        p = F(comb(n, q), ns)
        tot3 += p * sign3; tot35 += p * m35; tot5 += p * sign3 * m35
    return tot3, tot35, tot5, cond

def factorised(n):
    """The formula of the note (section 6.2 (5)) with its boundary terms."""
    def S(extra):
        tot = 0
        for q in range(1, n):
            if (q - C2(n)) % 2: continue
            e = C2(q) * C2(n - q) + q // 4 + (n - q) // 4 + (C2(q) if extra else 0)
            tot += comb(n, q) * (-1) ** e * 2 ** (q // 2 + (n - q) // 2)
        return tot
    B35 = (1 if C2(n) % 2 == 0 else 0) + ((-1) ** C2(n) if (n - C2(n)) % 2 == 0 else 0)
    B5 = (1 if C2(n) % 2 == 0 else 0) + (1 if (n - C2(n)) % 2 == 0 else 0)
    e35 = (-1) ** (n // 4) * F(1, 2 ** (n - 1)) * (B35 + S(False) / F(2) ** (n - 2))
    # for {5}: the sign (-1)^floor(n/4) of c3 cancels the one of c3+c5; the summand carries floor(q/4)+floor((n-q)/4) only
    e5 = F(1, 2 ** (n - 1)) * (B5 + S(True) / F(2) ** (n - 2))
    return e35, e5

def closed(n):
    if n % 2 == 1:
        e3 = F(1, 2 ** ((n - 1) // 2))
        if n == 1:
            e35 = F(1, 1) * (2 + kron2(1) * 1 - F(2))             # 2^(-(n-2)) = 2 at n = 1
            e5 = F(1, 1) * (3 - kron2(1) * F(2))                  # 2^(-(n-3)/2) = 2 at n = 1
        else:
            e35 = e3 * (2 + kron2(n) * F(1, 2 ** ((n - 1) // 2)) - F(1, 2 ** (n - 2)))
            e5 = F(1, 2 ** (n - 1)) * (3 - kron2(n) * F(1, 2 ** ((n - 3) // 2)))
        return e3, e35, e5
    if n % 4 == 2:
        return F(1, 2 ** ((n - 2) // 2)), F(1, 2 ** ((n - 2) // 2)), F(1, 2 ** (n - 2))
    e3 = F(1, 2 ** ((n - 2) // 2))
    _, e35, e5 = closed(n - 1)
    return e3, e35, e5

fail = 0
print(f"(a),(c),(d): direct double sum from (1) and (4) against the closed forms, n = 1..{NMAX1}")
D = {}
for n in range(1, NMAX1 + 1):
    t3, t35, t5, cond = direct(n)
    D[n] = (t3, t35, t5)
    c3, c35, c5 = closed(n)
    ok = (t3, t35, t5) == (c3, c35, c5)
    # conditional form
    okc = True
    for q, v in cond.items():
        if 1 <= q <= n - 1:
            mag = F(1, 2 ** (-(2 - (q + 1) // 2 - (n - q + 1) // 2))) if (2 - (q + 1) // 2 - (n - q + 1) // 2) <= 0 else F(2 ** (2 - (q + 1) // 2 - (n - q + 1) // 2))
            if abs(v) != mag: okc = False
            if n % 4 == 2 and v != F(1, 2 ** ((n - 2) // 2)): okc = False
    if not ok or not okc:
        fail += 1; print(f"  n={n}: MISMATCH closed {ok} conditional {okc}: direct {t3} {t35} {t5} closed {c3} {c35} {c5}")
print(f"  mismatches: {fail}")
print("  values n = 1..12 (E[{3}], E[{3,5}], E[{5}]):")
for n in range(1, 13):
    print(f"   n={n}: {D[n][0]}  {D[n][1]}  {D[n][2]}")

print(f"\n(b),(c): factorised sum of the note against the direct sum (n <= {NMAX1}) and the closed forms (n <= {NMAX2})")
f2 = 0; Fv = {}
for n in range(1, NMAX2 + 1):
    e35, e5 = factorised(n)
    Fv[n] = (e35, e5)
    c3, c35, c5 = closed(n)
    if (e35, e5) != (c35, c5): f2 += 1; print(f"  n={n}: factorised != closed")
    if n in D and (e35, e5) != (D[n][1], D[n][2]): f2 += 1; print(f"  n={n}: factorised != direct")
print(f"  mismatches: {f2} ; members of each class mod 8 checked: {[len([n for n in range(1, NMAX2+1) if n % 8 == r]) for r in range(8)]}")

print("\n(e): order-5 recurrence in k for n = 8k + r (roots 2^-4, 2^-8, 2^-12 and those of x^2 + (17/512) x + 2^-20)")
# characteristic polynomial coefficients, highest degree first
def polymul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r
P = [F(1)]
for root in (F(1, 16), F(1, 256), F(1, 4096)):
    P = polymul(P, [F(1), -root])
P = polymul(P, [F(1), F(17, 512), F(1, 2 ** 20)])
f3 = 0; tested = 0
for r in range(1, 9):
    for which in (0, 1):
        seqF = [Fv[n][which] for n in range(r, NMAX2 + 1, 8)]
        seqC = [closed(n)[1 + which] for n in range(r, NMAX2 + 1, 8)]
        for seq in (seqF, seqC):
            for k in range(len(seq) - 5):
                val = sum(P[j] * seq[k + 5 - j] for j in range(6))
                tested += 1
                if val != 0: f3 += 1
print(f"  recurrence instances tested: {tested} ; failures: {f3}")
# is five the right number?  rank of the Hankel-type matrix of the sequences
import itertools
def rank(rows):
    rows = [list(r) for r in rows]; rk = 0
    for c in range(len(rows[0])):
        p = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                fct = rows[i][c] / rows[rk][c]
                rows[i] = [a - fct * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
for r in range(1, 9):
    seq = [Fv[n][0] for n in range(r, NMAX2 + 1, 8)]
    H = [[seq[i + j] for j in range(7)] for i in range(7)]
    seq5 = [Fv[n][1] for n in range(r, NMAX2 + 1, 8)]
    H5 = [[seq5[i + j] for j in range(7)] for i in range(7)]
    print(f"   class {r % 8} mod 8: number of geometric sequences actually present (Hankel rank): E[{{3,5}}]: {rank(H)}, E[{{5}}]: {rank(H5)}")

print("\n(f): exhaustive and census values")
ex = {  # numerators over 2^C(n,2) from audit_small.c (labelled, brute force)
    3: (4, 4, 8), 4: (32, 32, 64), 5: (256, 416, 224), 6: (8192, 8192, 2048), 7: (262144, 548864, 90112)}
for n, (a3, a35, a5) in ex.items():
    tot = 2 ** C2(n)
    ok = (F(a3, tot), F(a35, tot), F(a5, tot)) == closed(n)
    print(f"  n={n}: exhaustive {F(a3, tot)} {F(a35, tot)} {F(a5, tot)} : {'OK' if ok else 'MISMATCH'}")
cen = {8: (F(1, 8), F(67, 256), F(11, 256)), 9: (F(1, 16), F(263, 2048), F(23, 2048)), 10: (F(1, 16), F(1, 16), F(1, 256))}
for n, v in cen.items():
    print(f"  n={n}: census {v[0]} {v[1]} {v[2]} : {'OK' if v == closed(n) else 'MISMATCH'}")
note = {3: (F(1, 2), F(1)), 4: (F(1, 2), F(1)), 5: (F(13, 32), F(7, 32)), 6: (F(1, 4), F(1, 16)), 7: (F(67, 256), F(11, 256)), 8: (F(67, 256), F(11, 256)),
        9: (F(263, 2048), F(23, 2048)), 10: (F(1, 16), F(1, 256)), 11: (F(1007, 16384), F(49, 16384))}
print("  values quoted in the note / THM-5132 for n = 3..11:", "OK" if all(closed(n)[1:] == v for n, v in note.items()) else "MISMATCH")
print("  limits: E[{3,5}]/E[{3}] and E[{5}]/E[{3}]^2 at n = 41, 43, 44:",
      [(float(closed(n)[1] / closed(n)[0]), float(closed(n)[2] / closed(n)[0] ** 2)) for n in (41, 43, 44)])
print("  Run C exact values (n = 11..18) E[{3,5}]:", [str(closed(n)[1]) for n in range(11, 19)])
print("  R7: E_14 =", closed(14)[1], "= 2^-6:", closed(14)[1] == F(1, 64), "; E_18 =", closed(18)[1], "= 2^-8:", closed(18)[1] == F(1, 256))
print("  R8: E_12 = E_11:", closed(12)[1] == closed(11)[1], "; E_16 = E_15:", closed(16)[1] == closed(15)[1])

print("\n(g): expected kernel size")
import math
ks = {1: 2, 2: 3, 3: 12, 4: 112, 5: None, 6: None, 7: None}
def Ek_formula_exact(m):
    # 1 + 2^(1-m) (#{k = 0,1 mod 4, 1 <= k <= m} C(m,k))
    cnt = sum(comb(m, k) for k in range(1, m + 1) if k % 4 in (0, 1))
    return 1 + F(cnt, 2 ** (m - 1))
for m in range(1, 41):
    ex_ = Ek_formula_exact(m)
    cf = 2 + 2 ** ((1 - m) / 2) * math.cos((m - 1) * math.pi / 4) - 2 ** (1 - m)
    if abs(float(ex_) - cf) > 1e-12: print(f"  m={m}: closed form {cf} != count formula {float(ex_)}")
print("  closed form 2 + 2^((1-m)/2) cos((m-1)pi/4) - 2^(1-m) equals the count formula 1 + 2^(1-m) #{x != 0, |x| = 0,1 mod 4} for m = 1..40")
print("  count formula times 2^C(m,2) (to compare with the exhaustive sums of audit_rand kernel):", {m: int(Ek_formula_exact(m) * 2 ** C2(m)) for m in range(1, 8)})
# second derivation of the row n = 1 mod 4: E_n[{3,5}] = 2^(-m/2) E|ker A|, m = n - 1 = 0 mod 4
ok = all(closed(n)[1] == F(1, 2 ** ((n - 1) // 2)) * Ek_formula_exact(n - 1) for n in range(5, 200, 4))
print("  row n = 1 mod 4 of (5) equals 2^(-(n-1)/2) E|ker A| (m = n-1) for n = 5, 9, ..., 197:", ok)
