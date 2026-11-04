"""aud_gf_check.py - independent audit (opus-2026-11-04-S75 audit) of the generating-function algebra of THM-4869 (C4),
the core-count formulas, the lineage (orphan) polynomials of the note's section 3.1, and the census values.

Exact integer power-series arithmetic only (no sympy). Census values are the audit's own (aud_census_*.out);
klein's order-5..10 values are from HYP-9459 / two_branch_rearming_sigma_sweep_20261104/README.md section 8.
"""
N = 40


def fib(k):
    if k < 0:
        # F(-k) = (-1)^(k+1) F(k)
        return (-1) ** (k + 1) * fib(-k) if k != 0 else 0
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def luc(k):
    return fib(k - 1) + fib(k + 1)


def mul(a, b):
    c = [0] * N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < N:
                    c[i + j] += x * y
    return c


def inv(a):
    # inverse of a power series with a[0] = +-1
    assert a[0] in (1, -1)
    b = [0] * N
    b[0] = a[0]
    for n in range(1, N):
        s = sum(a[k] * b[n - k] for k in range(1, n + 1) if k < len(a))
        b[n] = -s * a[0]
    return b


def poly(coeffs):
    p = [0] * N
    for e, c in coeffs.items():
        p[e] += c
    return p


den = poly({0: 1, 1: -3, 2: 1})
inv_den = inv(den)

# rooted prefixes: R(2) = 1, R(m) = F(2m - 5) (m >= 3)
R_series = poly({m: (1 if m == 2 else fib(2 * m - 5)) for m in range(2, N)})
R_claim = mul(poly({2: 1, 3: -2}), inv_den)
assert R_series == R_claim, "R(x) closed form"
print("R(x) = x^2 (1 - 2x)/(1 - 3x + x^2): OK (R(2) = 1, R(m) = F(2m-5))")


def c(k, d):
    if k == 3:
        return 2 if d == 5 else (3 * 2 ** (d - 6) if d >= 6 else 0)
    if k == 4:
        return 2 if d == 5 else (2 ** (d - 5) if d >= 6 else 0)
    if k == 5:
        return 4 if d == 6 else (5 * 2 ** (d - 7) if d >= 7 else 0)


C = {k: poly({d: c(k, d) for d in range(N)}) for k in (3, 4, 5)}
inv_1m2x = inv(poly({0: 1, 1: -2}))
C_claim = {3: mul(poly({5: 2, 6: -1}), inv_1m2x), 4: mul(poly({5: 2, 6: -2}), inv_1m2x), 5: mul(poly({6: 4, 7: -3}), inv_1m2x)}
for k in (3, 4, 5):
    assert C[k] == C_claim[k], k
print("C_3 = x^5(2-x)/(1-2x), C_4 = 2x^5(1-x)/(1-2x), C_5 = x^6(4-3x)/(1-2x): OK")

# U_k = C_k R / x, i.e. u_k(n) = sum_d c_k(d) R(n - d + 1)
U = {}
for k in (3, 4, 5):
    conv = [0] * N
    for n in range(N):
        conv[n] = sum(c(k, d) * (1 if n - d + 1 == 2 else (fib(2 * (n - d + 1) - 5) if n - d + 1 >= 3 else 0)) for d in range(n + 1))
    U[k] = conv
U_claim = {3: mul(poly({6: 2, 7: -1}), inv_den), 4: mul(poly({6: 2, 7: -2}), inv_den), 5: mul(poly({7: 4, 8: -3}), inv_den)}
for k in (3, 4, 5):
    assert U[k] == U_claim[k], k
print("U_3 = x^6(2-x)/den, U_4 = 2x^6(1-x)/den, U_5 = x^7(4-3x)/den: OK (as convolutions)")
for n in range(6, N):
    assert U[3][n] == fib(2 * n - 9), n
for n in range(7, N):
    assert U[4][n] == 2 * fib(2 * n - 11), n
    assert U[5][n] == 4 * fib(2 * n - 12) - 3 * fib(2 * n - 14) == luc(2 * n - 12) + fib(2 * n - 13), n
print("u_3(n) = F(2n-9) (n>=6), u_4(n) = 2F(2n-11) (n>=7), u_5(n) = 4F(2n-12) - 3F(2n-14) = L(2n-12) + F(2n-13) (n>=7): OK to n = %d" % (N - 1))
print("   at n = 6 the convolution gives u_4 = %d, u_5 = %d (census: 3 and 4, the extra being 1 and 4 2-strong classes)" % (U[4][6], U[5][6]))

# census values (audit's own runs; orders 5-10 exhaustive, 11-14 Moon extension) against the laws
census = {  # n: (u_3, u_4, u_5)
    6: (2, 3, 4), 7: (5, 4, 4), 8: (13, 10, 9), 9: (34, 26, 23), 10: (89, 68, 60),
    11: (233, 178, 157), 12: (610, 466, 411), 13: (1597, 1220, 1076), 14: (4181, 3194, 2817)}
for n, (a, b, cc) in census.items():
    two_strong = {6: (0, 1, 4)}.get(n, (0, 0, 0))
    assert a == U[3][n] + two_strong[0] and b == U[4][n] + two_strong[1] and cc == U[5][n] + two_strong[2], n
print("census u_3, u_4, u_5 at n = 6..14 = convolution + 2-strong classes (only at n = 6): OK")

# lineage: orphan cores by d+_D(u), from the explicit chain-plus-block description
def orphan_poly(k, d):
    """orphans = cores with p >= 2 and w_2 -> u, or p = 1 with A = W, or p = 0; returns {out-degree of u in D: count}"""
    blocks = {3: [(3, 3, 0), (4, 2, 2)], 4: [(4, 3, 1), (4, 3, 1)], 5: [(4, 4, 0), (5, 2, 3), (5, 2, 3), (5, 2, 3)]}
    # (block order w, |A|, |W - A|); u beats W - A, u -> w_1, free arcs u - w_3..w_p for orphans
    from math import comb
    res = {}
    for (w, a, wa) in blocks[k]:
        p = d - 1 - w
        if p >= 2:
            free = p - 2
            for j in range(free + 1):
                deg = 1 + wa + j
                res[deg] = res.get(deg, 0) + comb(free, j)
        elif p == 1 and a == w:
            deg = 1 + wa
            res[deg] = res.get(deg, 0) + 1
        elif p == 0 and a != w:
            res[wa] = res.get(wa, 0) + 1
    return res


def poly_from(coeffs_low, shift):
    return {shift + i: x for i, x in enumerate(coeffs_low) if x}


from math import comb
lineage_out = {  # from the session's fewcycle_core_lineage_20261104.out, orphans by d+_D(u), cores of order d = n - 1
    (3, 12): {1: 1, 2: 6, 3: 16, 4: 25, 5: 25, 6: 16, 7: 6, 8: 1},
    (4, 12): {2: 2, 3: 10, 4: 20, 5: 20, 6: 10, 7: 2},
    (4, 9): {2: 2, 3: 4, 4: 2},
    (3, 9): {1: 1, 2: 3, 3: 4, 4: 3, 5: 1},
    (5, 12): {1: 1, 2: 5, 3: 10, 4: 13, 5: 17, 6: 19, 7: 12, 8: 3},
}
for (k, d), obs in lineage_out.items():
    assert orphan_poly(k, d) == obs, (k, d, orphan_poly(k, d), obs)
print("orphan distributions recomputed from Lemma B4 match the session's lineage output at (k,d) =", sorted(lineage_out))
for d in range(7, 16):
    # note's k = 3 polynomial z(1+z)^(d-7)(1+z+z^2)
    k3 = {}
    for i in range(d - 7 + 1):
        for j, cf in enumerate((1, 1, 1)):
            k3[1 + i + j] = k3.get(1 + i + j, 0) + comb(d - 7, i) * cf
    assert k3 == orphan_poly(3, d), d
    k4_note = {2 + i: 2 * comb(d - 6, i) for i in range(d - 6 + 1)}
    k4_fix = {2 + i: 2 * comb(d - 7, i) for i in range(d - 7 + 1)}
    assert k4_fix == orphan_poly(4, d), d
    if d == 12:
        print("k=4, d=12: note's 2z^2(1+z)^(d-6) gives", [k4_note[x] for x in sorted(k4_note)],
              "; Lemma B4 and the lineage output give", [k4_fix[x] for x in sorted(k4_fix)], "= 2z^2(1+z)^(d-7)")
print("k = 3 orphan polynomial z(1+z)^(d-7)(1+z+z^2): OK for d = 7..15; k = 4: correct exponent is d - 7, not d - 6")

# attached share: exactly half only once every core has p >= 2
for k in (3, 4, 5):
    row = []
    for d in range(5, 10):
        tot = c(k, d)
        orph = sum(orphan_poly(k, d).values())
        row.append("d=%d: %d attached of %d" % (d, tot - orph, tot))
    print("k=%d:" % k, "; ".join(row))
