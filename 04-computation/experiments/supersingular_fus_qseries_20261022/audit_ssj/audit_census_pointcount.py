#!/usr/bin/env python3
"""Independent audit of section 2 (supersingular census) by brute-force point counting (audit_ssj).

No Hasse invariant, no Eisenstein series, no PARI:
  * n_p   = #{j in F_p : a_p(E_j) = 0}, a_p = -sum_x chi(x^3 + A x + B) over F_p (numpy), all 5 <= p < 400;
  * total = #{j in F_{p^2} : #E_j(F_{p^2}) == 1 (mod p)} by brute-force point counting over F_{p^2}
            for 5 <= p <= PMAX2 (numpy; chi on F_{p^2} via chi_p(norm));
  * class numbers h(D) by counting reduced primitive forms (own code);
  * genus of X_0(p) by the classical formula (own code), number of conjugate pairs vs
    (2g + 2 - h(-4p) - [p=3 mod 4] h(-p))/4 and vs Ogg's list.
The pair counts for larger p are done separately in audit_census_ffp2.gp (PARI ellcard over F_{p^2}).
Run: python3 -B audit_census_pointcount.py > audit_census_pointcount.out
"""
import sys
import numpy as np
from sympy import primerange, legendre_symbol

PMAX1 = 400
PMAX2 = 101
bad = 0


def fail(msg):
    global bad
    bad += 1
    print("FAIL:", msg)


def chi_table(p):
    t = np.full(p, -1, dtype=np.int64)
    t[0] = 0
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    t[sq] = 1
    return t


def model(j, p):
    """(A, B) in F_p with j(y^2 = x^3 + A x + B) = j"""
    j %= p
    if j == 0:
        return 0, 1
    if j == 1728 % p:
        return 1, 0
    k = j * pow((1728 - j) % p, -1, p) % p
    return 3 * k % p, 2 * k % p


def jinv(A, B, p):
    num = 1728 * 4 * pow(A, 3, p) % p
    den = (4 * pow(A, 3, p) + 27 * B * B) % p
    return num * pow(den, -1, p) % p


def count_Fp(p):
    chi = chi_table(p)
    x = np.arange(p, dtype=np.int64)
    x3 = x ** 3 % p
    ss = []
    for j in range(p):
        A, B = model(j, p)
        assert jinv(A, B, p) == j
        ap = -int(chi[(x3 + A * x + B) % p].sum())
        if ap % p == 0:
            ss.append(j)
    return ss


def count_Fp2_pairs(p):
    """number of supersingular j in F_{p^2} \\ F_p, by point counting over F_{p^2}"""
    chi = chi_table(p)
    nu = next(a for a in range(2, p) if legendre_symbol(a, p) == -1)
    # F_{p^2} = F_p[s]/(s^2 - nu); element u + v s; all elements as arrays
    U, V = np.meshgrid(np.arange(p, dtype=np.int64), np.arange(p, dtype=np.int64), indexing="ij")
    U = U.ravel(); V = V.ravel()

    def mul(a, b):
        return ((a[0] * b[0] + nu * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    X = (U, V)
    X2 = mul(X, X)
    X3 = mul(X2, X)

    def inv_scalar(a):
        n = (a[0] * a[0] - nu * a[1] * a[1]) % p
        ni = pow(int(n), -1, p)
        return (a[0] * ni % p, (-a[1]) * ni % p)

    def smul(c, a):  # scalar c in F_{p^2} times array a
        return ((c[0] * a[0] + nu * c[1] * a[1]) % p, (c[0] * a[1] + c[1] * a[0]) % p)

    nss = 0
    q = p * p
    for ju in range(p):
        for jv in range(1, p):          # j not in F_p (j = 0, 1728 are in F_p)
            j = (ju, jv)
            den = ((1728 - ju) % p, (-jv) % p)
            k = (j[0] * inv_scalar(den)[0] + nu * j[1] * inv_scalar(den)[1]) % p, \
                (j[0] * inv_scalar(den)[1] + j[1] * inv_scalar(den)[0]) % p
            A = (3 * k[0] % p, 3 * k[1] % p)
            B = (2 * k[0] % p, 2 * k[1] % p)
            AX = smul(A, X)
            fu = (X3[0] + AX[0] + B[0]) % p
            fv = (X3[1] + AX[1] + B[1]) % p
            nrm = (fu * fu - nu * fv * fv) % p
            s = int(chi[nrm].sum())
            card = q + 1 + s
            t = q + 1 - card
            if t % p == 0:
                nss += 1
    assert nss % 2 == 0
    return nss // 2


def h_forms(D):
    """class number of primitive positive definite forms of discriminant D < 0 (reduced forms)"""
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            from math import gcd
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


def genus_X0p(p):
    nu2 = 1 + legendre_symbol(-1 % p, p)          # 1 + (-4/p)
    nu3 = 1 + (legendre_symbol(-3 % p, p) if p != 3 else 0)
    num = 12 + (p + 1) - 3 * nu2 - 4 * nu3 - 12   # 12 g = 12 + psi - 3 nu2 - 4 nu3 - 6 * 2 cusps
    assert num % 12 == 0
    return num // 12


print("p, n_p (F_p count), DG prediction, h(-4p), h(-p), [pairs, g(X0(p)), (2g+2-h4-h1)/4] for p <= %d" % PMAX2)
ogg = []
rows = []
for p in primerange(5, PMAX1):
    ss = count_Fp(p)
    n_p = len(ss)
    h4 = h_forms(-4 * p)
    h1 = h_forms(-p) if p % 4 == 3 else 0
    pred = (h4 + h1) // 2
    if (h4 + h1) % 2 or n_p != pred:
        fail(f"DG at p={p}: {n_p} vs {(h4 + h1) / 2}")
    # the three-case form
    alt = h4 // 2 if p % 4 == 1 else (h1 if p % 8 == 7 else 2 * h1)
    if alt != pred or (p % 4 == 1 and h4 % 2):
        fail(f"DG three-case form at p={p}")
    g = genus_X0p(p)
    deg = p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
    if deg != g + 1:
        fail(f"degree law vs genus at p={p}")
    row = [p, n_p, pred, h4, h1]
    if p <= PMAX2:
        pairs = count_Fp2_pairs(p)
        total = n_p + 2 * pairs
        if total != deg:
            fail(f"total ss count at p={p}: {total} vs {deg}")
        gplus = (2 * g + 2 - h4 - h1)
        if gplus % 4 or pairs != gplus // 4:
            fail(f"pairs vs g(X0(p)^+) formula at p={p}")
        row += [pairs, g, gplus // 4, total]
        if pairs == 0:
            ogg.append(p)
    rows.append(row)
    print("  ", row)
    sys.stdout.flush()
print("primes 5 <= p <= %d with no conjugate pair (all ss j in F_p), by F_{p^2} point counting:" % PMAX2, ogg)
if ogg != [5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]:
    fail("Ogg list (p <= 101)")
print("Ogg's list restricted to 5 <= p <= 101: [5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]")
print("\nfailures:", bad)
print("ALL AUDIT CHECKS PASSED" if bad == 0 else "AUDIT FOUND FAILURES")
