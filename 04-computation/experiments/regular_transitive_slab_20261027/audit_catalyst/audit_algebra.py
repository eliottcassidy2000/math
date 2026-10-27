"""Independent audit of THM-4820 (A1), (B), (C), (E)-lower-bound algebra.

Everything here is recomputed from scratch; no producer code is imported.

Two independent value models are used for an Eisenstein integer:
  * coordinates (a, b) for a + b*w, Q = a^2 - a b + b^2;
  * the 'sqrt(-3) model'  value = (X + Y sqrt(-3)) / 2 with X = 2a - b, Y = b, so |value|^2 = (X^2 + 3 Y^2) / 4.
A count vector (p, q, r) (p ones, q omegas, r omega^2's) has a = p - r, b = q - r, i.e. X = 2p - q - r, Y = q - r.
"""
import itertools
import random
import sys
from fractions import Fraction

import numpy as np
import sympy

random.seed(20261027)
FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg)


def Qab(a, b):
    return a * a - a * b + b * b


def Qsqrt(X, Y):  # exact norm in the sqrt(-3) model
    v = Fraction(X * X + 3 * Y * Y, 4)
    assert v.denominator == 1
    return int(v)


def counts_to_XY(u):
    p, q, r = u
    return (2 * p - q - r, q - r)


def counts_to_ab(u):
    p, q, r = u
    return (p - r, q - r)


# ---------------------------------------------------------------------------------------------------------------
# (A1) count formula beta = 3 u.u' - |u||u'|, against the norm computed in the sqrt(-3) model
# ---------------------------------------------------------------------------------------------------------------
def beta_counts(u, v):
    return 3 * sum(x * y for x, y in zip(u, v)) - sum(u) * sum(v)


def beta_sqrt(u, v):
    X1, Y1 = counts_to_XY(u)
    X2, Y2 = counts_to_XY(v)
    return Qsqrt(X1 + X2, Y1 + Y2) - Qsqrt(X1, Y1) - Qsqrt(X2, Y2)


def beta_complex(u, v):  # floating complex model with omega = exp(2 pi i / 3), rounded
    w = complex(-0.5, 3 ** 0.5 / 2)
    y = u[0] + u[1] * w + u[2] * w * w
    z = v[0] + v[1] * w + v[2] * w * w
    return abs(y + z) ** 2 - abs(y) ** 2 - abs(z) ** 2


n_a1 = 0
for _ in range(20000):
    u = tuple(random.randint(0, 25) for _ in range(3))
    v = tuple(random.randint(0, 25) for _ in range(3))
    b1, b2, b3 = beta_counts(u, v), beta_sqrt(u, v), beta_complex(u, v)
    check(b1 == b2 and abs(b1 - b3) < 1e-6, f"A1 count formula {u} {v}: {b1} {b2} {b3}")
    n_a1 += 1
# exhaustive on small vectors as well
for u in itertools.product(range(5), repeat=3):
    for v in itertools.product(range(5), repeat=3):
        check(beta_counts(u, v) == beta_sqrt(u, v), f"A1 exhaustive {u} {v}")
        n_a1 += 1
print(f"(A1) count formula: {n_a1} pairs checked against the sqrt(-3) value model and the complex embedding")

# inner products of the units
w = complex(-0.5, 3 ** 0.5 / 2)
for i in range(3):
    for j in range(3):
        ip = (w ** i * (w ** j).conjugate()).real
        check(abs(ip - (1 if i == j else -0.5)) < 1e-12, f"unit inner product {i},{j}")

# (A2) neutral-triple merge sequence costs
units6 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]  # 1, -1, w, -w, -w^2 = 1 + w, w^2 = -1 - w
one, om, om2 = (1, 0), (0, 1), (-1, -1)


def add(y, z):
    return (y[0] + z[0], y[1] + z[1])


def beta(y, z):
    return Qab(*add(y, z)) - Qab(*y) - Qab(*z)


check(beta(one, om) == -1, "A2 beta(1,w) = -1")
check(add(add(one, om), om2) == (0, 0) and beta(add(one, om), om2) == -2, "A2 (1+w)+w^2 = 0 with beta -2")
check(all(beta(x, (0, 0)) == 0 for x in [(5, 3), (-7, 2), (11, 0)]), "A2 merging 0 costs 0")
# A4: negative units are free merges of two positive units with beta = -1
check(add(om, om2) == (-1, 0) and beta(om, om2) == -1, "A4 -1 = w + w^2")
check(add(one, om2) == (0, -1) and beta(one, om2) == -1, "A4 -w = 1 + w^2")
check(add(one, om) == (1, 1) and beta(one, om) == -1, "A4 -w^2 = 1 + w")
check((1, 1) == (-om2[0], -om2[1]), "A4 (1,1) is -w^2")
print("(A2)/(A4) free merges verified: 1+w (beta -1), (1+w)+w^2 = 0 (beta -2), x+0 (beta 0); -1, -w, -w^2 as beta -1 merges")

# (A3) all-ones trees: random trees on n ones always cost n(n-1)
for n in range(1, 40):
    for trial in range(20):
        piles = [1] * n
        cost = 0
        while len(piles) > 1:
            i, j = random.sample(range(len(piles)), 2)
            a, b = piles[i], piles[j]
            bt = (a + b) ** 2 - a * a - b * b
            check(bt == 2 * a * b > 0, "A3 beta = 2ab > 0")
            cost += max(0, bt)
            piles = [p for k, p in enumerate(piles) if k not in (i, j)] + [a + b]
        check(cost == n * (n - 1), f"A3 n={n} cost={cost}")
print("(A3) every random tree on n ones (n < 40, 20 trees each) costs n(n-1)")

# ---------------------------------------------------------------------------------------------------------------
# (B) potential bound, slack formula, tight characterisation, confinement of tight merges -- exhaustive in a box
# ---------------------------------------------------------------------------------------------------------------
R = 14
pts = [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1)]
Phi = lambda y: Qab(*y) // 3
n_b = 0
max_child_excess = {0: -10 ** 9, 1: -10 ** 9}
tight_beta_values = set()
for y in pts:
    qy = Qab(*y)
    check(qy % 3 in (0, 1), f"norm residue {y}")
    for z in pts:
        x = add(y, z)
        qx, qz = Qab(*x), Qab(*z)
        bt = qx - qy - qz
        r, s = qy % 3, qz % 3
        d = Phi(x) - Phi(y) - Phi(z)
        check(d == (r + s + bt) // 3, f"B floor identity {y} {z}")
        check(d <= max(bt, 0), f"B potential inequality {y} {z}")
        slack = max(bt, 0) - d
        check(slack >= 0, "B slack nonneg")
        tight_pred = (bt == 1 and r == 1 and s == 1) or (-(r + s) <= bt <= 0)
        check((slack == 0) == tight_pred, f"B tight characterisation {y} {z} beta={bt} r={r} s={s}")
        if slack == 0:
            tight_beta_values.add(bt)
            if y != (0, 0) and z != (0, 0):
                # confinement: children norms <= Q(x) + 1, and <= Q(x) when Q(x) = 1 mod 3
                ex = max(qy, qz) - qx
                max_child_excess[qx % 3] = max(max_child_excess[qx % 3], ex)
        if qx % 3 == 1 and y != (0, 0) and z != (0, 0):
            check((slack == 0) == (bt in (0, -1)), f"B Q=1 tight splits are beta in {{0,-1}} {y} {z}")
        n_b += 1
print(f"(B) {n_b} ordered pairs in the box |a|,|b| <= {R}: floor identity, Phi inequality, slack >= 0, tight "
      f"characterisation, and 'Q(x) = 1 mod 3 => tight split iff beta in {{0,-1}}' all verified")
print(f"    beta values of tight merges seen: {sorted(tight_beta_values)}")
print(f"    max over tight merges with both children nonzero of max(Q(y),Q(z)) - Q(x): "
      f"Q(x)=0 mod 3: {max_child_excess[0]}, Q(x)=1 mod 3: {max_child_excess[1]}  (confinement lemma predicts 1 and 0)")
check(max_child_excess[0] <= 1 and max_child_excess[1] <= 0, "confinement of tight merges")


# random trees: cost = Phi(root) + sum of slacks (telescoping)
def random_tree_check(nleaves):
    piles = [random.choice(units6) for _ in range(nleaves)]
    cost = 0
    slack_sum = 0
    while len(piles) > 1:
        i, j = random.sample(range(len(piles)), 2)
        y, z = piles[i], piles[j]
        x = add(y, z)
        bt = beta(y, z)
        cost += max(bt, 0)
        slack_sum += max(bt, 0) - (Phi(x) - Phi(y) - Phi(z))
        piles = [p for k, p in enumerate(piles) if k not in (i, j)] + [x]
    return cost, Phi(piles[0]) + slack_sum, Phi(piles[0])


for trial in range(3000):
    c, c2, ph = random_tree_check(random.randint(1, 60))
    check(c == c2 and c >= ph, f"telescoping {c} {c2}")
print("(B) telescoping cost = Phi(root) + sum slack verified on 3000 random six-unit trees")

# ---------------------------------------------------------------------------------------------------------------
# (C) norm-form criteria for rational n -- brute force over all y = a + b w with Q(y) <= n^2 + 2
# ---------------------------------------------------------------------------------------------------------------
def has_split_prime(n):
    return any(p % 3 == 1 for p in sympy.factorint(n))


def rep_b2_3t2_brute(M):
    t = 0
    while 3 * t * t <= M:
        r = M - 3 * t * t
        s = int(round(r ** 0.5))
        for c in (s - 1, s, s + 1):
            if c >= 0 and c * c == r:
                return (c, t)
        t += 1
    return None


def rep_b2_3t2_arith(M):  # classical: every prime = 2 mod 3 (incl. 2) to an even power
    return all(e % 2 == 0 for p, e in sympy.factorint(M).items() if p % 3 == 2)


NMAX = 241
rows = []
for n in range(1, NMAX + 1, 2):
    if n % 3 == 0:
        continue
    L = 2 * n + 2
    A, B = np.meshgrid(np.arange(-L, L + 1), np.arange(-L, L + 1), indexing="ij")
    A = A.astype(np.int64)
    B = B.astype(np.int64)
    qy = A * A - A * B + B * B
    za, zb = n - A, -B
    qz = za * za - za * zb + zb * zb
    bt = n * n - qy - qz
    # translation identity (u - n)^2 + 3 b^2 = n^2 - 2 beta, u = 2a - b
    U = 2 * A - B
    check(bool(np.all((U - n) ** 2 + 3 * B * B == n * n - 2 * bt)), f"C translation identity n={n}")
    nontriv = ~(((A == 0) & (B == 0)) | ((A == n) & (B == 0)))
    perp = bool(np.any((bt == 0) & nontriv))
    obt = bool(np.any(bt == -1))
    # all tight split parts have norm <= n^2 (box large enough: Q(y) <= n^2 forces |a|,|b| <= 2n/sqrt(3) < L)
    check(not np.any(((bt == 0) | (bt == -1)) & nontriv & (qy > n * n)), f"C tight parts have norm <= n^2, n={n}")
    M = (n * n + 2) // 3
    check(3 * M == n * n + 2 and M % 8 == 1, f"C M integral and = 1 mod 8, n={n}")
    perp_pred = has_split_prime(n)
    obt_pred = rep_b2_3t2_brute(M) is not None
    check(perp == perp_pred, f"C perpendicular criterion n={n}: brute {perp} predicted {perp_pred}")
    check(obt == obt_pred, f"C obtuse criterion n={n}: brute {obt} predicted {obt_pred}")
    check(rep_b2_3t2_arith(M) == obt_pred, f"C inert-odd-exponent criterion n={n}")
    if not perp and not obt:
        rows.append(n)
print(f"(C) brute-force tight splits for every odd n <= {NMAX}, 3 not dividing n: perpendicular-split criterion "
      f"(prime = 1 mod 3), obtuse criterion (M = b^2 + 3t^2, brute force) and the inert-odd-exponent form all agree")
print(f"    odd n <= {NMAX}, 3 not dividing n, with NO tight split: {rows}")

print("(C) table rows:")
for n in (11, 23, 25, 29, 41, 47):
    M = (n * n + 2) // 3
    fac = sympy.factorint(M)
    nfac = sympy.factorint(n)
    inert_odd = [p for p, e in fac.items() if p % 3 == 2 and e % 2 == 1]
    print(f"    n={n}: prime factors of n {dict(nfac)} (all = 2 mod 3: {all(p % 3 == 2 for p in nfac)}), "
          f"M={M} = {dict(fac)}, inert primes to odd power: {inert_odd}, brute b^2+3t^2: {rep_b2_3t2_brute(M)}")
    check(all(p % 3 == 2 for p in nfac) and inert_odd and rep_b2_3t2_brute(M) is None, f"C row n={n}")
# the 50 < n <= 86 predictions mentioned in (D)
for n in (53, 55, 59, 71, 83, 85):
    M = (n * n + 2) // 3
    print(f"    n={n}: n factors {dict(sympy.factorint(n))}, M={M} = {dict(sympy.factorint(M))}, "
          f"perp split {has_split_prime(n)}, obtuse rep {rep_b2_3t2_brute(M)}")

# (E) sum of beta over any tree = Q(root) - #leaves (three-unit leaves): random check
for trial in range(2000):
    n = random.randint(1, 15)
    k = random.randint(0, 10)
    piles = [one] * (n + k) + [om] * k + [om2] * k
    tot = 0
    while len(piles) > 1:
        i, j = random.sample(range(len(piles)), 2)
        tot += beta(piles[i], piles[j])
        x = add(piles[i], piles[j])
        piles = [p for kk, p in enumerate(piles) if kk not in (i, j)] + [x]
    check(piles[0] == (n, 0) and tot == n * n - (n + 3 * k), "E sum beta")
print("(E) sum of beta over a tree on (n+k,k,k) = n^2 - (n+3k): 2000 random trees")
print("(E) n=11 lower bound (121-11-41)/3 =", Fraction(121 - 11 - 41, 3))

print()
print("TOTAL FAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
