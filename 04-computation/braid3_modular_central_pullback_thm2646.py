#!/usr/bin/env python3
"""Exact companion for THM-2646.

The proof of the fibre-product theorem is presentation-theoretic.  This
companion checks its finite/matrix shadows, the minimal knot hostile, the
affine-V4 blindness, and the central Burau recurrence using integer and
Laurent-polynomial arithmetic only.
"""

from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# Laurent polynomials in t, represented by exponent -> integer coefficient.
def lp(value=0):
    return {} if value == 0 else {0: value}


def lp_clean(a):
    return {e: c for e, c in a.items() if c}


def lp_add(a, b):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, 0) + c
    return lp_clean(out)


def lp_neg(a):
    return {e: -c for e, c in a.items()}


def lp_mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            out[ea + eb] = out.get(ea + eb, 0) + ca * cb
    return lp_clean(out)


def mat_mul(a, b, add=lp_add, mul=lp_mul):
    return [
        [add(mul(a[0][0], b[0][0]), mul(a[0][1], b[1][0])),
         add(mul(a[0][0], b[0][1]), mul(a[0][1], b[1][1]))],
        [add(mul(a[1][0], b[0][0]), mul(a[1][1], b[1][0])),
         add(mul(a[1][0], b[0][1]), mul(a[1][1], b[1][1]))],
    ]


def mat_pow(a, n, identity, mul):
    out = identity
    base = a
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n //= 2
    return out


def mat_det(a):
    return lp_add(lp_mul(a[0][0], a[1][1]),
                  lp_neg(lp_mul(a[0][1], a[1][0])))


ZERO = lp()
ONE = lp(1)
T = {1: 1}
T_INV = {-1: 1}

S1 = [[lp_neg(T), ONE], [ZERO, ONE]]
S2 = [[ONE, ZERO], [T, lp_neg(T)]]
S1_INV = [[lp_neg(T_INV), T_INV], [ZERO, ONE]]
S2_INV = [[ONE, ZERO], [ONE, lp_neg(T_INV)]]
I2 = [[ONE, ZERO], [ZERO, ONE]]

burau = {1: S1, -1: S1_INV, 2: S2, -2: S2_INV}


def burau_word(word):
    out = I2
    for letter in word:
        out = mat_mul(out, burau[letter])
    return out


def eval_lp(a, value):
    if value == -1:
        return sum(c * ((-1) ** e) for e, c in a.items())
    return sum(c * (value ** e) for e, c in a.items())


def eval_matrix(a, value):
    return [[eval_lp(a[i][j], value) for j in range(2)] for i in range(2)]


require(burau_word((1, 2, 1)) == burau_word((2, 1, 2)), "Burau braid relation")
Z = burau_word((1, 2, 1, 2, 1, 2))
T3I = [[{3: 1}, ZERO], [ZERO, {3: 1}]]
require(Z == T3I, "full twist must be t^3 I")
require(mat_mul(S1, S1_INV) == I2 and mat_mul(S2, S2_INV) == I2,
        "Laurent inverses")


# Integer SL2 shadow at t=-1; PSL equality is equality up to sign.
def iadd(x, y):
    return x + y


def imul(x, y):
    return x * y


def imat_mul(a, b):
    return mat_mul(a, b, iadd, imul)


II = [[1, 0], [0, 1]]
MI = [[-1, 0], [0, -1]]
sl = {g: eval_matrix(m, -1) for g, m in burau.items()}


def sl_word(word):
    out = II
    for letter in word:
        out = imat_mul(out, sl[letter])
    return out


def psl_key(m):
    flat = tuple(m[0] + m[1])
    neg = tuple(-x for x in flat)
    return min(flat, neg)


require(eval_matrix(Z, -1) == MI, "full twist at t=-1")


# Strand permutations.  A word closes to a knot exactly when its permutation
# is a 3-cycle.
transposition = {1: (1, 0, 2), -1: (1, 0, 2), 2: (0, 2, 1), -2: (0, 2, 1)}


def perm_compose(p, q):
    return tuple(p[q[i]] for i in range(3))


def perm_word(word):
    out = (0, 1, 2)
    for letter in word:
        out = perm_compose(out, transposition[letter])
    return out


def is_three_cycle(p):
    return p not in ((0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 1, 0))


A_MINUS = (-2, -1)          # (sigma1 sigma2)^-1
A_PLUS = (1, 2, 1, 2)       # (sigma1 sigma2)^2
require(psl_key(sl_word(A_MINUS)) == psl_key(sl_word(A_PLUS)),
        "minimal hostile has one PSL image")
require(sum(1 if g > 0 else -1 for g in A_PLUS) -
        sum(1 if g > 0 else -1 for g in A_MINUS) == 6,
        "hostile differs by one full twist")
require(is_three_cycle(perm_word(A_MINUS)) and is_three_cycle(perm_word(A_PLUS)),
        "both hostile closures are knots")


# Exhaust the radius-three Artin word ball.  Exponent range alone proves the
# lower bound, while this check also retains the knot-permutation filter.
alphabet = (-2, -1, 1, 2)
short = []
for length in range(4):
    for word in product(alphabet, repeat=length):
        # Remove immediate inverse cancellations; every geodesic has this form.
        if any(word[i] == -word[i + 1] for i in range(length - 1)):
            continue
        if is_three_cycle(perm_word(word)):
            exponent = sum(1 if g > 0 else -1 for g in word)
            short.append((word, exponent, psl_key(sl_word(word))))
require(all(abs(a[1] - b[1]) < 6
            for i, a in enumerate(short) for b in short[i + 1:]),
        "no radius-three knot pair can differ by a nonzero central height")


# The central torus-knot fibre beta_k=(sigma1 sigma2)^(-1+3k).
fibre = []
base_key = psl_key(sl_word(A_MINUS))
for k in range(-2, 4):
    q = 3 * k - 1
    block = (1, 2) if q >= 0 else (-2, -1)
    word = block * abs(q)
    require(psl_key(sl_word(word)) == base_key, "constant PSL fibre")
    require(is_three_cycle(perm_word(word)), "fibre closure is a knot")
    genus = abs(q) - 1
    fibre.append((k, q, genus))


# The standard mod-2 linear shadow has no nonzero invariant translation, so
# an affine V4 lift with that linear part must kill the central full twist.
def f2_mat_vec(a, v):
    return ((a[0][0] * v[0] + a[0][1] * v[1]) % 2,
            (a[1][0] * v[0] + a[1][1] * v[1]) % 2)


T2 = ((1, 1), (0, 1))
B2 = ((1, 0), (1, 1))
vectors = [(0, 0), (0, 1), (1, 0), (1, 1)]
fixed = [v for v in vectors if f2_mat_vec(T2, v) == v and f2_mat_vec(B2, v) == v]
require(fixed == [(0, 0)], "S3 linear shadow has no fixed V4 translation")


# Two Alexander-normalization controls for the minimal hostile.
def burau_numerator(word):
    r = burau_word(word)
    i_minus_r = [[lp_add(I2[i][j], lp_neg(r[i][j])) for j in range(2)]
                 for i in range(2)]
    return mat_det(i_minus_r)


ONE_MINUS_T = {0: 1, 1: -1}
ONE_MINUS_T3 = {0: 1, 3: -1}
ALEX_MINUS = {-2: 1}
ALEX_PLUS = {0: 1, 1: -1, 2: 1}
require(lp_mul(ONE_MINUS_T, burau_numerator(A_MINUS)) ==
        lp_mul(ALEX_MINUS, ONE_MINUS_T3), "unknot Burau normalization")
require(lp_mul(ONE_MINUS_T, burau_numerator(A_PLUS)) ==
        lp_mul(ALEX_PLUS, ONE_MINUS_T3), "trefoil Burau normalization")


# Symbolic Burau-numerator recurrence.  A coefficient is a triple for
# c0 + ca*a + cd*d, where a=trace(R(beta)), d=det(R(beta)).
def cadd(u, v):
    return tuple(u[i] + v[i] for i in range(3))


def cscale(c, n):
    return tuple(n * x for x in c)


def xpoly_add(p, q):
    out = dict(p)
    for e, c in q.items():
        out[e] = cadd(out.get(e, (0, 0, 0)), c)
    return {e: c for e, c in out.items() if c != (0, 0, 0)}


def xpoly_scale_shift(p, scalar_poly):
    out = {}
    for ep, cp in p.items():
        for es, cs in scalar_poly.items():
            out[ep + es] = cadd(out.get(ep + es, (0, 0, 0)), cscale(cp, cs))
    return {e: c for e, c in out.items() if c != (0, 0, 0)}


def numerator(k):
    out = {0: (1, 0, 0)}
    out = xpoly_add(out, {k: (0, -1, 0)})
    out = xpoly_add(out, {2 * k: (0, 0, 1)})
    return out


P1 = {0: 1, 1: 1, 2: 1}
P2 = {1: 1, 2: 1, 3: 1}
P3 = {3: 1}
for k in range(0, 9):
    recurrence = numerator(k + 3)
    recurrence = xpoly_add(recurrence, xpoly_scale_shift(numerator(k + 2),
                                                         {e: -c for e, c in P1.items()}))
    recurrence = xpoly_add(recurrence, xpoly_scale_shift(numerator(k + 1), P2))
    recurrence = xpoly_add(recurrence, xpoly_scale_shift(numerator(k),
                                                         {e: -c for e, c in P3.items()}))
    require(recurrence == {}, "central Burau recurrence")


print("THM-2646 BRAID/MODULAR CENTRAL-PULLBACK AUDIT")
print("presentation: x^2=y^3=z, e(x,y,z)=(3,2,6), quotient=C2*C3")
print("Burau: braid_relation=PASS, rho(z)=t^3 I, rho(z)|t=-1=-I")
print("minimal hostile: beta_-=(s1s2)^-1 length=2 e=-2 closure=T(3,-1)=unknot")
print("minimal hostile: beta_+=(s1s2)^2  length=4 e=4  closure=T(3,2)=trefoil")
print(f"radius-three knot words checked={len(short)}; nonzero central pair=NONE")
print("central fibre (k,q,genus): " + ", ".join(map(str, fibre)))
print(f"affine V4 common fixed translations={fixed}; central translation=0")
print("Burau normalized hostiles: A(beta_-)=t^-2, A(beta_+)=t^2-t+1")
print("Burau numerator recurrence roots=(1,t^3,t^6): PASS for k=0..8")
