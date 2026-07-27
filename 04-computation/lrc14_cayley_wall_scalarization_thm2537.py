#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2537.

All linear algebra is over Z or Fraction.  The exhaustive Boolean census is
over every nonconstant mask on F_13 and every nonzero oriented slope.
"""

from collections import Counter
from fractions import Fraction


Q = 13
ROOTS = range(Q)
UNITS = range(1, Q)


def eye():
    return [[int(i == j) for j in ROOTS] for i in ROOTS]


def mat_add(a, b):
    return [[a[i][j] + b[i][j] for j in ROOTS] for i in ROOTS]


def mat_sub(a, b):
    return [[a[i][j] - b[i][j] for j in ROOTS] for i in ROOTS]


def mat_mul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in ROOTS) for j in ROOTS]
        for i in ROOTS
    ]


def mat_vec(a, x):
    return [sum(a[i][j] * x[j] for j in ROOTS) for i in ROOTS]


def scale(a, c):
    return [[c * a[i][j] for j in ROOTS] for i in ROOTS]


def shift(tau):
    # (P_tau x)_v=x_(v+tau).
    return [[int(j == (i + tau) % Q) for j in ROOTS] for i in ROOTS]


def power(a, n):
    ans = eye()
    for _ in range(n):
        ans = mat_mul(ans, a)
    return ans


def cayley(tau):
    p = shift(tau)
    ans = [[0 for _ in ROOTS] for _ in ROOTS]
    pk = eye()
    for d in range(1, Q):
        pk = mat_mul(pk, p)
        ans = mat_add(ans, scale(pk, 1 if d % 2 else -1))
    return ans


def cayley_vec(x, tau):
    return [
        sum((1 if d % 2 else -1) * x[(v + d * tau) % Q] for d in range(1, Q))
        for v in ROOTS
    ]


def green(tau):
    p = shift(tau)
    ans = [[Fraction(0) for _ in ROOTS] for _ in ROOTS]
    pk = eye()
    for d in range(1, Q):
        pk = mat_mul(pk, p)
        ans = mat_add(ans, scale(pk, Fraction(2 * d - Q, Q)))
    return ans


def bits(mask):
    return [(mask >> r) & 1 for r in ROOTS]


def word(e, tau, a):
    return tuple(e[(a + j * tau) % Q] for j in ROOTS)


def selector(e, tau):
    a = max(ROOTS, key=lambda r: word(e, tau, r))
    q = next(j for j in range(1, Q) if not e[(a + j * tau) % Q])
    s = (a + (q - 1) * tau) % Q
    t = (s + tau) % Q
    return a, q, s, t


def push(e, u, h):
    inv = pow(u, -1, Q)
    return [e[(inv * (r - h)) % Q] for r in ROOTS]


def correlation(e, d):
    return sum(e[r] * e[(r + d) % Q] for r in ROOTS)


def psi(e, tau):
    coeff = (7, -12, 8, -6, 7, -6, 2)
    return sum(coeff[j] * correlation(e, j * tau) for j in range(7))


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


I = eye()
J13 = [[Fraction(1, Q) for _ in ROOTS] for _ in ROOTS]

# Matrix identities, Green inverse, and the sharp centred denominator.
for tau in UNITS:
    P = shift(tau)
    C = cayley(tau)
    B = green(tau)
    assert mat_mul(mat_add(I, P), C) == mat_sub(P, I)
    assert mat_mul(C, mat_add(I, P)) == mat_sub(P, I)
    assert mat_mul(C, B) == mat_sub(I, J13)
    assert mat_mul(B, C) == mat_sub(I, J13)
    assert cayley((-tau) % Q) == scale(C, -1)
    for s in ROOTS:
        t = (s + tau) % Q
        rho = [int(r == s) - int(r == t) for r in ROOTS]
        brho = mat_vec(B, rho)
        expected = [Fraction(int(r == s) + int(r == t), 1) - Fraction(2, Q) for r in ROOTS]
        assert brho == expected
        assert all((Q * x).denominator == 1 for x in brho)
        assert any(x.denominator == Q for x in brho)
        assert sum(r * rho[r] for r in ROOTS) % Q == (-tau) % Q
        assert mat_vec(C, expected) == rho

# Deep-comb global boundary identity on every ray and on one independent
# positive integer combination of all 23 rays.
deep_masks = []
for j in range(1, Q):
    e = [0] * Q
    e[j] = 1
    deep_masks.append(("s", j, e))
for j in range(1, Q - 1):
    e = [0] * Q
    e[j] = e[j + 1] = 1
    deep_masks.append(("p", j, e))

def deep_profiles(weighted_masks):
    H = [0] * Q
    gp = [0] * Q
    gm = [0] * Q
    for w, e in weighted_masks:
        for a in ROOTS:
            H[a] += w * e[a]
            gp[a] += w * e[a] * (1 - e[(a - 1) % Q])
            gm[a] += w * e[a] * (1 - e[(a + 1) % Q])
    return H, gp, gm

deep_tests = [[(1, e)] for _, _, e in deep_masks]
deep_tests.append([(i + 1, e) for i, (_, _, e) in enumerate(deep_masks)])
P1, C1 = shift(1), cayley(1)
for test in deep_tests:
    H, gp, gm = deep_profiles(test)
    lhs = mat_vec(mat_sub(P1, I), H)
    rhs = [x - y for x, y in zip(mat_vec(P1, gp), gm)]
    assert lhs == rhs == mat_vec(mat_mul(mat_add(I, P1), C1), H)

# Exhaustive wall/global-boundary census.
pair_census = Counter()
fixed_edge_zero = 0
positive_masks = 0
positive_threshold_layers = 0
psi_total = 0
threshold_total = 0
scalarized_total = 0
affine_generator_checks = 0
for tau in UNITS:
    C = cayley(tau)
    P = shift(tau)
    for mask in range(1, 2**Q - 1):
        e = bits(mask)
        a, run, s, t = selector(e, tau)
        assert e[s] == 1 and e[t] == 0 and t == (s + tau) % Q
        ce = cayley_vec(e, tau)

        # Selected endpoint identity and its integral dual representative.
        assert ce[s] + ce[t] == -1
        # The centred Green dual was checked above for every slope/wall.
        assert sum(ce) == 0

        # Positive global boundary transform and Crofton/energy identity.
        bp = [e[v] * (1 - e[(v + tau) % Q]) for v in ROOTS]
        bm = [(1 - e[v]) * e[(v + tau) % Q] for v in ROOTS]
        grad = mat_vec(mat_sub(P, I), e)
        assert [bm[v] - bp[v] for v in ROOTS] == grad
        assert mat_vec(mat_add(I, P), ce) == grad
        assert sum(bp) == sum(bm) > 0
        assert dot(grad, grad) == 2 * sum(bp)

        # Root-independent owner weights commute with every identity.
        for g in (Fraction(0), Fraction(2, 3), Fraction(1)):
            assert g * (ce[s] + ce[t]) == -g

        # The THM-2527 score and its finite threshold scalarization.
        ps = psi(e, tau)
        assert 1 <= ps <= 98
        layers = sum(int(ps >= j) for j in range(1, 99))
        assert layers == ps
        positive_masks += 1
        positive_threshold_layers += layers
        psi_total += ps
        threshold_total += layers
        scalarized_total += -layers * (ce[s] + ce[t])
        pair_census[(ce[s], ce[t])] += 1

        # A fixed edge has no force without the selected-wall sidecar.
        if ce[0] + ce[tau] == 0:
            fixed_edge_zero += 1

        # Translation and multiplicative generators imply full affine
        # covariance; reflection is also checked explicitly.
        for u, h in ((1, 1), (2, 0), (-1 % Q, 0)):
            ep = push(e, u, h)
            ap, qp, sp, tp = selector(ep, (u * tau) % Q)
            assert (ap, qp, sp, tp) == (
                (u * a + h) % Q,
                run,
                (u * s + h) % Q,
                (u * t + h) % Q,
            )
            cp = cayley_vec(ep, (u * tau) % Q)
            assert cp[sp] + cp[tp] == ce[s] + ce[t] == -1
            affine_generator_checks += 1

        # Complement reverses this particular wall and slope.  Its own
        # canonical wall may be different, but is covered elsewhere in the
        # exhaustive loop.
        ec = [1 - x for x in e]
        cc = cayley_vec(ec, (-tau) % Q)
        assert ec[t] == 1 and ec[s] == 0
        assert cc[t] + cc[s] == -1

assert psi_total == threshold_total == scalarized_total

# Sharp small hostiles.
tau = 1
C = cayley(tau)
singleton = [1] + [0] * 12
adjacent = [1, 1] + [0] * 11
fixed_miss = [0, 0, 1] + [0] * 10
cs = mat_vec(C, singleton)
ca = mat_vec(C, adjacent)
cf = mat_vec(C, fixed_miss)
assert sum(cs) == sum(ca) == sum(cf) == 0
assert (cs[0], cs[1]) == (0, -1)       # source-only can vanish
assert (ca[1], ca[2]) == (-1, 0)       # head-only can vanish
assert cf[0] + cf[1] == 0              # fixed unselected edge cancels
assert cf[2] + cf[3] == -1             # selected edge does not

print("THM-2537 exact referee")
print("q=13 slopes=12 nonconstant_masks=8190 slope_mask_states=98280")
print("cayley_endpoint_identity=PASS green_dual=PASS denominator=13_sharp_on_augmentation")
print("global_boundary=Bminus-Bplus=(P-I)e=(I+P)Ce")
print("crofton_energy=PASS equal_boundary_masses=PASS")
print("deep_comb=(P-I)H=Pgamma_plus-gamma_minus=(I+P)C(H): PASS")
print(f"positive_mask_states={positive_masks} positive_threshold_layers={positive_threshold_layers}")
print(f"psi_sum=threshold_sum=selected_scalar_sum={psi_total}")
print(f"affine_generator_reflection_checks={affine_generator_checks}")
print(f"fixed_unselected_edge_zero_states={fixed_edge_zero}")
print("tau1_selected_endpoint_pair_census=" + repr(sorted((k, v // 12) for k, v in pair_census.items())))
print("hostiles=untwisted_average_zero,source_only_zero,head_only_zero,fixed_edge_zero: PASS")
print("THM-2537 PASS")
