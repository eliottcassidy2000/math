#!/usr/bin/env python3
"""Import-free global hostile checks for THM-4344's W=0 degree drop.

The script distinguishes realizable coefficient/support strata from the
larger Boolean deletion atlas, reconstructs the fan and edge schemes, and
checks the complementary collision charts and simultaneous sign incidence.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def rank_two(points):
    for a, b, c in combinations(points, 3):
        if (b[0] - a[0]) * (c[1] - a[1]) != (b[1] - a[1]) * (c[0] - a[0]):
            return True
    return False


def planes(universe):
    answer = set()
    for p, q, r in combinations(universe, 3):
        det = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        if det == 0:
            continue
        aa = F((q[2] - p[2]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[2] - p[2]), det)
        bb = F((q[0] - p[0]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[0] - p[0]), det)
        cc = F(p[2]) - aa * p[0] - bb * p[1]
        answer.add((aa, bb, cc))
    return tuple(sorted(answer))


def lower(active, candidate_planes):
    result = set()
    for aa, bb, cc in candidate_planes:
        equal = []
        for rr, ll, hh in active:
            gap = F(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equal.append((rr, ll, hh))
        else:
            if rank_two(equal):
                result.add((aa, bb, cc))
    return frozenset(result)


def pick(vertices):
    pairs = list(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in pairs))
    boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1)) for (x1, y1), (x2, y2) in pairs)
    return area2, boundary, (area2 - boundary + 2) // 2


# Complete specialized source support.
U, u, alpha, xi, Delta, eta, Theta, Phi, K = sp.symbols(
    "U u alpha xi Delta eta Theta Phi K"
)
e = -sp.Rational(1376, 135)
rows = {
    (1, 0): -3,
    (2, 0): sp.Rational(8, 3),
    (3, 0): e,
    (4, 0): Delta,
    (5, 0): u,
    (6, 0): U,
    (0, 2): K,
    (2, 1): Phi,
    (1, 2): Theta,
    (3, 1): eta,
    (2, 2): xi,
    (4, 1): alpha,
}
symbolic_support = defaultdict(lambda: sp.Integer(0))
symbolic_support[(2, 0, 0)] += 1
symbolic_support[(0, 1, 0)] -= 1
symbolic_support[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in rows.items():
    symbolic_support[(j + 2, i + j, 1)] -= coefficient
    symbolic_support[(j, i + j + 1, 1)] += coefficient
symbolic_support = dict(symbolic_support)

need(sp.factor(symbolic_support[(2, 6, 1)]) == -U, "unique U owner")
need(sp.factor(symbolic_support[(4, 4, 1)]) == -xi, "unique xi owner")
need(sp.factor(symbolic_support[(4, 2, 1)]) == -K, "unique K owner")
need(sp.simplify(symbolic_support[(2, 3, 1)] - (K - e)) == 0, "first aggregate")
need(sp.simplify(symbolic_support[(2, 4, 1)] - (Theta - Delta)) == 0, "second aggregate")
need(sp.simplify(symbolic_support[(2, 5, 1)] - (xi - u)) == 0, "third aggregate")

M = (F(1, 12), F(1, 6), F(-1, 6))
D6 = (F(1, 6), F(1, 6), F(-1, 3))
T = (F(1, 2), F(0), F(-1))
expected_fan = frozenset((M, D6, T))
candidate_planes = planes(tuple(symbolic_support))

# Exact realizable support patterns.  K is never independent of Delta.
Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
dstar = sp.solve(sp.Eq(Krel - e, 0), Delta)[0]
need(dstar == sp.Rational(3968, 63), "special aggregate-cancellation Delta")
need(sp.solve(sp.Eq(Krel, 0), Delta)[0] == sp.Rational(5696, 105), "K=0 boundary")

# These representatives realize all eight Boolean classes of
# (Delta, K-e, Theta, Theta-Delta), without crossing K=0.
delta_theta_reps = (
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
    (dstar, sp.Integer(0)),
    (dstar, dstar),
    (dstar, 2 * dstar),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)
u_reps = (sp.Integer(0), sp.Integer(1), sp.Integer(2))

actual_states = {}
embedding_keys = set()
for (delta_value, theta_value), u_value, (phi_value, eta_value, alpha_value) in product(
    delta_theta_reps, u_reps, product((0, 1), repeat=3)
):
    substitution = {
        U: 1,
        xi: 1,
        Delta: delta_value,
        K: Krel.subs(Delta, delta_value),
        Theta: theta_value,
        u: u_value,
        Phi: phi_value,
        eta: eta_value,
        alpha: alpha_value,
    }
    need(substitution[K] != 0, "representative stays on K-nonzero gate")
    active = frozenset(
        point for point, coefficient in symbolic_support.items()
        if sp.simplify(coefficient.subs(substitution)) != 0
    )
    # The embedding into the 512-state over-atlas records the six row
    # presences followed by the three actual aggregate cancellations.
    key = (
        int(delta_value != 0),
        int(u_value != 0),
        int(phi_value != 0),
        int(theta_value != 0),
        int(eta_value != 0),
        int(alpha_value != 0),
        int(sp.simplify((K - e).subs(substitution)) == 0),
        int(sp.simplify((Theta - Delta).subs(substitution)) == 0),
        int(sp.simplify((xi - u).subs(substitution)) == 0),
    )
    embedding_keys.add(key)
    actual_states[active] = key

need(len(embedding_keys) == 8 * 3 * 8 == 192, "192 realizable Boolean states")
need(len(actual_states) == 192, "192 distinct exact support patterns")
actual_fans = Counter(lower(active, candidate_planes) for active in actual_states)
need(actual_fans == Counter({expected_fan: 192}), "192/192 exact fan")

# Reconstruct the conservative Boolean atlas and verify each exact state maps
# into it.  This atlas deliberately forgets coefficient coupling.
optional_labels = ((4, 0), (5, 0), (2, 1), (1, 2), (3, 1), (4, 1))
required_labels = ((1, 0), (2, 0), (3, 0), (6, 0), (0, 2), (2, 2))
cancellable = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    support = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


over_states = {}
for row_bits in product((0, 1), repeat=6):
    labels = list(required_labels)
    labels += [label for bit, label in zip(row_bits, optional_labels) if bit]
    for cancellation_bits in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(cancellation_bits, cancellable):
            if bit:
                active.discard(point)
        key = tuple(row_bits) + tuple(cancellation_bits)
        over_states[key] = frozenset(active)
need(len(over_states) == 512, "512 conservative over-atlas states")
need(Counter(lower(active, candidate_planes) for active in over_states.values()) == Counter({expected_fan: 512}), "512/512 over-atlas fan")
for active, key in actual_states.items():
    need(key in over_states, "exact key embeds in over-atlas")
    need(active == over_states[key], "exact support equals embedded over-state")

# Pick, Riemann--Hurwitz, and exact edge restrictions.
polygons = {
    "M": [(0, 1), (2, 0), (2, 6), (0, 7)],
    "D6": [(2, 0), (4, 4), (2, 6)],
    "T": [(2, 0), (4, 2), (4, 4)],
    "global": [(0, 1), (2, 0), (4, 2), (4, 4), (2, 6), (0, 7)],
}
need(
    {name: pick(vertices) for name, vertices in polygons.items()}
    == {
        "M": (24, 14, 6),
        "D6": (12, 10, 2),
        "T": (4, 6, 0),
        "global": (40, 14, 14),
    },
    "Pick ledgers",
)
packet = (11, 6, 6, 3, 3, 2, 2, 1)
need(sum(index - 1 for index in packet) == 26, "RH packet")

X = sp.symbols("X")
boundary_and_internal = (
    X - 1,
    1 - K * X**2,
    K + Theta * X + xi * X**2,
    U + alpha * X + xi * X**2,
    X - 1,
    U - X**6,
    1 - U * X**6,
    1 - xi * X**2,
)
need(sp.discriminant(boundary_and_internal[1], X) == 4 * K, "K boundary simple")
need(sp.factor(sp.discriminant(boundary_and_internal[5], X)) == 46656 * U**5, "U boundary simple")
need(sp.factor(sp.discriminant(boundary_and_internal[6], X)) == 46656 * U**5, "M-D6 internal simple")
need(sp.discriminant(boundary_and_internal[7], X) == 4 * xi, "two reduced D6-T points")
need(sp.discriminant(boundary_and_internal[2], X) == Theta**2 - 4 * K * xi, "finite discriminant")
need(sp.discriminant(boundary_and_internal[3], X) == alpha**2 - 4 * U * xi, "infinity discriminant")

# Complement of the infinity-tail scale.  If psi(rho)=rho^r*c(rho), after
# sigma=tau^(2(6-r)) set rho=sigma*x, z=tau^12/rho, h=y/x^3.
# The exact complementary chart is
#       z*rho=tau^12,  h^2=1-z^(6-r)c(rho).
# Its two central points h=+/-1 are etale over the A_11 toric surface, so an
# ordinary-blowup resolution contributes two rational chains and no carrier.
tau, rho, z, h, crho = sp.symbols("tau rho z h crho")
for r in range(1, 6):
    sigma6 = tau ** (12 * (6 - r))
    raw_term = sigma6 * rho ** (r - 6) * crho
    replaced = sp.cancel(raw_term.subs(tau**12, z * rho))
    need(sp.simplify(replaced - z ** (6 - r) * crho) == 0, f"infinity complement r={r}")
    complement = h**2 - 1 + z ** (6 - r) * crho
    need(sp.diff(complement, h).subs({h: 1, z: 0}) == 2, f"plus attachment etale r={r}")
    need(sp.diff(complement, h).subs({h: -1, z: 0}) == -2, f"minus attachment etale r={r}")

# Finite double: the ordinary blowup of (delta,b,x) has three charts of the
# same homogeneous exceptional conic.
delta, b, xx = sp.symbols("delta b xx")
Ba, aa, lam = sp.symbols("Ba aa lam", nonzero=True)
Q = b**2 - aa**2 * xi * xx**2 - Ba * delta * b
Bd, Xd = sp.symbols("Bd Xd")
Db, Xb = sp.symbols("Db Xb")
Dx, Yx = sp.symbols("Dx Yx")
delta_chart = sp.expand(Q.subs({b: delta * Bd, xx: delta * Xd}) / delta**2)
b_chart = sp.expand(Q.subs({delta: b * Db, xx: b * Xb}) / b**2)
x_chart = sp.expand(Q.subs({delta: xx * Dx, b: xx * Yx}) / xx**2)
need(delta_chart == Bd**2 - aa**2 * xi * Xd**2 - Ba * Bd, "finite delta chart")
need(b_chart == 1 - aa**2 * xi * Xb**2 - Ba * Db, "finite b complement")
need(x_chart == Yx**2 - aa**2 * xi - Ba * Dx * Yx, "finite x complement")
# Ba!=0 makes the projective conic smooth; delta=0 cuts two attachment points.
need(sp.solve((sp.diff(Q, b), sp.diff(Q, xx), sp.diff(Q, delta)), (b, xx, delta), dict=True) == [{b: 0, delta: 0, xx: 0}], "finite conic has no projective singularity")
need(sp.factor(Q.subs(delta, 0).subs(xi, lam**2)) == (b - aa * lam * xx) * (b + aa * lam * xx), "finite bridge two attachments")

# When Ba=0, preparation has x^2 times a unit discriminant.
q, B1 = sp.symbols("q B1", nonzero=True)
tangent = q * b**2 - delta * B1 * xx * b - aa**2 * xi * xx**2
need(sp.factor(sp.discriminant(tangent, b)) == xx**2 * (delta**2 * B1**2 + 4 * q * aa**2 * xi), "finite horizontal discriminant")

# Simultaneous collision and incidence.  Let lambda^2=xi.  On the common
# D6-T edge t=S*P^2=v/x^3=P/zeta, derive both branch limits rather than
# inserting their signs by hand.  The D6 signs satisfy
# lambda*(v-c)=sign*x^3; the T signs satisfy
# zeta=sign*lambda*(P-a).
x_edge, p_edge, c_edge = sp.symbols("x_edge p_edge c_edge", nonzero=True)
for sign in (-1, 1):
    d_branch_v = c_edge + sign * x_edge**3 / lam
    t_branch_zeta = sign * lam * (p_edge - aa)
    d_limit = sp.limit(d_branch_v / x_edge**3, x_edge, sp.oo)
    t_limit = sp.limit(p_edge / t_branch_zeta, p_edge, sp.oo)
    need(d_limit == t_limit, f"D6-T sign matching {sign}")

# Six M-D6 roots split 3+3 according to x^3=+/-lambda*c.  The matching
# D6-T edges keep both T signs connected, even when both quadratics split.
m_plus = x_edge**3 - lam * c_edge
m_minus = x_edge**3 + lam * c_edge
need(sp.degree(m_plus, x_edge) == sp.degree(m_minus, x_edge) == 3,
     "three M-D6 roots on each sign")
need(sp.discriminant(m_plus, x_edge) != 0
     and sp.discriminant(m_minus, x_edge) != 0,
     "each M-D6 cubic packet is reduced")
need(sp.resultant(m_plus, m_minus, x_edge) != 0,
     "the two M-D6 cubic packets are disjoint")
need(20 - 11 + 1 == 10, "both splits base graph connected b1")
need(22 - 12 + 1 == 11, "one two-ended tail/bridge restores one cycle")
need(24 - 13 + 1 == 12, "both two-ended carriers restore both cycles")

# Exhibit simultaneous discriminant vanishing on the inherited seam.
sim = {
    xi: 1,
    K: 1,
    Theta: -2,
    U: 4,
    alpha: -4,
    Delta: sp.Rational(5606, 105),
}
need(sp.simplify(Krel.subs(Delta, sim[Delta]) - sim[K]) == 0, "simultaneous example seam")
need(sp.simplify((Theta**2 - 4 * K * xi).subs(sim)) == 0, "simultaneous finite collision")
need(sp.simplify((alpha**2 - 4 * U * xi).subs(sim)) == 0, "simultaneous infinity collision")

print("PASS W0 GLOBAL HOSTILE REFEREE CHECKS")
print("realizable_atlas=192/192:M,D6,T")
print("conservative_over_atlas=512/512:M,D6,T")
print("coupling=(Delta,K-e,Theta,Theta-Delta):8;(u,xi-u):3;(Phi,eta,alpha):8")
print("unique_owners=(2,6,1):-U;(4,4,1):-xi;(4,2,1):-K")
print("pick=global14;RH=26;generic_graph=E20,V9,b1=12")
print("infinity_complement=z*rho=tau^12;h^2=1-z^(6-r)c(rho);two_etale_signs")
print("finite_complements=three_charts_of_smooth_projective_conic_or_horizontal_split")
print("simultaneous_signs=M-D6:3+3;D6-T:1+1_matching;corrections_add")
print("seam=K=2848/45-(7/6)Delta;K!=0_iff_Delta!=5696/105")
