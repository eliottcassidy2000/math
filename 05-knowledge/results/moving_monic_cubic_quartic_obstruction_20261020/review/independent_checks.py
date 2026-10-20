#!/usr/bin/env python3
"""Independent exact checks for the moving cubic/quartic obstruction.

No producer module is imported.  All failures use explicit exceptions, so
python -O preserves the checks.  The universal proof is audited separately.
"""

import json
import sympy as S

y, t, r, s = S.symbols("y t r s")
a, b, c, u, v = S.symbols("a b c u v")
ap, bp, cp, j, lam, x, B = S.symbols("ap bp cp j lam x B")
h = y**3 + a*y
k = y**4 + b*y**2 + c*y
H = s**3 - (b - 2*a)*s - c
I = 3*c**2 + 4*a*(b - S.Rational(2, 3)*a)**2
conductor = (
    y**6 + (a+b)*y**4 + 2*c*y**3
    + (a*a-a*b+b*b)*y*y + (a+b)*c*y
    + a*(a-b)**2 + c*c
)
checks = []


def eq(name, lhs, rhs=0):
    diff = S.cancel(S.expand(lhs-rhs))
    if diff != 0:
        raise RuntimeError(f"{name}: nonzero difference {diff}")
    checks.append(name)


def require(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)


def reduce_gaps(poly, hv=h, kv=k):
    """Reduce using the unique monic h^i k^j, 0 <= j < 3, of each nongap degree."""
    rem = S.expand(poly)
    gaps = 0
    representation = 0
    while rem != 0:
        d = S.degree(rem, y)
        lc = S.Poly(rem, y).LC()
        candidates = [
            (int((d-4*q)//3), q)
            for q in range(3)
            if d >= 4*q and (d-4*q) % 3 == 0
        ]
        if candidates:
            if len(candidates) != 1:
                raise RuntimeError("nonunique degree basis")
            p, q = candidates[0]
            basis = hv**p * kv**q
            representation += lc*basis
            rem = S.expand(rem-lc*basis)
        else:
            gaps += lc*y**d
            rem = S.expand(rem-lc*y**d)
    return S.expand(gaps), S.expand(representation)


# 1. Two unrelated constructions of the conductor.
F = -S.resultant(
    t**3+a*t-u, t**4+b*t*t+c*t-v, t
)
eq("implicit equation vanishes", F.subs({u: h, v: k}))
eq("implicit equation monic in v", S.Poly(F, v).LC(), 1)
require("implicit equation degree three in v", S.degree(F, v) == 3)
eq("trace-dual conductor identity",
   S.diff(F, v).subs({u: h, v: k}), S.diff(h, y)*conductor)
eq("collision-resultant conductor",
   S.resultant(H, y*y-s*y+s*s+a, s), conductor)
eq("derivative resultant equals 9I",
   S.resultant(S.diff(h, y), S.diff(k, y), y), 9*I)
eq("separation resultant equals 9I",
   S.resultant(H, 3*s*s+4*a, s), 9*I)
eq("conductor times one belongs to E",
   conductor, h*h+(b-a)*k+2*c*h+a*(a-b)**2+c*c)
eq("conductor times y belongs to E",
   conductor*y, h*k+c*k+(a-b)**2*h)
eq("conductor times y squared belongs to E",
   conductor*y*y, k*k+(a-b)*h*h-a*(a-b)*k+c*(a-b)*h)
for exponent in range(7):
    remainder, _ = reduce_gaps(conductor*y**exponent)
    eq(f"gap reduction conductor y^{exponent}", remainder)

# 2. Collision equations and actual transverse-node criterion.
pair = y*y-s*y+s*s+a
eq("cubic common-value identity",
   S.rem(h+c+(b-a)*s, pair, y), -H)
eq("quartic common-value identity",
   S.rem(k, pair, y),
   S.rem(k.subs(y, s-y), pair, y) + S.rem(-H*(2*y-s), pair, y))
eq("pair separation identity", s*s-4*(s*s+a), -3*s*s-4*a)
hp_r = 3*r*r+a
kp_r = 4*r**3+2*b*r+c
hp_t = 3*t*t+a
kp_t = 4*t**3+2*b*t+c
tangent = S.cancel((hp_r*kp_t-hp_t*kp_r)/(t-r))
tangent = S.expand(tangent.subs(t, s-r))
tangent = S.rem(tangent, r*r-s*r+s*s+a, r)
eq("collision tangent determinant identity",
   tangent.subs(c, s**3-(b-2*a)*s),
   (3*s*s+4*a)*S.diff(H, s))
eq("repeated cubic root factorization",
   H.subs({b: 2*a+3*t*t, c: -2*t**3}),
   (s-t)**2*(s+2*t))
eq("repeated root unit factorization",
   I.subs({b: 2*a+3*t*t, c: -2*t**3}),
   S.Rational(4, 9)*(a+3*t*t)*(4*a+3*t*t)**2)
eq("triple-image locus unit equation",
   I.subs(b, a), 3*c*c+S.Rational(4, 9)*a**3)

# 3. Independent coefficient/gap computation for both reference classes.
J = S.expand(ap*y*S.diff(k, y) - (bp*y*y+cp*y)*S.diff(h, y))
J_display = (
    (4*ap-3*bp)*y**4 - 3*cp*y**3
    + (2*b*ap-a*bp)*y*y + (c*ap-a*cp)*y
)
eq("standard reference Jacobian", J, J_display)
W = j*y-lam*S.integrate(J, y)
eq("standard leading gap five", S.expand(W).coeff(y, 5),
   -lam*(4*ap-3*bp)/5)
standard_rem, _ = reduce_gaps(W.subs(bp, S.Rational(4, 3)*ap))
eq("standard gap two",
   standard_rem.coeff(y, 2),
   -lam*((3*b-2*a)*cp+2*c*ap)/4)
require("only standard gaps one and two remain",
        set(S.Poly(standard_rem, y).monoms()).issubset({(1,), (2,)}))
eq("second differential equation is product derivative",
   ((3*b-2*a)*cp+2*c*ap).subs(b, S.Rational(4, 3)*a+B),
   (2*a+3*B)*cp+2*ap*c)
second_k = k.subs(b, S.Rational(2, 3)*a)
second_q = h*h+S.Rational(3, 2)*c*h
J_second = (
    S.diff(second_k, a)*ap*S.diff(second_q, y)
    - S.diff(second_k, y)*S.diff(second_q, a)*ap
)
T = (2*h+S.Rational(3, 2)*c)*(
    2*y**4+S.Rational(2, 3)*a*y*y+c*y
)
eq("second reference Jacobian", J_second, -ap*T)
second_rem, second_rep = reduce_gaps(
    S.integrate(T, y), hv=h, kv=second_k
)
eq("second primitive gap remainder", second_rem, c*c*y*y/4)
required_rem, _ = reduce_gaps(
    j*y-lam*S.integrate(J_second, y), hv=h, kv=second_k
)
eq("second required primitive impossible remainder",
   required_rem, j*y+lam*ap*c*c*y*y/4)
eq("second-class unit condition",
   I.subs(b, S.Rational(2, 3)*a), 3*c*c)
recovery = (
    second_q**2-second_k**3-2*a*second_k*second_q
    +S.Rational(2, 3)*a*a*second_k**2
    +(S.Rational(2, 27)*a**3+S.Rational(3, 4)*c*c)*second_q
    -(a**4/9+a*c*c)*second_k
)
eq("second-class entire recovery", recovery, c**3*h/8)

# 3b. Recheck the finite terminal-pair calculation used after the cited bound.
beta, gamma = S.symbols("beta gamma")
candidate_p = k+beta*h
candidate_q = h*h+gamma*h
cancelled = S.expand(candidate_q**2-candidate_p**3)
require("wrong-beta cancellation degree at most eleven",
        S.degree(cancelled, y) == 11)
eq("wrong-beta degree eleven coefficient", cancelled.coeff(y, 11), -3*beta)
rho = 4*a-3*b
cancelled = S.expand(candidate_q**2-k**3-rho*k*candidate_q)
require("wrong-gamma cancellation degree at most nine",
        S.degree(cancelled, y) == 9)
eq("wrong-gamma degree nine coefficient", cancelled.coeff(y, 9), 2*gamma-3*c)
eta = -2*a*(a-b)
cancelled = S.expand(cancelled.subs(gamma, S.Rational(3, 2)*c)-eta*k*k)
require("off-locus cancellation degree at most seven",
        S.degree(cancelled, y) == 7)
eq("off-locus degree seven coefficient", cancelled.coeff(y, 7),
   -c*(2*a-3*b)/2)
for label, degrees, codimension in [
    ("wrong beta", [0, 6, 11, 17], 7),
    ("wrong gamma", [0, 6, 9, 15], 6),
    ("off second locus", [0, 6, 7, 13], 5),
]:
    require(f"{label} basis residues distinct",
            sorted(d % 4 for d in degrees) == [0, 1, 2, 3])
    require(f"{label} exact gap count",
            sum((d-d % 4)//4 for d in degrees) == codimension)
nondividing_pairs = [
    (d, e) for d in [3, 4, 6, 7] for e in [3, 4, 6, 7]
    if d > e and d % e != 0
]
require("all bounded nondividing terminal degree pairs",
        nondividing_pairs == [(4, 3), (6, 4), (7, 3), (7, 4), (7, 6)])

# 4. Cotangent and primitive identities, including a generic-only hostile.
# Representatives differing by a multiple of F agree to first order on
# F=F_u=F_v=0.  This is checked modulo the formal singular-point ideal.
z, fu, fv, A, Bu, Av, Bv, Au = S.symbols("z fu fv A Bu Av Bv Au")
D11, D12, D21, D22 = S.symbols("D11 D12 D21 D22")
matrix = S.Matrix([
    [D11+Au*z+A*fu, D12+Av*z+A*fv],
    [D21+Bu*z+B*fu, D22+Bv*z+B*fv],
])
eq("cotangent determinant ignores implicit-relation representatives",
   matrix.det().subs({z: 0, fu: 0, fv: 0}),
   D11*D22-D12*D21)
P_ambient = u+x*v*v
Q_ambient = v+x*P_ambient*P_ambient
eq("nontrivial entire ambient family determinant",
   S.diff(P_ambient, u)*S.diff(Q_ambient, v)
   - S.diff(P_ambient, v)*S.diff(Q_ambient, u), 1)
alpha_u = S.expand(P_ambient*S.diff(Q_ambient, u))
alpha_v = S.expand(P_ambient*S.diff(Q_ambient, v)-u)
eq("difference Liouville form is closed",
   S.diff(alpha_v, u)-S.diff(alpha_u, v))
potential = S.integrate(alpha_u, u)
potential += S.integrate(S.expand(alpha_v-S.diff(potential, v)), v)
eq("difference Liouville primitive u", S.diff(potential, u), alpha_u)
eq("difference Liouville primitive v", S.diff(potential, v), alpha_v)
eq("generic-only hostile determinant", S.diff(x*u, u), x)
require("generic-only hostile loses rank at x=0", (x*u).subs(x, 0) == 0)

# 5. Hostiles: immersion is not enough to make the coefficients fixed.
moving_one = {a: x, b: S.Rational(2, 3)*x, c: 1}
eq("moving second-class family is everywhere immersed",
   I.subs(moving_one), 3)
moving_a = (1-3*x*x)/4
moving_two = {a: moving_a, b: S.Rational(2, 3)*moving_a+1, c: x}
eq("moving standard-only family is everywhere immersed",
   I.subs(moving_two), 1)
require("second moving family is not on second-class locus",
        S.expand((2*a-3*b).subs(moving_two)) == -3)
eq("moving standard-only leading obstruction is nonzero",
   4*S.diff(moving_a, x)-3*S.diff(moving_two[b], x), -3*x)

# 6. Boundary controls: node, triple point, cusp, repeated-pair-sum locus.
eq("three-node conductor", conductor.subs({a: -1, b: 0, c: 0}),
   (y*y-1)*(y**4+1))
eq("triple-point conductor", conductor.subs({a: 0, b: 0, c: 1}),
   (y**3+1)**2)
eq("cusp conductor", conductor.subs({a: 0, b: 0, c: 0}), y**6)
eq("cusp fails immersion", I.subs({a: 0, b: 0, c: 0}), 0)
repeated_values = {a: 1, b: 5, c: -2}
eq("repeated-pair-sum immersed control",
   I.subs(repeated_values), S.Rational(784, 9))
eq("repeated-pair-sum cubic",
   H.subs(repeated_values), (s-1)**2*(s+2))

print(json.dumps({
    "status": "PASS",
    "check_count": len(checks),
    "checks": checks,
    "scope": "Universal symbolic identities plus explicit hostile controls; "
             "the all-generating-pair argument and classical degree bound "
             "are mathematical proof inputs, not inferred from samples.",
}, indent=2))
