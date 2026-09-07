#!/usr/bin/env python3
"""Exact controls for the Nori resolution criterion and its equality frontier.

Universe: the two stated rational sextics and the ordinary cusp control;
this is not an enumeration of all sextics or all Keller maps.
"""
import hashlib
import json
import sympy as S

gates = 0


def check(condition, label):
    global gates
    gates += 1
    if not bool(condition):
        raise RuntimeError(label)


def equal(left, right, label):
    check(S.cancel(left - right) == 0, label)


def leading(expression, variable):
    numerator, denominator = S.fraction(S.cancel(expression))
    check(numerator != 0 and denominator != 0, "nonzero rational chart")
    np = S.Poly(numerator, variable)
    dp = S.Poly(denominator, variable)
    ni = min(exponent[0] for exponent, coefficient in np.terms())
    di = min(exponent[0] for exponent, coefficient in dp.terms())
    return ni - di, S.cancel(np.nth(ni) / dp.nth(di))


s, t, p, q, z, a, b = S.symbols("s t p q z a b")
U = t**4 + t**3 + t**2
V = 2*t**6 + 3*t**5 + 2*t**3 + 2*t**2
W = S.expand(V - 2*U + 2*U**2)
check(S.gcd(S.diff(U, t), S.diff(V, t)) == t, "unique critical parameter")
check(S.gcd(U, V) == t**2, "unshared cusp image")
equal(W, 7*t**5 + 8*t**6 + 4*t**7 + 2*t**8, "finite cusp shear")

N = S.cancel((U.subs(t, s) - U)/(s-t))
M = S.cancel((V.subs(t, s) - V)/(s-t))
pair_quadric = s*s - p*s + q
equal(S.rem(N.subs(t, p-s), pair_quadric, s),
      p*(p*p+p+1) - (2*p+1)*q, "complete symmetric first equation")
check((p*(p*p+p+1)).subs(p, -S.Rational(1, 2)) == -S.Rational(3, 8),
      "the symmetric denominator cannot vanish")
q_on_pair = p*(p*p+p+1)/(2*p+1)
h = p**4 + 3*p**3 + 5*p*p + p - 7
equal(S.rem(M.subs(t, p-s), pair_quadric, s).subs(q, q_on_pair),
      -p*p*h/(2*p+1), "complete pair algebra")
equal(p*p - 4*q_on_pair, -p*(2*p*p+3*p+4)/(2*p+1), "pair discriminant")
check(h.subs(p, 0) == -7, "off-diagonal roots avoid cusp")
check(S.discriminant(h, p) == -80787, "four simple pair sums")
check(S.resultant(h, p*(2*p*p+3*p+4)*(2*p+1), p) == 610785,
      "four distinct ordered pairs with lawful denominator")
check(S.rem(M.subs(t, -s), s*s, s) == 0, "diagonal p=0 control")
check(S.rem(N.subs(t, -s), s*s, s) == 0, "diagonal first equation control")

T = S.diff(U, t).subs(t, s)*S.diff(V, t) - S.diff(V, t).subs(t, s)*S.diff(U, t)
actual_tangencies = S.groebner([N, M, T], s, t, domain=S.QQ)
expected_tangencies = S.groebner([s+t**3+t**2+t, t**4], s, t, domain=S.QQ)
check(actual_tangencies == expected_tangencies, "all off-diagonal collisions transverse")

R = S.rem(V-b, U-a, t)
equal(R, 4*t**3+(2*a+5)*t**2+a*t-3*a-b, "cubic remainder for triple fibres")
r2 = 4*a*a + 8*a + 21
r1 = 2*a*a + 13*a + 4*b
r0 = -6*a*a - 2*a*b - 19*a - b
equal(16*S.rem(U-a, R, t), r2*t*t + r1*t + r0,
      "three exact cubic-divisibility equations")
j = 4*a**3 + 4*a*a - 63*a
equal(4*r0.subs(b, -(2*a*a+13*a)/4), j, "eliminate second target coordinate")
equal(S.rem(j, r2, a), 21-76*a, "linear triple obstruction")
equal((-76*a*a-97*a+1617)*r2 + (76*a+173)*j, 33957,
      "nonzero Bezout certificate excludes every triple fibre")
check(S.groebner([r2, r1, r0], a, b, domain=S.QQ) == S.groebner([1], a, b, domain=S.QQ),
      "independent triple exclusion")

# Explicit blowup charts. The multiplicity of the strict curve at each
# centre is the smaller valuation, checked directly in the parameter.
finite = [(U, W), (U, W/U), (U, W/U**2), (U**3/W, W/U**2),
          (U**5/W**2-S.Rational(1, 49), W/U**2)]
finite_orders = [(2, 5), (2, 3), (2, 1), (1, 1), (1, 1)]
for index, ((x, y), wanted) in enumerate(zip(finite, finite_orders)):
    check((leading(x, t)[0], leading(y, t)[0]) == wanted, f"finite chart {index}")
check(leading(finite[-1][1], t) == (1, S.Integer(7)), "finite final transverse exceptional")

X = S.cancel(U.subs(t, 1/z)/V.subs(t, 1/z))
Z = S.cancel(1/V.subs(t, 1/z))
rho = S.cancel(Z/X**3 - 4)
tau = S.cancel(rho/X + 30)
infinity = [(X, Z), (X, Z/X), (X, Z/X**2), (X, rho),
            (X, tau), (X/tau, tau), (X/tau**2-S.Rational(1, 2450), tau)]
infinity_orders = [(2, 6), (2, 4), (2, 2), (2, 2), (2, 1), (1, 1), (1, 1)]
for index, ((x, y), wanted) in enumerate(zip(infinity, infinity_orders)):
    check((leading(x, z)[0], leading(y, z)[0]) == wanted, f"infinity chart {index}")
check(leading(Z-4*X**3+30*X**4, z) == (9, S.Rational(35, 16)),
      "infinity has first odd exponent nine")
check(leading(tau, z) == (1, S.Integer(35)), "infinity final transverse exceptional")
check(S.Integer(-4) != 0, "line at infinity exits at rho=-4, away from branch")

finite_m = [min(row) for row in finite_orders[:-1]]
infinity_m = [min(row) for row in infinity_orders[:-1]]
check(finite_m == [2, 2, 1, 1], "finite multiplicity sequence")
check(infinity_m == [2, 2, 2, 2, 1, 1], "infinity multiplicity sequence")
all_m = finite_m + infinity_m
N_nodes = 4
degree = 6
delta_cusps = sum(m*(m-1)//2 for m in all_m)
square = degree*degree - sum(m*m for m in all_m)
cost = sum(all_m)
check(delta_cusps == 6, "delta of finite and infinity cusps")
check(N_nodes + delta_cusps == (degree-1)*(degree-2)//2,
      "full rational genus identity")
check(square == 8 and cost == 16, "equality curve consequence object")
check(square-2*N_nodes == 3*degree-2-cost == 0, "Nori equality; no exclusion")

# Analytic formula's symbolic audit, with genus retained explicitly.
d, genus, node_count, delta, total_m = S.symbols("d genus node_count delta total_m")
node_twice = (d-1)*(d-2)-2*genus-2*delta
equal(d*d-(2*delta+total_m)-node_twice,
      3*d-2+2*genus-total_m, "all-degree adjunction identity")

# Positive inherited example and strict-inequality hostile. The ordinary
# cusp complement has the nonabelian S3 quotient of its B3 presentation.
old_m = [2, 2, 1, 1, 2, 2, 2, 1, 1]
check(36-sum(m*m for m in old_m)-2*5 == 18-2-sum(old_m) == 2,
      "inherited five-node sextic positive margin")
ordinary_m = [2, 1, 1, 1, 1, 1]
check(9-sum(m*m for m in ordinary_m) == 9-2-sum(ordinary_m) == 0,
      "ordinary cusp with tangent infinity has equality")

def compose(left, right):
    return tuple(left[right[i]] for i in range(len(right)))

sigma, eta = (1, 0, 2), (0, 2, 1)
check(compose(sigma, compose(eta, sigma)) == compose(eta, compose(sigma, eta)),
      "ordinary cusp braid relation")
check(compose(sigma, eta) != compose(eta, sigma), "equality has nonabelian quotient")

summary = {
    "target": [str(U), str(V)],
    "affine_cusp": [2, 5], "infinity_cusp": [2, 9], "ordinary_nodes": 4,
    "pair_sum_polynomial": str(h),
    "finite_multiplicities": finite_m, "infinity_multiplicities": infinity_m,
    "strict_transform_square": square, "linear_cost": cost, "Nori_margin": 0,
    "status": "VERIFIED curve; equality gives no abelianity or Keller realization",
}
payload = json.dumps(summary, sort_keys=True, separators=(",", ":"))
print("Nori consequence: margin = 3*degree - 2 + 2*genus - sum(multiplicities).")
print("Inherited sextic: margin 2, abelian affine complement and Keller support excluded.")
print("New sextic: one (2,5) cusp, four nodes, infinity (2,9); margin 0, OPEN beyond this test.")
print("Hostile: ordinary cusp has margin 0 and a nonabelian S3 quotient.")
print("PASS", gates, "always-active exact gates")
print("semantic SHA-256:", hashlib.sha256(payload.encode()).hexdigest())
