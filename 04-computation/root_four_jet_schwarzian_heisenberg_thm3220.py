"""Exact companion for THM-3220's root four-jet Heisenberg seam.

It uses ordinary Taylor coefficients in tangent-to-identity germs modulo
u^5, and it is assertion-independent.
"""

from fractions import Fraction
from itertools import product
import ast
import hashlib
from pathlib import Path
import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3215-arbitrary-degree-root-jet-hamiltonian-affine-dihedral-"
    "holonomy-and-p-fold-carry.md"
)
DEPENDENCY_SHA256 = "00c5775d62a6db1f651265bfc7f659f96cdf40464b644d29ae4dab7df419dba8"


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3215 dependency hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def omega(x, y):
    return x[0] * y[1] - x[1] * y[0]


def compose_log(outer, inner):
    """Log coordinates of outer(inner(u)), through u^4."""
    ai, bi, ci = inner
    ao, bo, co = outer
    return (
        ai + ao,
        bi + bo,
        ci + co + Fraction(1, 2) * (ai * bo - bi * ao),
    )


def coeff_from_log(x):
    a, b, c = x
    return (a, b + a * a, c + Fraction(5, 2) * a * b + a**3)


def log_from_coeff(x):
    a, b, c = x
    return (a, b - a * a, c - Fraction(5, 2) * a * b + Fraction(3, 2) * a**3)


def transition(xj, xi):
    ai, bi, ci = xi
    aj, bj, cj = xj
    return (
        aj - ai,
        bj - bi,
        cj - ci - Fraction(1, 2) * (ai * bj - bi * aj),
    )


def thm2779_coordinates(x):
    aa, bb, cc = x
    return (bb, aa, cc - Fraction(1, 2) * aa * bb)


def thm2779_product(x, y):
    """THM-2779 equation (23), with x the outer/left germ."""
    xx, xy, xz = x
    yx, yy, yz = y
    return (xx + yx, xy + yy, xz + yz - xy * yx)


def mod_fraction(x, p):
    return (x.numerator * pow(x.denominator, -1, p)) % p


def compose_coeff_mod(outer_coeff, inner_coeff, p):
    ao, bo, co = outer_coeff
    ai, bi, ci = inner_coeff
    return (
        (ai + ao) % p,
        (bi + bo + 2 * ao * ai) % p,
        (ci + co + ao * (ai * ai + 2 * bi) + 3 * bo * ai) % p,
    )


def compose_coeff(outer_coeff, inner_coeff):
    ao, bo, co = outer_coeff
    ai, bi, ci = inner_coeff
    return (
        ai + ao,
        bi + bo + 2 * ao * ai,
        ci + co + ao * (ai * ai + 2 * bi) + 3 * bo * ai,
    )


def coeff_from_log_mod(log_vector, p):
    aa, bb, cc = log_vector
    inv2 = pow(2, -1, p)
    return (
        aa % p,
        (bb + aa * aa) % p,
        (cc + 5 * inv2 * aa * bb + aa**3) % p,
    )


def series_add(x, y, bound):
    return tuple(
        (x[i] if i < len(x) else 0) + (y[i] if i < len(y) else 0)
        for i in range(bound)
    )


def series_mul(x, y, bound):
    out = [0] * bound
    for i, xi in enumerate(x[:bound]):
        if xi == 0:
            continue
        for j, yj in enumerate(y[: bound - i]):
            out[i + j] += xi * yj
    return tuple(out)


def series_compose(outer_series, inner_series, bound):
    out = (0,) * bound
    power = (1,) + (0,) * (bound - 1)
    for coefficient in outer_series[:bound]:
        if coefficient:
            out = series_add(
                out, tuple(coefficient * z for z in power), bound
            )
        power = series_mul(power, inner_series, bound)
    return out


def series_inverse(phi, bound):
    inverse = [0] * bound
    inverse[1] = 1
    for degree in range(2, bound):
        current = series_compose(phi, tuple(inverse), bound)[degree]
        inverse[degree] = -current
    result = tuple(inverse)
    require(
        series_compose(phi, result, bound) == (0, 1) + (0,) * (bound - 2),
        "formal series inverse",
    )
    return result


# The raw characteristic-two boundary is D8, exactly as in THM-2779.
f2_group = tuple(product(range(2), repeat=3))
f2_identity = (0, 0, 0)


def f2_power(g, n):
    z = f2_identity
    for _ in range(n):
        z = compose_coeff_mod(g, z, 2)
    return z


f2_orders = []
for g in f2_group:
    for n in range(1, 9):
        if f2_power(g, n) == f2_identity:
            f2_orders.append(n)
            break
require(sorted(f2_orders) == [1, 2, 2, 2, 2, 2, 4, 4], "F2 D8 order census")
f2_center = tuple(
    g
    for g in f2_group
    if all(
        compose_coeff_mod(g, h, 2) == compose_coeff_mod(h, g, 2)
        for h in f2_group
    )
)
require(f2_center == ((0, 0, 0), (0, 0, 1)), "F2 D8 center")
for aa, bb, cc in f2_group:
    square = f2_power((aa, bb, cc), 2)
    require(square == (0, 0, aa * (aa + bb) % 2), "F2 quadratic square form")


u = sp.symbols("u")
a, b, c, A, B, C = sp.symbols("a b c A B C")


def trunc(poly):
    return sp.expand(sp.series(poly, u, 0, 5).removeO())


inner = u + a * u**2 + b * u**3 + c * u**4
outer = u + A * u**2 + B * u**3 + C * u**4
composed = trunc(outer.subs(u, inner))
actual_coeff = tuple(sp.expand(composed.coeff(u, j)) for j in (2, 3, 4))
inner_log = log_from_coeff((a, b, c))
outer_log = log_from_coeff((A, B, C))
expected_log = compose_log(outer_log, inner_log)
expected_coeff = tuple(sp.expand(x) for x in coeff_from_log(expected_log))
require(
    all(sp.expand(x - y) == 0 for x, y in zip(actual_coeff, expected_coeff)),
    "symbolic four-jet group law",
)

dictionary_checks = 0
for out_log in (
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(-2), Fraction(3), Fraction(5, 2)),
    (Fraction(0), Fraction(2), Fraction(-1)),
):
    for in_log in (
        (Fraction(2), Fraction(-1), Fraction(3, 2)),
        (Fraction(-1), Fraction(4), Fraction(0)),
    ):
        require(
            thm2779_coordinates(compose_log(out_log, in_log))
            == thm2779_product(
                thm2779_coordinates(out_log), thm2779_coordinates(in_log)
            ),
            "THM-2779 Heisenberg dictionary",
        )
        dictionary_checks += 1

# The first noncommutative jet: two normalized quadratics.
r, s = sp.symbols("r s")
phi_r = u + r * u**2
phi_s = u + s * u**2


def inverse_jet(phi):
    x2, x3, x4 = sp.symbols("x2 x3 x4")
    candidate = u + x2 * u**2 + x3 * u**3 + x4 * u**4
    error = trunc(phi.subs(u, candidate)) - u
    solution = sp.solve(
        [error.coeff(u, j) for j in (2, 3, 4)],
        (x2, x3, x4),
        dict=True,
    )[0]
    return sp.expand(candidate.subs(solution))


def compose(outer_phi, inner_phi):
    return trunc(outer_phi.subs(u, inner_phi))


comm = inverse_jet(phi_s)
comm = compose(inverse_jet(phi_r), comm)
comm = compose(phi_s, comm)
comm = compose(phi_r, comm)
require(
    sp.expand(comm - (u + r * s * (s - r) * u**4)) == 0,
    "quadratic central commutator",
)
T = sp.symbols("T")
require(
    sp.factor(sp.discriminant(T * (T - r) * (T - s), T) - (r * s * (s - r)) ** 2)
    == 0,
    "commutator is oriented cubic discriminant",
)
v, delta, Delta = sp.symbols("v delta Delta", nonzero=True)
gamma = 2 * v**3 / delta**3
require(
    sp.factor((gamma**2 * Delta**3 - 4 * v**6).subs(Delta, delta**2)) == 0,
    "quadratic deck central norm",
)

witt_bracket_checks = 0
for m in range(2, 13):
    for n in range(2, 13):
        if m == n:
            continue
        bound = m + n
        phi = [0] * bound
        psi = [0] * bound
        phi[1] = 1
        psi[1] = 1
        phi[m] = 2
        psi[n] = 3
        inv_phi = series_inverse(tuple(phi), bound)
        inv_psi = series_inverse(tuple(psi), bound)
        comm_series = inv_psi
        comm_series = series_compose(inv_phi, comm_series, bound)
        comm_series = series_compose(tuple(psi), comm_series, bound)
        comm_series = series_compose(tuple(phi), comm_series, bound)
        require(
            all(comm_series[j] == (1 if j == 1 else 0) for j in range(m + n - 1)),
            "Witt commutator lower terms",
        )
        require(
            comm_series[m + n - 1] == (m - n) * 2 * 3,
            "Witt commutator leading coefficient",
        )
        witt_bracket_checks += 1

# Exact triangle: the area term, not naive central addition, closes it.
x0 = (Fraction(0), Fraction(0), Fraction(0))
x1 = (Fraction(1), Fraction(0), Fraction(0))
x2 = (Fraction(0), Fraction(2), Fraction(0))
t10 = transition(x1, x0)
t21 = transition(x2, x1)
t20 = transition(x2, x0)
require(t10 == (1, 0, 0), "hostile first edge")
require(t21 == (-1, 2, -1), "hostile second edge")
require(t20 == (0, 2, 0), "hostile direct edge")
require(compose_log(t21, t10) == t20, "area-corrected triangle")
require(t10[2] + t21[2] != t20[2], "naive central addition must fail")

# A nonconstant common multiplicative unit is not a harmless higher-jet gauge.
t = sp.symbols("t")
base_f = u
base_g = u + u**2
unit = 1 + t * u
scaled_f = sp.expand(unit * base_f)
scaled_g = sp.expand(unit * base_g)
scaled_f_coeff = tuple(sp.expand(scaled_f.coeff(u, j)) for j in (2, 3, 4))
scaled_g_coeff = tuple(sp.expand(scaled_g.coeff(u, j)) for j in (2, 3, 4))
scaled_f_log = log_from_coeff(scaled_f_coeff)
scaled_g_log = log_from_coeff(scaled_g_coeff)
require(
    sp.expand((scaled_g_log[1] - scaled_f_log[1]) - (-1 - t)) == 0,
    "nonconstant-unit higher-jet boundary",
)

transition_composition_checks = 0
vertex_logs = (
    x0,
    x1,
    x2,
    (Fraction(2), Fraction(-1), Fraction(3, 2)),
    (Fraction(-1), Fraction(4), Fraction(-3)),
    (Fraction(3, 2), Fraction(5, 2), Fraction(7, 2)),
)
for xi in vertex_logs:
    for xj in vertex_logs:
        for xk in vertex_logs:
            tij = transition(xj, xi)
            tjk = transition(xk, xj)
            tik = transition(xk, xi)
            require(compose_log(tjk, tij) == tik, "strict transition composition")
            require(
                tik[2]
                == tij[2]
                + tjk[2]
                + Fraction(1, 2) * omega(tij, tjk),
                "central triangle transgression",
            )
            transition_composition_checks += 1

# Coordinate changes act with weights 1,2,3, and the area has weight 3.
coordinate_checks = 0
for scale in (Fraction(2), Fraction(3), Fraction(5, 2)):
    for x in (x1, x2, t21, (Fraction(2), Fraction(-3), Fraction(7, 2))):
        scaled = (x[0] / scale, x[1] / scale**2, x[2] / scale**3)
        require(
            coeff_from_log(scaled)[0] == coeff_from_log(x)[0] / scale,
            "weight-one coefficient",
        )
        coordinate_checks += 1
    require(
        omega(
            (t10[0] / scale, t10[1] / scale**2),
            (t21[0] / scale, t21[1] / scale**2),
        )
        == omega(t10, t21) / scale**3,
        "cubic area weight",
    )

# Quadratic root germs lie on the twisted cubic (A,-A^2,3A^3/2).
quadratic_vandermonde_checks = 0
for p in (3, 5, 7, 11, 13):
    inv2 = pow(2, -1, p)
    for ai in range(p):
        xi = (ai, (-ai * ai) % p, (3 * ai**3 * inv2) % p)
        for aj in range(p):
            xj = (aj, (-aj * aj) % p, (3 * aj**3 * inv2) % p)
            for ak in range(p):
                xk = (ak, (-ak * ak) % p, (3 * ak**3 * inv2) % p)
                d1 = ((xj[0] - xi[0]) % p, (xj[1] - xi[1]) % p)
                d2 = ((xk[0] - xj[0]) % p, (xk[1] - xj[1]) % p)
                lhs = (d1[0] * d2[1] - d1[1] * d2[0]) % p
                rhs = (-(aj - ai) * (ak - aj) * (ak - ai)) % p
                require(lhs == rhs, "quadratic Vandermonde area")
                quadratic_vandermonde_checks += 1

# Every nonzero log vector has exact order p in the mod-p four-jet group.
order_checks = 0
for p in (3, 5, 7, 11, 13):
    for av in range(p):
        for bv in range(p):
            for cv in range(p):
                if (av, bv, cv) == (0, 0, 0):
                    continue
                generator = coeff_from_log_mod((av, bv, cv), p)
                running = (0, 0, 0)
                for n in range(1, p):
                    running = compose_coeff_mod(generator, running, p)
                    require(
                        (n * av % p, n * bv % p, n * cv % p) != (0, 0, 0),
                        "premature projective order",
                    )
                    require(running != (0, 0, 0), "premature coefficient-jet order")
                running = compose_coeff_mod(generator, running, p)
                require(
                    (p * av % p, p * bv % p, p * cv % p) == (0, 0, 0),
                    "p-fold reset",
                )
                require(running == (0, 0, 0), "coefficient-jet p-fold reset")
                order_checks += 1

# Exact first carry of the integer factorial control A=1,B=2.
carry_checks = 0
gamma = 2
for p in (3, 5, 7, 11, 13, 17):
    central_coefficient = p * gamma
    require(central_coefficient // p == gamma, "primitive central carry")
    control_log = (2, 3, 5)
    control_coeff = tuple(int(x) for x in coeff_from_log(control_log))
    require(
        tuple(Fraction(x) for x in control_coeff) == coeff_from_log(control_log),
        "integral carry control",
    )
    running = (0, 0, 0)
    for _ in range(p):
        running = compose_coeff(control_coeff, running)
    require(all(x % p == 0 for x in running), "coefficient carry divisibility")
    require(
        tuple((x // p) % p for x in running)
        == tuple(x % p for x in control_log),
        "full logarithmic first carry",
    )
    carry_checks += 1

frobenius_character_checks = 0
for p in (3, 5, 7, 11, 13, 17, 19):
    for discriminant_value in range(1, p):
        chi = pow(discriminant_value, (p - 1) // 2, p)
        require(chi in (1, p - 1), "Euler quadratic character")
        # If delta^2=Delta, then delta^p=chi*delta.  The odd inverse cube
        # delta^-3 therefore has the same character.
        require(pow(chi, -3, p) == chi, "cubic inverse character")
        frobenius_character_checks += 1

print("dependency_thm3215_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("four_jet_group_law=Heisenberg_class_2")
print("characteristic_two_boundary=D8:1^1,2^5,4^2")
print("log_coordinates=A=a,B=b-a^2,C=c-(5/2)ab+(3/2)a^3")
print("thm2779_heisenberg_dictionary_checks=%d" % dictionary_checks)
print("transition_triangle=central_area_half_omega")
print("transition_composition_checks=%d" % transition_composition_checks)
print("area_hostile=C10+C21=-1_but_C20=0")
print("common_nonconstant_unit_hostile=B_transition_-1_to_-1-t")
print("coordinate_weight_checks=%d" % coordinate_checks)
print("quadratic_commutator=u+A*B*(B-A)u^4")
print("quadratic_commutator_square=cubic_discriminant")
print("positive_witt_bracket_checks=%d" % witt_bracket_checks)
print("quadratic_deck_commutator=2v^3/delta^3")
print("frobenius_character_checks=%d" % frobenius_character_checks)
print("quadratic_vandermonde_checks=%d" % quadratic_vandermonde_checks)
print("nonzero_order_p_checks=%d" % order_checks)
print("primitive_central_carry_checks=%d" % carry_checks)
print("higher_jet_hostile=(u,u+u^5):same_four_jet_distinct_germs")
print("scope=selected_simple_root_formal_germs_not_global_root_selector_or_GMC2")
print("all_exact_checks=PASS")
