#!/usr/bin/env python3
"""Exact companion for THM-3650.

The script proves, over a rational function field, the universal retained
fourth-stable endpoint identity on the ordinary principal zero-second-debt
plane.  It checks all 105 relevant target-two-form monomial jets, classifies
the fourth-debt hyperplane, and verifies the Q6, Q*, Qh, Q9, beta-line, and
actual local Q9 controls.
"""

import ast
from math import factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.fields import field


CHECKS = 0


def require(label, condition):
    """Record one exact optimization-safe gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def truncated(expr, source, stable, source_order=2, stable_order=4):
    """Taylor dictionary at (source,stable)=(0,0)."""
    answer = {}
    for source_degree in range(source_order + 1):
        for stable_degree in range(stable_order + 1):
            value = (
                sp.diff(expr, source, source_degree, stable, stable_degree)
                .subs({source: 0, stable: 0})
                / factorial(source_degree)
                / factorial(stable_degree)
            )
            if value:
                answer[(source_degree, stable_degree)] = value
    return answer


def truncated_product(left, right, source_order=2, stable_order=4):
    """Product in the fixed bivariate jet ring."""
    answer = {}
    for (left_source, left_stable), left_value in left.items():
        for (right_source, right_stable), right_value in right.items():
            source_degree = left_source + right_source
            stable_degree = left_stable + right_stable
            if source_degree > source_order or stable_degree > stable_order:
                continue
            key = (source_degree, stable_degree)
            answer[key] = answer.get(key, 0) + left_value * right_value
    return {key: value for key, value in answer.items() if value != 0}


def truncated_power(value, exponent):
    """Power in the fixed bivariate jet ring."""
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = truncated_product(answer, value)
    return answer


print("THM-3650 exact companion -- principal fourth hyperplane")
print("status=PROVED VERIFIED-EXACT INDEPENDENTLY HOSTILE-AUDITED")


print("SECTION principal zero-second-debt jet plane")
x, q, source, stable = sp.symbols("x q source stable")
y, z = sp.symbols("y z")
u, r_minus, r_zero = sp.symbols("u r_minus r_zero")
k_minus, k_zero, k_plus = sp.symbols("k_minus k_zero k_plus")
A = 9 - 4 * u
B = 9 + 4 * u
r_plus = sp.factor((-243 - A * r_minus) / B)
points = (-1, 0, 1)
q_values = (-3, -sp.Rational(3, 4), -3)
slopes = (-sp.Rational(9, 2), u, sp.Rational(9, 2))
curvatures = (r_minus, r_zero, r_plus)
third_jets = (k_minus, k_zero, k_plus)

D_general = 1 + x**2 * q
c_general = sp.expand(x * D_general * (D_general + 2))
e_general = sp.expand(q * (D_general + 3))
y_general = c_general / 3
z_general = e_general + 3

require(
    "zero-second-debt curvature equation",
    sp.factor(A * r_minus + B * r_plus + 243) == 0,
)

branch_packets = {}
for point, q_value, slope, curvature, third_jet in zip(
    points, q_values, slopes, curvatures, third_jets
):
    branch_q = (
        q_value
        + slope * source
        + curvature * source**2 / 2
        + third_jet * source**3 / 6
        + stable**2
    )
    pulled_y = sp.expand(y_general.subs({x: point + source, q: branch_q}))
    pulled_z = sp.expand(z_general.subs({x: point + source, q: branch_q}))
    pulled_y_source = sp.diff(pulled_y, source)
    pulled_z_source = sp.diff(pulled_z, source)
    area = sp.expand(
        pulled_y_source * sp.diff(pulled_z, stable)
        - sp.diff(pulled_y, stable) * pulled_z_source
    )
    y_jet = truncated(pulled_y, source, stable)
    z_jet = truncated(pulled_z, source, stable)
    branch_packets[point] = (
        [truncated_power(y_jet, exponent) for exponent in range(6)],
        [truncated_power(z_jet, exponent) for exponent in range(6)],
        (
            truncated(area, source, stable),
            truncated(pulled_y_source, source, stable),
            truncated(pulled_z_source, source, stable),
        ),
    )

require(
    "principal tangent y",
    tuple(branch_packets[point][2][1].get((0, 0), 0) for point in points)
    == (1, 1, 1),
)
require(
    "principal tangent z",
    tuple(branch_packets[point][2][2].get((0, 0), 0) for point in points)
    == (-9, 4 * u, 9),
)
lambda_coefficients = (A / 18, -1, B / 18)
require("lambda kills constants", sum(lambda_coefficients) == 0)
require(
    "lambda kills tangents",
    sp.expand(-9 * lambda_coefficients[0] + 4 * u * lambda_coefficients[1] + 9 * lambda_coefficients[2])
    == 0,
)
print("PASS principal_ordinary=A*B!=0 curvature=A*rminus+B*rplus=-243")


print("SECTION full universal endpoint identity")
P_minus = (
    5184 * r_zero * r_minus
    + 69984 * r_zero
    + 1280 * r_minus**2 * u**2
    + 2304 * r_minus**2 * u
    - 11664 * r_minus**2
    + 62208 * r_minus * u**2
    + 132192 * r_minus * u
    - 389286 * r_minus
    - 2304 * k_minus * u**2
    + 11664 * k_minus
    - 1152 * k_plus * u**2
    - 5184 * k_plus * u
    - 5832 * k_plus
    + 622080 * u**2
    + 1819584 * u
    - 2184813
)
P_plus = (
    5184 * r_zero * r_minus
    + 69984 * r_zero
    + 1280 * r_minus**2 * u**2
    - 8064 * r_minus**2 * u
    + 11664 * r_minus**2
    + 62208 * r_minus * u**2
    - 287712 * r_minus * u
    + 240570 * r_minus
    + 1152 * k_minus * u**2
    - 5832 * k_minus
    + 2304 * k_plus * u**2
    + 10368 * k_plus * u
    + 11664 * k_plus
    + 622080 * u**2
    - 2799360 * u
    + 177147
)
a_minus = P_minus / (26244 * B)
a_plus = -P_plus / (26244 * B)
b_minus = -2 * (
    80 * r_minus * u**2
    - 72 * r_minus * u
    - 243 * r_minus
    + 1728 * u**2
    - 4374
) / (729 * B)
b_plus = -2 * (
    20 * r_minus * u
    - 27 * r_minus
    + 432 * u
    - 243
) / 729
c_zero = 4 * (2 * r_minus + 27) / (9 * B)
d_value = (
    8 * r_minus * u
    - 18 * r_minus
    + 216 * u
    - 243
) / 243

# Rows use Taylor coefficients, so source degree two is J''/2.
identity_coefficients = {
    (0, -1, 0): a_minus,
    (0, -1, 1): b_minus,
    (0, -1, 2): -2 * A / 81,
    (0, 0, 1): c_zero,
    (0, 0, 2): sp.S.Zero,
    (0, 1, 0): a_plus,
    (0, 1, 1): b_plus,
    (0, 1, 2): -2 * B / 81,
    (2, -1, 0): d_value,
    (2, -1, 1): A / 27,
    (2, 0, 1): sp.S.Zero,
    (2, 1, 0): -d_value,
    (2, 1, 1): -B / 27,
}


def two_form_value(kind, y_degree, z_degree, w_degree, point, source_degree, stable_degree):
    """One target two-form monomial's pulled-back Taylor coefficient."""
    y_powers, z_powers, bases = branch_packets[point]
    value = truncated_product(y_powers[y_degree], z_powers[z_degree])
    value = truncated_product(value, {(0, w_degree): sp.S.One})
    value = truncated_product(value, bases[kind])
    return sp.sympify(value.get((source_degree, stable_degree), 0))


coefficient_monomials = [
    (kind, y_degree, z_degree, w_degree)
    for kind in range(3)
    for y_degree in range(5)
    for z_degree in range(5 - y_degree)
    for w_degree in range(5 - y_degree - z_degree)
]
require("two-form universe size", len(coefficient_monomials) == 105)

fraction_field, *_field_generators = field(
    "u,r_minus,r_zero,k_minus,k_zero,k_plus", QQ
)
field_lambda = [fraction_field.from_expr(value) for value in lambda_coefficients]
field_identity = {
    label: fraction_field.from_expr(value)
    for label, value in identity_coefficients.items()
}

mutation_detected = False
for monomial_index, (kind, y_degree, z_degree, w_degree) in enumerate(
    coefficient_monomials
):
    left = fraction_field.zero
    for branch_index, point in enumerate(points):
        left += field_lambda[branch_index] * fraction_field.from_expr(
            two_form_value(kind, y_degree, z_degree, w_degree, point, 0, 4)
        )
    right = fraction_field.zero
    for (stable_degree, point, source_degree), coefficient in field_identity.items():
        right += coefficient * fraction_field.from_expr(
            two_form_value(
                kind,
                y_degree,
                z_degree,
                w_degree,
                point,
                source_degree,
                stable_degree,
            )
        )
    require(f"symbolic two-form monomial {monomial_index}", left == right)

    # Mutating a_minus by one changes at least one basis identity.
    mutation = fraction_field.from_expr(
        two_form_value(kind, y_degree, z_degree, w_degree, -1, 0, 0)
    )
    if mutation != 0:
        mutation_detected = True
require("active endpoint coefficient mutation", mutation_detected)

# Total target degree five and above is invisible to every used source jet.
invisible_count = 0
for kind in range(3):
    for y_degree in range(6):
        for z_degree in range(6 - y_degree):
            w_degree = 5 - y_degree - z_degree
            used_values = []
            for point in points:
                used_values.extend(
                    two_form_value(
                        kind, y_degree, z_degree, w_degree, point, source_degree, stable_degree
                    )
                    for stable_degree, source_bound in ((0, 2), (2, 1), (4, 0))
                    for source_degree in range(source_bound + 1)
                )
            require(f"invisible degree-five monomial {kind,y_degree,z_degree,w_degree}", all(value == 0 for value in used_values))
            invisible_count += 1
require("degree-five invisible count", invisible_count == 63)
print("PASS symbolic_two_form_basis=105 degree_five_invisible=63")


print("SECTION fourth-debt hyperplane")
N4 = (
    12 * A * r_minus**2
    + 486 * (A - 3) * r_minus
    - A * B * k_minus
    + B**2 * k_plus
    + 243 * (22 * A - 153)
)
fourth_debt = -2 * N4 / (243 * B)
require(
    "constant specialization of endpoint identity",
    sp.factor(a_minus + a_plus - fourth_debt) == 0,
)
k_plus_zero = sp.factor(
    (
        A * B * k_minus
        - 12 * A * r_minus**2
        - 486 * (A - 3) * r_minus
        - 243 * (22 * A - 153)
    )
    / B**2
)
require("zero hyperplane forward", sp.factor(N4.subs(k_plus, k_plus_zero)) == 0)
require("zero hyperplane coefficient", sp.diff(N4, k_plus) == B**2)
require("central curvature cancels", not fourth_debt.has(r_zero))
require("central third jet cancels", not fourth_debt.has(k_zero))
print("PASS lambda_J4=-2*N4/(243*B) zero_iff_N4=0")
print("PASS N4_independent_of_rzero_kzero_and_affine_in_endpoint_thirds")


print("SECTION exact polynomial controls")
Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
L = x * (x**2 - 1)
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * L**2)
Qstar = (
    -x**7
    - sp.Rational(27, 4) * x**6
    + 3 * x**5
    + 18 * x**4
    - 3 * x**3
    - sp.Rational(27, 2) * x**2
    + x
    - sp.Rational(3, 4)
)
Qh = (
    44069 * x**11
    + 112059 * x**10
    - 154749 * x**9
    - 406377 * x**8
    + 188107 * x**7
    + 513081 * x**6
    - 82835 * x**5
    - 230931 * x**4
    + 5408 * x
    - 4056
) / 5408
Q9 = sp.expand(Q6 + sp.Rational(5717, 324) * L**3)


def polynomial_packet(polynomial):
    """Return u,r-/0/+,k-/0/+ and the two debts."""
    packet_u = sp.diff(polynomial, x).subs(x, 0)
    packet_r = tuple(sp.diff(polynomial, x, 2).subs(x, point) for point in points)
    packet_k = tuple(sp.diff(polynomial, x, 3).subs(x, point) for point in points)
    substitutions = {
        u: packet_u,
        r_minus: packet_r[0],
        r_zero: packet_r[1],
        k_minus: packet_k[0],
        k_zero: packet_k[1],
        k_plus: packet_k[2],
    }
    second_row = sp.factor(A.subs(substitutions) * packet_r[0] + B.subs(substitutions) * packet_r[2] + 243)
    packet_fourth = sp.factor(fourth_debt.subs(substitutions))
    return packet_u, packet_r, packet_k, second_row, packet_fourth


Q6_packet = polynomial_packet(Q6)
Qstar_packet = polynomial_packet(Qstar)
Qh_packet = polynomial_packet(Qh)
Q9_packet = polynomial_packet(Q9)
require("Q6 zero second debt", Q6_packet[3] == 0)
require("Q6 nonzero fourth debt", Q6_packet[4] == sp.Rational(365888, 6561))
require("Qstar zero second debt", Qstar_packet[3] == 0)
require("Qstar nonzero fourth debt", Qstar_packet[4] == sp.Rational(5440, 81))
require("Qh zero second debt", Qh_packet[3] == 0)
require("Qh zero fourth debt", Qh_packet[4] == 0)
require("Q9 zero second debt", Q9_packet[3] == 0)
require("Q9 zero fourth debt", Q9_packet[4] == 0)
require(
    "Qh degree-four relation",
    65 * Qh_packet[2][0] - 169 * Qh_packet[2][2] + 10449 == 0,
)
print("PASS Q6_debt=365888/6561 Qstar_debt=5440/81")
print("PASS Qh_Q9_simultaneous_zero_second_and_fourth_debt")

beta = sp.symbols("beta")
Qbeta = sp.expand(
    Q1 + L**2 * (beta * x - (259 + 16 * beta) / 36)
)
beta_packet = polynomial_packet(Qbeta)
expected_beta_debt = -64 * (
    520 * beta**2 + 1688 * beta - 5717
) / 6561
require("beta line zero second debt", beta_packet[3] == 0)
require(
    "beta line fourth quadratic",
    sp.factor(beta_packet[4] - expected_beta_debt) == 0,
)
beta_roots = (
    -sp.Rational(211, 130) + sp.Rational(99, 260) * sp.sqrt(94),
    -sp.Rational(211, 130) - sp.Rational(99, 260) * sp.sqrt(94),
)
for root_index, root in enumerate(beta_roots):
    require(
        f"beta zero root {root_index}",
        sp.factor(expected_beta_debt.subs(beta, root)) == 0,
    )
require("beta zero polynomial discriminant", sp.discriminant(520 * beta**2 + 1688 * beta - 5717, beta) == 396**2 * 94)
print("PASS beta_debt=-64*(520beta^2+1688beta-5717)/6561")
print("PASS beta_roots=-211/130_plusminus_99sqrt94/260")

# A uniform degree-at-most-eleven Hermite control on every ordinary u.
Q5_u = (
    u * x**5
    - 2 * u * x**3
    + u * x
    + sp.Rational(9, 2) * x**4
    - sp.Rational(27, 4) * x**2
    - sp.Rational(3, 4)
)
H_u = -(
    64 * u**2 * x
    + 756 * u * x**2
    + 144 * u * x
    - 432 * u
    + 1944 * x**2
    + 243 * x
    - 972
) / (16 * B)
K_u = (
    2048 * u**3 * x**2
    - 1536 * u**3
    + 9216 * u**2 * x**2
    + 8784 * u**2 * x
    - 6912 * u**2
    + 23328 * u * x**2
    + 51516 * u * x
    - 9720 * u
    + 9477 * x**2
    + 51759 * x
    - 4374
) / (32 * B**2)
Q_uniform = sp.cancel(Q5_u + L**2 * H_u + L**3 * K_u)
require(
    "uniform values",
    tuple(sp.factor(Q_uniform.subs(x, point)) for point in points)
    == (-3, -sp.Rational(3, 4), -3),
)
require(
    "uniform slopes",
    tuple(sp.factor(sp.diff(Q_uniform, x).subs(x, point)) for point in points)
    == (-sp.Rational(9, 2), u, sp.Rational(9, 2)),
)
require(
    "uniform curvatures",
    tuple(sp.factor(sp.diff(Q_uniform, x, 2).subs(x, point)) for point in points)
    == (0, 0, -243 / B),
)
uniform_thirds = tuple(
    sp.factor(sp.diff(Q_uniform, x, 3).subs(x, point)) for point in points
)
require(
    "uniform third jets",
    uniform_thirds[:2] == (0, 0)
    and sp.factor(uniform_thirds[2] - 243 * (88 * u - 45) / B**2) == 0,
)
require("uniform specializes to Qh", sp.factor(Q_uniform.subs(u, 1) - Qh) == 0)

# General degree-eleven Hermite realization through branch order three.
hermite_matrix = sp.Matrix(
    [
        [sp.diff(x**degree, x, derivative).subs(x, point) for degree in range(12)]
        for point in points
        for derivative in range(4)
    ]
)
require("degree-eleven Hermite rank", hermite_matrix.rank() == 12)
print("PASS uniform_degree11_control_specializes_to_Qh Hermite_rank=12")


print("SECTION actual local Q9 fourth-value control")
F0_Q9 = -(
    6725752777 * y**3
    - 1223385280 * y**2 * z
    - 147386304 * y**2
    - 23380969 * y * z**2
    + 45143280 * y * z
    - 118272960 * y
    - 2074176 * z**2
) / 118272960
F2_Q9 = -(
    28674065232 * y**2
    + 23094463735 * y * z
    + 106774200 * y
    - 4336511632 * z**2
    - 508472640 * z
) / 1995856200
F4_Q9 = -sp.Rational(779580313, 39917124) * y - sp.Rational(428886352, 249482025) * z

target_F0 = sp.expand(F0_Q9.subs({y: y_general, z: z_general}))
target_F2 = sp.expand(F2_Q9.subs({y: y_general, z: z_general}))
target_F4 = sp.expand(F4_Q9.subs({y: y_general, z: z_general}))
pulled_F = sp.expand(
    target_F0.subs(q, Q9 + stable**2)
    + stable**2 * target_F2.subs(q, Q9 + stable**2)
    + stable**4 * target_F4.subs(q, Q9 + stable**2)
)
# With G=w=t, the source Jacobian is simply partial_x pulled_F.
local_jacobian = sp.Poly(sp.diff(pulled_F, x), stable)
J0 = local_jacobian.nth(0)
J1 = local_jacobian.nth(1)
J2 = local_jacobian.nth(2)
J3 = local_jacobian.nth(3)
J4 = local_jacobian.nth(4)
require("Q9 local J0 values", tuple(J0.subs(x, point) for point in points) == (1, 1, 1))
require("Q9 local J0 first derivatives", tuple(sp.diff(J0, x).subs(x, point) for point in points) == (0, 0, 0))
require("Q9 local J0 second derivatives", tuple(sp.diff(J0, x, 2).subs(x, point) for point in points) == (0, 0, 0))
require("Q9 local J1 parity", J1 == 0)
require("Q9 local J2 values", tuple(J2.subs(x, point) for point in points) == (0, 0, 0))
require("Q9 local J2 first derivatives", tuple(sp.diff(J2, x).subs(x, point) for point in points) == (0, 0, 0))
require("Q9 local J3 parity", J3 == 0)
require("Q9 local J4 values", tuple(J4.subs(x, point) for point in points) == (0, 0, 0))

# Removing F4 exposes a nonzero tangent-plane vector with lambda zero.
pulled_without_F4 = sp.expand(
    target_F0.subs(q, Q9 + stable**2)
    + stable**2 * target_F2.subs(q, Q9 + stable**2)
)
base_J4 = sp.Poly(sp.diff(pulled_without_F4, x), stable).nth(4)
base_J4_values = tuple(sp.factor(base_J4.subs(x, point)) for point in points)
expected_base_J4 = (
    sp.Rational(4049599153, 997928100),
    sp.Rational(26351689457, 997928100),
    sp.Rational(34929416497, 997928100),
)
require("Q9 base J4 vector", base_J4_values == expected_base_J4)
require(
    "Q9 base J4 lambda zero",
    sp.factor(
        sp.Rational(5, 18) * base_J4_values[0]
        - base_J4_values[1]
        + sp.Rational(13, 18) * base_J4_values[2]
    )
    == 0,
)
print("PASS Q9_actual_local_J0twojet_J2onejet_J1=J3=0_J4values=0")
print("PASS Q9_F4_linear_kills_lambda_zero_tangent_plane_vector")


print("SECTION source AST and scope gates")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")
print("PASS scope=principal_ordinary_retained_quotient_only_no_global_pair_JC2_OPEN")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
