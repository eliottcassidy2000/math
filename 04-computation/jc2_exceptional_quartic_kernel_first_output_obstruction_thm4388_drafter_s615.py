"""Exact scout for the exceptional-quartic kernel first-output obstruction.

This is a session artifact, not a canonical theorem companion.  It extracts
the rank-18 potential module behind ker(L_0), freezes an explicit K[Z]-basis,
and checks the retained-triple obstruction to hitting the seminormal
derivative class.  The differentiated map is the fixed-x, chosen-target-
representative map; it is not identified with a descended moving-family
transgression in the sense of THM-4067.
"""

import contextlib
import hashlib
import io
import subprocess
import types
from pathlib import Path


CONDUCTOR_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py"
)
CONDUCTOR_CORE_MARKER = (
    "# The three divided-difference resultants give an intrinsic conductor formula."
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_conductor_core(path, name):
    if path.is_file():
        source = path.read_text()
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
    require(source.count(CONDUCTOR_CORE_MARKER) == 1, "unique conductor marker")
    source = source.split(CONDUCTOR_CORE_MARKER, 1)[0]
    module = types.ModuleType(name)
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    return module


certificate = load_conductor_core(CONDUCTOR_PATH, "quartic_kernel_s615")
sagbi = certificate.sagbi
field = certificate.field
ring = certificate.ring
X = certificate.X
zero = certificate.zero
one = certificate.one
Z = sagbi.Z
B = sagbi.B
C = sagbi.C
E = sagbi.E
Q = sagbi.Q
C_prime = C.diff(X)
E_prime = E.diff(X)

L = X * (X**2 - 1)
c = certificate.conductor
h172 = certificate.quotient
r = L * h172
require(c == L * r, "conductor factorization")
require((L.degree(), r.degree(), c.degree()) == (3, 175, 178), "factor degrees")

# At an ordinary node, distinct full target tangents imply distinct projected
# (C,E)-tangents unless the projection of the Russell surface ramifies.  For
# c^2 e=b(b+4), that ramification divisor is exactly b=-2.  A positive
# resultant certificate at THM-4034's good fibre shows that none of the 86
# node branches lies above it.
h172_mod = certificate.polynomial_mod(h172)
B_plus_2_mod = certificate.polynomial_mod(B + 2)
require(h172_mod.degree() == 172, "h172 degree survives good reduction")
require(B_plus_2_mod.degree() == 30, "B+2 degree survives good reduction")
node_projection_gcd_mod = h172_mod.gcd(B_plus_2_mod)
require(node_projection_gcd_mod.degree() == 0, "no node over B=-2")


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


def field_coordinates(value):
    return certificate.field_coordinates(value)


def serialize_field(value):
    return ",".join(str(item) for item in field_coordinates(value))


def polynomial_hash(polynomial):
    payload = ";".join(
        f"{degree}:{serialize_field(coefficient(polynomial, degree))}"
        for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


points = (-1, 0, 1)
weights = (
    field.convert(5) / field.convert(18),
    -one,
    field.convert(13) / field.convert(18),
)


def evaluate(polynomial, point):
    return polynomial.evaluate(X, field.convert(point))


def retained_lambda(values):
    return sum(
        (weight * value for weight, value in zip(weights, values)),
        zero,
    )


def polynomial_lambda(polynomial):
    return retained_lambda(tuple(evaluate(polynomial, point) for point in points))


# The retained tangent plane is ker Lambda.  These are the two tangent rows
# that span it; B has zero tangent row at the triple.
C_tangent = tuple(evaluate(C_prime, point) for point in points)
E_tangent = tuple(evaluate(E_prime, point) for point in points)
require(C_tangent == (3 * one, 3 * one, 3 * one), "C tangent row")
require(E_tangent == (-9 * one, 4 * one, 9 * one), "E tangent row")
require(retained_lambda(C_tangent) == zero, "Lambda C tangent")
require(retained_lambda(E_tangent) == zero, "Lambda E tangent")


# If U is a kernel potential, the conductor-fibre structure of THM-4381
# first forces U=rH.  Membership of C'rH and E'rH in S is then equivalent
# to the following two retained first-jet equations on H modulo L.  Solve
# them with a monic quadratic H_0=x^2+u*x+v.
r_prime_values = tuple(evaluate(r.diff(X), point) for point in points)


def condition_row(multiplier_values):
    return tuple(
        retained_lambda(
            tuple(
                multiplier * r_prime * field.convert(point) ** power
                for multiplier, r_prime, point in zip(
                    multiplier_values, r_prime_values, points
                )
            )
        )
        for power in range(3)
    )


condition_C = condition_row(C_tangent)
condition_E = condition_row(E_tangent)
condition_minor = condition_C[0] * condition_E[1] - condition_C[1] * condition_E[0]
if not condition_minor:
    condition_minor = condition_C[0] * condition_E[2] - condition_C[2] * condition_E[0]
if not condition_minor:
    condition_minor = condition_C[1] * condition_E[2] - condition_C[2] * condition_E[1]
require(bool(condition_minor), "two independent retained conditions")

# Columns are the coefficients of 1,x,x^2.  Set the x^2 coefficient to one.
determinant = condition_C[1] * condition_E[0] - condition_C[0] * condition_E[1]
require(bool(determinant), "constant-linear solve pivot")
u = (-condition_C[2] * condition_E[0] + condition_C[0] * condition_E[2]) / determinant
v = (-condition_C[1] * condition_E[2] + condition_C[2] * condition_E[1]) / determinant
H0 = X**2 + u * X + v
H0_values = tuple(evaluate(H0, point) for point in points)
require(
    retained_lambda(
        tuple(c_tangent * r_prime * value for c_tangent, r_prime, value in zip(C_tangent, r_prime_values, H0_values))
    )
    == zero,
    "H0 C condition",
)
require(
    retained_lambda(
        tuple(e_tangent * r_prime * value for e_tangent, r_prime, value in zip(E_tangent, r_prime_values, H0_values))
    )
    == zero,
    "H0 E condition",
)

U0 = r * H0
potential_basis = (U0,) + tuple(c * X**power for power in range(17))
require(tuple(item.degree() for item in potential_basis) == tuple(range(177, 195)), "potential degrees")
require(
    tuple(item.degree() % 18 for item in potential_basis)
    == (15, 16, 17) + tuple(range(15)),
    "one leading degree per Z residue",
)

# Z vanishes at all three retained points.  The factor theorem therefore
# gives Z=L*W with W monic of degree 15, and Z*U0=c*A with A=W*H0 monic of
# degree 17.  The triangular basis 1,x,...,x^16,A of K[x] over K[Z] proves
# that the displayed 18 potentials generate K*U0+cK[x], while distinct
# leading residues prove freeness.  We check the three exact values rather
# than invoking PolyElement.exquo: SymPy 1.13 retains an explicit zero ANP
# coefficient in this polynomial, which breaks its exact-quotient predicate.
retained_fibre_gcd = B.gcd(C).gcd(E + 3).monic()
require(retained_fibre_gcd == L, "retained fibre gcd")
require(
    tuple(evaluate(Z, point) for point in points) == (zero, zero, zero),
    "Z vanishes on retained triple",
)


# Extend the inherited monic Apéry reduction just far enough to certify that
# both entries of all 18 restriction channels really lie in S.
restriction_pairs = tuple((C_prime * potential, E_prime * potential) for potential in potential_basis)
max_restriction_degree = max(
    polynomial.degree() for pair in restriction_pairs for polynomial in pair
)
canonical = {}
for basis in sagbi.module_basis:
    polynomial = basis
    while polynomial.degree() <= max_restriction_degree:
        require(polynomial.degree() not in canonical, "unique filtered S degree")
        canonical[polynomial.degree()] = polynomial
        polynomial *= Z
require(
    all(degree in canonical for degree in range(170, max_restriction_degree + 1)),
    "filtered S tail coverage",
)


def normal_form(polynomial):
    remainder = polynomial
    for degree in range(remainder.degree(), -1, -1):
        basis = canonical.get(degree)
        if basis is None:
            continue
        value = coefficient(remainder, degree)
        if value:
            remainder -= value * basis
    require(not remainder or remainder.degree() <= 169, "gap normal-form bound")
    return remainder


for index, ((P, R), potential) in enumerate(zip(restriction_pairs, potential_basis)):
    require(not normal_form(P), f"channel {index} first entry in S")
    require(not normal_form(R), f"channel {index} second entry in S")
    require(C_prime * R - E_prime * P == ring.zero, f"channel {index} in ker L0")
    require(
        potential == sagbi.field.convert(1) * potential,
        f"channel {index} typed potential",
    )


# For every kernel channel, U vanishes at the retained triple.  Since P' and
# Q' are tangent rows, Lambda(U')=Lambda(E'U')=0; the intersection is the
# constant line.  Check this exact conclusion on the extracted basis.
potential_derivative_rows = tuple(
    tuple(evaluate(potential.diff(X), point) for point in points)
    for potential in potential_basis
)
potential_lambdas = []
for index, (potential, derivative_row) in enumerate(
    zip(potential_basis, potential_derivative_rows)
):
    require(
        tuple(evaluate(potential, point) for point in points) == (zero, zero, zero),
        f"channel {index} retained vanishing",
    )
    require(
        derivative_row == (derivative_row[0],) * 3,
        f"channel {index} constant derivative row",
    )
    potential_lambdas.append(derivative_row[0])
require(bool(potential_lambdas[0]), "exceptional channel has nonzero first jet")
require(all(value == zero for value in potential_lambdas[1:]), "conductor channels have zero first jet")


# Fixed-x source-normal sidecars for the Russell generators.  At the triple,
# the sidecar of a target representative is determined by its ordinary
# derivative row.  Thus a channel with U'=lambda*(1,1,1) has sidecar values
# gamma_1(P~)=lambda*n_C and gamma_1(Q~)=lambda*n_E for every chosen actual
# target lift P~,Q~ of its restrictions.
D = 1 + X**2 * Q
normal_B = 3 * X**2 * D * (D + 2)
normal_C = 2 * X**3 * (D + 1)
normal_E = 2 * (D + 1)
normal_rows = tuple(
    tuple(evaluate(polynomial, point) for point in points)
    for polynomial in (normal_B, normal_C, normal_E)
)
require(normal_rows[0] == (zero, zero, zero), "normal B row")
require(normal_rows[1] == (2 * one, zero, -2 * one), "normal C row")
require(normal_rows[2] == (-2 * one, 4 * one, -2 * one), "normal E row")

# The full first output is the s-derivative of
# L_s(P~,Q~)=C_s' gamma_s(Q~)-E_s' gamma_s(P~).  Terms containing P,Q vanish
# at the retained triple, leaving the contraction below.  Its transverse
# sidecar components cancel to a constant row, annihilated by Lambda.
first_output_rows = []
for lam in potential_lambdas:
    sidecar_P = tuple(lam * value for value in normal_rows[1])
    sidecar_Q = tuple(lam * value for value in normal_rows[2])
    first_output = tuple(
        c_tangent * q_value - e_tangent * p_value
        for c_tangent, e_tangent, p_value, q_value in zip(
            C_tangent, E_tangent, sidecar_P, sidecar_Q
        )
    )
    require(first_output == (12 * lam,) * 3, "contracted first-output row")
    require(retained_lambda(first_output) == zero, "first-output Lambda obstruction")
    first_output_rows.append(first_output)


# Lambda descends to N/dS: every S-derivative has retained row in the tangent
# plane.  The seminormal generator r supplies the missing derivative row.
require(polynomial_lambda(B.diff(X)) == zero, "dB in retained plane")
require(polynomial_lambda(C.diff(X)) == zero, "dC in retained plane")
require(polynomial_lambda(E.diff(X)) == zero, "dE in retained plane")
r_prime_lambda = polynomial_lambda(r.diff(X))
require(bool(r_prime_lambda), "seminormal derivative is transverse")


potential_hash = hashlib.sha256(
    ";".join(polynomial_hash(item) for item in potential_basis).encode("ascii")
).hexdigest()
pair_hash = hashlib.sha256(
    ";".join(
        f"{polynomial_hash(P)}:{polynomial_hash(R)}" for P, R in restriction_pairs
    ).encode("ascii")
).hexdigest()
H0_coefficient_payload = f"u={serialize_field(u)};v={serialize_field(v)}"

print("scope=exceptional_quartic_fixed_x_chosen_representative_first_output;not_THM4067_moving_family")
print("rings=S=K[B,C,E]_subset_N=K[x];L0(P,Q)=Cprime*Q-Eprime*P;Z_degree=18")
print("node_projection_certificate=gcd_mod137_alpha44(h172,B+2)=1;positive_resultant_direction")
print("kernel_potential_module=M={U_in_N:Cprime*U,Eprime*U_in_S}")
print("local_characterization=M={r*H:Lambda(rprime*H)=Lambda(Eprime*rprime*H)=0}")
print("H0=x^2+u*x+v;H0_coefficients=" + H0_coefficient_payload)
print("H0_coefficient_hash=" + hashlib.sha256(H0_coefficient_payload.encode("ascii")).hexdigest())
print("M=K*U0+c*N;U0=r*H0;Z*U0=c*(Z/L)*H0")
print("KZ_potential_basis=U0,c,x*c,...,x^16*c")
print("potential_degrees=" + ",".join(str(item.degree()) for item in potential_basis))
print("potential_residues_mod18=" + ",".join(str(item.degree() % 18) for item in potential_basis))
print("restriction_channel_degrees=" + ";".join(f"{P.degree()},{R.degree()}" for P, R in restriction_pairs))
print("potential_basis_hash=" + potential_hash)
print("restriction_pair_hash=" + pair_hash)
print("retained_Uprime_pattern=nonzero_constant,0x17")
print("normal_rows_B_C_E=" + ";".join(",".join(serialize_field(value) for value in row) for row in normal_rows))
print("raw_sidecar_exceptional_channel_is_locally_transverse=True")
print("contracted_first_output_rows=12*lambda*(1,1,1);Lambda_zero=True")
print("Lambda_descends_to_N_mod_dS=True;rprime_Lambda_nonzero=True")
print("rprime_Lambda_hash=" + hashlib.sha256(serialize_field(r_prime_lambda).encode("ascii")).hexdigest())
print("CONCLUSION=[rprime]_not_in_kernel_constrained_first_output_mod_dS")
print("TYPE_FIREWALL=no_descended_gamma1_on_S;no_moving_endpoint_or_graph_family_identification")
print("NO_CLAIM=J8_payment_or_Keller_pair_or_JC2_or_DC2")
print("RESULT=PASS")
