"""Exact restriction-ring SAGBI/module certificate for THM-3703."""

import hashlib
import importlib.util
from pathlib import Path
import subprocess
import types

import sympy as sp
from sympy.polys.rings import ring


PREDECESSOR = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_modular_lift_thm3687.py"
)
EXPECTED_APERY = (
    0, 145, 92, 21, 166, 113, 42, 187, 134,
    63, 154, 83, 30, 175, 104, 51, 124, 71,
)
EXPECTED_GAPS = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
    37, 38, 40, 41, 43, 44, 45, 46, 47, 49, 50, 52, 53, 55, 56,
    58, 59, 61, 62, 64, 65, 67, 68, 70, 73, 74, 76, 77, 79, 80,
    82, 85, 86, 88, 91, 94, 95, 97, 98, 100, 103, 106, 109, 112,
    115, 116, 118, 121, 127, 130, 133, 136, 139, 148, 151, 157, 169,
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_predecessor():
    if PREDECESSOR.is_file():
        spec = importlib.util.spec_from_file_location("thm3687_modular", PREDECESSOR)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    source = subprocess.check_output(
        ["git", "show", f"HEAD:{PREDECESSOR.as_posix()}"], text=True
    )
    module = types.ModuleType("thm3687_modular")
    exec(compile(source, PREDECESSOR.as_posix(), "exec"), module.__dict__)
    return module


probe = load_predecessor()
alpha_symbol = sp.symbols("alpha")
field = sp.QQ.alg_field_from_poly(
    sp.Poly(probe.F_QUARTIC.subs(probe.r, alpha_symbol), alpha_symbol),
    "alpha",
)
alpha = field.from_sympy(field.ext.as_expr())
polynomial_ring, X = ring("X", field)


def field_coefficient(expression):
    return field.from_sympy(
        sp.cancel(expression).subs(probe.r, field.ext.as_expr())
    )


Q = polynomial_ring.zero
for (degree,), coefficient in probe.Q_R.terms():
    Q += field_coefficient(coefficient.as_expr()) * X**degree

D = 1 + X**2 * Q
B = (D - 1) * (D + 2) ** 2
C = X * D * (D + 2)
E = Q * (D + 3)
require(C * C * E == B * (B + 4), "compiler relation")
require((B.degree(), C.degree(), E.degree()) == (30, 21, 18), "compiler degrees")
require(field_coefficient(probe.F_QUARTIC) == field.zero, "quartic field relation")


def monic(polynomial, label):
    require(bool(polynomial), f"nonzero {label}")
    answer = polynomial.monic()
    require(answer.LC == field.one, f"monic {label}")
    return answer


# The collision point has E=-3, B=C=0.  Scaling by field units and shifting
# E by 3 do not change the generated restriction algebra.
Z = monic(E + 3, "Z")
C0 = monic(C, "C0")
B0 = monic(B, "B0")

# Three successive leading-term collisions produce the missing residue
# classes.  Both summands in each parenthesis are monic of the same degree.
G71 = monic(Z**4 - C0**2 * B0, "G71")
G83 = monic(C0**4 - Z**3 * B0, "G83")
G124 = monic(G71 * Z**3 - G83 * C0**2, "G124")
generators = (Z, C0, B0, G71, G83, G124)
require(tuple(item.degree() for item in generators) == (18, 21, 30, 71, 83, 124), "SAGBI degrees")

# Apéry representatives in residue order modulo deg Z=18.
module_basis = (
    polynomial_ring.one,
    C0 * G124,
    C0 * G71,
    C0,
    G83**2,
    B0 * G83,
    C0**2,
    C0 * G83**2,
    C0 * B0 * G83,
    C0**3,
    G71 * G83,
    G83,
    B0,
    C0 * G71 * G83,
    C0 * G83,
    C0 * B0,
    G124,
    G71,
)
require(all(item.LC == field.one for item in module_basis), "monic module basis")
apery = tuple(item.degree() for item in module_basis)
require(apery == EXPECTED_APERY, "Apéry degrees")
require(tuple(degree % 18 for degree in apery) == tuple(range(18)), "Apéry residues")


def field_coordinates(value):
    entries = list(value.to_list())
    entries = [sp.Rational(0)] * (4 - len(entries)) + [sp.Rational(entry) for entry in entries]
    return tuple(reversed(entries))


def serialize_field(value):
    return ",".join(str(item) for item in field_coordinates(value))


reduction_steps = []


def reduce_module(polynomial, label):
    remainder = polynomial
    local_steps = 0
    while remainder:
        degree = remainder.degree()
        residue = degree % 18
        base = module_basis[residue]
        difference = degree - base.degree()
        require(difference >= 0 and difference % 18 == 0, f"module leading degree {label}")
        coefficient = remainder.LC
        power = difference // 18
        reduction_steps.append((label, degree, residue, power, coefficient))
        remainder -= coefficient * Z**power * base
        local_steps += 1
    return local_steps


step_counts = []
for generator_label, generator in (("Z", Z), ("C0", C0), ("B0", B0)):
    for index, basis_element in enumerate(module_basis):
        step_counts.append(
            reduce_module(generator * basis_element, f"{generator_label}:{index}")
        )
require(len(step_counts) == 54, "module product count")
require((max(step_counts), sum(step_counts)) == (120, 522), "module reduction anatomy")

# The module presentation makes the leading-degree semigroup and filtered
# codimension literal.  This is the numerical-semigroup conductor, not the
# global conductor ideal inside K[X].
degree_semigroup = {
    base + 18 * power
    for base in apery
    for power in range(30)
}
gaps = tuple(degree for degree in range(376) if degree not in degree_semigroup)
require(gaps == EXPECTED_GAPS, "degree gaps")
require(len(gaps) == 89, "degree genus")
require(all(degree in degree_semigroup for degree in range(170, 376)), "semigroup conductor upper")
require(169 not in degree_semigroup, "semigroup conductor sharp")


def polynomial_hash(polynomial):
    payload = ";".join(
        f"{degree}:{serialize_field(polynomial.get((degree,), field.zero))}"
        for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


generator_hashes = tuple(polynomial_hash(item) for item in generators)
module_hash = hashlib.sha256(
    ";".join(polynomial_hash(item) for item in module_basis).encode("ascii")
).hexdigest()
reduction_hash = hashlib.sha256(
    ";".join(
        f"{label}:{degree}:{residue}:{power}:{serialize_field(coefficient)}"
        for label, degree, residue, power, coefficient in reduction_steps
    ).encode("ascii")
).hexdigest()

print("field=THM3683_exceptional_quartic")
print("sagbi_degrees=18,21,30,71,83,124")
print("apery_mod18=" + ",".join(map(str, apery)))
print("degree_genus=89;degree_semigroup_conductor=170;global_conductor_ideal=OPEN")
print("module=direct_sum_r=0^17_K[Z]p_r;products_checked=54;max_steps=120;total_steps=522")
print("generator_hashes=" + ",".join(generator_hashes))
print(f"module_hash={module_hash}")
print(f"reduction_hash={reduction_hash}")
print("normalization=K[X];fraction_field=K(X);collision_retained=True")
print("RESULT=PASS")
