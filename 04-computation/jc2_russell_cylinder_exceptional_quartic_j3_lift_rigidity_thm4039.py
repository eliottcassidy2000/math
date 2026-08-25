"""Exact retained-point J3/Lambda proof for the frozen THM-3688 certificate.

Unlike the full-polynomial scout, this propagates only the four local data
(gamma H, d_x gamma H, gamma partial_q H, d_x gamma partial_q H) at the
three retained points through the canonical target expressions.  It also
gates the retained tangent and normal rows used by the one- and two-step
solution-rigidity proofs.
"""

import contextlib
import hashlib
import io
import subprocess
import types

from flint import fmpq


PATH = "04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.py"
source = subprocess.check_output(["git", "show", f"HEAD:{PATH}"], text=True)
source_sha256 = hashlib.sha256(source.encode("utf-8")).hexdigest()
module = types.ModuleType("thm3688_frozen")
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(source, PATH, "exec"), module.__dict__)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


POINTS = (-1, 0, 1)
V, DX, DQ, DXDQ = range(4)


def escale(value, scalar):
    scalar = fmpq(scalar)
    return tuple(scalar * coordinate for coordinate in value)


def eadd(left, right):
    return module.elem_add(left, right)


def esub(left, right):
    return module.elem_sub(left, right)


def emul(left, right):
    return module.elem_mul(left, right)


def eval_poly(polynomial, point):
    point = fmpq(point)
    return tuple(component(point) for component in polynomial)


def raw_node(polynomial, normal):
    dx = module.kpoly_diff(polynomial)
    dxdq = module.kpoly_diff(normal)
    jets = tuple(
        (
            eval_poly(polynomial, point),
            eval_poly(dx, point),
            eval_poly(normal, point),
            eval_poly(dxdq, point),
        )
        for point in POINTS
    )
    return polynomial, jets


def nadd(left, right):
    return (
        module.kpoly_add(left[0], right[0]),
        tuple(tuple(eadd(a, b) for a, b in zip(ljet, rjet)) for ljet, rjet in zip(left[1], right[1])),
    )


def nsub(left, right):
    return (
        module.kpoly_sub(left[0], right[0]),
        tuple(tuple(esub(a, b) for a, b in zip(ljet, rjet)) for ljet, rjet in zip(left[1], right[1])),
    )


def nmul(left, right):
    jets = []
    for u, v in zip(left[1], right[1]):
        jets.append(
            (
                emul(u[V], v[V]),
                eadd(emul(u[DX], v[V]), emul(u[V], v[DX])),
                eadd(emul(u[DQ], v[V]), emul(u[V], v[DQ])),
                eadd(
                    eadd(emul(u[DXDQ], v[V]), emul(u[DQ], v[DX])),
                    eadd(emul(u[DX], v[DQ]), emul(u[V], v[DXDQ])),
                ),
            )
        )
    return module.kpoly_mul(left[0], right[0]), tuple(jets)


def npow(value, exponent):
    answer = raw_node(module.kpoly_constant(module.elem_one()), module.kpoly_zero())
    for _ in range(exponent):
        answer = nmul(answer, value)
    return answer


def nscale(value, scalar):
    return (
        module.kpoly_mul_elem(value[0], scalar),
        tuple(tuple(emul(component, scalar) for component in jet) for jet in value[1]),
    )


def nmonic(value, label):
    degree = module.kpoly_degree(value[0])
    scalar = module.elem_inv(module.kpoly_coeff(value[0], degree))
    answer = nscale(value, scalar)
    require(answer[0] == module.kpoly_monic(value[0], label), f"{label} monic value")
    return answer


B = raw_node(module.B, module.from_sympy_poly(module.exact.dB))
C = raw_node(module.C, module.dC)
E = raw_node(module.E, module.dE)
three = raw_node(
    module.kpoly_constant((fmpq(3), fmpq(0), fmpq(0), fmpq(0))),
    module.kpoly_zero(),
)
Z = nmonic(nadd(E, three), "Z-jet")
C0 = nmonic(C, "C0-jet")
B0 = nmonic(B, "B0-jet")
G71 = nmonic(nsub(npow(Z, 4), nmul(npow(C0, 2), B0)), "G71-jet")
G83 = nmonic(nsub(npow(C0, 4), nmul(npow(Z, 3), B0)), "G83-jet")
G124 = nmonic(nsub(nmul(G71, npow(Z, 3)), nmul(G83, npow(C0, 2))), "G124-jet")

basis = (
    raw_node(module.kpoly_constant(module.elem_one()), module.kpoly_zero()),
    nmul(C0, G124), nmul(C0, G71), C0, npow(G83, 2), nmul(B0, G83),
    npow(C0, 2), nmul(C0, npow(G83, 2)), nmul(nmul(C0, B0), G83),
    npow(C0, 3), nmul(G71, G83), G83, B0, nmul(nmul(C0, G71), G83),
    nmul(C0, G83), nmul(C0, B0), G124, G71,
)
require(tuple(node[0] for node in basis) == module.module_basis, "module-basis restrictions")

canonical = []
for residue, node in enumerate(basis):
    power = 0
    value = node
    while module.kpoly_degree(value[0]) <= 375:
        canonical.append(((residue, power), value))
        power += 1
        value = nmul(value, Z)
canonical.sort(key=lambda item: module.kpoly_degree(item[1][0]))
require([address for address, _node in canonical] == module.monomials, "canonical address order")
require([node[0] for _address, node in canonical] == [item[1] for item in module.canonical_items], "canonical restrictions")


def combine(coefficients, component):
    result = [module.elem_zero() for _point in POINTS]
    for coefficient, (_address, node) in zip(coefficients, canonical):
        if not any(coefficient):
            continue
        for index in range(len(POINTS)):
            result[index] = eadd(result[index], emul(coefficient, node[1][index][component]))
    return tuple(result)


F2v, F2x, F2q = (combine(module.F2_values, component) for component in (V, DX, DQ))
G2v, G2x, G2q = (combine(module.G2_values, component) for component in (V, DX, DQ))
F3v, F3x = (combine(module.F3_values, component) for component in (V, DX))
G3v, G3x = (combine(module.G3_values, component) for component in (V, DX))


def values(polynomial):
    return tuple(eval_poly(polynomial, point) for point in POINTS)


Cv, Cx = values(module.C), values(module.C_prime)
Ev, Ex = values(module.E), values(module.E_prime)
F1v, F1x = values(module.F1), values(module.F1_prime)
G1v, G1x = values(module.G1), values(module.G1_prime)
dBv = values(module.from_sympy_poly(module.exact.dB))
dCv, dCx = values(module.dC), values(module.kpoly_diff(module.dC))
dEv, dEx = values(module.dE), values(module.kpoly_diff(module.dE))
dF1v, dF1x = values(module.delta_F1), values(module.kpoly_diff(module.delta_F1))
dG1v, dG1x = values(module.delta_G1), values(module.kpoly_diff(module.delta_G1))

point_terms = []
point_direct = []
for index, point in enumerate(POINTS):
    a2v, a2x = eadd(F2v[index], dCv[index]), eadd(F2x[index], dCx[index])
    b2v, b2x = eadd(G2v[index], dEv[index]), eadd(G2x[index], dEx[index])
    a3v, a3x = eadd(F3v[index], dF1v[index]), eadd(F3x[index], dF1x[index])
    b3v, b3x = eadd(G3v[index], dG1v[index]), eadd(G3x[index], dG1x[index])
    a4v = eadd(F2q[index], (fmpq(point**5), fmpq(0), fmpq(0), fmpq(0)))
    b4v = eadd(G2q[index], (fmpq(point**2), fmpq(0), fmpq(0), fmpq(0)))

    term04 = escale(esub(emul(Cx[index], b4v), emul(a4v, Ex[index])), 4)
    term13 = esub(escale(emul(F1x[index], b3v), 3), emul(F1v[index], b3x))
    term22 = escale(esub(emul(a2x, b2v), emul(a2v, b2x)), 2)
    term31 = esub(emul(a3x, G1v[index]), escale(emul(a3v, G1x[index]), 3))
    grouped = (term04, term13, term22, term31)
    point_terms.append(grouped)

    av = (Cv[index], F1v[index], a2v, a3v, a4v)
    ax = (Cx[index], F1x[index], a2x, a3x, module.elem_zero())
    bv = (Ev[index], G1v[index], b2v, b3v, b4v)
    bx = (Ex[index], G1x[index], b2x, b3x, module.elem_zero())
    direct = module.elem_zero()
    for i in range(5):
        j = 4 - i
        direct = eadd(direct, esub(escale(emul(ax[i], bv[j]), j), escale(emul(av[i], bx[j]), i)))
    require(direct == eadd(eadd(term04, term13), eadd(term22, term31)), f"direct D3 at {point}")
    point_direct.append(direct)


def Lambda(point_values):
    return escale(
        eadd(esub(escale(point_values[0], 5), escale(point_values[1], 18)), escale(point_values[2], 13)),
        fmpq(1, 18),
    )


contributions = tuple(Lambda(tuple(point_terms[p][term] for p in range(3))) for term in range(4))
gate = Lambda(tuple(point_direct))
require(gate == eadd(eadd(contributions[0], contributions[1]), eadd(contributions[2], contributions[3])), "Lambda additivity")
require(source_sha256 == "02cd67446b18b3863bc3665d48a6c5cccda81c394f94b754d2b90b1597c53ba6", "frozen THM3688 source hash")
require(not any(gate), "exact J3 scalar gate")

# The proof of adjacent-solution rigidity uses only these exact tangent rows.
# They are derived here from the same frozen quartic-field restrictions rather
# than trusted as copied constants.  For any target restriction P=p(B,C,E),
# the chain rule puts its retained derivative row in their span.
B1 = values(module.kpoly_diff(module.B))
C1 = values(module.C_prime)
E1 = values(module.E_prime)
zero = module.elem_zero()


def rational_element(value):
    return (fmpq(value), fmpq(0), fmpq(0), fmpq(0))


require(B1 == (zero, zero, zero), "retained B-prime row")
require(C1 == (rational_element(3),) * 3, "retained C-prime row")
require(
    E1 == (rational_element(-9), rational_element(4), rational_element(9)),
    "retained E-prime row",
)
require(not any(Lambda(B1)), "Lambda annihilates B-prime")
require(not any(Lambda(C1)), "Lambda annihilates C-prime")
require(not any(Lambda(E1)), "Lambda annihilates E-prime")
require(fmpq(3) * fmpq(4) - fmpq(3) * fmpq(-9) != 0, "retained tangent rank two")

# If an admissible L0-kernel has potential U, then U vanishes at the three
# retained points.  Admissibility of the first derivatives says both U' and
# E'U' lie in the tangent hyperplane.  These two rows have independent
# annihilators, and their intersection is exactly the constant line.
lambda_numerator = (fmpq(5), fmpq(-18), fmpq(13))
eprime_rational = (fmpq(-9), fmpq(4), fmpq(9))
weighted_row = tuple(a * b for a, b in zip(lambda_numerator, eprime_rational))
require(
    lambda_numerator[0] * weighted_row[1]
    - lambda_numerator[1] * weighted_row[0]
    == fmpq(-1170),
    "two-step tangent intersection rank",
)
require(not any(Lambda((rational_element(1),) * 3)), "constant tangent line")

# The vertical compiler rows turn the surviving constant U' gauge into a
# constant normal response, again killed by Lambda.  The retained J1 identity
# supplies the other cancellation used in the J5 two-step response proof.
require(dBv == (zero, zero, zero), "retained delta B row")
require(
    dCv == (rational_element(2), rational_element(0), rational_element(-2)),
    "retained delta C row",
)
require(
    dEv == (rational_element(-2), rational_element(4), rational_element(-2)),
    "retained delta E row",
)
normal_response = tuple(
    esub(escale(dEv[index], 3), emul(E1[index], dCv[index]))
    for index in range(len(POINTS))
)
require(normal_response == (rational_element(12),) * 3, "two-step normal response")
require(not any(Lambda(normal_response)), "Lambda kills two-step normal response")
require(F1v == (zero, zero, zero), "retained F1 values")
require(G1v == (rational_element(fmpq(1, 3)),) * 3, "retained G1 values")
require(not any(Lambda(F1x)), "Lambda kills F1 derivative")
a2_retained = tuple(eadd(F2v[index], dCv[index]) for index in range(len(POINTS)))
b2_retained = tuple(eadd(G2v[index], dEv[index]) for index in range(len(POINTS)))
j1_tangent_relation = tuple(
    eadd(
        esub(escale(b2_retained[index], 3), emul(a2_retained[index], E1[index])),
        escale(F1x[index], fmpq(1, 6)),
    )
    for index in range(len(POINTS))
)
require(j1_tangent_relation == (zero, zero, zero), "retained J1 tangent relation")


def etext(value):
    return ",".join(str(coordinate) for coordinate in value)


print("field=THM3683_exceptional_quartic")
print("representative=THM3688_frozen_canonical_target_certificate")
print("observer=retained_value_dx_dq_dxdq_jets_at_-1_0_1")
print("equation=J3=D3+4*L0(F4,G4)")
print(f"thm3688_source_sha256={source_sha256}")
for label, contribution in zip(("term04", "term13", "term22", "term31"), contributions):
    payload = etext(contribution)
    print(f"Lambda_{label}_nonzero={int(any(contribution))};chars={len(payload)};sha256={hashlib.sha256(payload.encode('ascii')).hexdigest()}")
print(f"Lambda_D3={etext(gate)}")
print(f"J3_stagewise_solvable={int(not any(gate))}")
print("retained_tangent_rows=Bprime:0,0,0;Cprime:3,3,3;Eprime:-9,4,9")
print("retained_tangent_rank=2;annihilator=Lambda")
print("retained_normal_rows=deltaB:0,0,0;deltaC:2,0,-2;deltaE:-2,4,-2")
print("adjacent_kernel_cokernel_response=0_for_all_m_ge_2")
print("two_step_J5_solution_choice_response=0")
print("all_four_embeddings_uniform=1")
print("NO_CLAIM=J5_value_or_coherent_all_order_or_global_pair_or_Keller_map_or_JC2")
print("RESULT=" + ("PASS" if not any(gate) else "FAIL"))
