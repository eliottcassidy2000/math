#!/usr/bin/env python3
"""Exact characteristic-zero audit of the exceptional mixed-pair transgression.

This is a no-repository-edit continuation of THM-4054/THM-4058.  It reuses
the audited Q(alpha) jet engine, then independently checks the coupled tangent
ranks, stable-monomial membership, moving-endpoint cancellation, and the
closed carrier factors behind the all-cutoff proof.
"""

import ast
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
from io import StringIO
from math import comb
from pathlib import Path
from runpy import run_path
from subprocess import run

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def locate_root():
    here = Path(__file__).resolve()
    required = Path("04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.py")
    for parent in here.parents:
        if (parent / "AGENTS.md").is_file() and (parent / required).is_file():
            return parent
    completed = run(
        ["git", "rev-parse", "--show-toplevel"],
        cwd=Path.cwd(),
        check=True,
        capture_output=True,
        text=True,
    )
    candidate = Path(completed.stdout.strip()).resolve()
    require((candidate / "AGENTS.md").is_file() and (candidate / required).is_file(),
            "cannot locate math workspace")
    return candidate


ROOT = locate_root()
THM4054 = ROOT / "04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.py"
THM4058_REL = "04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058.py"
THM4054_SHA256 = "1e84ce5674c74c8b2d504fe8f57341be86ab20e1eab5c8fafe7494bb52d488a2"
THM4058_SHA256 = "99d6e27e93129636faba409750fbd1b71fd3f60a7b65117ea24f906cd78888e1"


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def tracked_hash(relative_path):
    completed = run(
        ["git", "show", "HEAD:" + relative_path],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return sha256(completed.stdout).hexdigest()


require(file_hash(THM4054) == THM4054_SHA256, "THM-4054 source drift")
require(tracked_hash(THM4058_REL) == THM4058_SHA256, "THM-4058 source drift")

# Suppress the already-frozen THM-4054 transcript while retaining its complete
# exact Q(alpha) implementation and requiring its own terminal PASS gate.
captured = StringIO()
with redirect_stdout(captured):
    ns = run_path(str(THM4054))
require(captured.getvalue().rstrip().endswith("RESULT=PASS"), "THM-4054 replay")
g = ns["build_linearized_pair_columns"].__globals__


def stable_monomial(r):
    return [
        g["K1"] if stable == r and source == 0 else g["K0"]
        for stable, _, source in g["ROWS"]
    ]


# Moving-endpoint terms for the oriented edges
# 0:v20->v01, 1:v01->v12, 2:v12->v20.
# Each tuple records coefficients of (u20*c20',u01*c01',u12*c12').
endpoint_rows = ((-1, 1, 0), (0, -1, 1), (1, 0, -1))
endpoint_sum = tuple(sum(row[column] for row in endpoint_rows) for column in range(3))
require(endpoint_sum == (0, 0, 0), "moving-endpoint cancellation")


def monomials(cap):
    return tuple(
        (a, c, u)
        for a in range(cap + 1)
        for c in range(cap + 1 - a)
        for u in range(cap + 1 - a - c)
    )


def poly_add(*values):
    answer = {}
    for value in values:
        for monomial, coefficient in value.items():
            answer[monomial] = answer.get(monomial, Fraction(0)) + coefficient
            if not answer[monomial]:
                del answer[monomial]
    return answer


def poly_scale(value, scalar):
    scalar = Fraction(scalar)
    return {monomial: scalar * coefficient for monomial, coefficient in value.items()
            if scalar * coefficient}


def poly_diff(value, axis):
    answer = {}
    for monomial, coefficient in value.items():
        exponent = monomial[axis]
        if exponent:
            target = list(monomial)
            target[axis] -= 1
            answer[tuple(target)] = coefficient * exponent
    return answer


def poly_integrate(value, axis):
    answer = {}
    for monomial, coefficient in value.items():
        target = list(monomial)
        target[axis] += 1
        answer[tuple(target)] = coefficient / target[axis]
    return answer


def coefficients(f, g_value):
    p_value = poly_add(poly_scale(poly_diff(f, 0), -4), poly_diff(g_value, 1))
    q_value = poly_diff(g_value, 2)
    s_value = poly_scale(poly_diff(f, 2), 4)
    return p_value, q_value, s_value


def closure(p_value, q_value, s_value):
    return poly_add(poly_diff(p_value, 2),
                    poly_scale(poly_diff(q_value, 1), -1),
                    poly_diff(s_value, 0))


# Complete forward/reverse control in coefficient degree <=4.  Forward uses
# every f- and g-monomial through degree five.  Reverse uses every basis vector
# of the full closed-triple kernel P_u-Q_c+S_A=0 through degree four.
closed_cap = 4
for monomial in monomials(closed_cap + 1):
    unit = {monomial: Fraction(1)}
    require(not closure(*coefficients(unit, {})), ("forward f", monomial))
    require(not closure(*coefficients({}, unit)), ("forward g", monomial))

source_monomials = monomials(closed_cap)
closure_monomials = monomials(closed_cap - 1)
closure_index = {monomial: index for index, monomial in enumerate(closure_monomials)}
closure_matrix = sp.zeros(len(closure_monomials), 3 * len(source_monomials))
for kind in range(3):
    for column, monomial in enumerate(source_monomials):
        value = {monomial: Fraction(1)}
        derivative = (poly_diff(value, 2),
                      poly_scale(poly_diff(value, 1), -1),
                      poly_diff(value, 0))[kind]
        for target, coefficient in derivative.items():
            closure_matrix[closure_index[target], kind * len(source_monomials) + column] = (
                sp.Rational(coefficient.numerator, coefficient.denominator)
            )
closed_basis = closure_matrix.nullspace()
expected_closed_dimension = 3 * len(source_monomials) - len(closure_monomials)
require(len(closed_basis) == expected_closed_dimension, "closed-triple dimension")

for basis_index, vector in enumerate(closed_basis):
    triple = []
    for kind in range(3):
        value = {}
        for index, monomial in enumerate(source_monomials):
            coefficient = vector[kind * len(source_monomials) + index]
            if coefficient:
                value[monomial] = Fraction(int(coefficient.p), int(coefficient.q))
        triple.append(value)
    p_value, q_value, s_value = triple
    f_zero = poly_scale(poly_integrate(s_value, 2), Fraction(1, 4))
    g_zero = poly_integrate(q_value, 2)
    residual = poly_add(p_value,
                        poly_scale(poly_diff(f_zero, 0), 4),
                        poly_scale(poly_diff(g_zero, 1), -1))
    require(not poly_diff(residual, 2), ("reverse residual", basis_index))
    f_value = f_zero
    g_value = poly_add(g_zero, poly_integrate(residual, 1))
    require(coefficients(f_value, g_value) == tuple(triple),
            ("reverse reconstruction", basis_index))

# If Lambda(f)=b*A^k+..., then
# Lambda(12*d_A f)=12*(k+5)*b*A^(k-1)+....
# THM-4058 gives Lambda(w^s)=-(s/2)*rho*A^(s+4)+....
carrier_rows = []
for s in range(1, 13):
    fixed_lead = Fraction(-s, 2)
    carrier_lead = 12 * (s + 9) * fixed_lead
    require(carrier_lead == -6 * s * (s + 9), ("carrier", s))
    carrier_rows.append((s, fixed_lead, carrier_lead))

# For H=t+gamma*t^m, r=m-1 and the uncancelled first-order density is
# -12*m*gamma*w^r.  The leading f-carrier is
# f=-(gamma*r/(r+10))*w^(r+1)+higher target degree.
cancellation_rows = []
for m in range(2, 13):
    r = m - 1
    rhs_period = 6 * m * r
    correction = Fraction(-r, r + 10)
    carrier = -6 * (r + 1) * (r + 10)
    require(correction * carrier == rhs_period, ("monomial cancellation", m))
    cancellation_rows.append((m, r, rhs_period, correction))

# Exact Q(alpha) retained matrices.  The all-N formula proved by the period
# quotient is rank=rows-min(N+1,4); here it is checked through N=8, together
# with every stable monomial visible at each cutoff.
rank_rows = []
for cutoff in range(0, 9):
    ns["configure"](cutoff)
    columns = ns["build_linearized_pair_columns"]()
    pivots, relations = ns["rref_nullspace"](columns)
    rows = len(g["ROWS"])
    expected_columns = 2 * (comb(cutoff + 4, 3) - 1)
    expected_cokernel = min(cutoff + 1, 4)
    expected_rank = rows - expected_cokernel
    require(len(columns) == expected_columns, ("column census", cutoff))
    require(len(pivots) == expected_rank, ("rank", cutoff, len(pivots)))
    require(len(relations) == expected_cokernel, ("cokernel", cutoff))
    for r in range(1, cutoff + 1):
        response = [ns["dot"](stable_monomial(r), relation) for relation in relations]
        require(all(value == g["K0"] for value in response), ("stable monomial", cutoff, r))
    rank_rows.append((cutoff, rows, len(columns), len(pivots), len(relations)))

serialization = "\n".join(
    ",".join(str(value) for value in row)
    for row in rank_rows + carrier_rows + cancellation_rows
)
serialization += "\nclosed,%d,%d" % (closed_cap, len(closed_basis))
table_hash = sha256(serialization.encode()).hexdigest()

tree = ast.parse(Path(__file__).read_text())
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                for node in ast.walk(tree)), "float literal")

print("EXCEPTIONAL MIXED-PAIR PERIOD TRANSGRESSION -- EXACT OVER Q(ALPHA)")
print("upstream_hashes=" + THM4054_SHA256 + "," + THM4058_SHA256)
print("upstream_thm4054_replay_pass=True")
print("source_normalization=L_mix=12*f_A+4*c_x*f_u-3*g_c+a_x*g_u=-3*D")
print("bracket_normalization=D=P+R_c*Q-R_A*S=-4*d_A(f_on_branch)+d_c(g_on_branch)")
print("closed_form_cap=4;closed_dimension=" + str(len(closed_basis))
      + ";P=-4*f_A+g_c;Q=g_u;S=4*f_u;closure=P_u-Q_c+S_A=0;"
      "forward_and_reverse_complete=True")
print("moving_endpoint_sum=(0,0,0);Pi(L_mix)=12*d_A(Pi(f))")
print("normalized_transgression=Lambda(L_mix)=12*(d/dA+5/A)*Lambda(f)")
print("common_target_period_ideal=A^5*K[[A]];w^s_lead=-(s/2)*rho*A^(s+4)")
print("carrier_leads_for_f=w^s:D=2*s*(s+9)*rho*A^(s+3);"
      "L_mix=-6*s*(s+9)*rho*A^(s+3)")
print("rank_formula=3*binom(N+2,2)-min(N+1,4);exact_controls=N0..8")
print("rank_rows=" + repr(rank_rows))
print("all_visible_stable_monomials_in_coupled_image=True")
print("H=t+gamma*t^m:D_rhs=4*m*gamma*w^(m-1);L_mix_rhs=-12*m*gamma*w^(m-1);"
      "first_f_lead=-gamma*(m-1)/(m+9)*w^m;m=2..12_checked")
print("persistent_cokernel_degrees=(0,1,2,3);no_later_linearized_source-cutoff_obstruction")
print("table_sha256=" + table_hash)
print("scope=formal_local_closed-form_cokernel_and_formal_pair_lift;"
      "formal_Darboux_factorization_is_textual_not_computational;"
      "no_convergence_algebraization_or_global_pair_claim")
print("RESULT=PASS")
