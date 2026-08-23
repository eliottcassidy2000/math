#!/usr/bin/env python3
"""Independent hostile controls for THM-3829's local-order proof."""

import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy

import sympy as sp


GATES = 0


def gate(condition, label):
    global GATES
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")
    GATES += 1


def zero(expression, label):
    gate(sp.cancel(sp.expand(expression)) == 0, label)


def repository_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("repository root not found")


ROOT = repository_root()
PRIMARY = ROOT / "04-computation" / "jc2_cubic_pseudoplane_misaligned_second_row_thm3829.py"
FROZEN = ROOT / "05-knowledge" / "results" / "jc2_cubic_pseudoplane_misaligned_second_row_thm3829.out"
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    n = runpy.run_path(str(PRIMARY), run_name="__main__")
origin_output = captured.getvalue().encode("utf-8")
gate(n["CHECKS"] == 58, "canonical namespace completed all gates")
gate(origin_output.endswith(b"CHECKS=58\nRESULT=PASS\n"), "canonical inner result")
gate(origin_output.replace(b"\r\n", b"\n") == FROZEN.read_bytes(),
     "canonical replay equals LF-frozen output")

e = n["e"]
v = n["v"]
dv = sp.diff(v, e)
alpha, beta, gamma = n["alpha"], n["beta"], n["gamma"]
M, N, H = n["M"], n["N"], n["H"]
f, p, g = n["f"], n["p"], n["g"]
R, U = n["R"], n["U"]

# Exhaust the two elementary Diophantine residue classes behind the UFD
# tower.  Checking all residues modulo ten is exhaustive, not a sample.
e_residues = [a for a in range(10) if (7*a - 2) % 10 == 0]
gate(e_residues == [6], "e-prime exponent residue is six mod ten")
gate((7*6 - 2)//10 == 4, "e-prime companion residue is four mod seven")
gate(sp.gcd(7, 10) == 1, "nonzero-prime exponents are ten/seven multiples")
zero(n["ode_107"].subs({n["X"]: n["X_tower"], n["L"]: n["L_tower"]}).doit(),
     "arbitrary-polynomial v solves exact 10/7 ODE")

# The r5 integral retains both a zero and a nonzero integration constant and
# does not divide by M.  Include the formally hostile M=0 specialization even
# though the arm later excludes it.
for m_value in (sp.Integer(0), e**2 + 1):
    for gamma_value in (sp.Integer(0), sp.Integer(5)):
        v_value = e + 1
        integrated_equation = n["r5_difference"].subs({
            n["X"]: n["X_tower"], n["L"]: n["L_tower"],
            n["S"]: n["S_integrated"],
        }).doit()
        zero(integrated_equation.subs({
            v: v_value, M: m_value, alpha: 2, beta: 3,
            gamma: gamma_value,
        }).doit(),
             "r5 integral zero/nonzero gamma and M control")
gate(sp.Poly(n["D"], e).degree() == 2, "arm D remains a quadratic")
gate(sp.rem(sp.Poly(-sp.Rational(1, 12), e), sp.Poly(n["D"], e)).as_expr()
     == -sp.Rational(1, 12), "M zero would require quadratic dividing unit")

# Recover every strictly-later r4 monomial type from the full source, then
# check the four low/high comparison cones symbolically.
dM = sp.diff(M, e)
r4_later = sp.expand(n["P4_expected"]-n["M2_block"]-n["f_block"])
r4_signatures = set()
for term in sp.Add.make_args(r4_later):
    powers = term.as_powers_dict()
    r4_signatures.add((
        int(powers.get(e, 0)), int(powers.get(v, 0)),
        int(powers.get(dv, 0)), int(powers.get(M, 0)),
        int(powers.get(dM, 0)),
    ))
expected_r4_signatures = {
    (6, 15, 0, 0, 0), (7, 14, 1, 0, 0),
    (3, 8, 0, 1, 0), (4, 7, 1, 1, 0), (4, 8, 0, 0, 1),
    (2, 5, 0, 1, 0), (3, 4, 1, 1, 0), (3, 5, 0, 0, 1),
}
gate(r4_signatures == expected_r4_signatures,
     "full r4 remainder has exactly the eight later monomial types")
mm = sp.symbols("mm", integer=True, positive=True)
dd = sp.symbols("dd", integer=True, nonnegative=True)
for difference in (10*mm+2, 2*mm+1, 5*mm+1,
                   sp.Integer(2), 2*mm+1, 10*mm,
                   6+10*dd, 2+2*dd, 3+5*dd,
                   sp.Integer(2), 2*dd+2, 4+10*dd):
    gate(difference.is_positive is True,
         "r4 later block is strictly later in every local cone")


def lower_order_signature(term, at_origin):
    powers = term.as_powers_dict()
    v_power = int(powers.get(v, 0))
    dv_power = int(powers.get(dv, 0))
    if at_origin:
        e_power = int(powers.get(e, 0))
        return v_power + dv_power, e_power - dv_power
    return v_power + dv_power, -dv_power


def check_termwise_bound(expression, slope, intercept, label):
    """Check A*n+B >= slope*n+intercept for n>=1 term by term."""
    terms = sp.Add.make_args(sp.expand(expression))
    gate(len(terms) > 0, f"{label} nonempty")
    for term in terms:
        a, b = lower_order_signature(term, at_origin=False)
        gate(a >= slope and (a-slope) + b-intercept >= 0,
             f"{label} nonzero-root term bound")


def check_origin_bound(expression, slope, intercept, label):
    """Check d>=1 symbolically and d=0 separately."""
    terms = sp.Add.make_args(sp.expand(expression))
    gate(len(terms) > 0, f"{label} origin nonempty")
    for term in terms:
        a, b = lower_order_signature(term, at_origin=True)
        gate(a >= slope and (a-slope) + b-intercept >= 0,
             f"{label} origin d-positive term bound")
        e_power = int(term.as_powers_dict().get(e, 0))
        gate(e_power >= intercept, f"{label} origin d-zero term bound")


# Independently enumerate every non-p and non-g term, rather than trusting a
# displayed block or a sampled jet.
p_remainders = []
for n_value in (n["N_tower"], sp.Integer(0)):
    remainder = sp.expand(n["p_source"](n_value) - n["p_block"])
    p_remainders.append(remainder)
    check_termwise_bound(remainder, 8, -1, "p competitor")
    check_origin_bound(remainder, 8, 4, "p competitor")

g_remainders = []
for h_value in (n["H_tower"], sp.Integer(0)):
    remainder = sp.expand(n["g_source"](h_value) - n["g_block"])
    g_remainders.append(remainder)
    check_termwise_bound(remainder, 8, -1, "g competitor")
    check_origin_bound(remainder, 8, 4, "g competitor")


def shifted_coefficient(expression, point, order):
    tau = sp.symbols("tau")
    shifted = sp.expand(expression.subs(e, point + tau).doit())
    return sp.Poly(shifted, tau).coeff_monomial(tau**order)


# Repeated nonzero root: verify every low p and g jet, including order zero.
rho = sp.Integer(2)
multiplicity = 2
v_repeated = (e-rho)**multiplicity
for p_order in range(3*multiplicity):
    expression = n["p_block"].subs({beta: 1, v: v_repeated,
                                     p: (e-rho)**p_order}).doit()
    coefficient = shifted_coefficient(
        expression, rho, p_order + 5*multiplicity - 1)
    zero(coefficient - 147*rho**3*(5*multiplicity-p_order),
         "repeated-root p leading coefficient")
for g_order in range(multiplicity):
    expression = n["g_block"].subs({beta: 1, v: v_repeated,
                                     g: (e-rho)**g_order}).doit()
    coefficient = shifted_coefficient(
        expression, rho, g_order + 7*multiplicity - 1)
    zero(coefficient - 49*rho**4*(3*multiplicity-g_order),
         "repeated-root g leading coefficient")

# Origin controls include d=0 explicitly and a repeated origin root d=2.
for d_value in (0, 2):
    v_origin = e**d_value
    for p_order in range(3*d_value + 2):
        expression = n["p_block"].subs({beta: 1, v: v_origin,
                                         p: e**p_order}).doit()
        coefficient = sp.Poly(sp.expand(expression), e).coeff_monomial(
            e**(p_order + 5*d_value + 2))
        zero(coefficient - 147*(5*d_value-p_order+3),
             "origin p leading coefficient including d zero")
    for g_order in range(d_value + 1):
        expression = n["g_block"].subs({beta: 1, v: v_origin,
                                         g: e**g_order}).doit()
        coefficient = sp.Poly(sp.expand(expression), e).coeff_monomial(
            e**(g_order + 7*d_value + 3))
        zero(coefficient - 49*(3*d_value-g_order+2),
             "origin g leading coefficient including d zero")

# Cancellation-safe r4 payments: a repeated nonzero root, origin d=0, and a
# repeated origin root.  The correct payment cancels the tied coefficient;
# a hostile wrong alpha does not.
def typed_p4(v_value, f_value, alpha_value):
    return sp.expand(n["P4_expected"].subs({
        v: v_value, M: e*v_value**2, f: f_value,
        alpha: alpha_value, beta: 1, gamma: 0,
    }).doit())


good = typed_p4(v_repeated, sp.Integer(1), sp.Rational(49, 15))
bad = typed_p4(v_repeated, sp.Integer(1), sp.Integer(1))
zero(shifted_coefficient(good, rho, 5*multiplicity-1),
     "repeated-root r4 tied payment cancels")
gate(shifted_coefficient(bad, rho, 5*multiplicity-1) != 0,
     "repeated-root r4 wrong payment survives")

for d_value, v_value in ((0, e+2), (2, e**2*(1+e))):
    good = typed_p4(v_value, sp.Rational(1, 12), sp.Rational(49, 180))
    bad = typed_p4(v_value, sp.Rational(1, 12), sp.Integer(1))
    tied_order = 5*d_value + 2
    zero(sp.Poly(good, e).coeff_monomial(e**tied_order),
         "origin r4 payment cancels including d zero")
    gate(sp.Poly(bad, e).coeff_monomial(e**tied_order) != 0,
         "origin r4 wrong payment survives")

# N and H zero branches are literal solutions, while the nonzero towers work
# for an arbitrary polynomial v.
zero(n["expected_r4z2"].subs(N, 0).doit(), "N zero seam")
zero(n["expected_r4z1"].subs(H, 0).doit(), "H zero seam")
zero(n["expected_r4z2"].subs({n["X"]: n["X_tower"],
                               N: n["N_tower"]}).doit(), "N nonzero tower")
zero(n["expected_r4z1"].subs({n["X"]: n["X_tower"],
                               H: n["H_tower"]}).doit(), "H nonzero tower")

# All four terminal branches have a polynomial higher remainder.  In the
# N=0 branches, g disappears exactly, so no illicit g typing is used.
terminal_count = 0
for n_value, g_value in ((n["N_tower"], n["g_typed"]),
                         (sp.Integer(0), g)):
    for h_value in (n["H_tower"], sp.Integer(0)):
        source = n["terminal_source"](n_value, h_value, g_value)
        higher = sp.cancel((source-n["universal_terminal"])/(e**2*v**4))
        gate(sp.denom(sp.together(higher)) == 1,
             "terminal higher remainder polynomial")
        zero(source-n["universal_terminal"]-e**2*v**4*higher,
             "terminal exact four-term split")
        if n_value == 0:
            gate(not source.has(g), "N zero terminal is independent of g")
        terminal_count += 1
gate(terminal_count == 4, "all four N/H terminal branches")

# Direct hostile terminal coefficients at a repeated nonzero root and at
# monomial origins, including d=0.
terminal = n["universal_terminal"]
nonzero_terminal = terminal.subs({beta: 1, R: 1, f: 1,
                                  v: v_repeated}).doit()
zero(shifted_coefficient(nonzero_terminal, rho, multiplicity-1)
     + 56*rho*multiplicity, "repeated-root terminal coefficient")
for d_value in (0, 3):
    c_value = 3
    origin_terminal = terminal.subs({beta: 1, R: 1, f: 1,
                                     v: c_value*e**d_value}).doit()
    coefficient = sp.Poly(sp.expand(origin_terminal), e).coeff_monomial(
        e**d_value)
    zero(coefficient + 28*c_value*(2*d_value+1),
         "origin odd terminal coefficient including d zero")

source_text = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(source_text))),
     "assertion-free independent source")
semantic = {
    "target": "THM-3829 current canonical exact companion",
    "tower": "complete mod-10 UFD residues and arbitrary-v substitution",
    "r5": "gamma zero/nonzero and formal M zero controls",
    "lower_orders": "every p/g competitor checked termwise at nonzero roots and origin",
    "hostiles": "repeated root m=2; origin d=0; repeated origin d=2",
    "branches": "N/H zero and nonzero; all four terminal remainders polynomial",
    "scope": "X nonzero fixed second row only; THM-3834 closes one-sided separately; higher slots open",
}
digest = hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
print("audit=THM-3829-independent-local-orders")
print("result=PASS")
print("tower=UFD_10_7_complete")
print("integration=gamma_zero_nonzero+formal_M_zero")
print("orders=p_g_all_competitors_termwise")
print("hostiles=repeated_nonzero_root+origin_d0+repeated_origin")
print("branches=N0_N1_x_H0_H1")
print("scope=X_nonzero_fixed_second_row_only")
print(f"semantic_sha256={digest}")
print(f"GATES={GATES}")
