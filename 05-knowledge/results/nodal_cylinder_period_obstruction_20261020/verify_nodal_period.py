#!/usr/bin/env python3
"""Exact finite controls for the nodal generating-pair period obstruction.

The mathematical input about *all* embeddings is cited in the proof, not
established by this finite test. No Python assert is used.
"""
import argparse
import json
from pathlib import Path
import sympy as s

x, y, W, Z = s.symbols("x y W Z")
w, z = y**2 - 1, y * (y**2 - 1)
checks = 0
records = []


def equal(label, actual, expected=0):
    global checks
    checks += 1
    if s.expand(actual - expected) != 0:
        raise RuntimeError(f"{label}: {actual!s} != {expected!s}")


def truth(label, value):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)


def operation(F, G, op, reverse=False):
    kind, coefficient, power = op
    sign = -1 if reverse else 1
    if kind == "shear_first":
        F = F + sign * coefficient * G**power
    elif kind == "shear_second":
        G = G + sign * coefficient * F**power
    elif kind == "translate_first":
        F = F + sign * coefficient
    elif kind == "translate_second":
        G = G + sign * coefficient
    elif kind == "scale_first":
        F = F / coefficient if reverse else coefficient * F
    elif kind == "scale_second":
        G = G / coefficient if reverse else coefficient * G
    elif kind == "swap":
        F, G = G, F
    else:
        raise RuntimeError(f"Unknown operation: {kind}")
    return s.expand(F), s.expand(G)


def period(P, Q):
    return s.integrate(s.expand(P * s.diff(Q, y)), (y, -1, 1))


def jac(P, Q, a=x, b=y):
    return s.expand(s.diff(P, a) * s.diff(Q, b)
                    - s.diff(P, b) * s.diff(Q, a))


def run_family(name, operations):
    F, G = W, Z
    expected_det = s.Integer(1)
    for op in operations:
        F, G = operation(F, G, op)
        if op[0] in ("scale_first", "scale_second"):
            expected_det *= op[1]
        elif op[0] == "swap":
            expected_det *= -1
    inverse_F, inverse_G = F, G
    for op in reversed(operations):
        inverse_F, inverse_G = operation(inverse_F, inverse_G, op, True)
    equal(f"{name}: inverse first", inverse_F, W)
    equal(f"{name}: inverse second", inverse_G, Z)
    equal(f"{name}: ambient determinant", jac(F, G, W, Z), expected_det)
    P = s.expand(F.subs({W: w, Z: z}, simultaneous=True))
    Q = s.expand(G.subs({W: w, Z: z}, simultaneous=True))
    equal(f"{name}: P endpoint", P.subs(y, 1), P.subs(y, -1))
    equal(f"{name}: Q endpoint", Q.subs(y, 1), Q.subs(y, -1))
    equal(f"{name}: first node value", P.subs(y, 1),
          F.subs({W: 0, Z: 0}))
    equal(f"{name}: second node value", Q.subs(y, 1),
          G.subs({W: 0, Z: 0}))
    L = s.Matrix([[s.diff(F, W), s.diff(F, Z)],
                  [s.diff(G, W), s.diff(G, Z)]]).subs({W: 0, Z: 0})
    equal(f"{name}: conormal determinant", L.det(), expected_det)
    A = period(P, Q)
    equal(f"{name}: area", A, s.Rational(8, 15) * expected_det)
    J = jac(P, Q)
    equal(f"{name}: integral of J", s.integrate(J, (y, -1, 1)))
    equal(f"{name}: integration by parts",
          s.diff(A, x) - s.integrate(J, (y, -1, 1))
          - (P * s.diff(Q, x)).subs(y, 1)
          + (P * s.diff(Q, x)).subs(y, -1))
    truth(f"{name}: not nonzero constant J", J == 0 or J.has(x, y))
    # Exactness of F dG - det(L) W dZ by explicit antidifferentiation.
    coeff_W = s.expand(F * s.diff(G, W))
    coeff_Z = s.expand(F * s.diff(G, Z) - expected_det * W)
    H0 = s.integrate(coeff_W, W)
    remainder = s.expand(coeff_Z - s.diff(H0, Z))
    truth(f"{name}: closed remainder independent of W", not remainder.has(W))
    H = H0 + s.integrate(remainder, Z)
    equal(f"{name}: exact W coefficient", s.diff(H, W), coeff_W)
    equal(f"{name}: exact Z coefficient", s.diff(H, Z), coeff_Z)
    pull_H = s.expand(H.subs({W: w, Z: z}, simultaneous=True))
    equal(f"{name}: exact loop zero", pull_H.subs(y, 1), pull_H.subs(y, -1))
    records.append({
        "name": name, "P_y_degree": s.degree(P, y),
        "Q_y_degree": s.degree(Q, y), "conormal_det": str(expected_det),
        "area": str(A), "jacobian": str(J),
        "operations": [(k, str(c), p) for k, c, p in operations],
    })


families = [("standard", [])]
for d in (1, 2, 3):
    families.extend([
        (f"first_shear_{d}", [("shear_first", x**d, 2)]),
        (f"second_shear_{d}", [("shear_second", x**d + 1, 3)]),
        (f"translated_{d}", [("translate_first", x**d, 0),
                             ("translate_second", x**(d+1), 0)]),
        (f"reversed_scaled_{d}", [("swap", 0, 0),
                                  ("scale_first", s.Integer(-d), 0),
                                  ("scale_second", s.Integer(d+1), 0),
                                  ("translate_second", x**d, 0)]),
    ])
families.extend([
    ("two_quadratic_shears", [("shear_first", x, 2),
                              ("shear_second", x + 1, 2)]),
    ("shift_then_shear", [("translate_first", x, 0),
                          ("shear_second", x - 1, 2),
                          ("translate_second", x*x, 0)]),
    ("shear_swap_shear", [("shear_first", x, 2),
                          ("swap", 0, 0),
                          ("shear_first", x*x + 1, 2)]),
    ("hostile_known_non_keller", [("translate_second", x, 0)]),
])
for name, ops in families:
    run_family(name, ops)

# Strongest-survivor tests: distinguish global generation, generic generation,
# and the conormal test. These are not positive generating families.
equal("standard relation", z*z - w*w*(w+1))
equal("standard period", period(w, z), s.Rational(8, 15))
equal("known non-Keller J", jac(w, x+z), -2*y)
equal("generic-only period", period(w, x*z), s.Rational(8, 15)*x)
equal("generic-only conormal determinant",
      s.det(s.Matrix([[1, 0], [0, x]])), x)
equal("generic-only fails generation at zero", (x*z).subs(x, 0), 0)

P_bad, Q_bad = w, z + x*w*w*z
equal("unit-conormal nongenerator area",
      period(P_bad, Q_bad), s.Rational(8, 15)+s.Rational(64, 315)*x)
equal("unit-conormal nongenerator determinant",
      jac(W, Z+x*W*W*Z, W, Z).subs({W: 0, Z: 0}), 1)
equal("extra pair same P at x=-1",
      s.rem(P_bad-1, y*y-2, y))
equal("extra pair same Q at x=-1",
      s.rem(Q_bad.subs(x, -1), y*y-2, y))
equal("extra pair distinct standard z",
      z.subs(y, s.sqrt(2))-z.subs(y, -s.sqrt(2)), 2*s.sqrt(2))
truth("extra pair is not standard node", w.subs(y, s.sqrt(2)) != 0)
equal("nongenerator integration by parts",
      s.diff(period(P_bad, Q_bad), x)
      - s.integrate(jac(P_bad, Q_bad), (y, -1, 1)))

parser = argparse.ArgumentParser()
parser.add_argument("--output", default="results.json")
args = parser.parse_args()
result = {
    "status": "PASS",
    "checks": checks,
    "generating_families": len(families),
    "universe": "17 declared compositions of elementary C[x]-plane automorphisms; "
                "two nongenerator controls and the standard relation",
    "not_claimed": "Finite checks do not prove all-embedding uniqueness or the "
                   "all-generating-pairs theorem.",
    "families": records,
    "hostile_controls": {
        "generic_only": {"P": str(w), "Q": str(x*z), "area": "8*x/15"},
        "unit_conormal_nongenerator": {
            "P": str(P_bad), "Q": str(Q_bad),
            "area": "8/15 + 64*x/315",
            "collision": "x=-1, y=+sqrt(2) and y=-sqrt(2)",
        },
    },
}
Path(args.output).write_text(
    json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
print(f"PASS: {checks} exact checks; {len(families)} generating families.")
