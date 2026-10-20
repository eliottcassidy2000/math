#!/usr/bin/env python3
"""Finite exact controls for the fixed and moving cylinder arguments.

Uses only explicit polynomial/rational identities. It does not test all
embeddings or supply a classification of source algebras. No Python assert.
"""
import argparse
import json
from pathlib import Path
import sympy as s

x, y, t, W, Z, c, kappa = s.symbols("x y t W Z c kappa", nonzero=True)
checks = 0
records = []


def eq(label, got, want=0):
    global checks
    checks += 1
    if s.cancel(got - want) != 0:
        raise RuntimeError(f"{label}: {s.factor(got-want)}")


def true(label, value):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)


def jac(F, G, u, v):
    return s.expand(s.diff(F, u)*s.diff(G, v)
                    - s.diff(F, v)*s.diff(G, u))


def moment(poly, endpoint):
    out = 0
    for (power,), coeff in s.Poly(s.expand(poly), y).terms():
        if power % 2 == 0:
            out += 2*coeff*endpoint**(power+1)/s.Integer(power+1)
    return s.expand(out)


# Symbolic coefficient and moving-endpoint identities for degrees 1..6.
# These test the algebra independently of any all-fiber embedding theorem.
for d in range(1, 7):
    aa = [c] + list(s.symbols(f"a1:{d+1}"))
    ap = [s.Integer(0)] + list(s.symbols(f"p1:{d+1}"))
    g = sum(aa[i]*W**i for i in range(d+1))
    gx = sum(ap[i]*W**i for i in range(1, d+1))
    gy = g.subs(W, y*y)
    z = y*gy
    J = -2*y*y*gx.subs(W, y*y)
    eq(f"d{d}: standard Jacobian", J,
       -2*sum(ap[i]*y**(2*i+2) for i in range(1, d+1)))
    I = moment(y*y*s.diff(z, y), t)
    I_formula = sum(2*s.Rational(2*i+1, 2*i+3)*aa[i]*t**(2*i+3)
                    for i in range(d+1))
    eq(f"d{d}: standard period", I, I_formula)
    tp = -gx.subs(W, t*t)/(2*t*s.diff(g, W).subs(W, t*t))
    total_I = sum(s.diff(I, aa[i])*ap[i] for i in range(1, d+1))
    total_I += s.diff(I, t)*tp
    expected = -4*sum(ap[i]*t**(2*i+3)/s.Integer(2*i+3)
                      for i in range(1, d+1))
    eq(f"d{d}: endpoint derivative defect",
       total_I-expected,
       -t*g.subs(W, t*t)*gx.subs(W, t*t)/s.diff(g, W).subs(W, t*t))
    eq(f"d{d}: integrated standard Jacobian", moment(J, t), expected)
    H = sum(ap[i]*W**(i+1)/s.Integer(2*i+3)
            for i in range(1, d+1))
    A, B = kappa*aa[1]/c**2, -kappa/c
    residual = s.Poly(s.expand(H-kappa-(A*W+B)*g), W)
    eq(f"d{d}: constant coefficient", residual.nth(0))
    eq(f"d{d}: linear coefficient", residual.nth(1))
    eq(f"d{d}: leading coefficient", residual.nth(d+1),
       ap[d]/s.Integer(2*d+3)-kappa*aa[1]*aa[d]/c**2)
    for i in range(2, d+1):
        eq(f"d{d}: coefficient {i}", residual.nth(i),
           ap[i-1]/s.Integer(2*i+1)-A*aa[i-1]-B*aa[i])
    records.append({"kind": "symbolic_degree", "degree": d,
                    "leading_exponent_prime": bool(s.isprime(2*d+1)),
                    "all_fiber_embedding_claimed_by_test": False})

# Independently eliminate the quadratic node values.
a, b, ap, bp = s.symbols("a b ap bp")
g = a*W*W+b*W+c
H = ap*W**3/7+bp*W*W/5
rem = s.rem(H, g, W)
eq("quadratic root difference coefficient", s.expand(rem).coeff(W),
   (5*ap*(b*b-a*c)-7*a*b*bp)/(35*a*a))
eq("quadratic common value",
   rem.subs(bp, 5*ap*(b*b-a*c)/(7*a*b)), ap*c*c/(7*a*b))

# Explicit fibers, global polynomial representatives, and relation changes.
# These are generating pairs because they are elementary ambient maps of
# (w,z), with added multiples of the defining equation.
families = [
    ("node_fixed", W+1),
    ("node_degree_drop", x*W+1),  # empty collision fiber x=0: boundary control
    ("two_nodes_moving", W*W-x*W+1),
    ("quadratic_leading_drop", x*W*W+W+1),
    ("cubic_moving", W**3+x*W*W+(x*x+1)*W+1),
    ("cubic_leading_drop", x*W**3+W*W+x*W+1),
    ("double_root_moving", (W-1)**2*(x*W-1)),
    ("constant_tacnode", (W-1)**2),
    ("constant_prime11", W**5+x*W+1),
    ("prime11_to_composite9_drop", x*W**5+W**4+1),
]
for name, g in families:
    h = Z*Z-W*g*g
    F0, G0 = W+x*Z*Z, Z+x*x
    eq(f"{name}: ambient determinant", jac(F0, G0, W, Z), 1)
    F = F0+(x+W)*h
    G = G0+(Z-x)*h
    z = y*g.subs(W, y*y)
    subs = {W: y*y, Z: z}
    P, Q = s.expand(F.subs(subs, simultaneous=True)), s.expand(G.subs(subs, simultaneous=True))
    eq(f"{name}: relation", h.subs(subs, simultaneous=True))
    eq(f"{name}: representative P", P, F0.subs(subs, simultaneous=True))
    eq(f"{name}: representative Q", Q, G0.subs(subs, simultaneous=True))
    D = jac(F, G, W, Z).subs(Z, 0)
    eq(f"{name}: node scalar descent", s.rem(D-1, g, W))
    eq(f"{name}: singular derivative w divisible g",
       s.rem(s.diff(h, W).subs(Z, 0), g, W))
    eq(f"{name}: singular derivative z", s.diff(h, Z).subs(Z, 0))
    # The normalization stays immersed at zero exactly because g(x,0) is a unit.
    eq(f"{name}: normalization derivative at zero",
       s.diff(z, y).subs(y, 0), g.subs(W, 0))
    true(f"{name}: constant term nonzero constant",
         not g.subs(W, 0).has(x) and g.subs(W, 0) != 0)
    records.append({"kind": "explicit_family", "name": name, "g": str(g)})

# Repeated-factor collision: an actual generating pair has opposite
# nonzero Jacobians at the two generic preimages.
g = (W-1)**2*(x*W-1)
z = y*g.subs(W, y*y)
P, Q = x+z, y*y
J = jac(P, Q, x, y)
for sign in (-1, 1):
    eq(f"repeated root g {sign}", g.subs(W, 1))
    eq(f"repeated root gx {sign}", s.diff(g, x).subs(W, 1))
    eq(f"repeated root gw {sign}", s.diff(g, W).subs(W, 1))
    eq(f"repeated root Jacobian {sign}", J.subs(y, sign), 2*sign)

# A root-cover unit is not necessarily a descended scalar.
g = W*W-x*W+1
eq("cover unit inverse", s.rem(W*(x-W)-1, g, W))
eq("cover unit nonconstant remainder", s.rem(W, g, W), W)
true("cover unit fails joint scalar condition", s.degree(s.rem(W, g, W), W) == 1)

# If a fiber loses *all* collisions, rational scalar poles need not disappear.
g = x*W+1
eq("empty fiber pole identity", x*(-W)-1, -g)
eq("empty fiber constant", g.subs(x, 0), 1)
eq("generic root pole", (-W).subs(W, -1/x), 1/x)

# Moving quadratic family: exact canonical periods and their derivatives.
r, rp = s.symbols("r rp", nonzero=True)
g = W*W-x*W+c
z = y*g.subs(W, y*y)
I = moment(y*y*s.diff(z, y), t)
formula = s.Rational(8, 105)*t**5*(3*t*t-7*(x-t*t))
eq("moving quadratic loop period",
   (I-formula).subs(c, x*t*t-t**4))
tp = t/(2*(2*t*t-x))  # differentiate r^2-x*r+c=0, r=t^2
derivative = s.diff(I, x)+s.diff(I, t)*tp
eq("moving quadratic loop derivative",
   (derivative-s.Rational(4, 5)*t**5).subs(c, x*t*t-t**4))
eq("moving quadratic Jacobian", jac(y*y, z, x, y), 2*y**4)

# Common-zero fiber supplies C[y], but the generic model is proper.
g = 1+x*W+x*x*W*W+x**3*W**3
z = y*g.subs(W, y*y)
eq("common-zero second generator", z.subs(x, 0), y)
true("generic algebra is proper criterion", s.degree(g, W) > 0)

# Fixed-cylinder period controls, independent of node normalization.
for power in (1, 2, 3):
    for scale in (-2, 1, 3):
        g = (W-1)**power
        z = y*g.subs(W, y*y)
        P, Q = y*y+x*z*z, scale*z+x*x
        A_period = moment(P*s.diff(Q, y), s.Integer(1))
        base_period = moment(y*y*s.diff(z, y), s.Integer(1))
        eq(f"fixed{power},{scale}: area transport", A_period, scale*base_period)
        eq(f"fixed{power},{scale}: zero area derivative", s.diff(A_period, x))
        eq(f"fixed{power},{scale}: zero Jacobian integral",
           moment(jac(P, Q, x, y), s.Integer(1)))
        records.append({"kind": "fixed_period", "multiplicity": power,
                        "scale": scale, "period": str(A_period)})

# Quadratic-coordinate recognition controls and the principal-conductor boundary.
g = W**3+x*W+1
even_multiplier, odd_multiplier = x+W*W, W+x*x
f = g.subs(W, y*y)*(even_multiplier.subs(W, y*y)
                    + y*odd_multiplier.subs(W, y*y))
f_even = s.expand((f+f.subs(y, -y))/2)
f_odd_coefficient = s.expand((f-f.subs(y, -y))/(2*y))
eq("recognition conductor even component",
   f_even, (g*even_multiplier).subs(W, y*y))
eq("recognition conductor odd component",
   f_odd_coefficient, (g*odd_multiplier).subs(W, y*y))
eq("recognition conductor times y odd component",
   s.expand((y*f-(y*f).subs(y, -y))/(2*y)), f_even)
b_shear, d_shear = x*x+x, x**3-1
U = y*y+b_shear*y+d_shear
Y = y+b_shear/2
eq("recognition completing square", U-d_shear+b_shear*b_shear/4, Y*Y)
eq("recognition source shear Jacobian", jac(x, Y, x, y), 1)
# I=(x,W) gives S=C[x,W] + y I, conductor (x,y²)B: not principal.
eq("nonprincipal conductor gcd", s.gcd(x, y*y), 1)
gb = s.groebner([x, y*y], x, y)
eq("nonprincipal conductor remains proper", gb.reduce(s.Integer(1))[1], 1)
true("nonprincipal conductor has two generators", len(gb.polys) == 2)
# An odd proposed conductor generator fails the even/odd module condition.
odd_C = y*(y*y-1)
even_intersection_generator = y*y*(y*y-1)
true("odd conductor cannot contain its coefficient",
     s.rem(y*y-1, odd_C, y) != 0)
true("odd conductor fails extension from even part",
     s.rem(odd_C, even_intersection_generator, y) != 0)
records.append({"kind": "recognition_controls",
                "nonprincipal_boundary": "I=(x,W), conductor=(x,y^2)B"})

parser = argparse.ArgumentParser()
parser.add_argument("--output", default="results.json")
args = parser.parse_args()
result = {
    "status": "PASS", "checks": checks, "families": records,
    "family_count": len(records),
    "universe": "symbolic degrees1..6; ten explicit moving/fixed source models; "
                "nine fixed-period transforms; distinct hostile boundary controls",
    "not_claimed": "No finite test proves uniqueness of all plane embeddings, "
                   "the all-pair theorem, or recognition of arbitrary source algebras.",
}
Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(f"PASS: {checks} exact checks; {len(records)} declared families.")
