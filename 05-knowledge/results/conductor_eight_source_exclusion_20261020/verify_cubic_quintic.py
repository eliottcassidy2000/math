#!/usr/bin/env python3
"""Universal symbolic checks of the direct terminal-(3,5) proof."""
import json
import sympy as S

checks = []


def check(label, okay):
    if not okay:
        raise RuntimeError("FAILED: " + label)
    checks.append(label)


def zero(label, expr):
    check(label, S.expand(expr) == 0)


t, x = S.symbols("t x")
a, b, c, d, e, f, g = S.symbols("a b c d e f g")
ap, bp, cp, dp, ep, fp, gp = S.symbols("ap bp cp dp ep fp gp")
E, F, G = S.symbols("E F G")
variables = [a, b, c, d, e, f, g]
derivatives = [ap, bp, cp, dp, ep, fp, gp]
Dx = lambda h: sum(S.diff(h, v)*vp for v, vp in zip(variables, derivatives))
p = t**3+a*t+b
q = t**5+c*t**4+d*t**3+e*t*t+f*t+g
J = S.Poly(Dx(p)*S.diff(q, t)-S.diff(p, t)*Dx(q), t)
raw = [-a*gp+b*0+bp*f, -a*fp+ap*f+2*bp*e,
       -a*ep+2*ap*e+3*bp*d-3*gp,
       -a*dp+3*ap*d+4*bp*c-3*fp,
       -a*cp+4*ap*c+5*bp-3*ep,
       5*ap-3*dp, -3*cp]
for n, expected in enumerate(raw):
    zero("raw coefficient t^" + str(n), J.nth(n)-expected)

normalized = {
    c: c,
    d: 5*a/S.Integer(3),
    e: 4*c*a/S.Integer(3)+5*b/S.Integer(3)+E,
    f: 5*a*a/S.Integer(9)+4*c*b/S.Integer(3)+F,
    g: 2*c*a*a/S.Integer(9)+10*a*b/S.Integer(9)+2*E*a/S.Integer(3)+G,
}
replacements = dict(normalized)
for v, vp in zip([c, d, e, f, g], [cp, dp, ep, fp, gp]):
    value = normalized[v]
    replacements[vp] = S.diff(value, a)*ap+S.diff(value, b)*bp
Jn = S.Poly(J.as_expr().subs(replacements, simultaneous=True), t)
for n in range(2, 7):
    zero("all integrated coefficients t^" + str(n), Jn.nth(n))
U = -5*a*a/S.Integer(9)+4*c*b/S.Integer(3)+F
V = 4*c*a/S.Integer(3)+10*b/S.Integer(3)+2*E
H = -5*a**3/S.Integer(27)+4*c*a*b/S.Integer(3)+F*a+5*b*b/S.Integer(3)+2*E*b
zero("full reduced Jacobian",
     Jn.as_expr()-((U*ap+V*bp)*t+U*bp-a*V*ap/3))
zero("first integral derivative",
     S.diff(H, a)*ap+S.diff(H, b)*bp-U*ap-V*bp)

# Symbolic coefficient for arbitrary positive k; no bounded-test inference.
k, A, B = S.symbols("k A B", nonzero=True)
zero("universal residual leading coefficient",
     (-5*A*A/S.Integer(9))*(3*k*B)
     -(A*(10*B/S.Integer(3))/3)*(2*k*A)
     +35*k*A*A*B/S.Integer(9))

# Hostile exact family: constant H with nonconstant a,b.
for degree in [1, 2, 3, 5, 9]:
    vals = {a: x**(2*degree), b: x**(3*degree)/3, c: 0, E: 0, F: 0, G: 0}
    P = p.subs(vals)
    Q = q.subs(normalized, simultaneous=True).subs(vals)
    zero("constant-H hostile k=" + str(degree), H.subs(vals))
    J_actual = S.diff(P, x)*S.diff(Q, t)-S.diff(P, t)*S.diff(Q, x)
    zero("hostile full Jacobian k=" + str(degree),
         J_actual+S.Rational(35, 27)*degree*x**(7*degree-1))
    check("hostile actual y-degrees k=" + str(degree),
          (S.degree(P, t), S.degree(Q, t)) == (3, 5))

# The trivial stationary family has Jacobian zero, not a nonzero constant.
zero("stationary family", Jn.as_expr().subs({ap: 0, bp: 0}))

print(json.dumps({"status": "PASS", "count": len(checks), "checks": checks,
                  "scope": "Monic constant-leading y-degree (3,5) Keller obstruction",
                  "not_claimed": "No conclusion about other terminal degree patterns"},
                 indent=2, sort_keys=True))
