#!/usr/bin/env python3
"""Independent exact controls for the quadratic period-three index budget.

Run: python audit_period_three_index.py
Requires: Python 3 and SymPy. No network, checkout, or Git operations.
"""

import sympy as s

x, c, v, m, t, lam = s.symbols("x c v m t lam")
checks = []


def check(label, actual, expected=0):
    delta = s.cancel(s.expand(actual - expected))
    if delta != 0:
        raise RuntimeError(f"FAIL {label}: {actual} != {expected}; difference={delta}")
    checks.append(label)


f = x**2 + c
F = f
for _ in range(2):
    F = s.expand(F.subs(x, f))
q = f - x
Phi, remainder = s.div(F - x, q, x)
check("exact dynatomic quotient", remainder)
D = 16*c**2 + 4*c + 7
check("fixed-dynatomic resultant", s.resultant(q, Phi, x), D)
check("dynatomic discriminant", s.discriminant(Phi, x), -(4*c + 7)**3 * D**2)
check("full fixed-locus discriminant", s.discriminant(F-x, x), (4*c-1)*(4*c+7)**3*D**4)

# The fixed-point multipliers are lam and 2-lam, with 4c=lam(2-lam).
S = 1/(1-lam**3) + 1/(1-(2-lam)**3)
check("simple fixed-point index sum", S, (6/D).subs(c, lam*(2-lam)/4))
check("removable index value at c=1/4", (6/D).subs(c, s.Rational(1, 4)), s.Rational(2, 3))
check("actual merged fixed-point residue at c=1/4", s.residue(1/(x-F.subs(c, s.Rational(1, 4))), x, s.Rational(1, 2)), s.Rational(2, 3))

# All parameter calculations below are in Q(v)[x].
cv = -(v**2 + 7)/4
sigma = (v-1)/2
Pv = s.expand(x**3 - sigma*x**2 - (sigma**2 + 2*sigma + 3)*x + sigma**3 + 2*sigma**2 + 3*sigma + 1)
Pm = Pv.subs(v, -v)
A_v = v**2 + v + 7
mu_v = 1 - v**2 - v*(7+v**2)
check("two trace-cubic factorization", Phi.subs(c, cv), Pv*Pm)
check("trace-cubic root multiplier", -8*Pv.subs(x, 0), mu_v)
check("oriented-gap resultant", s.resultant(q.subs(c, cv), Pv, x), A_v)
check("trace-cubic discriminant", s.discriminant(Pv, x), A_v**2)
check("multiplier-gap factorization", 1-mu_v, v*A_v)
check("resultant pair factorization", D.subs(c, cv), A_v*A_v.subs(v, -v))
check("two-cycle index sum", 1/(1-mu_v) + 1/(1-mu_v.subs(v, -v)), -2/(v**4+13*v**2+49))
check("generic finite index budget", (6/D).subs(c, cv) + 3*(1/(1-mu_v)+1/(1-mu_v.subs(v, -v))))
check("two-cycle coalescence limit", s.limit(1/(1-mu_v) + 1/(1-mu_v.subs(v, -v)), v, 0), -s.Rational(2, 49))
cycle_multiplier_polynomial = m**2 - 8*(c+2)*m + 64*(c**3+2*c**2+c+1)
check("multiplier resultant", s.resultant(Phi, m-s.diff(F, x), x), cycle_multiplier_polynomial**3)

# The distinguished cycle has three distinct roots, none fixed or critical.
star = -s.Rational(7, 4)
P = x**3 + x**2/2 - s.Rational(9, 4)*x - s.Rational(1, 8)
qs = q.subs(c, star)
check("parabolic dynatomic square", Phi.subs(c, star), P**2)
check("parabolic resultant sign", s.resultant(qs, P, x), 7)
check("parabolic cycle discriminant", s.discriminant(P, x), 49)
check("parabolic cycle multiplier", -8*P.subs(x, 0), 1)
check("no fixed roots in parabolic cubic", s.gcd(qs, P), 1)
check("cubic invariance under f", s.rem(P.subs(x, f.subs(c, star)), P, x))
check("global parabolic cycle index", -(6/D).subs(c, star)/3, -s.Rational(2, 49))

# Direct independent local residue check. F(r+y)=r+y+a*y^2+b*y^3+...
# Thus Res dy/(y-(F(r+y)-r))=b/a^2, not -b/a^2.
a = qs*s.diff(P, x)**2
b = s.diff(qs, x)*s.diff(P, x)**2 + qs*s.diff(P, x)*s.diff(P, x, 2)
check("nonzero local quadratic coefficients", s.gcd(a, P), 1)
check("direct local index -2/49", s.rem(49*b+2*a**2, P, x))
check("iterative residue convention", 1-(-s.Rational(2, 49)), s.Rational(51, 49))

# Fixed-chart real saddle-node unfolding, with epsilon=c-c_star>0.
alpha = s.rem(s.diff(F, c).subs(c, star), P, x)
check("parameter unfolding derivative", alpha, 2*x**2-s.Rational(11, 2))
check("normalized positive unfolding constant", s.rem(a*alpha, P, x), 49)

# Inherited cyclotomic trace and its exact orientation.
M = t**3+t**2-2*t-1
check("seventh-root affine identification", P.subs(x, -t-s.Rational(1, 2)), -M)
check("conjugated quadratic", -s.Rational(1, 2)-f.subs(c, star).subs(x, -t-s.Rational(1, 2)), 1-t-t**2)
check("tripling equals conjugated quadratic on trace cubic", s.rem(t**3-3*t-(1-t-t**2), M, t))

print(f"PASS: {len(checks)} exact controls (SymPy {s.__version__}).")
for label in checks:
    print(f"  PASS {label}")
print("Scope: rational identities extend meromorphically; separate simple-cycle formulas")
print("exclude c=1/4, c=-7/4, and D(c)=0 as appropriate. The global distinct-point")
print("residue budget itself is valid for every complex c.")
