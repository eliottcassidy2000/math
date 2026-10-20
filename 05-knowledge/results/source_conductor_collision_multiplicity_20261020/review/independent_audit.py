#!/usr/bin/env python3
"""Independent finite exact probes, not a search for noninvertible Keller maps.

Collision lengths use standard monomials of a Groebner basis, rather than
degrees of the resultant. Norms and Fitting determinants use multiplication
matrices. For coprime-degree curves, delta is counted from missing leading
degrees of the free power basis, rather than from the index determinant.
"""
from pathlib import Path
import json
import math
import sympy as s

x, y, t, u, v = s.symbols("x y t u v")
checks = []
records = []


def check(value, label):
    if not bool(value):
        raise RuntimeError(label)
    checks.append(label)


def zero(value):
    return s.cancel(value) == 0


def mult_matrix(f, modulus, var):
    n = int(s.degree(modulus, var))
    cols = []
    for j in range(n):
        rem = s.Poly(s.rem(s.expand(f * var**j), modulus, var), var)
        cols.append(s.Matrix([rem.nth(i) for i in range(n)]))
    return s.Matrix.hstack(*cols) if n else s.zeros(0, 0)


def quotient_length(gens):
    gb = s.groebner(gens, y, t, order="lex")
    lm = [p.LM(order=gb.order).exponents for p in gb.polys]
    if (0, 0) in lm:
        return 0
    y_bounds = [a for a, b in lm if b == 0 and a > 0]
    t_bounds = [b for a, b in lm if a == 0 and b > 0]
    if not y_bounds or not t_bounds:
        return None
    return sum(
        not any(i >= a and j >= b for a, b in lm)
        for i in range(min(y_bounds))
        for j in range(min(t_bounds))
    )


def collision(p, q):
    U = s.cancel((p.subs(y, t) - p) / (t - y)).expand()
    V = s.cancel((q.subs(y, t) - q) / (t - y)).expand()
    return U, V


# In C[u] + y*h(u)*C[u], membership and multiplication by y give the
# conductor h(u)*C[y] directly. B/S has basis y*u^j for 0 <= j < deg h.
for h in [u - 1, (u - 1)**2, (u - 1)**4, u, u**2,
          (u - 1)*(u - 4), (u - 1)**2*(u - 4)]:
    p, q = y**2, y*h.subs(u, y**2)
    c = h.subs(u, y**2)
    U, V = collision(p, q)
    deg_h = int(s.degree(h, u))
    tag = f"rank-two h={s.expand(h)}"
    check(zero(s.rem(V, U, t) - c), tag + ": collision quotient is C[y]/(c)")
    check(quotient_length([U, V]) == 2*deg_h, tag + ": direct collision length")
    check(zero(mult_matrix(c, p-u, y).det()-h**2), tag + ": matrix norm square")
    check(quotient_length([U, V, t-y]) == (1 if h.subs(u, 0) == 0 else 0),
          tag + ": diagonal contamination exactly when ramification is present")
    check(zero(s.rem(c, p-u, y)-h), tag + ": conductor even component")
    check(zero(s.rem(y*c, p-u, y)-y*h), tag + ": conductor odd component")
    records.append({"family": tag, "delta_direct": deg_h,
                    "collision_length": quotient_length([U, V]),
                    "conductor": str(s.factor(c))})


# Coprime leading degrees give distinct leading degrees n*i + m*j
# for i >= 0 and 0 <= j < n. This basis counts B/S without using trace duals.
for n, m in [(2, 3), (3, 2), (3, 4), (3, 5), (4, 3), (5, 2), (5, 3)]:
    check(math.gcd(n, m) == 1, f"coprime degrees {n},{m}")
    p, q = y**n, y**m+y
    U, V = collision(p, q)
    D = s.factor(s.resultant(U, V, t))
    C = mult_matrix(q, p-u, y)
    R = s.expand(C.charpoly(v).as_expr())
    c = s.cancel(s.diff(R, v).subs({u:p, v:q}, simultaneous=True)/s.diff(p, y))
    bound = n*m+n+m
    lead_degrees = {n*i+m*j for j in range(n)
                    for i in range(bound//n+1) if n*i+m*j <= bound}
    gaps = [k for k in range(bound+1) if k not in lead_degrees]
    delta = len(gaps)
    index_cols = []
    for j in range(n):
        rem = s.Poly(s.rem(q**j, p-u, y), y)
        index_cols.append(s.Matrix([rem.nth(i) for i in range(n)]))
    index = s.Matrix.hstack(*index_cols).det()
    tag = f"coprime degree curve {n},{m}"
    check(s.gcd(s.diff(p, y), s.diff(q, y)) == 1, tag + ": immersion")
    check(zero(R.subs({u:p, v:q}, simultaneous=True)), tag + ": matrix graph equation")
    check(zero(c-D), tag + ": derivative quotient equals D-resultant")
    check(zero(mult_matrix(V, U, t).det()-D), tag + ": actual Fitting determinant")
    check(zero(mult_matrix(c, p-u, y).det()-index**2), tag + ": matrix norm square")
    check(quotient_length([U, V]) == 2*delta, tag + ": Groebner length equals twice gap delta")
    check(int(s.degree(D, y)) == 2*delta, tag + ": D-resultant degree")
    check(quotient_length([U, V, t-y]) == 0, tag + ": no diagonal points")
    check(delta == (n-1)*(m-1)//2, tag + ": full coprime leading-degree gap count")
    records.append({"family": tag, "delta_direct": delta, "gaps": gaps,
                    "collision_length": quotient_length([U, V]),
                    "conductor": str(D)})


# The requested node control uses P=y^2-1, not y^2.
p, q = y*y-1, x+y*(y*y-1)
R = (v-x)**2-u*u*(u+1)
J = s.diff(p,x)*s.diff(q,y)-s.diff(p,y)*s.diff(q,x)
c = y*y-1
check(J == -2*y, "requested node: explicitly not Keller")
check(zero(R.subs({u:p,v:q}, simultaneous=True)), "requested node: relation")
check(zero(s.diff(R,x).subs({u:p,v:q}, simultaneous=True)-J*c),
      "requested node: raw graph derivative contains nonconstant Jacobian")
check(s.gcd(s.diff(p,y),s.diff(q,y)) == 1, "requested node: all slices immersed")
check(quotient_length(list(collision(p,q))) == 2, "requested node: collision length two")
check(s.det(s.Matrix([[s.diff(p,y).subs(y,1),s.diff(p,y).subs(y,-1)],
                      [s.diff(q,y).subs(y,1),s.diff(q,y).subs(y,-1)]])) == 8,
      "requested node: transverse branch tangents")


# Explicit hostile to unjustified conductor/image base change:
# S=C[x,y^2,xy], conductor xB; at x=0 the image is C[y^2],
# but the abstract tensor product is C[u,v]/(v^2), and c specializes to zero.
p, q = y*y, x*y
R = v*v-x*x*u
check(zero(R.subs({u:p,v:q}, simultaneous=True)), "base-change hostile: graph equation")
check(s.cancel(s.diff(R,v).subs({u:p,v:q}, simultaneous=True)/s.diff(p,y)) == x,
      "base-change hostile: global conductor generator x")
check(R.subs(x,0) == v*v, "base-change hostile: nonreduced abstract special fiber")
check(quotient_length(list(collision(p,q.subs(x,0)))) is None,
      "base-change hostile: nonfinite collision length after loss of primitivity")


# Scalar bookkeeping is exact only in the prescribed monic-v normalization.
p, q = y**3, y**4+y
C = mult_matrix(q, p-u, y)
R = s.expand(C.charpoly(v).as_expr())
c = s.cancel(s.diff(R,v).subs({u:p,v:q}, simultaneous=True)/s.diff(p,y))
for scalar in [-1, 2, 3]:
    check(zero(mult_matrix(scalar*c,p-u,y).det()
               -scalar**3*mult_matrix(c,p-u,y).det()),
          f"graph rescaling {scalar}: norm scales by scalar^n")

output = {
    "status": "PASS",
    "checks": len(checks),
    "labels": checks,
    "records": records,
    "scope": "Finite exact controls only. No noninvertible Keller map or conjecture closure.",
}
Path(__file__).with_name("independent_results.json").write_text(
    json.dumps(output, indent=2, sort_keys=True)+"\n"
)
print(json.dumps({"status":"PASS", "checks":len(checks), "families":len(records)}, sort_keys=True))
