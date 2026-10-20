#!/usr/bin/env python3
"""Finite exact controls for the cubic/quartic boundary packet.

These controls supplement, and do not replace, the proofs in the note.
Run normally and with python -O; no correctness test uses assert.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import sympy as sp

x, y, t, u, v = sp.symbols("x y t u v")
checks: dict[str, bool] = {}


def check(name: str, condition) -> None:
    if not bool(condition):
        raise RuntimeError(f"Failed exact check: {name}")
    checks[name] = True


def zero(expr) -> bool:
    return sp.cancel(sp.expand(expr)) == 0


def rem(expr, divisor, var=y):
    return sp.rem(sp.Poly(expr, var), sp.Poly(divisor, var)).as_expr()


def column(expr, n):
    p = sp.Poly(expr, y)
    return sp.Matrix([p.nth(i) for i in range(n)])


def multiplication(expr, relation, n):
    return sp.Matrix.hstack(
        *(column(rem(expr * y**j, relation), n) for j in range(n))
    )


def algebra_mod_conductor(f, g, conductor, name):
    """Generate A/(c) to closure in B/(c), with exact row-rank tests."""
    conductor = sp.Poly(conductor, y).monic().as_expr()
    n = sp.degree(conductor, y)
    basis = [sp.Integer(1)]
    vectors = column(1, n)
    queue = [sp.Integer(1)]
    while queue:
        b = queue.pop(0)
        for generator in (f, g):
            candidate = rem(b * generator, conductor)
            new_vectors = vectors.row_join(column(candidate, n))
            if new_vectors.rank() > vectors.cols:
                basis.append(candidate)
                queue.append(candidate)
                vectors = new_vectors
    descending = vectors.T[:, ::-1]
    _, pivots = descending.rref()
    leading_degrees = {n - 1 - p for p in pivots}
    gaps = [i for i in range(n) if i not in leading_degrees]
    check(name + ":closed_under_f", all(
        vectors.row_join(column(rem(f*b, conductor), n)).rank() == len(basis)
        for b in basis
    ))
    check(name + ":closed_under_g", all(
        vectors.row_join(column(rem(g*b, conductor), n)).rank() == len(basis)
        for b in basis
    ))
    return {"delta": int(n - len(basis)), "gaps": gaps}


def characteristic_patterns():
    semigroup = lambda n: any(
        3*a + 4*b == n for a in range(n//3 + 1) for b in range(n//4 + 1)
    )
    allowed = [i for i in range(1, 8) if semigroup(i)]
    found = []
    # The initial gcd is at most 7 and strictly falls after the first
    # nonprincipal pair, so length <=4 is an exhaustive bound here.
    for length in range(2, 5):
        for ds in itertools.product(allowed, repeat=length):
            if ds[0] <= ds[1] or ds[0] % ds[1] == 0:
                continue
            if any(a > 5 for a in ds[2:]) or not {3, 4}.issubset(ds):
                continue
            gcds = [ds[0]]
            for a in ds[1:]:
                gcds.append(math.gcd(gcds[-1], a))
            if gcds[-1] != 1:
                continue
            ns = [gcds[i-1] // gcds[i] for i in range(1, length)]
            if any(n <= 1 for n in ns):
                continue
            if any(ds[i] >= ds[i-1] * ns[i-2] for i in range(2, length)):
                continue

            def generated(n, generators):
                possible = {0}
                for j in range(n+1):
                    if j in possible:
                        possible.update(j+a for a in generators if j+a <= n)
                return n in possible

            if any(not generated(ns[i-1]*ds[i], ds[:i])
                   for i in range(1, length)):
                continue
            conductor = 1-ds[0] + sum(ns[i-1]*ds[i]-ds[i]
                                      for i in range(1, length))
            if conductor != 6:
                continue
            found.append(list(ds))
    return sorted(found)


def run():
    f, g = y**4, x*y**3+y**2+y
    relation = y**4-u
    Mg = multiplication(g, relation, 4)
    H = (v*sp.eye(4)-Mg).det().expand()
    expected_H = (
        v**4-(4*x+2)*u*v**2-4*u*v-4*x**2*u**2*v-x**4*u**3
        +(2*x**2-4*x+1)*u**2-u
    )
    check("hostile:characteristic_equation", zero(H-expected_H))
    check("hostile:equation_on_parametrization", zero(H.subs({u:f, v:g})))
    c = (x*y**2+1)*(x**2*y**4+2*x*y**3+2*y**2+2*y+1)
    check("hostile:trace_dual_conductor",
          zero(sp.diff(H, v).subs({u:f,v:g}) - sp.diff(f,y)*c))
    U = sp.cancel((t**4-y**4)/(t-y)).expand()
    V = sp.cancel(((x*t**3+t**2+t)-g)/(t-y)).expand()
    V_reduced = sp.rem(sp.Poly(V,t),sp.Poly(U,t)).as_expr()
    divided_resultant = sp.resultant(U,V_reduced,t)
    check("hostile:divided_difference_conductor", zero(divided_resultant-c))
    check("hostile:jacobian_not_keller",
          zero(sp.diff(f,x)*sp.diff(g,y)-sp.diff(f,y)*sp.diff(g,x)+4*y**6))
    check("hostile:derivative_at_only_possible_ramification",
          sp.diff(g,y).subs(y,0) == 1)
    check("hostile:derivative_resultant",
          sp.resultant(sp.diff(f,y),sp.diff(g,y),y) == 16)
    check("hostile:generic_conductor_degree", sp.degree(c,y) == 6)
    check("hostile:special_conductor",
          zero(c.subs(x,0)-(2*y**2+2*y+1)))
    check("hostile:special_conductor_squarefree",
          sp.gcd(c.subs(x,0),sp.diff(c.subs(x,0),y)) == 1)
    check("hostile:special_birational_identity",
          zero(f - ((y**2+y)**2+(y**2+y)-(2*(y**2+y)+1)*y)))
    q0 = y**2+y
    check("hostile:special_collision_Q", zero(q0-q0.subs(y,-1-y)))
    check("hostile:special_collision_P",
          zero(sp.rem(sp.Poly(f-f.subs(y,-1-y),y),
                      sp.Poly(2*y**2+2*y+1,y)).as_expr()))
    check("hostile:generic_opposite_pair_P", zero(f-f.subs(y,-y)))
    check("hostile:generic_opposite_pair_Q",
          zero(g-g.subs(y,-y)-2*y*(x*y**2+1)))
    check("hostile:translated_quadratic",
          zero(q0.subs(y,t-sp.Rational(1,2))-(t**2-sp.Rational(1,4))))
    even_part = t**4+sp.Rational(3,2)*t**2+sp.Rational(1,16)
    check("hostile:translated_odd_generator",
          zero(f.subs(y,t-sp.Rational(1,2))-even_part
               +2*t*(t**2+sp.Rational(1,4))))
    check("hostile:primitive_cubic_coefficient_gcd", sp.gcd(x,1) == 1)
    check("hostile:cubic_leading_coefficient_is_nonunit", sp.degree(x,x) == 1)
    check("hostile:generic_monic_cubic_pole",
          sp.Poly(g/x,y).nth(2) == 1/x and sp.Poly(g/x,y).nth(1) == 1/x)

    # Check c*y^j lies in S over R[u], j=0,...,3; higher powers follow
    # because y^4=u. This independently checks conductor containment.
    Qpowers = sp.Matrix.hstack(*(
        column(rem(g**j,relation),4) for j in range(4)
    ))
    for j in range(4):
        target = column(rem(c*y**j,relation),4)
        coefficients = Qpowers.inv()*target
        coefficients = coefficients.applyfunc(sp.cancel)
        check(f"hostile:conductor_module_{j}",
              all(sp.denom(a) == 1 for a in coefficients))
        check(f"hostile:conductor_module_identity_{j}",
              all(zero(a) for a in Qpowers*coefficients-target))
    index = sp.factor(Qpowers.det())
    norm_c = sp.factor(multiplication(c,relation,4).det())
    check("hostile:conductor_norm_index_square", zero(norm_c-index**2))

    fiber_records = {}
    for a in (sp.Integer(0),sp.Integer(1),sp.Integer(-1),sp.Integer(2),
              sp.Rational(1,2)):
        record = algebra_mod_conductor(f,g.subs(x,a),c.subs(x,a),f"fiber:{a}")
        check(f"fiber:{a}:delta", record["delta"] == (1 if a == 0 else 3))
        check(f"fiber:{a}:gaps", record["gaps"] == ([1] if a == 0 else [1,2,5]))
        fiber_records[str(a)] = record

    h, k = y**3, y**4+y
    z = 8*h**2+12*h+3
    check("embedding:entire_generation_identity", zero(z**2-64*k**3-9-8*h))
    F = v**3-u*(u+1)**3
    G = 4096*u**6-128*u**3*v**2+384*u**3+v**4-6*v**2-8*v-3
    check("embedding:standard_equation", zero(F.subs({u:h,v:k})))
    check("embedding:alternate_equation", zero(G.subs({u:k,v:z})))
    check("embedding:alternate_characteristic_equation",
          zero(sp.resultant(t**4+t-u,8*t**6+12*t**3+3-v,t)-G))
    check("embedding:immersion", sp.gcd(sp.diff(h,y),sp.diff(k,y)) == 1)
    check("embedding:standard_conductor",
          zero(sp.cancel(sp.diff(F,v).subs({u:h,v:k})/sp.diff(h,y))
               -(y**3+1)**2))
    for expr, name in ((F,"equation"),(sp.diff(F,u),"u_gradient"),
                       (sp.diff(F,v),"v_gradient")):
        check(f"embedding:singular_point_{name}", expr.subs({u:-1,v:0}) == 0)
    D = sp.Matrix([[0,1],[16*u+12,0]]).det()
    check("embedding:forced_local_determinant", D.subs(u,-1) == 4)
    L1 = sp.Rational(4,7)*y**7+sp.Rational(1,4)*y**4
    L2 = sp.Rational(24,5)*y**10+12*y**7+9*y**4
    check("embedding:standard_period_primitive", zero(sp.diff(L1,y)-h*sp.diff(k,y)))
    check("embedding:alternate_period_primitive", zero(sp.diff(L2,y)-k*sp.diff(z,y)))
    check("embedding:standard_period_on_triple",
          zero(rem(L1,y**3+1)-sp.Rational(9,28)*y))
    check("embedding:alternate_period_on_triple",
          zero(rem(L2,y**3+1)+sp.Rational(9,5)*y))
    ratio = (-sp.Rational(9,5))/sp.Rational(9,28)
    check("embedding:period_ratio", ratio == -sp.Rational(28,5))
    check("embedding:period_obstruction", ratio != 4)
    check("embedding:three_distinct_conductor_roots", sp.discriminant(y**3+1,y) != 0)
    alpha = sp.symbols("alpha")
    check("embedding:normalization_symmetry_h",
          sp.rem(sp.Poly(h.subs(y,alpha*y)-h,alpha),
                 sp.Poly(alpha**3-1,alpha)).as_expr() == 0)
    check("embedding:normalization_symmetry_k",
          sp.rem(sp.Poly(k.subs(y,alpha*y)-alpha*k,alpha),
                 sp.Poly(alpha**3-1,alpha)).as_expr() == 0)
    triple = algebra_mod_conductor(h,k,(y**3+1)**2,"triple_standard")
    alternate = algebra_mod_conductor(k,z,(y**3+1)**2,"triple_alternate")
    check("embedding:equal_finite_algebras",
          triple == alternate == {"delta":3,"gaps":[1,2,5]})
    patterns = characteristic_patterns()
    check("semigroup:bounded_characteristic_patterns", patterns == [[4,3],[6,4,3]])
    return {
        "status": "PASS: finite exact controls, not a substitute for the proofs",
        "sympy_version": sp.__version__,
        "check_count":len(checks),
        "checks":checks,
        "hostile":{
            "P":str(f),"Q":str(g),"relation":str(H),"conductor":str(c),
            "index":str(index),"fibers":fiber_records,"jacobian":"-4*y**6"
        },
        "embedding":{
            "generation_identity":"z^2-64*k^3-9=8*h",
            "forced_ambient_determinant":4,
            "period_ratio":str(ratio),
            "standard_degrees":[3,4],"alternate_degrees":[4,6]
        },
        "characteristic_patterns":patterns,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(result)
    print(result, end="")
