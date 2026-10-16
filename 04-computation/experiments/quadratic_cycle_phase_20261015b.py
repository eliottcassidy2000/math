#!/usr/bin/env python3
"""Exact controls for rational-cycle phase and changing-generator Pell growth."""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as s

def need(p, message):
    if not p:
        raise AssertionError(message)

x, u = s.symbols('x u')
B = s.Matrix([[0, -1], [1, 1]])
U = s.Matrix([[1, u], [0, 1]])
C = U.inv()*B.inv()*U
K = s.Matrix([[2, 1], [1, 2]])
need(B**3 == -s.eye(2) and B**6 == s.eye(2), 'order-six lift')
need(s.simplify(C**3) == -s.eye(2) and s.simplify(C**6) == s.eye(2), 'second order-six lift')
need(B.T*K*B == K, 'positive invariant quadratic form')
KC = U.T*K*U
need(s.simplify(C.T*KC*C-KC) == s.zeros(2), 'second invariant form with changed coordinate frame')
product = s.simplify(B*C)
need(product == s.Matrix([[1, u], [u, u*u+1]]), 'changing-generator composition')
need(product.det() == 1 and product.trace() == u*u+2, 'hyperbolic trace outside u=0')
change = s.Matrix([[1, -1], [0, 1]])
Pell = change.inv()*product.subs(u, 2)*change
need(Pell == s.Matrix([[3, 4], [2, 3]]), 'exact Pell conjugacy')
# Residual covariance, rather than stationary covariance, detects pointwise constraints.
A = s.Matrix([[0, -1, 0], [0, 0, -1], [-1, 0, 0]])
need(A**3 == -s.eye(3) and A.T*A == s.eye(3), 'signed three-cycle moment hostile')
need((s.eye(3)-A).det() == 2, 'pointwise anti-label solution only zero')
need((s.eye(3)-A)*(s.eye(3)-A).T != s.zeros(3), 'nonzero residual covariance')

cases = []
for c, denom, nums in ((F(-7,4), 2, [-3,-1,1,3]), (F(-29,16),4,[-7,-5,-3,-1,1,3,5,7])):
    points = [F(a, denom) for a in nums]
    edges = {z: z*z+c for z in points}
    need(all(y in points for y in edges.values()), 'all denominator-and-escape candidates stay in declared set')
    cycles = set()
    for z in points:
        path, where = [], {}
        while z not in where:
            where[z] = len(path)
            path.append(z)
            z = edges[z]
        cycle = path[where[z]:]
        rotations = [tuple(cycle[i:]+cycle[:i]) for i in range(len(cycle))]
        cycles.add(min(rotations))
    need(all(len(cycle) != 6 for cycle in cycles), 'complete fixed-parameter rational preperiodic census excludes six')
    cycle_data = []
    for cycle in sorted(cycles):
        multiplier = F(1)
        for z in cycle:
            multiplier *= 2*z
        cycle_data.append({'points':[str(z) for z in cycle], 'multiplier':str(multiplier)})
    # Independent polynomial gcd: all rational candidates fixed by f^6 have period 2 or 3 here.
    f = x*x+s.Rational(c.numerator,c.denominator)
    f2 = s.expand(f.subs(x,f))
    f3 = s.expand(f.subs(x,f2))
    f6 = s.expand(f3.subs(x,f3))
    numer = s.Poly((f6-x)*(f-x),x)
    denom_poly = s.Poly((f3-x)*(f2-x),x)
    phi6, rem = s.div(numer,denom_poly)
    need(rem.is_zero and phi6.degree() == 54, 'exact sixth dynatomic polynomial')
    for z in points:
        zz = s.Rational(z.numerator,z.denominator)
        need(phi6.eval(zz) != 0, 'no candidate creates a false formal six-period')
    cases.append({'parameter':str(c), 'common_denominator':denom,
        'all_rational_preperiodic_points':[str(z) for z in points],
        'edges':[[str(z),str(y)] for z,y in edges.items()], 'cycles':cycle_data,
        'dynatomic_six_degree':phi6.degree()})

# General fixed-point formal-six factor: the geometric sums cancel to Phi_6(lambda).
lam = s.symbols('lambda')
need(s.cancel((lam**6-1)*(lam-1)/((lam**3-1)*(lam**2-1))) == lam*lam-lam+1,
     'fixed-point formal-six multiplier factor')
need(s.discriminant(lam*lam-lam+1,lam) == -3, 'no rational primitive-six multiplier')
need(s.discriminant(lam*lam+lam+1,lam) == -3, 'no rational primitive-three multiplier')
# At c=-7/4 the algebraic parabolic period-three polynomial is squared.
f = x*x-s.Rational(7,4)
f3 = s.expand(f.subs(x,f).subs(x,f))
P = x**3+x*x/2-s.Rational(9,4)*x-s.Rational(1,8)
need(s.expand(f3-x-(f-x)*P*P) == 0, 'parabolic third dynatomic double factor')
need(s.Poly(P,x).is_irreducible, 'parabolic three-cycle is not rational')

out = {'status':'PROVED algebraic identities; fixed-parameter rational completeness uses the valuation/escape proof in the note',
       'order_six_lift':B.tolist(), 'Pell_conjugate':Pell.tolist(),
       'changing_product':[[str(v) for v in row] for row in product.tolist()], 'quadratic_cases':cases,
       'rational_neutral_cycle_lemma':'Only fixed points x=1/2,c=1/4 and x=-1/2,c=-3/4 have multiplier a rational root of unity.',
       'general_rational_six_cycle_nonexistence':'NOT PROVED; cited Stoll theorem is conditional.'}
target = Path(__file__).resolve().parents[2]/'05-knowledge/results/quadratic_cycle_phase_20261015b.json'
target.write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps(out,indent=2,default=int))
