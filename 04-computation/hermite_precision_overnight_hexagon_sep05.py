#!/usr/bin/env python3
"""Independent full-matrix audit of the universal two-jet precision formula.

No code is imported from a Smith-profile or minor-envelope experiment.
The formula uses only F'(node), F''(node); the control uses literal integer
Smith reduction, with direct rational inversion on a separate small bank.
"""
from itertools import combinations
from math import gcd, lcm, prod
from hashlib import sha256
import json
import random
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

GATES = 0


def need(test, detail):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(detail)


def matrix(nodes):
    n = 2*len(nodes)
    return s.Matrix([[x**j if r==0 else j*x**(j-1) if j else 0
                      for j in range(n)] for x in nodes for r in (0,1)])


def formula(nodes):
    factors = []
    data = []
    for i,x in enumerate(nodes):
        differences = [x-y for j,y in enumerate(nodes) if j!=i]
        first = prod(differences)
        second = 2*sum(prod(d for j,d in enumerate(differences) if j!=k)
                       for k in range(len(differences)))
        factor = abs(first)**3//gcd(abs(first),second)
        factors.append(factor)
        data.append((first,second,factor))
    return lcm(*factors),data


def main():
    universe = [(x,) for x in range(-10,11)]
    universe += list(combinations(range(-5,11),2))
    universe += list(combinations(range(-3,9),3))
    for n,height in ((4,12),(5,10),(6,10)):
        universe += [(0,)+row for row in combinations(range(1,height+1),n-1)]
    rng = random.Random(240300737)
    for n,count in ((3,40),(4,40),(5,40),(6,40),(7,50),(8,30)):
        for _ in range(count):
            base = rng.sample(range(-100,101),n)
            scale = rng.choice((1,2,3,8,9,12,25,72))
            shift = rng.randrange(-1000,1001)
            universe.append(tuple(shift+scale*x for x in base))
    rows = []
    for nodes in universe:
        actual = smith_normal_form(matrix(nodes),domain=s.ZZ)
        factors = tuple(abs(int(actual[i,i])) for i in range(2*len(nodes)))
        expected,data = formula(nodes)
        need(factors[-1] == expected, ('largest integer Smith factor',nodes,factors,expected,data))
        need(prod(factors) == prod(abs(x-y)**4 for x,y in combinations(nodes,2)),
             ('full determinant control',nodes))
        rows.append((nodes,expected))
    inverse_cases = universe[:21] + universe[21:41] + universe[-30:]
    for nodes in inverse_cases:
        inverse = matrix(nodes).inv()
        denominator = lcm(*(int(s.denom(entry)) for entry in inverse))
        need(denominator == formula(nodes)[0], ('literal inverse denominator',nodes))
    # Directly verify the cardinal polynomial proof, not just its denominator.
    X = s.symbols('X')
    for nodes in ((-3,),(-2,7),(0,1,2),(0,8,16,40),(0,8,24,32),(0,9,18,27)):
        for i,x in enumerate(nodes):
            Q = s.sympify(s.prod(X-y for j,y in enumerate(nodes) if i!=j))
            first = Q.subs(X,x)
            tau = s.diff(Q,X).subs(X,x)/first
            value = s.Poly((1-2*tau*(X-x))*Q**2/first**2,X)
            derivative = s.Poly((X-x)*Q**2/first**2,X)
            for j,y in enumerate(nodes):
                need(value.eval(y)==int(i==j) and value.diff().eval(y)==0,
                     ('value cardinal polynomial',nodes,i,j))
                need(derivative.eval(y)==0 and derivative.diff().eval(y)==int(i==j),
                     ('derivative cardinal polynomial',nodes,i,j))
    # These isometric four-node controls have the same largest dyadic loss,
    # but their complete Smith lists differ: precision is not full torsion.
    for base in ((0,1,2,5),(0,1,3,4)):
        print('CONTROL',tuple(8*x for x in base),formula(tuple(8*x for x in base)))
    print('UNIVERSE',len(universe),'node_sizes=1..8',
          'literal_inverse_controls',len(inverse_cases),'EXPLICIT_GATES',GATES)
    print('SEMANTIC_SHA256',sha256(json.dumps(rows).encode()).hexdigest())
    print('RESULT=PASS universal two-jet largest-factor formula; no full Smith extrapolation')


if __name__ == '__main__':
    main()
