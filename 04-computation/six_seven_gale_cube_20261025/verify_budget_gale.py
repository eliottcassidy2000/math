#!/usr/bin/env python3
"""Exact controls for BG1/BG2; no assertions disabled by optimized Python."""
import itertools
import json
import math
import sympy as sp


def require(test, message):
    if not test:
        raise RuntimeError(message)


def hilbert_coefficient(n, c, r, degree):
    if degree < 0:
        return 0
    return sum((-1)**j * math.comb(c, j) * math.comb(n + degree-r*j-1, n-1)
               for j in range(min(c, degree//r)+1))


def polynomial_coefficient(n, c, r, degree):
    if degree < 0:
        return 0
    # Independent monic-basis Hilbert series:
    # (1+...+t^(r-1))^c /(1-t)^(n-c).
    coefficients = [1]
    for _ in range(c):
        out = [0] * (len(coefficients)+r-1)
        for i, value in enumerate(coefficients):
            for j in range(r):
                out[i+j] += value
        coefficients = out
    return sum(value * math.comb(n-c+degree-i-1, n-c-1)
               for i, value in enumerate(coefficients) if i <= degree)


def main():
    gates = 0
    cases = []
    for d in range(1, 6):
        for r in range(2, 6):
            for n in range(d+2, 2*d+6):
                c = n-d-1
                a = r*c-n
                require(a == (r-1)*n-r*(d+1), 'canonical coefficient')
                p = hilbert_coefficient(n,c,r,a)
                require(p == polynomial_coefficient(n,c,r,a), 'two Hilbert paths')
                gates += 2
                if r == 2 and n in (2*d+2, 2*d+3):
                    require(p == (1 if n == 2*d+2 else n), 'threshold sections')
                    cases.append(dict(d=d,n=n,r=r,canonical_coefficient=a,
                                      canonical_sections=p,degree=r**c))
                    gates += 1
    require(4*(4-2-1)-4 == 0, 'quartic K3 hostile')
    u,v = sp.symbols('u v')
    A = sp.Matrix([[1,1,1,1,1,1],[0,1,-1,0,0,u],[0,0,0,1,-1,v]])
    C = sp.Matrix([[-2,1,1,0,0,0],[-2,0,0,1,1,0],[u+v-1,-u,0,-v,0,1]])
    M = sp.Matrix([[u,v,2],[u,0,0],[0,v,0]])
    permutation = (5,2,1,4,3,0)
    diagonal = sp.diag(-2,-u,u,-v,v,2)
    require(A*C.T == sp.zeros(3), 'Gale orthogonality')
    require(sp.simplify(M*C-A[:,permutation]*diagonal) == sp.zeros(3,6), 'Gale return')
    require(M.det() == 2*u*v, 'Gale determinant')
    require(diagonal.det() == -4*u**2*v**2, 'double-plane minus-one twist')
    require(M.subs(u,0).det() == 0 and M.subs(v,0).det() == 0, 'boundary hostile')
    gates += 6
    a,c = A.subs({u:1,v:1}), C.subs({u:1,v:1})
    def triples(matrix):
        return {frozenset(t) for t in itertools.combinations(range(6),3)
                if matrix[:,t].det() == 0}
    ta,tc = triples(a),triples(c)
    variables = sp.symbols('m:9')
    T = sp.Matrix(3,3,variables)
    compatible,matches = 0,[]
    for p in itertools.permutations(range(6)):
        if {frozenset(p[i] for i in t) for t in tc} != ta:
            continue
        compatible += 1
        equations = [e for i in range(6) for e in (T*c[:,i]).cross(a[:,p[i]])]
        kernel = sp.linear_eq_to_matrix(equations,variables)[0].nullspace()
        # For this finite universe every compatible solution space has dim<=1.
        require(len(kernel) <= 1, 'unexpected higher-dimensional projective maps')
        if kernel and sp.Matrix(3,3,kernel[0]).det() != 0:
            matches.append(list(p))
        gates += 1
    require(compatible == 8 and len(matches) == 2, 'complete projective census')
    gates += 1
    print(json.dumps(dict(status='PASS',gates=gates,
                         universe='d=1..5,r=2..5,d+2<=n<=2d+5; symbolic u,v; all720 permutations at u=v=1',
                         square_root_thresholds=cases,matroid_compatible=compatible,
                         projective_matches=matches),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
