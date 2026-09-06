#!/usr/bin/env python3
"""Independent symbolic and literal referee; no primary-checker import.

Build entries by differentiating monomials, use SymPy's independent
polynomial determinant algorithm on all 886 minors, reconstruct intrinsic
closest-pair unit coordinates, and compute literal 9x9 integer Smith forms.
"""
from fractions import Fraction
from itertools import combinations
from math import factorial
import hashlib
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

GATES = 0
DIGEST = hashlib.sha256()


def require(ok, data):
    global GATES
    if not ok:
        raise RuntimeError(repr(data))
    GATES += 1
    DIGEST.update(repr(data).encode()+b'\n')


def v2(a):
    a = abs(int(a))
    if not a:
        raise ValueError("zero valuation")
    return (a & -a).bit_length()-1


def main():
    a, x = sp.symbols('a x')
    degrees = tuple(range(3,9))
    orders = (0,1,2,0,1,2)
    residual = sp.Matrix([[sp.diff(x**j,x,r).subs(x,node)/factorial(r)
                           for j in degrees]
                          for node in (1,a) for r in range(3)])
    expected = {1:{(1,0,0)}, 2:{(5,0,0),(4,0,1),(3,1,2)},
                3:{(9,0,0),(7,1,2)},4:{(13,1,1),(12,4,5)}}
    count = 0
    for k in range(1,5):
        costs = set()
        for rr in combinations(range(6),k):
            for cc in combinations(range(6),k):
                det = residual.extract(rr,cc).det(method='domain-ge')
                polynomial = sp.Poly(det,a,domain=sp.ZZ)
                require(not polynomial.is_zero, ('nonzero',rr,cc))
                weight = sum(degrees[c] for c in cc)-sum(orders[r] for r in rr)
                for (power,),coefficient in polynomial.terms():
                    costs.add((weight,power,power+v2(coefficient)))
                count += 1
        front = {c for c in costs if not any(d!=c and all(t<=u for t,u in zip(d,c)) for d in costs)}
        require(front==expected[k], ('independent-front',k,sorted(front)))
        print('INDEPENDENT DOMAIN DETERMINANTS rank',k,'FRONT',sorted(front))
    witnesses = [((2,),(0,),3), ((0,2),(0,1),3), ((1,2),(0,1),6),
                 ((2,5),(0,1),18*a*(a-1)), ((0,1,2),(0,1,2),1),
                 ((1,2,5),(0,1,2),30*a*(a-1)*(2*a-1)),
                 ((0,1,2,5),(0,1,2,3),3*a*(a-1)*(5*a*a-5*a+1)),
                 ((1,2,4,5),(0,1,2,3),90*a**4*(a-1)**4)]
    for rr,cc,factor in witnesses:
        require(sp.expand(residual.extract(rr,cc).det(method='domain-ge')-factor)==0,
                ('independent-attainment',rr,cc))
    controls = 0
    for e0 in range(7):
        for d0 in range(1,6):
            for outside,close in ((1,1),(3,1),(5,3),(9,7),(15,5)):
                shift = e0-3*d0
                sign = -1 if (e0+d0+outside)%2 else 1
                nodes = tuple(shift+sign*2**e0*z for z in (0,outside,2**d0*close))
                matrix = sp.Matrix([[sp.diff(x**j,x,r).subs(x,node)/factorial(r)
                                     for j in range(9)] for node in nodes for r in range(3)])
                snf = smith_normal_form(matrix,domain=sp.ZZ)
                actual = tuple(v2(snf[i,i]) for i in range(9))
                depths = {(i,j):v2(nodes[i]-nodes[j]) for i,j in combinations(range(3),2)}
                e = min(depths.values())
                pair = max(depths,key=depths.get)
                d = depths[pair]-e
                outsider = next(i for i in range(3) if i not in pair)
                t = Fraction(2*(nodes[outsider]-nodes[pair[0]]),nodes[pair[1]]-nodes[pair[0]])
                # t is odd only in the shallow d=1 unit-class case.
                if d==1:
                    require(t.numerator%2==t.denominator%2==1, ('intrinsic-unit',nodes))
                    residue = (t.numerator*pow(t.denominator,-1,16))%16
                    gamma = 3 if residue%4==3 else 2 if residue%8==5 else 1 if residue==1 else 0
                    largest = max(7*e+4,8*e+gamma)
                else:
                    largest = 8*e+5*d-1
                D = (0,0,0,0,e,min(5*e,4*e+1,3*e+d+1),
                     min(9*e,7*e+d+1),min(13*e+d,12*e+4*d+1),
                     27*e+9*d-largest,27*e+9*d)
                predicted = tuple(D[j+1]-D[j] for j in range(9))
                require(actual==predicted, ('independent-Smith',nodes,e,d,actual,predicted))
                controls += 1
    require(count==886, ('complete-symbolic-universe',count))
    require(min(13*5+1,12*5+4+1)==65, 'missing-slope-hostile')
    print('SYMBOLIC MINORS',count,'ATTAINING FACTORS',len(witnesses))
    print('ARBITRARY SIGN/TRANSLATION/UNIT FULL9x9 SMITH CONTROLS',controls)
    print('GATES',GATES,'SEMANTIC SHA256',DIGEST.hexdigest())
    print('PASS independent full dyadic three-node three-jet partition audit')


if __name__=='__main__':
    main()
