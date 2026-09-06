#!/usr/bin/env python3
"""Independent literal integer Smith audit; imports no primary probe code.

Frozen universe: all triples -8<=x0<x1<x2<=16; 300 seeded translated,
nonunit-dilated triples; complete p-distance depth controls at p=2,3,5,7,11
through depth six, with several unit lifts. Mutation controls deliberately
remove the exceptional factors two/three and must disagree.
"""
from itertools import combinations, permutations
from math import gcd, prod
import hashlib
import json
import random
import sympy as s
from sympy.matrices.normalforms import smith_normal_form


def require(test, detail):
    if not test:
        raise RuntimeError(detail)


def literal(nodes):
    matrix = s.Matrix([[x**j if r == 0 else j*x**(j-1) if j else 0
                        for j in range(6)] for x in nodes for r in (0,1)])
    normal = smith_normal_form(matrix, domain=s.ZZ)
    return tuple(abs(int(normal[i,i])) for i in range(6))


def formula(nodes):
    z,x,y = sorted(nodes)
    a,b = x-z,y-z
    g = gcd(a,b)
    product = a*b*(b-a)
    divisors = (1,g*gcd(g,2),gcd(g**4,6*product),
                2*g**4*product*gcd(gcd(3,g),(a+b)//g),product**4)
    require(all(divisors[i]%divisors[i-1] == 0 for i in range(1,5)),
            ('integral quotients',nodes))
    answer = (1,1)+tuple(divisors[i]//divisors[i-1] for i in range(1,5))
    require(all(answer[i]%answer[i-1] == 0 for i in range(1,6)),
            ('Smith divisibility',nodes,answer))
    return answer


def vp(n,p):
    require(n!=0, 'finite valuation required')
    count = 0
    while n%p == 0:
        n//=p
        count+=1
    return count


def main():
    rows = []
    universe = list(combinations(range(-8,17),3))
    rng = random.Random(240300737)
    for _ in range(300):
        base = rng.sample(range(-200,201),3)
        scale = rng.choice((2,3,5,6,12,30,72))
        shift = rng.randrange(-1000,1001)
        universe.append(tuple(shift+scale*x for x in base))
    depth_cases = []
    for p in (2,3,5,7,11):
        for e in range(7):
            for f in range(e,7):
                if p==2 and f==e:
                    continue
                for lift in (1,p+1,2*p+1):
                    nodes = (0,p**e,(2 if f==e else 1)*p**f*lift)
                    depth_cases.append((p,e,f,nodes))
                    universe.append(nodes)
    for nodes in universe:
        actual, expected = literal(nodes),formula(nodes)
        require(actual==expected, ('Smith formula',nodes,actual,expected))
        require(prod(actual)==prod(abs(x-y)**4 for x,y in combinations(nodes,2)),
                ('Vandermonde determinant',nodes))
        rows.append((nodes,actual))
    for p,e,f,nodes in depth_cases:
        got = tuple(vp(x,p) for x in literal(nodes))
        correction = int(p==3 and e==f and e>0)
        divisors = (0,min(2*e,e+int(p==2)),
                    min(4*e,2*e+f+int(p in (2,3))),
                    6*e+f+int(p==2)+correction,8*e+4*f)
        predicted = (0,0)+tuple(divisors[i]-divisors[i-1] for i in range(1,5))
        require(got==predicted, ('metric formula',p,e,f,nodes,got,predicted))
        require(got[-1]==2*e+3*f-int(p==2)-correction,'precision law')
    for nodes in ((0,1,2),(0,2,4),(0,3,6),(0,8,16),(-7,3,41)):
        for order in permutations(nodes):
            require(literal(order)==formula(nodes),('node permutation',order))
    require(literal((0,1,2)) != (1,1,1,1,2,8), 'factor-two mutation detected')
    require(literal((0,3,6)) != (1,1,3,27,108,972), 'factor-three mutation detected')
    print('PASS independent literal integer Smith audit')
    print('UNIVERSE triples',len(universe),'depth_controls',len(depth_cases),
          'permutations',30,'mutation_controls',2)
    print('CONTROLS', [(x,literal(x)) for x in ((0,1,2),(0,2,4),(0,3,6),(0,8,16))])
    print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(rows).encode()).hexdigest())


if __name__=='__main__':
    main()
