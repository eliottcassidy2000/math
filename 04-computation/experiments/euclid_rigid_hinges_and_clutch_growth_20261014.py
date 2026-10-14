#!/usr/bin/env python3
"""Exact finite tests for the Euclid/geometry construction. Requires SymPy.
No floating-point tests, network calls, assertions, or external graph packages.
"""
import json
from itertools import combinations
import sympy as s


def require(test, message):
    if not test:
        raise RuntimeError(message)


def norm(v):
    i, j = v
    return i*i+i*j+j*j


def patch_test(n):
    vertices = [(i,j) for i in range(n+1) for j in range(n+1)]
    edges = [(a,b) for a,b in combinations(vertices,2)
             if norm((a[0]-b[0],a[1]-b[1])) == 1]
    require(len(edges)==3*n*n+2*n, 'patch edge count')
    require(all((a[0]-a[1])%3 != (b[0]-b[1])%3 for a,b in edges), 'patch colouring')
    cross=[]
    for u in vertices:
        for v in vertices:
            if u == (0,0) or v == (0,0):
                continue  # common origin, edges already counted inside a patch
            i,j=u; k,l=v
            if i*l-j*k:
                continue  # nonzero irrational coefficient, an exact rejection
            twice_dot=2*i*k+i*l+j*k+2*j*l
            if 6*n*n*(norm(u)+norm(v)-1)==(6*n*n-1)*twice_dot:
                cross.append((u,v))
    require(cross==[((n,n),(n,n))], 'unexpected hinge seam')
    # Explicit four-colouring: the second copy tip alone receives colour 3.
    colour=lambda side,v: 3 if side==1 and v==(n,n) else (v[0]-v[1])%3
    require(all(colour(side,a)!=colour(side,b) for side in [0,1] for a,b in edges), 'four-colour certificate')
    require(colour(0,(n,n))!=colour(1,(n,n)), 'tip-edge certificate')
    return {'n':n, 'delta':12*n*n-1, 'vertices':2*(n+1)**2-1,
            'edges':2*len(edges)+1, 'cross_edges':len(cross),
            'four_colouring_checked':True, 'three_colour_obstruction':'triangle propagation forces both tips to equal the common origin'}

x,y=s.symbols('x y')

def critical_polynomial(V,A,B):
    J=2*V*s.diff(B,x)-B*s.diff(V,x)
    K=(-A**3*J**3+12*A**2*J**2*V**2-4*A*B**3*J**2*V
       +4*A*B**2*J*V**2+24*A*B*J*V**3-48*A*J*V**4
       -16*A*V**4-8*B**4*J*V**2+8*B**3*J*V**3
       +32*B**2*V**4-96*B*V**5+64*V**6)
    return s.Poly(s.expand(K),x),s.expand(J)


def clutch_test(N):
    V=x; A=2*x; B=1 if N==0 else 1+x**N
    K,J=critical_polynomial(V,A,B)
    order=min(m[0] for m,c in K.terms())
    # A separate resultant calculation, rather than a second use of K's formula.
    L=y*y+B*y
    R1=2*L*(2*y+B)+V*A
    R2=V**3+V**2*y+J*y*L
    result=s.Poly(s.resultant(R1,R2,y),x)
    require(result==s.Poly(V**3*K.as_expr(),x), 'independent resultant disagreement')
    require(order==3, 'boundary multiplicity changed')
    if N:
        require(K.degree()-order==5*N-1, 'interior degree law')
    if N>=2:
        require(K.LC()==-16*N*(2*N-1), 'leading coefficient law')
    return {'N':N,'degree':K.degree(),'boundary_order':order,
            'interior_resultant_multiplicity':K.degree()-order,
            'leading_coefficient':int(K.LC()),'independent_resultant_checked':True}


def tower_test(count=3):
    n=1; old=[]; rows=[]
    for j in range(count):
        D=12*n*n-1
        require(D%3==2 and D%8==3,'local class / nonsquare check')
        require(all(s.gcd(D,e)==1 for e in old),'new divisor is not fresh')
        rows.append({'stage':j+1,'n':n,'delta':D,'pairwise_coprime':True,
                     'global_field_degree':2**(j+2),'local_2adic_field':'Q_2(sqrt(3))'})
        old.append(D); n*=D
    return rows


def main():
    output={'scope':'finite exact checks; general statements are proved in the companion note; hypotheses remain open',
            'rigid_hinges':[patch_test(n) for n in range(1,9)],
            'odd_euclid_tower':tower_test(),
            'clutch_tests':[clutch_test(N) for N in range(9)],
            'hostiles':{'fresh_factor_does_not_force_quadratic_extension':'3+1=2^2; sqrt(4) is rational',
                        'even_n_local_class':{'n':2,'delta':47,'mod8':7},
                        'cover_defect':'unramified-at-escape Kummer base change multiplies escape mass by its degree, so normalized defect is unchanged'}}
    print(json.dumps(output,indent=2,sort_keys=True))

if __name__=='__main__': main()
