#!/usr/bin/env python3
"""Finite-exact controls for the polynomial-bit signed-contact construction.

The proof is in the companion note. No floating-point values or assert
statements are used. Universal coefficient bounds are checked as integers.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd
import json
import sympy as sp

COUNTS = {'symbolic': 0, 'integer': 0, 'rational': 0}


def check(ok, why, kind='rational'):
    COUNTS[kind] += 1
    if not ok:
        raise RuntimeError(why)


def bits(v):
    v = F(v)
    return max(abs(v.numerator).bit_length(), v.denominator.bit_length())


def point_bits(P):
    return 0 if P is None else max(map(bits, P))


def neg(P):
    return None if P is None else (P[0], -P[1])


def add(P, Q, a2, a4):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    xx, yy = Q
    if x == xx:
        if y == -yy:
            return None
        if y != yy:
            raise RuntimeError('same x with incompatible y')
        slope = (3*x*x + 2*a2*x + a4)/(2*y)
    else:
        slope = (yy-y)/(xx-x)
    out_x = slope*slope-a2-x-xx
    return out_x, -y+slope*(x-out_x)


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def on_curve(P, k):
    return P is None or P[1]**2 == (P[0]-k*k)*(P[0]**2-1)


def forward(u, z, k):
    X = (k*z+k*k+2*u+k*k*u*u)/(2*u*u)
    return X, (X-k*k)*(X*u-1)/k


def inverse(P, k):
    if P is None or P[0] == 0 or P[0] == k*k:
        return None
    X, Y = P
    u = (k*Y+X-k*k)/(X*(X-k*k))
    z = (2*u*u*X-k*k-2*u-k*k*u*u)/k
    return u, z


def primitive(t):
    n, m = t.numerator, t.denominator
    d = 2 if n % 2 and m % 2 else 1
    a, b, c = (n*n-m*m)//d, 2*n*m//d, (n*n+m*m)//d
    if a % 2 == 0:
        a, b = b, a
    check(a>0 and b>0 and a*a+b*b==c*c, 'positive Pythagorean output')
    check(gcd(a,b)==gcd(a,c)==gcd(b,c)==1, 'primitive output')
    check(F(2*a*b,c*c)==spread(t), 'spread retained')
    return a,b,c,d


def normalized_division_polynomials(x, aa, bb, expand=lambda p:p):
    G = x**3+aa*x+bb
    base = {
        1: 1, 2: 1,
        3: 3*x**4+6*aa*x*x+12*bb*x-aa*aa,
        4: 2*(x**6+5*aa*x**4+20*bb*x**3-5*aa*aa*x*x
              -4*aa*bb*x-aa**3-8*bb*bb),
    }
    @lru_cache(None)
    def fn(n):
        if n in base:
            return expand(base[n])
        m = n//2
        if n%2:
            first = fn(m+2)*fn(m)**3
            second = fn(m-1)*fn(m+1)**3
            return expand(16*G*G*first-second if m%2==0
                          else first-16*G*G*second)
        return expand(fn(m)*(fn(m+2)*fn(m-1)**2
                             -fn(m-2)*fn(m+1)**2))
    return fn


def symbolic_checks():
    x,a2,a4,a6 = sp.symbols('x a2 a4 a6')
    xi=9*x+3*a2
    aa=81*a4-27*a2*a2
    bb=54*a2**3-243*a2*a4+729*a6
    check(sp.expand(xi**3+aa*xi+bb-729*(x**3+a2*x*x+a4*x+a6))==0,
          'short model scaling', 'symbolic')
    x,aa,bb = sp.symbols('x alpha beta')
    fn=normalized_division_polynomials(x,aa,bb,sp.expand)
    rows=[]
    for n in range(1,13):
        p=sp.Poly(fn(n),x,aa,bb)
        expected=n*n-(1 if n%2 else 4)
        check(all(2*i+4*j+6*k==expected for i,j,k in p.monoms()),
              'weighted degree n=%d'%n, 'symbolic')
        norm=sum(abs(int(c)) for c in p.coeffs())
        check(norm <= (n+1)**(4*n*n), 'universal norm n=%d'%n, 'integer')
        rows.append({'n':n,'weighted_degree':expected,
                     'monomials':len(p.monoms()),'coefficient_norm_bits':norm.bit_length()})
    return rows


def test_anchor(t, epsilon, scalar_checks=16):
    a,b,C,d0=primitive(t)
    f0=spread(t)
    check(t>1 and 0<f0<F(1,2), 'admissible input')
    # Use the actual input triple's inverse parameter, as in the theorem.
    t0=F(C+a,b)
    check(spread(t0)==f0, 'input parameter reconstruction')
    A=abs(a+epsilon*b)
    k=F(A,C)
    check(k*k==1+epsilon*f0 and gcd(A,C)==1 and A%2==C%2==1,
          'controller parameter')
    P=(F(0),k)
    Q=forward(epsilon*t0,t0*t0+1,k)
    check(on_curve(P,k) and on_curve(Q,k), 'anchor cubic points')
    check(inverse(Q,k)==(epsilon*t0,t0*t0+1), 'anchor inverse chart')
    a2=-A*A; a4=-C**4; a6=A*A*C**4
    aa=81*a4-27*a2*a2
    bb=54*a2**3-243*a2*a4+729*a6
    xi=3*a2; eta=27*A*C*C
    Rsize=6*C
    check(abs(xi)<=Rsize**2 and abs(2*eta)<=Rsize**3
          and abs(aa)<=Rsize**4 and abs(bb)<=Rsize**6,
          'weighted input size', 'integer')
    check(eta*eta==xi**3+aa*xi+bb, 'short image of P', 'integer')
    fn=normalized_division_polynomials(xi,aa,bb)
    def psi(n):
        return 0 if n==0 else fn(n)*(1 if n%2 else 2*eta)
    Ps=(F(xi),F(eta))
    S=None
    scalar_max_bits=0
    for n in range(1,scalar_checks+1):
        S=add(S,Ps,0,aa)
        check(S is not None, 'nonzero multiple')
        check(psi(n)!=0, 'nonzero division denominator', 'integer')
        check(S[0] == F(xi)-F(psi(n-1)*psi(n+1),psi(n)**2),
              'division polynomial vs exact additions')
        check(abs(psi(n)) <= Rsize**(n*n-1)*(n+1)**(4*n*n),
              'explicit division evaluation bound', 'integer')
        check(S[1]**2==S[0]**3+aa*S[0]+bb, 'multiple short curve')
        scalar_max_bits=max(scalar_max_bits,point_bits(S))
    # Every exceptional group-law and skipped inverse-chart case.
    check(add(P,None,-k*k,-1)==P and add(None,P,-k*k,-1)==P,
          'identity cases')
    check(add(P,neg(P),-k*k,-1) is None, 'opposite-point case')
    torsion=(F(-1),F(0))
    check(add(torsion,torsion,-k*k,-1) is None, 'vertical tangent case')
    for excluded in [None,P,neg(P),(k*k,F(0))]:
        check(on_curve(excluded,k) and inverse(excluded,k) is None,
              'excluded inverse chart exact test')
    D=2+f0 if epsilon==1 else F(2)
    J=0
    while (1<<J)*f0 < D:
        J+=1
    N=100*(4+J)//7
    check((1<<J)*f0>=D and (J==0 or (1<<(J-1))*f0<D), 'least J')
    check(J<=2*(C-1).bit_length()+1, 'logarithmic stopping budget')
    B=add(P,P,-k*k,-1)
    R=None
    max_bits=max(point_bits(Q),point_bits(B))
    for i in range(1,N+1):
        R=add(R,B,-k*k,-1)
        for sign in (1,-1):
            S=add(Q,R if sign==1 else neg(R),-k*k,-1)
            check(on_curve(S,k), 'search candidate on curve')
            max_bits=max(max_bits,point_bits(R),point_bits(S))
            inv=inverse(S,k)
            if inv is None:
                continue
            u,z=inv
            check(z*z==k*k*(u*u+1)**2-4*u*(u*u-1), 'inverse quartic')
            t1=epsilon*u
            if t1<=1 or z<=0:
                continue
            f1=spread(t1)
            if not (0<f1<F(1,2)) or f1==f0:
                continue
            aa1,bb1,CC,d=primitive(t1)
            root=z/(u*u+1)
            check(root>0 and root*root==1+epsilon*f0-epsilon*f1,
                  'exact contact certificate')
            check((aa1,bb1,CC)!=(a,b,C), 'distinct primitive output')
            check(CC*CC>4*aa1*bb1, 'admissible output')
            q=2*sign*i
            check(0<abs(q)<=2*N and q%2==0, 'certified index')
            return {'input_t':str(t), 'input_C':str(C),
                    'input_C_bits':C.bit_length(), 'port':epsilon,
                    'q':q, 'budget':2*N, 'J':J,
                    'neighbor_C_bits':CC.bit_length(),
                    'search_coordinate_max_bits':max_bits,
                    'scalar_control_max_bits':scalar_max_bits,
                    'certificate_bits':bits(root),
                    'normalization_divisor':d}
    raise RuntimeError('certified search budget exhausted')


def main():
    symbolic=symbolic_checks()
    params=[F(6,5),F(9),F(11,9),F(101),F(257,256),F(256),
            F(65537,65536),F(65536)]
    rows=[test_anchor(t,e) for t in params for e in (1,-1)]
    print(json.dumps({'status':'FINITE-EXACT','checks':COUNTS,
                      'symbolic_polynomial_controls':symbolic,
                      'witnesses':rows},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
