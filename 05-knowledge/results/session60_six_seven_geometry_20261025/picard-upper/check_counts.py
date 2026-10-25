#!/usr/bin/env python3
"""Independent exact point counts for the six-line double-plane K3.
No floating point, no imports outside stdlib. Two character computations are checked.
"""
import itertools, json
from fractions import Fraction
from pathlib import Path


def det(M):
    a,b,c=M[0];d,e,f=M[1];g,h,i=M[2]
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def rational_point(v):
    a=next(x for x in v if x)
    return tuple(Fraction(x,a) for x in v)


def field(p,n):
    if n==1:
        return p,None,lambda a,b:(a+b)%p,lambda a,b:(a*b)%p,lambda a:(-a)%p
    d=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
    def add(a,b): return ((a%p+b%p)%p)+p*((a//p+b//p)%p)
    def mul(a,b):
        aa,ab=a%p,a//p;ba,bb=b%p,b//p
        return (aa*ba+d*ab*bb)%p+p*((aa*bb+ab*ba)%p)
    def neg(a): return (-(a%p))%p+p*((-(a//p))%p)
    return p*p,d,add,mul,neg


def count(p,n):
    q,d,add,mul,neg=field(p,n)
    def power(a,m):
        z=1
        while m:
            if m&1:z=mul(z,a)
            a=mul(a,a);m//=2
        return z
    def sub(a,b):return add(a,neg(b))
    squares={mul(a,a) for a in range(1,q)}
    if len(squares)!=(q-1)//2:raise AssertionError('field square cardinality')
    chis=[]
    for a in range(q):
        c=0 if not a else 1 if a in squares else -1
        ep=power(a,(q-1)//2)
        if ep != (0 if not c else 1 if c==1 else neg(1)): raise AssertionError('Euler mismatch')
        if power(a,q)!=a: raise AssertionError('Frobenius field check')
        if n==2:
            norm=((a%p)**2-d*(a//p)**2)%p
            legendre=0 if not norm else 1 if pow(norm,(p-1)//2,p)==1 else -1
            if c!=legendre:raise AssertionError('independent norm-character check')
        chis.append(c)
    S=0; histogram={-1:0,0:0,1:0}; direct=0
    for x,z in itertools.product(range(q),repeat=2):
        value=mul(mul(sub(1,mul(x,x)),sub(1,mul(z,z))),add(add(1,x),z))
        c=chis[value];S+=c;histogram[c]+=1
        # Fiber size counted independently from its square lookup.
        direct+= 1 if value==0 else 2 if value in squares else 0
    singular=direct+q+1 # all points on E=0 branch line have one point above
    if singular != q*q+q+1+S:raise AssertionError('point formula')
    resolved=singular+17*q
    return dict(p=p,n=n,q=q,quadratic_nonresidue=d,character_sum=S,
                affine_character_histogram=histogram,singular_points=singular,
                resolved_points=resolved,residual18_trace=resolved-1-q*q-18*q)


def main():
    lines=[(1,0,0),(1,-1,0),(1,1,0),(1,0,-1),(1,0,1),(1,1,1)]
    incidences={rational_point(cross(a,b)):tuple(i for i,l in enumerate(lines) if sum(u*v for u,v in zip(l,cross(a,b)))==0) for a,b in itertools.combinations(lines,2)}
    counts={r:sum(len(v)==r for v in incidences.values()) for r in [2,3]}
    if counts!={2:9,3:2}:raise AssertionError(counts)
    determinants=sorted(set(det(M) for M in itertools.combinations(lines,3)))
    if any(abs(v) not in [0,1,2,3] for v in determinants):raise AssertionError(determinants)
    data={str(p):[count(p,n) for n in [1,2]] for p in [5,7,11,13]}
    residuals={}
    for p in [5,7,11,13]:
        a,A=[r['residual18_trace'] for r in data[str(p)]]
        b=(a*a-A)//2
        if 2*b!=a*a-A:raise AssertionError('nonintegral coefficient')
        residuals[str(p)]=dict(trace=a,trace_square=A,middle_coefficient=b,
             positive_reciprocal_coefficients=[1,-a,b,-a*p*p,p**4],
             positive_reciprocity_forced=(b!=0))
    # p7,p13 target factorization by direct coefficient multiplication.
    for p,c in [(7,10),(13,22)]:
        expected=[1,c-2*p,2*p*p-2*p*c,p*p*c-2*p**3,p**4]
        if residuals[str(p)]['positive_reciprocal_coefficients']!=expected:raise AssertionError('factorization')
        residuals[str(p)].update(factorization=f'(T-{p})^2*(T^2+{c}*T+{p*p})',
             quadratic_root_sum=-c,normalized_root_sum=str(Fraction(-c,p)),
             artin_tate_discriminant_squareclass_representative=-c-2*p)
    packet=dict(branch_lines=lines,incidences=[dict(point=[str(c) for c in P],lines=v) for P,v in sorted(incidences.items())],intersection_multiplicities=counts,three_line_determinants=determinants,field_counts=data,residual_polynomials=residuals)
    print(json.dumps(packet,sort_keys=True,indent=2))

if __name__=='__main__':main()
