#!/usr/bin/env python3
"""Complete earlier-memory probe, valuation>=12 through row15.

FINITE-EXACT exploration until the result and its global boundary consumer
are independently audited. No repository mathematical imports.
"""
import sympy as S
from sympy.polys.matrices import DomainMatrix
from math import comb

def main():
    x=S.Symbol('x')
    A=[1+x*x/4,S.Rational(4,3)+2*x*x,
       -S.Rational(32,9)-S.Rational(4,5)*x*x,
       S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4,
       S.Rational(224,9)*x**4-S.Rational(256,75)*x*x-S.Rational(37376,405)]
    C=[-3*x/4-x**3/8,-4*x-S.Rational(3,2)*x**3,
       S.Rational(88,15)*x-S.Rational(12,5)*x**3,
       -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5,
       -S.Rational(184,15)*x**5-S.Rational(2048,25)*x**3+S.Rational(98944,675)*x]
    da,dc,unknown={},{},[]
    for n in range(11,16):
        va=S.symbols(f'a{n}_0:{n+2}');vc=S.symbols(f'c{n}_0:{n+3}')
        da[n]=sum(v*x**j for j,v in enumerate(va))
        dc[n]=sum(v*x**j for j,v in enumerate(vc))
        unknown.extend(va+vc)
    low=[(v-2*b,b) for v in range(12,16) for b in range(v//2+1)
         if 2*v-b<=23]
    monomials=low+[(3,6)]
    rs=S.symbols(f'r0:{len(monomials)}')
    source={n:S.Integer(0) for n in range(11,16)}
    for r,(a,b) in zip(rs,monomials):
        for j in range(a+b+1):
            n=a+2*b+j
            if n in source: source[n]+=r*comb(a+b,j)*x**(b+2*j)
    ga=[S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))
                +(S.Rational(3,4) if k==0 else 0)) for k in range(5)]
    eq=[]
    for n in range(11,16):
        pe=-source[n];je=0
        for k in range(5):
            j=n-k
            if j not in da:continue
            pe+=ga[k]*da[j]+2*C[k]*dc[j]
            je+=(j*S.diff(A[k],x)*dc[j]+k*S.diff(da[j],x)*C[k]
                 -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        eq.extend(S.Poly(S.expand(pe),x).all_coeffs())
        eq.extend(S.Poly(S.expand(je),x).all_coeffs())
    for d,v in [(2,da),(3,dc)]:
        for ell in range(1,31):
            start=(ell+1)//2;stop=min(15,ell+d);rho=(ell+2)//3
            if start>stop:continue
            for q in range(max(0,ell+d-stop),rho):
                eq.append(sum((-1)**(n-start)*comb(stop+q-n,q)
                              *S.expand(v[n]).coeff(x,2*n-ell)
                              for n in range(start,stop+1) if n in v))
    matrix,rhs=S.linear_eq_to_matrix(eq,unknown+list(rs))
    if rhs!=S.zeros(len(eq),1):raise RuntimeError('inhomogeneous probe')
    dm,piv=DomainMatrix.from_Matrix(matrix).convert_to(S.QQ).rref()
    red=dm.to_Matrix();nu=len(unknown)
    raw=[list(red[i,nu:]) for i in range(red.rows)
         if all(red[i,j]==0 for j in range(nu)) and any(red[i,nu:])]
    dual=S.Matrix(raw)
    print('VALUATION12 RAW',len(eq),'UNKNOWN',nu,'SOURCE',len(rs))
    print('MONOMIALS',[(a,b,a+2*b,2*a+3*b) for a,b in monomials])
    print('TANGENT_RANK',sum(p<nu for p in piv),'SOURCE_RANK',len(raw))
    print('RELATIONS',dual.tolist())
    for w in range(18,24):
        idx=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=w]
        lower=dual[:,idx]
        before=lower.rank();after=lower.row_join(dual[:,-1]).rank()
        print('WEIGHT',w,'rank',before,'with_high',after)
        if before==after:
            ans=S.linsolve((lower,dual[:,-1]))
            print('REPLACEMENTS r_high=-1',[(monomials[i]) for i in idx],ans)
    print('FINITE-EXACT; lower valuations and later rows untested')

if __name__=='__main__':main()
