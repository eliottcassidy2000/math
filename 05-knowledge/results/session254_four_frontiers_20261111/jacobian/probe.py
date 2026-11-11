#!/usr/bin/env python3
"""Exact row16 response certificate; not an inherited debt reconstruction."""
import sympy as s
from math import gcd
from functools import reduce
x=s.Symbol('x')
def depth_matrix(d,n):
    coords=[(i,j) for i in range(n+1) for j in range(i+d+1)]
    ix={a:i for i,a in enumerate(coords)}
    cols=[(a,b,c,e) for b in range(d+1) for a in range(d-b+1) for e in range(n//2+1) for c in range(n+1) if b+c+2*e<=n]
    M=s.zeros(len(coords),len(cols))
    for k,(a,b,c,e) in enumerate(cols):
        for j in range(c+e+1):
            v=b+c+2*e+j
            if v<=n:M[ix[v,a+2*b+e+2*j],k]+=s.binomial(c+e,j)
    return coords,M
def L(m,f):return s.expand(-x*f+(x*x+6)*s.diff(f,x)/(2*m))
def coeffs(f,n):return s.Matrix([s.expand(f).coeff(x,i) for i in range(n+1)])
def moment(m,j):
    return s.S(0) if j%2 else s.prod(s.Rational(6*(2*i-1),2*m-2*i+1) for i in range(1,j//2+1))
def main():
    selectors=[]
    for d,derivative in [(2,x/2),(3,-s.Rational(3,4)-3*x*x/8)]:
        coords,M=depth_matrix(d,15);null=M.T.nullspace()
        N=s.Matrix.vstack(*(v.T for v in null))
        V=s.Matrix.hstack(*(s.Matrix([s.expand(derivative*x**j).coeff(x,r) if n==15 else 0 for n,r in coords]) for j in range(16)))
        selectors.append(N*V)
        print('DEPTH',d,'shape',M.shape,'rank',M.rows-len(null),'cokernel',len(null))
    D=s.Matrix.vstack(*selectors)
    assert D[:,:10]==s.zeros(D.rows,10) and D[:,10:].rank()==6
    rp=list(D[:,10:].T.rref()[1])
    print('TERMINAL15_KERNEL=span(1,x,...,x^9); selector shape',D.shape,'rank',D.rank())
    print('DEPTH_NONZERO_MINOR rows',rp,'cols10..15 det',D.extract(rp,list(range(10,16))).det())
    selectors16=[]
    for d,derivative in [(2,x/2),(3,-s.Rational(3,4)-3*x*x/8)]:
        coords,M=depth_matrix(d,16);null=M.T.nullspace()
        N=s.Matrix.vstack(*(v.T for v in null))
        V=s.Matrix.hstack(*(s.Matrix([s.expand(derivative*x**j).coeff(x,r) if n==16 else 0 for n,r in coords]) for j in range(17)))
        selectors16.append(N*V)
        print('DEPTH16',d,'shape',M.shape,'rank',M.rows-len(null),'cokernel',len(null))
    D16=s.Matrix.vstack(*selectors16)
    assert D16[:,:11]==s.zeros(D16.rows,11) and D16[:,11:].rank()==6
    rp16=list(D16[:,11:].T.rref()[1])
    print('TERMINAL16_KERNEL=span(1,x,...,x^10); selector shape',D16.shape,'rank',D16.rank())
    print('DEPTH16_NONZERO_MINOR rows',rp16,'cols11..16 det',D16.extract(rp16,list(range(11,17))).det())
    T=s.Matrix.hstack(*(coeffs(L(16,x**j),16) for j in range(10)))
    S=s.eye(17)[:,:9]
    Q=s.Matrix([[moment(16,j) if j<=10 else 0 for j in range(17)]]).col_join(s.eye(17)[11:,:])
    assert T.rank()==10 and Q.rank()==7 and Q*T==s.zeros(7,10)
    assert (Q*S).rank()==1 and T.row_join(S).rank()==11
    den=s.ilcm(*(s.denom(moment(16,j)) for j in range(9)))
    obs=[int(den*moment(16,j)) for j in range(9)];g=reduce(gcd,obs)
    print('ROW16_QUOTIENT_DIM=7; RESPONSE_RANK=1; JOINT_RANK=11')
    print('ROW16_RESPONSE_PRIMITIVE=',[v//g for v in obs])
    print('ROW16_MOMENTS_TO_10=',[moment(16,j) for j in range(11)])
    print('WEIGHT24_PAYER=y^8; moment=',moment(16,8))
    assert (Q*S)[:,1]==s.zeros(7,1) and (Q*S)[:,8]!=s.zeros(7,1)
    for k in range(11,17):assert T.row_join(S).row_join(coeffs(x**k,16)).rank()==12
    print('HOSTILES=x^11,...,x^16 are each unpayable; x^7 source invisible')
    target=sum((j+1)*x**j for j in range(11));debt=sum((j+1)*moment(16,j) for j in range(11))
    paid=target-debt/moment(16,8)*x**8
    theta=T.gauss_jordan_solve(coeffs(paid,16))[0]
    assert s.expand(sum(theta[j]*L(16,x**j) for j in range(10))-paid)==0
    print('POSITIVE_CONTROL=coefficients 1..11, unique tangent after y8 payment; PASS')
    for m in range(2,21):
        for d in range(m):
            T=s.Matrix.hstack(*(coeffs(L(m,x**j),d+1) for j in range(d+1)))
            ell=s.Matrix([[moment(m,j) for j in range(d+2)]])
            assert ell*T==s.zeros(1,d+1) and T.rank()==d+1
    print('GENERAL_QUOTIENT_CONTROLS=m2..20,d0..m-1; PASS')
if __name__=='__main__':main()
