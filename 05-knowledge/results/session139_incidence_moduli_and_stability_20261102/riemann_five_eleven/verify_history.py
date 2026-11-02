#!/usr/bin/env python3
"""Exact release cross-checks; SymPy is used only for polynomial arithmetic.
No claim of exhaustive literature review or external novelty is made.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations, permutations
import json
import sympy as S


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mutate',choices=['tetranacci','paley','incidence']);args=ap.parse_args()
    out={}
    def check(label,condition,value):
        if not condition: raise RuntimeError(label)
        out[label]=str(value)
    x,z=S.symbols('x z')
    M=S.Matrix([[1,2*x,0],[0,0,1],[1,x,0]])
    f=z**3-z**2-x*z-x
    check('transfer_characteristic',S.expand(M.charpoly(z).as_expr()-f)==0,f)
    d=S.discriminant(f,z)
    check('transfer_discriminant',S.expand(d-4*x*(x*x-11*x-1))==0,S.factor(d))
    check('quadratic_discriminant',S.discriminant(x*x-11*x-1,x)==125,125)
    d4=S.discriminant(z**4-z**3-z**2-z-1,z)
    check('tetranacci_correction',d4==(-283 if args.mutate=='tetranacci' else -563),d4)
    d5=S.discriminant(z**5-z**4-x*(z**3+z*z+z+1),z)
    expected=-16*x**3*(4*x**4-20*x**3-14*x*x-553*x-16)
    check('pentanacci_parameter',S.expand(d5-expected)==0,S.factor(d5))
    check('pentanacci_at_one',d5.subs(x,1)==9584,d5.subs(x,1))
    check('bernoulli_denominator_repeat',S.denom(S.bernoulli(6))==S.denom(S.bernoulli(114))==42,42)
    n=11;res={pow(i,2,n) for i in range(1,n)}
    A=[[int(i!=j and (j-i)%n in res) for j in range(n)] for i in range(n)]
    partial=[];v=0
    for t in range(1,n): v+=1 if t in res else -1;partial.append(v)
    natural=sum(A[i][j] for i in range(n) for j in range(i+1,n))
    dp=[0]*(1<<n)
    for mask in range(1,1<<n):
        dp[mask]=max(dp[mask^(1<<j)]+sum(A[i][j] for i in range(n) if i!=j and mask>>i&1) for j in range(n) if mask>>j&1)
    check('paley11_partial_sums',min(partial)>=0,partial)
    check('paley11_local_not_global',natural==33 and dp[-1]==(33 if args.mutate=='paley' else 35),(natural,dp[-1]))
    n=5;w=S.Matrix([1,1,-2,0,0]);Q=S.eye(n)-w*w.T/3;K=S.zeros(n)
    K[0,3]=S.Rational(1,4);K[3,0]=-K[0,3];K[1,4]=S.Rational(1,2);K[4,1]=-K[1,4]
    AA=(S.ones(n)-S.eye(n)+K)/2;BB=Q*AA*Q.T
    def frac(v):return F(int(S.numer(v)),int(S.denom(v)))
    def inv(M):
        a=[[frac(M[i,j]) for j in range(n)] for i in range(n)]
        dp={}
        for j in range(n): dp[1<<j,j]=F(1)
        for mask in range(1,1<<n):
            for j in range(n):
                if mask>>j&1 and mask!=(1<<j):
                    prev=mask^(1<<j);dp[mask,j]=sum(dp.get((prev,i),0)*a[i][j] for i in range(n) if prev>>i&1)
        H=sum(dp[(1<<n)-1,j] for j in range(n))
        tri=lambda t:a[t[0]][t[1]]*a[t[1]][t[2]]*a[t[2]][t[0]]+a[t[0]][t[2]]*a[t[2]][t[1]]*a[t[1]][t[0]]
        tau=sum(tri(t) for t in combinations(range(n),3))
        m23=sum(a[i][j]*a[j][i]*tri([k for k in range(n) if k not in (i,j)]) for i,j in combinations(range(n),2))
        u=frac((S.eye(n)+M).det());v=-frac(M.det())
        return H,tau,m23,H-u-abs(v)
    aa,bb=inv(AA),inv(BB)
    check('rational_pair_conjugacy',Q*Q.T==S.eye(n) and Q*S.ones(n,1)==S.ones(n,1) and AA.charpoly(z)==BB.charpoly(z),'orthogonal, fixes 1, cospectral')
    check('held_karp_hamiltonian_difference',bb[0]-aa[0]==F(1,3456),bb[0]-aa[0])
    check('delta_difference',bb[3]-aa[3]==F(1,3456),bb[3]-aa[3])
    check('incidence_difference',bb[2]-aa[2]==F(1,3456 if args.mutate=='incidence' else 6912),bb[2]-aa[2])
    for label,M,I in [('A',AA,aa),('B',BB,bb)]:
        c5=frac(S.Poly((S.eye(n)-z*M).det(),z).nth(5))
        check('enhanced_recovery_'+label,I[0]==1+2*I[1]+2*I[2]-2*c5,I[0])
    print(json.dumps({'status':'PASS','sympy_version':S.__version__,'checks':out},sort_keys=True,indent=2))

if __name__=='__main__':main()
