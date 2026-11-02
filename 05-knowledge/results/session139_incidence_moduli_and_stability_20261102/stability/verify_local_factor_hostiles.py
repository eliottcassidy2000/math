#!/usr/bin/env python3
"""Independent Fraction determinants; exact negative controls for local factors."""
from fractions import Fraction as F
import json

def check(value,msg):
    if not value:
        raise RuntimeError(msg)

def det(A):
    A=[list(map(F,row)) for row in A]
    n=len(A); ans=F(1)
    for j in range(n):
        k=next((k for k in range(j,n) if A[k][j]),None)
        if k is None:
            return F(0)
        if k!=j:
            A[j],A[k]=A[k],A[j];ans=-ans
        p=A[j][j];ans*=p
        for k in range(j+1,n):
            z=A[k][j]/p
            for l in range(j+1,n):
                A[k][l]-=z*A[j][l]
    return ans

def margins(A):
    n=len(A)
    u=det([[A[i][j]+(i==j) for j in range(n)] for i in range(n)])
    da=det(A);v=(-1)**n*da
    return {'u':u,'detA':da,'fplus':u+v,'fminus':u-v,'minimum':u-abs(da)}

def five(x):
    C=[[0,1,1],[-1,0,1],[-1,-1,0]]
    U=[[F(1,2),F(1,2),0],[-F(1,4),F(1,4),F(1,2)]]
    K=[[F(0) for _ in range(5)] for _ in range(5)]
    K[0][1]=x;K[1][0]=-x
    for i in range(2):
        for j in range(3):
            K[i][j+2]=U[i][j];K[j+2][i]=-U[i][j]
    for i in range(3):
        for j in range(3):K[i+2][j+2]=C[i][j]
    return [[F(0) if i==j else (1+K[i][j])/2 for j in range(5)] for i in range(5)]

def ordinal_pair_two_triangles(p):
    sizes=[2,3,3];n=sum(sizes)
    A=[[F(0) for _ in range(n)] for _ in range(n)]
    block=[0,0,1,1,1,2,2,2]
    for i in range(n):
        for j in range(n):
            if block[i]<block[j]:A[i][j]=F(1)
    A[0][1]=p;A[1][0]=1-p
    for s in [2,5]:
        for i,j in [(0,1),(1,2),(2,0)]:A[s+i][s+j]=F(1)
    return A

def valid(A):
    n=len(A)
    return all(A[i][i]==0 for i in range(n)) and all(0<=A[i][j]<=1 and A[i][j]+A[j][i]==1 for i in range(n) for j in range(n) if i!=j)

def stringify(obj):
    if isinstance(obj,F):return str(obj)
    if isinstance(obj,dict):return {k:stringify(v) for k,v in obj.items()}
    return obj

x=F(3,4);p=(1+x)/2;q=p*(1-p)
A=five(x);check(valid(A),'n5 validity')
five_result=margins(A)
lo=margins(five(-1));hi=margins(five(1))
L=(1-p)*lo['fplus']+p*hi['fplus']
check(five_result['u']==F(165,256),'n5 u')
check(five_result['detA']==F(33,128),'n5 detA')
check(five_result['fplus']==F(99,256),'n5 fplus')
check(lo['fplus']==F(169,256) and hi['fplus']==F(121,256),'n5 endpoints')
check(all(o['fplus']==o['minimum'] for o in [five_result,lo,hi]),'n5 minimum branch')
residual=five_result['fplus']-(1-2*q)*L
check(residual==-F(7,8192),'n5 negative factor residual')
check(residual<0,'n5 hostile sign')
five_result.update(q=q,L=L,proposed_factor=1-2*q,residual=residual)

p=F(1,2);q=p*(1-p)
A=ordinal_pair_two_triangles(p);check(valid(A),'n8 validity')
eight=margins(A)
lo=margins(ordinal_pair_two_triangles(F(0)));hi=margins(ordinal_pair_two_triangles(F(1)))
L=(1-p)*lo['fplus']+p*hi['fplus']
check(eight['u']==4*(1-q) and eight['detA']==-q,'n8 determinant formulas')
check(eight['fplus']==4-5*q,'n8 margin formula')
check(L==4,'n8 endpoint interpolation')
residual=eight['fplus']-(1-q)*L
check(residual==-q and residual<0,'n8 hostile sign')
check(eight['minimum']==eight['fplus'],'n8 minimum branch')
check(eight['minimum']>1-q,'n8 not an absolute-floor counterexample')
eight.update(q=q,L=L,proposed_factor=1-q,residual=residual)
print(json.dumps(stringify({'status':'PASS','n5_fractional_cross_factor_1_minus_2q':five_result,'n8_matching_support_factor_1_minus_q':eight}),sort_keys=True,indent=2))
