#!/usr/bin/env python3
"""Exact finite checks; Python standard library, no theorem inferred from sampling."""
import itertools, json, math
from fractions import Fraction
from pathlib import Path

def det(a):
    a=[r[:] for r in a]; n=len(a); prev=1; sign=1
    for k in range(n-1):
        if not a[k][k]:
            r=next(i for i in range(k+1,n) if a[i][k]);a[k],a[r]=a[r],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                if numerator%prev: raise ArithmeticError('Nonexact Bareiss division')
                a[i][j]=numerator//prev
        for i in range(k+1,n):a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

def paths_fixed_start(A):
    n=len(A); m=n-1; size=1<<m; dp=[0]*(size*m)
    incoming=[sum(1<<u for u in range(m) if A[u+1][v+1]) for v in range(m)]
    for mask in range(1,size):
        bits=mask
        while bits:
            b=bits&-bits;bits-=b;v=b.bit_length()-1;prior=mask^b
            if not prior:dp[mask*m+v]=A[0][v+1]
            else:
                ubs=prior&incoming[v];total=0
                while ubs:
                    ub=ubs&-ubs;ubs-=ub;u=ub.bit_length()-1
                    total+=dp[prior*m+u]
                dp[mask*m+v]=total
    return sum(dp[(size-1)*m:])

def cyclotomic_norm(p,S):
    # Multiplication by 1 + sum_{s in S} x^s on Z[x]/(1+...+x^(p-1)).
    n=p-1; matrix=[[0]*n for _ in range(n)]
    for j in range(n):
        for s in [0]+list(S):
            e=(j+s)%p
            if e==p-1:
                for i in range(n):matrix[i][j]-=1
            else:matrix[e][j]+=1
    return det(matrix)

def circulant_rows(p):
    m=(p-1)//2;seen=set();rows=[]
    for mask in range(1<<m):
        S=frozenset(i+1 if (mask>>i)&1 else p-i-1 for i in range(m))
        if S in seen:continue
        orbit={frozenset(u*s%p for s in S) for u in range(1,p)}
        assert len(orbit)==p-1
        seen.update(orbit)
        A=[[int((j-i)%p in S) for j in range(p)] for i in range(p)]
        H=p*paths_fixed_start(A)
        K=[[A[i][j]-A[j][i] for j in range(p)] for i in range(p)]
        M=[[int(i==j)+K[i][j] for j in range(p)] for i in range(p)]
        determinant=det(M);d=Fraction(determinant,2**(p-1))
        assert d.denominator==1 and H%p==0 and H%2==1
        norm=cyclotomic_norm(p,S); assert norm==d
        B=[[sum(K[i][k]*K[j][k] for k in range(p)) for j in range(p)] for i in range(p)]
        variance=sum((B[i][j]-p*int(i==j)+1)**2 for i in range(p) for j in range(p))
        if p%4==1:assert variance>=4*p*(p-1)
        if p<=5:
            direct=sum(all(A[v[i]][v[i+1]] for i in range(p-1)) for v in itertools.permutations(range(p)))
            assert direct==H
        rows.append({'p':p,'S':sorted(S),'H':H,'disc':int(d),'H_over_p_disc':str(Fraction(H,p*d)), 'spectral_variance':variance,'cyclotomic_norm':norm})
    assert len(seen)==1<<m
    assert len(rows)==(1<<m)//(p-1)
    return rows

def fadd(x,y):
    return sum(((x//3**k)%3+(y//3**k)%3)%3*3**k for k in range(3))
def fneg(x):return sum((-(x//3**k)%3)%3*3**k for k in range(3))
def fmul(x,y):
    a=[(x//3**k)%3 for k in range(3)];b=[(y//3**k)%3 for k in range(3)]
    c=[sum(a[i]*b[k-i] for i in range(3) if 0<=k-i<3)%3 for k in range(5)]
    for k in (4,3):
        c[k-2]=(c[k-2]+c[k])%3;c[k-3]=(c[k-3]+c[k])%3;c[k]=0
    return sum(c[k]*3**k for k in range(3))
def fpow(x,k):
    r=1
    while k:
        if k&1:r=fmul(r,x)
        x=fmul(x,x);k//=2
    return r

def paley27():
    q=27;Q={fmul(a,a) for a in range(1,q)};assert len(Q)==13
    assert all(fpow(x,26)==1 for x in range(1,q))
    assert Q.isdisjoint({fneg(x) for x in Q})
    A=[[int(fadd(y,fneg(x)) in Q) for y in range(q)] for x in range(q)]
    K=[[A[i][j]-A[j][i] for j in range(q)] for i in range(q)]
    for i in range(q):
        for j in range(q):assert sum(K[i][k]*K[j][k] for k in range(q))==q*int(i==j)-1
    perms=set()
    for a,b,j in itertools.product(Q,range(q),range(3)):
        perm=tuple(fadd(fmul(a,fpow(x,3**j)),b) for x in range(q));assert len(set(perm))==q
        assert all(A[x][y]==A[perm[x]][perm[y]] for x in range(q) for y in range(q))
        perms.add(perm)
    assert len(perms)==1053
    determinant=det([[int(i==j)+K[i][j] for j in range(q)] for i in range(q)])
    assert determinant==28**13
    return {'field':'F_3[t]/(t^3-t-1)','q':27,'certified_automorphism_subgroup_order':len(perms), 'disc':7**13,'s':27,'strong_base_target':27*7**13, 'spectral_variance':0, 'construction':'x -> a*x^(3^j)+b, square a, j=0,1,2'}

def main():
    rows=sum((circulant_rows(p) for p in (3,5,17)),[])
    p17=[x for x in rows if x['p']==17]
    result={'snapshot':'0f7b1b200b1b96a6ce61bbe304bb92c603515bbf','status':'FINITE-EXACT; classification theorem proved separately', 'fermat_class_counts':{'3':'1','5':'1','17':'16','257':'2^120','65537':'2^32752'},'circulant_rows':rows,'p17_min_ratio':min(p17,key=lambda r:Fraction(r['H_over_p_disc'])),'p17_max_disc':max(p17,key=lambda r:r['disc']),'p17_min_variance':min(r['spectral_variance'] for r in p17),'p17_all_strong_base_pass':all(r['H']>=r['p']*r['disc'] for r in p17),'paley27':paley27()}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
