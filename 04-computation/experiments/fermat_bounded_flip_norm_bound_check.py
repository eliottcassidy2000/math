#!/usr/bin/env python3
"""Exact checks for the bounded-carousel-flip theorem; standard library only."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json

CHECKS=0

def require(condition,message):
    global CHECKS
    CHECKS+=1
    if not condition:raise ArithmeticError(message)

def determinant(a):
    a=[row[:] for row in a]; n=len(a); previous=1; sign=1
    for k in range(n-1):
        if not a[k][k]:
            pivot_row=next((i for i in range(k+1,n) if a[i][k]),None)
            if pivot_row is None:return 0
            a[k],a[pivot_row]=a[pivot_row],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                if numerator%previous:raise ArithmeticError('Nonexact Bareiss division')
                a[i][j]=numerator//previous
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]

def cyclotomic_norm(p,coefficients):
    # Multiplication in Z[x]/Phi_p, independent of tournament determinant route.
    n=p-1; matrix=[[0]*n for _ in range(n)]
    for j in range(n):
        for exponent,coefficient in enumerate(coefficients):
            e=(j+exponent)%p
            if e==p-1:
                for i in range(n):matrix[i][j]-=coefficient
            else:matrix[e][j]+=coefficient
    return determinant(matrix)

def symbol(p,flips):return (set(range(1,(p+1)//2))-set(flips))|{p-j for j in flips}

def reduced_polynomial(p,flips):
    c=[0]*p;c[0]=1
    for j in flips:
        for e,sign in ((-2*j,1),(1-2*j,1),(2*j,-1),(2*j+1,-1)):
            c[e%p]+=sign
    require(all(x in (-1,0,1) for x in c),'Same-sign coefficient collision')
    require(sum(c)==1,'Polynomial at one changed')
    require(sum(x*x for x in c)<=1+4*len(flips),'Coefficient energy bound failed')
    return c

def tournament(p,S):return [[int((j-i)%p in S) for j in range(p)] for i in range(p)]
def disc(A):
    p=len(A); d=determinant([[int(i==j)+A[i][j]-A[j][i] for j in range(p)]for i in range(p)])
    require(d%2**(p-1)==0,'Nonintegral normalized determinant')
    return d//2**(p-1)

def triangle_packing_polynomial(A):
    n=len(A)
    triangles=[(1<<i)|(1<<j)|(1<<k) for i,j,k in combinations(range(n),3)
               if A[i][j]+A[j][k]+A[k][i] in (0,3)]
    containing=[[t for t in triangles if t&(1<<v)]for v in range(n)]
    require(len(triangles)==n*(n*n-1)//24,'Regular triangle total')
    require(all(len(ts)==(n*n-1)//8 for ts in containing),'Regular rooted triangle total')
    @lru_cache(None)
    def packs(mask):
        if not mask:return (1,)
        v=(mask&-mask).bit_length()-1
        out=list(packs(mask^(1<<v)))
        for t in containing[v]:
            if mask&t==t:
                child=packs(mask^t)
                while len(out)<len(child)+1:out.append(0)
                for k,count in enumerate(child):out[k+1]+=count
        return tuple(out)
    coefficients=packs((1<<n)-1)
    require(coefficients[1]==len(triangles),'Triangle packing coefficient one')
    if len(coefficients)>2:
        direct=sum(not(a&b) for a,b in combinations(triangles,2))
        require(coefficients[2]==direct,'Independent disjoint-triangle-pair count')
    return list(coefficients)

def main():
    symbol_controls=0
    for p in (3,5,7,11,17):
        m=(p-1)//2
        for mask in range(1<<m):
            flips=[j+1 for j in range(m) if (mask>>j)&1]
            c=reduced_polynomial(p,flips);S=symbol(p,flips);A=tournament(p,S)
            B={0}|{2*s%p for s in S};bits=[int(e in B) for e in range(p)]
            require(c==[bits[e]+bits[(e-1)%p]-1 for e in range(p)],'Cyclic run-defect coefficient identity')
            require(sum(x*x for x in c)==sum(bits[e]==bits[(e-1)%p] for e in range(p)),'Cyclic equal-neighbor energy identity')
            d=disc(A); norm=cyclotomic_norm(p,c)
            require(d==norm,'Norm carrier disagrees with skew determinant')
            require(d%p==1,'Prime cyclic norm congruence')
            if flips:
                # Square both sides to avoid nonintegral exponent notation.
                B=Fraction((1+4*len(flips))*p-1,p-1)
                require(d*d<=B**(p-1),'Parseval-AMGM norm bound')
            else:require(d==1,'Unflipped carousel must have unit norm')
            symbol_controls+=1
    small_rows=[]
    for p in (3,5,17):
        for flip in range(1,(p+1)//2):
            S=symbol(p,[flip]);A=tournament(p,S);d=disc(A)
            coeff=triangle_packing_polynomial(A)
            lower=sum((2**i)*c for i,c in enumerate(coeff))
            require(lower>=p*d,'Small Fermat triangle-packing lower bound')
            small_rows.append({'p':p,'flipped_step':flip,'disc':d,'triangle_packing_coefficients':coeff,'triangle_only_OCF_lower':lower,'target_p_disc':p*d})
    C=Fraction(321,64);threshold=[]
    for n in (257,259,261,263,265):
        lower=Fraction(n*n-1,12)**(n//10); upper=n*C**((n-1)//2)
        ratio=lower/upper
        require(ratio>1,'One-flip threshold base case')
        threshold.append({'n':n,'k':n//10,'exact_lower_exceeds_upper':True,'approximate_ratio':float(ratio)})
    step=Fraction(267**2-1,12)*Fraction(257,267)/C**5
    require(step>1,'Ten-step threshold induction')
    require(6**3>=22*(3+1),'General bounded-flip exponential base')
    print(json.dumps({'status':'PROVED in accompanying note; finite arithmetic controls below','symbol_controls':symbol_controls,'one_flip_small_fermat_controls':small_rows,'one_flip_threshold_controls':threshold,'ten_step_ratio_lower':str(step),'general_sufficient_threshold':'p >= max(257,4*(4*r+2)^3), r>=1','p65537_certified_flips_by_general_threshold':5,'checks':CHECKS},indent=2,sort_keys=True))
if __name__=='__main__':main()
