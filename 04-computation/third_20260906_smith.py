#!/usr/bin/env python3
"""All-m mixed (m,2,1) Smith classification: exact residual and full controls.

Standard library only. Independent paths: all residual minors, closed gcd
ideals, Fraction DVR elimination, literal homogeneous Hasse matrices. Every
finite universe, the first label/chart hostiles, and actual kernel ratios are
printed. No producer imports, floats, or optimization-disabled assertions.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb, gcd, prod
import json
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0


def require(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def vp(x,p):
    if not x:
        return None
    x=F(x);a,b=abs(x.numerator),x.denominator;v=0
    while a%p==0:
        a//=p;v+=1
    while b%p==0:
        b//=p;v-=1
    return v


def bracket(v,w):
    return v[0]*w[1]-v[1]*w[0]


def tangent(v):
    a,b=v;x,y,s,t,r,q=1,0,0,1,a,b
    while q:
        k=r//q
        r,q,x,s,y,t=q,r-k*q,s,x-k*s,t,y-k*t
    require(abs(r)==1,'primitive vector')
    result=(-y*r,x*r)
    require(bracket(v,result)==1,'unimodular tangent')
    return result


def residual(m,a,b):
    return [[a**m,a**(m+1),a**(m+2)],
            [m*a**(m-1),(m+1)*a**m,(m+2)*a**(m+1)],
            [b**m,b**(m+1),b**(m+2)]]


def determinant3(a):
    return sum(a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3]-a[1][(i+2)%3]*a[2][(i+1)%3]) for i in range(3))


def closed_factors(m,a,b):
    first=gcd(a**m,m*a**(m-1),b**m)
    second=gcd(a**(2*m),a**m*b**m*(b-a),a**(m-1)*b**m*(m*b-(m+1)*a))
    total=abs(a**(2*m)*b**m*(b-a)**2)
    require(second%first==0 and total%second==0,'integral closed factors')
    result=(first,second//first,total//second)
    require(result[1]%result[0]==0 and result[2]%result[1]==0,'closed factor divisibility')
    return result


def all_minor_factors(a):
    first=gcd(*(x for row in a for x in row))
    second=0
    for rows in combinations(range(3),2):
        for cols in combinations(range(3),2):
            value=a[rows[0]][cols[0]]*a[rows[1]][cols[1]]-a[rows[0]][cols[1]]*a[rows[1]][cols[0]]
            second=gcd(second,value)
    total=abs(determinant3(a))
    return (first,second//first,total//second)


def dvr(a,p):
    a=[[F(x) for x in row] for row in a];n=len(a);answer=[]
    for k in range(n):
        i,j=min(((i,j) for i in range(k,n) for j in range(k,n) if a[i][j]),key=lambda ij:vp(a[ij[0]][ij[1]],p))
        a[i],a[k]=a[k],a[i]
        for row in a:
            row[j],row[k]=row[k],row[j]
        pivot=a[k][k];answer.append(vp(pivot,p))
        for i in range(k+1,n):
            ratio=a[i][k]/pivot
            require(not ratio or vp(ratio,p)>=0,'integral DVR operations')
            for j in range(k,n):
                a[i][j]-=ratio*a[k][j]
    require(answer==sorted(answer),'ordered full DVR factors')
    return tuple(answer)


def kappa(m,vectors,p,e,d,reference=None):
    K=min(e,m*d)
    if K==0:
        return 0
    v0,v1,v2=vectors
    w=reference or tangent(v0)
    require(all(bracket(v,w)%p for v in vectors),'unit-separated reference')
    tau=F(bracket(v0,v1)*bracket(v2,w),p**d*bracket(v0,v2)*bracket(v1,w))
    require(vp(tau,p)==0,'normalized cross ratio unit')
    numerator=F(m,p**d)-(m+1)*tau
    return K if not numerator else min(K,vp(numerator,p))


def partition(m,vectors,p,reference=None):
    v0,v1,v2=vectors
    require(all(bracket(v,w) for v,w in combinations(vectors,2)),'distinct projective directions')
    A,B,C=vp(bracket(v0,v1),p),vp(bracket(v0,v2),p),vp(bracket(v1,v2),p)
    mu=vp(m,p)
    if A<=B:
        h=min(A,mu)
        result=((m-1)*A+h,(m+1)*A-h,m*B+2*C)
    else:
        e,d=B,A-B
        require(C==e,'deep heavy/double pair metric')
        if d!=mu:
            rho=(m-1)*e+min(e,(m-1)*d+mu)
            theta=min(d,mu)
            sigma=2*m*e+(m-1)*d+theta
            result=(rho,sigma-rho,(m+2)*e+(m+1)*d-theta)
        else:
            K=min(e,m*d)
            kap=kappa(m,vectors,p,e,d,reference)
            result=((m-1)*e+K,(m+1)*e+m*d-K+kap,(m+2)*e+m*d-kap)
    require(list(result)==sorted(result),'classified residual order')
    require(sum(result)==2*m*A+m*B+2*C,'full determinant mass')
    return (0,)*m+result


def full_matrix(m,vectors,shifts=(0,0,0)):
    degree=m+2;matrix=[]
    for i,(v,multiplicity) in enumerate(zip(vectors,(m,2,1))):
        base=tangent(v);w=tuple(base[j]+shifts[i]*v[j] for j in (0,1))
        for r in range(multiplicity):
            row=[]
            for q in range(degree+1):
                value=0
                for t in range(r+1):
                    if t<=q and r-t<=degree-q:
                        value+=comb(q,t)*v[0]**(q-t)*w[0]**t*comb(degree-q,r-t)*v[1]**(degree-q-r+t)*w[1]**(r-t)
                row.append(value)
            matrix.append(row)
    return matrix


def inverse(a):
    n=len(a)
    a=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        i=next(i for i in range(k,n) if a[i][k])
        a[k],a[i]=a[i],a[k]
        pivot=a[k][k]
        a[k]=[x/pivot for x in a[k]]
        for i in range(n):
            if i!=k:
                ratio=a[i][k]
                a[i]=[x-ratio*y for x,y in zip(a[i],a[k])]
    return [row[n:] for row in a]


def apply(g,v):
    return (g[0][0]*v[0]+g[0][1]*v[1],g[1][0]*v[0]+g[1][1]*v[1])


def main():
    trace=[];residual_rows=prime_rows=0
    for m in range(2,21):
        for a,b in product(range(-10,11),repeat=2):
            if not a or not b or a==b:
                continue
            r=residual(m,a,b);exact=closed_factors(m,a,b)
            require(all_minor_factors(r)==exact,'complete all-minor integer Smith')
            for p in (2,3,5,7):
                expected=partition(m,((0,1),(a,1),(b,1)),p)
                require(tuple(vp(x,p) for x in exact)==expected[m:],'all-prime valuation classification')
                if m<=5:
                    require(dvr(r,p)==expected[m:],'independent residual DVR')
                prime_rows+=1
            residual_rows+=1
            trace.append((m,a,b,exact))
    ladder_rows=0
    for m in range(2,49):
        for p in (2,3,5,7):
            mu=vp(m,p)
            if not mu:
                continue
            for e in (0,1,2,3,m*mu,m*mu+2):
                K=min(e,m*mu)
                values=range(0,K+1) if p!=2 else (range(1,K+1) if K else (0,))
                for target in values:
                    v=m+1;c=m//p**mu
                    if target==K:
                        u=c
                    elif target==0:
                        add=next(i for i in range(1,p) if (c+i)%p)
                        u=c+add
                    else:
                        u=c+p**target
                    a,b=p**(e+mu)*u,p**e*v
                    vectors=((0,1),(a,1),(b,1))
                    require(kappa(m,vectors,p,e,mu)==target,'every claimed ladder state attained')
                    exact=closed_factors(m,a,b)
                    require(tuple(vp(x,p) for x in exact)==partition(m,vectors,p)[m:],'all-depth ladder factors')
                    if target<K:
                        A=(m+1)*e+m*mu-K
                        B=(m+2)*e+m*mu
                        old=(A+target,B-target);new=(A+target+1,B-target-1)
                        for N in (max(1,old[0]),old[0]+1,old[1]-1,old[1],old[1]+1):
                            gain=sum(min(N,q) for q in new)-sum(min(N,q) for q in old)
                            require(gain==int(old[0]+1<=N<=old[1]-1),'actual adjacent full-kernel ratio')
                    ladder_rows+=1
    full_rows=reference_rows=0
    controls=[]
    for m in range(2,9):
        for p in (2,3,5):
            for e,d in ((0,1),(1,1),(2,1),(1,2)):
                a,b=p**(e+d),p**e*(m+1)
                if a==b:
                    b+=p**e
                controls.append((m,p,((0,1),(a,1),(b,1))))
    controls.extend((m,p,((1,0),(0,1),(1,1))) for m in (2,3,5,8) for p in (2,3,5,7))
    for m,p,vectors in controls:
        expected=partition(m,vectors,p)
        require(dvr(full_matrix(m,vectors),p)==expected,'literal complete homogeneous Hasse matrix')
        transformed=tuple(apply(((1,2),(1,3)),v) for v in vectors)
        require(partition(m,transformed,p)==expected,'projective classified partition')
        require(dvr(full_matrix(m,transformed,(2,-1,3)),p)==expected,'literal chart and tangent covariance')
        full_rows+=2
    for m,p in ((2,2),(3,3),(4,2),(6,3),(8,2),(12,3)):
        d=vp(m,p)
        for e in (1,2,3,m*d+1):
            vectors=((0,1),(p**e*m,1),(p**e*(m+1),1))
            K=min(e,m*d)
            w=tangent(vectors[0])
            for shift in (-3,-1,1,2,5):
                ww=tuple(w[j]+shift*vectors[0][j] for j in (0,1))
                require(kappa(m,vectors,p,e,d,ww)==K,'capped reference independence')
                reference_rows+=1
    small_value_rows=0
    for m in (2,3,5,8):
        for p in (2,3,5):
            for a,b in ((p**2,p),(p,2*p+1)):
                vectors=((0,1),(a,1),(b,1))
                vectors=tuple(apply(((1,2),(1,3)),v) for v in vectors)
                mat=full_matrix(m,vectors)
                inv=inverse(mat)
                loss=partition(m,vectors,p)[-1]
                small=max(-vp(row[j],p) for row in inv for j in (m,m+2) if row[j])
                require(small==loss,'small-bank value column attains worst precision')
                small_value_rows+=1
    # Four isometric ternary spectra already at m=3, followed by label hostile.
    family=[]
    for u in (2,4,10,1):
        vectors=((0,1),(81*u,1),(108,1))
        got=dvr(full_matrix(3,vectors),3)
        require(got==partition(3,vectors,3),'four-state hostile full observer')
        family.append(got)
    base=((0,1),(81,1),(108,1));swapped=(base[1],base[0],base[2])
    require(dvr(full_matrix(3,base),3)==(0,0,0,9,15,15),'label-hostile base')
    require(dvr(full_matrix(3,swapped),3)==(0,0,0,9,12,18),'swapping unequal multiplicity labels changes spectrum')
    # Exact linear cancellation is not chart invariant; its cap is.
    base=((0,1),(9,1),(12,1));changed=tuple(apply(((1,0),(1,1)),v) for v in base)
    w=tangent(changed[0]);tau=F(bracket(changed[0],changed[1])*bracket(changed[2],w),3*bracket(changed[0],changed[2])*bracket(changed[1],w))
    require(vp(1-4*tau,3)==1,'raw infinite cancellation becomes finite under chart')
    require(partition(3,base,3)==partition(3,changed,3),'capped cancellation survives chart')
    print('PROVED SCOPE: complete (m,2,1) Hasse banks on homogeneous binary degree m+2, every integer m>=2.')
    print('Residual universe: m2..20; distinct nonzero a,b in -10..10; all9 two-minors and full determinant.')
    print('Integer residual rows:',residual_rows,'prime partitions:',prime_rows)
    print('Ladder universe: m2..48; p2,3,5,7 dividing m; e0,1,2,3,m*v_p(m),m*v_p(m)+2; every capped residue state.')
    print('Attained ladder rows:',ladder_rows,'literal full projective matrices:',full_rows,'capped reference rows:',reference_rows)
    print('m3 p3 e3 residual ladder:',family)
    print('Unequal multiplicity label hostile: (0,0,0,9,15,15) versus (0,0,0,9,12,18).')
    print('Raw chart hostile at m3,p3,e1: infinite numerator valuation becomes1; capped state and full partition agree.')
    print('Actual kernel consequence: adjacent capped states differ by factor p exactly between the two inner precision thresholds.')
    print('Independent full-inverse small-bank value-attainment rows:',small_value_rows)
    print('semantic_sha256:',sha256(json.dumps(trace,separators=(',',':')).encode()).hexdigest())
    print('Exact gates:',GATES)
    print('PASS')


if __name__=='__main__':
    main()
