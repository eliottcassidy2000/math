"""Independent marked-graph coefficient compiler; no unmarked partition loop.

This is a proof prototype, not the practical interval backend. It implements
arbitrary fixed correction order, with factorial graph enumeration constants.
"""
from fractions import Fraction as F
from math import gcd, comb, factorial
from itertools import combinations
from functools import lru_cache
from pathlib import Path
import sys, json, time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'enumeration'))
from backbone import configs,bweight,sector,parts,elementary

def assignments(n,m,allowed):
    if not m:
        yield ();return
    for a in allowed:
        if a>n:break
        for tail in assignments(n-a,m-1,allowed):yield (a,)+tail

@lru_cache(None)
def graphs(m,k):
    edges=list(combinations(range(m),2));out=[]
    for chosen in combinations(edges,k):
        if set(v for e in chosen for v in e)==set(range(m)):out.append(chosen)
    return out

@lru_cache(None)
def residual_dp(N,B,ms,A):
    """Unmarked cycles in the external field of backbone and marked cycles."""
    v=[F(0)]*(N+1);v[0]=F(1)
    for b in range(1,N+1,2):
        if b in B:continue
        R=sum(c*gcd(a,b) for a,c in zip(B,ms))+sum(min(a,b) for a in A)
        w=[F(0)]*(N+1)
        for s,x in enumerate(v):
            if not x:continue
            for j in range((N-s)//b+1):
                e=j*s+j*(b-1)//2+b*comb(j,2)+j*R
                w[s+b*j]+=x*F(1<<e,b**j*factorial(j))
        v=w
    return v[N]

def marked_coefficient(n,k,B=()):
    B=tuple(sorted(B));answer=F(0);nonzero=0
    for ms in configs(n,B):
        N=n-sum(a*c for a,c in zip(B,ms));bw=bweight(B,ms)
        if k==0:
            answer+=bw*residual_dp(N,B,ms,());continue
        allowed=tuple(a for a in range(3,N+1,2) if a not in B)
        for m in range(2,2*k+1):
            if comb(m,2)<k:continue
            for A in assignments(N,m,allowed):
                Q={e:1-F(1,1<<(min(A[e[0]],A[e[1]])-gcd(A[e[0]],A[e[1]]))) for e in combinations(range(m),2)}
                mark=F(0)
                for G in graphs(m,k):
                    prod=F(1)
                    for e in G:prod*=Q[e]
                    mark+=prod
                if not mark:continue
                exponent=sum((a-1)//2 for a in A)+sum(min(a,b) for a,b in combinations(A,2))
                exponent+=sum(c*gcd(a,b) for a in A for b,c in zip(B,ms))
                denom=factorial(m)
                for a in A:denom*=a
                answer+=bw*mark*F(1<<exponent,denom)*residual_dp(N-sum(A),B,ms,tuple(sorted(A)))
                nonzero+=1
    return answer,nonzero

def main():
    start=time.time();rows=[]
    for n in range(19):
        for B in [(),(3,)]:
            ref=[F(0)]*3
            for p in parts(n):
                w,q=sector(p,B)
                for k,e in enumerate(elementary(q,2)):ref[k]+=w*e
            vals=[];marks=[]
            for k in range(3):
                c,num=marked_coefficient(n,k,B)
                if c!=ref[k]:raise ArithmeticError((n,B,k,c,ref[k]))
                vals.append(str(c));marks.append(num)
            rows.append({'n':n,'backbone':list(B),'C0_C1_C2':vals,'nonzero_marked_assignments':marks})
    # A three-edge star at total14 is the first k3 boundary in the plain model.
    for n,B in [(14,()),(18,()),(22,(3,))]:
        c,num=marked_coefficient(n,3,B)
        ref=sum((w*elementary(q,3)[3] for p in parts(n) for w,q in [sector(p,B)]),F(0))
        if c!=ref:raise ArithmeticError(('degree3',n,B,c,ref))
        rows.append({'n':n,'backbone':list(B),'C3':str(c),'nonzero_marked_assignments':num})
    print(json.dumps({'status':'PASS','seconds':time.time()-start,'rows':rows},indent=2))
def order64():
    # Complete producer: no reference partition enumeration, no known answer.
    n=64;B=(3,5,7,9,11);first=min(a for a in range(3,n+2,2) if a not in B)
    r=n//first;K=comb(r,2);rows=[];value=F(0)
    for k in range(K+1):
        start=time.perf_counter();c,num=marked_coefficient(n,k,B)
        rows.append({'k':k,'coefficient':str(c),'nonzero_marked_assignments':num,
                     'seconds':time.perf_counter()-start})
        value+=(-1)**k*c
    if value.denominator!=1:raise ArithmeticError('nonintegral completed count')
    print(json.dumps({'status':'PASS','n':n,'backbone':list(B),'max_residual_cycles':r,
        'max_defect_edges':K,'exact_count':str(value.numerator),'rows':rows},indent=2))

if __name__=='__main__':
    if '--order64' in sys.argv:order64()
    else:main()
