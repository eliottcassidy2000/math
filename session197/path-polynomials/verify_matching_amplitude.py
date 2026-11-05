#!/usr/bin/env python3
from verify_extension import *
from math import factorial
import random

def cross_matching(a,r):
    O=[i for i in range(len(a)) if a[r]>>i&1];I=[i for i in range(len(a)) if a[i]>>r&1];hist={}
    def visit(i,used,k):
        if i==len(O):hist[k]=hist.get(k,0)+1;return
        visit(i+1,used,k)
        for v in I:
            if v not in used and a[O[i]]>>v&1:visit(i+1,used|{v},k+1)
    visit(0,set(),0);nu=max(hist);return nu,hist[nu]

def balanced(k,crossbits,insidebits):
    # Root0, O1..k, I(k+1)..2k. Inside orientations are arbitrary.
    a=[0]*(2*k+1)
    for u in range(1,k+1):a[0]|=1<<u
    for v in range(k+1,2*k+1):a[v]|=1
    b=0
    for u in range(1,k+1):
        for v in range(k+1,2*k+1):
            if crossbits>>b&1:a[u]|=1<<v
            else:a[v]|=1<<u
            b+=1
    b=0
    for side in (range(1,k+1),range(k+1,2*k+1)):
        for u,v in combinations(side,2):
            x,y=(u,v) if insidebits>>b&1 else (v,u);a[x]|=1<<y;b+=1
    validate(a);return tuple(a)

def top_run_from_counts(a,r,nu):
    # Independent triangular inversion of ordered Stirling clone-run counts.
    st=[[0]*(nu+2) for _ in range(nu+2)];st[0][0]=1
    for m in range(1,nu+2):
        for k in range(1,m+1):st[m][k]=st[m-1][k-1]+k*st[m-1][k]
    runs=[0]*(nu+2)
    for m in range(1,nu+2):
        h=h_subsets(clone(a,r,m))[-1]
        remainder=h-sum(runs[k]*factorial(k)*st[m][k] for k in range(1,m))
        require(remainder%factorial(m)==0,'Stirling integral');runs[m]=remainder//factorial(m)
    return runs[-1]

def run():
    cases=[]
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask)
            for r in range(n):cases.append((a,r))
    for crossbits in range(16):
        for insidebits in range(4):cases.append((balanced(2,crossbits,insidebits),0))
    rng=random.Random(197901)
    for _ in range(12):cases.append((balanced(3,rng.getrandbits(9),rng.getrandbits(6)),0))
    for _ in range(20):
        n=rng.choice((6,7));cases.append((tournament(n,rng.getrandbits(n*(n-1)//2)),rng.randrange(n)))
    cases.append((balanced(3,0b100010001,0),0))
    cases.append((balanced(3,(1<<9)-1,0),0))
    checks=0;balanced_checks=0;high=[]
    for a,r in cases:
        nu,mu=cross_matching(a,r);op=skeleton_operator(a,r);N=len(a)-1
        coeff=[0]*(N+1)
        for (power,j),v in op.items():
            if j==nu:coeff[power]+=v
        coeff=trim(coeff);require(len(coeff)-1==N-nu and coeff[-1]==mu,('leading',a,r,nu,mu,coeff))
        rise=[1]
        for j in range(nu):rise=mul(rise,[j,1])
        lower=scale(mul(mon(N-2*nu),rise),mu)
        diff=add(coeff,scale(lower,-1));require(all(v>=0 for v in diff),('coefficient lower bound',a,r,coeff,lower))
        if N==2*nu:require(coeff==lower,'balanced equality');balanced_checks+=1
        top=top_run_from_counts(a,r,nu);require(top==value(coeff,1)>=mu*factorial(nu),'top amplitude')
        if nu>=2:high.append({'rows':a,'root':r,'nu':nu,'mu':mu,'forward_polynomial':coeff,'top_amplitude':top})
        checks+=1
    print(json.dumps({'rooted_templates':checks,'balanced_templates':balanced_checks,'higher_matching_controls':high},sort_keys=True,indent=2))
if __name__=='__main__':run()
