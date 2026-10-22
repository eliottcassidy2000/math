#!/usr/bin/env python3
"""Independent exact arithmetic controls; no imported producer code; assertions explicit."""
import math,json,hashlib
from pathlib import Path

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def strong(N,q):
    t=[pow(2,n*(n-1)//2,q) for n in range(N+1)]
    s=[0]*(N+1)
    for n in range(1,N+1):
        z=t[n]; c=1
        for k in range(1,n):
            c=c*(n-k+1)//k
            z-=c*s[k]*t[n-k]
        s[n]=z%q
    return s

rec=[]
for p,rs in [(3,range(1,6)),(5,range(1,5)),(7,range(1,4)),(11,range(1,3)),(13,range(1,3)),(17,range(1,3)),(19,range(1,3))]:
    for r in rs:
        q=p**r; P=q//p*(p-1); eps=1 if pow(2,(p-1)//2,p)==1 else -1
        N=P+max(4*r,2*p+20); s=strong(N,q)
        checks=0
        for n in range(2*r-1,N-P+1):
            need((s[n+P]-eps*s[n])%q==0,('recurrence',p,r,n))
            checks+=1
        rec.append(dict(p=p,r=r,checks=checks,lower_failures=[n for n in range(2*r-1) if (s[n+P]-eps*s[n])%q]))

# Integer-valued polynomial Newton coefficients, including degrees >= p.
newton=[]
for j in range(13):
    vals=[math.comb(n*(n-1)//2,j) if n*(n-1)//2>=j else 0 for n in range(2*j+2)]
    coeff=[]
    while vals:
        coeff.append(vals[0]); vals=[b-a for a,b in zip(vals,vals[1:])]
    need(coeff[-1]==0,('degree',j))
    for n in range(2*j+2,2*j+20):
        need(sum(a*math.comb(n,d) for d,a in enumerate(coeff))==math.comb(n*(n-1)//2,j),('Newton',j,n))
    newton.append(dict(j=j,coefficients=coeff[:-1]))

# R=(Z/p^r)[sqrt(D)]; exact roots/Fourier support and finite group-ring inverse.
ring=[]
for p,r in [(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2)]:
    q=p**r; L=2*(p-1); eps=1 if pow(2,(p-1)//2,p)==1 else -1
    D=next(d for d in range(2,p) if pow(d,(p-1)//2,p)==p-1)
    zero=(0,0); one=(1,0)
    def add(a,b): return ((a[0]+b[0])%q,(a[1]+b[1])%q)
    def mul(a,b): return ((a[0]*b[0]+D*a[1]*b[1])%q,(a[0]*b[1]+a[1]*b[0])%q)
    def scale(a,b): return (a[0]*b%q,a[1]*b%q)
    def power(a,n):
        ans=one
        while n:
            if n&1: ans=mul(ans,a)
            a=mul(a,a); n//=2
        return ans
    roots=[]
    for a in range(p):
        for b in range(p):
            z=power((a,b),p**(2*(r-1)))
            if power(z,L)==one and all(power(z,d)!=one for d in range(1,L) if L%d==0): roots.append(z)
    need(bool(roots),('root missing',p,r)); z=roots[0]
    omega=pow(2,p**(r-1),q)
    tt=[pow(omega,n*(n-1)//2,q) for n in range(L)]
    four=[]
    for k in range(L):
        alpha=power(z,k); c=zero
        for n in range(L): c=add(c,scale(power(alpha,(-n)%L),tt[n]))
        c=scale(c,pow(L,-1,q))
        if c!=zero:
            need(power(alpha,p-1)==(eps%q,0),('support',p,r,k,c))
            four.append((alpha,c))
    theta=one if eps==1 else z
    line={scale(theta,b):b for b in range(q)}
    A=[zero]*q
    for a,c in four: A[line[a]]=add(A[line[a]],c)
    aug=zero
    for a in A: aug=add(aug,a)
    need(aug==one,('augmentation',p,r))
    def conv(a,b):
        out=[zero]*q
        for i,x in enumerate(a):
            if x==zero: continue
            for j,y in enumerate(b):
                if y!=zero: out[(i+j)%q]=add(out[(i+j)%q],mul(x,y))
        return out
    unit=[one]+[zero]*(q-1)
    def gpower(a,n):
        out=unit
        while n:
            if n&1: out=conv(out,a)
            a=conv(a,a); n//=2
        return out
    need(gpower(A,p**(2*r-1))==unit,('finite group inverse exponent',p,r))
    for n in range(2*L+1):
        val=zero
        for alpha,c in four: val=add(val,mul(c,power(alpha,n)))
        need(val==(pow(omega,n*(n-1)//2,q),0),('Fourier expansion',p,r,n))
    u=((2*pow(omega,-1,q)%q)-1)//p if r>1 else 0
    for n in range(2*L+10):
        val=zero
        for j in range(r):
            for d,coeff in enumerate(newton[j]['coefficients']):
                if d>n: continue
                for alpha,c in four:
                    term=mul(mul(c,power(alpha,d)),power(alpha,n-d))
                    val=add(val,scale(term,p**j*u**j*coeff*math.comb(n,d)))
        need(val==(pow(2,n*(n-1)//2,q),0),('filtered T lift',p,r,n))
    ring.append(dict(p=p,r=r,D=D,root=z,theta=theta,fourier_terms=len(four),group_size=q,inverse_exponent=p**(2*r-1)))

receipt=dict(status='PASS',recurrence=rec,recurrence_checks=sum(x['checks'] for x in rec),newton=newton,rings=ring)
print(json.dumps(receipt,indent=2))
