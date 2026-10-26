#!/usr/bin/env python3
"""Exact aliquot inverses; finite universes are theorem bounded, not heuristic.
All arithmetic and primality checks are deterministic. No third-party packages.
"""
from math import isqrt, gcd

def icbrt(n):
    lo, hi = 0, 1 << ((n.bit_length()+2)//3)
    while lo+1 < hi:
        md=(lo+hi)//2
        if md**3 <= n: lo=md
        else: hi=md
    return hi if hi**3 <= n else lo

def spf_sieve(n):
    s=list(range(n+1))
    for p in range(2,isqrt(n)+1):
        if s[p]==p:
            for k in range(p*p,n+1,p):
                if s[k]==k: s[k]=p
    return s

def factors(n,spf=None):
    ans={}
    if spf is not None and n<len(spf):
        while n>1:
            p=spf[n]; ans[p]=ans.get(p,0)+1; n//=p
        return ans
    p=2
    while p*p<=n:
        while n%p==0: ans[p]=ans.get(p,0)+1; n//=p
        p=3 if p==2 else p+2
    if n>1: ans[n]=ans.get(n,0)+1
    return ans

def sigma_factorization(fs):
    ans=1
    for p,e in fs.items(): ans *= (p**(e+1)-1)//(p-1)
    return ans

def aliquot(n):
    return sigma_factorization(factors(n))-n

def sigma_sieve(n):
    # Ordinary divisor sieve; intentionally independent from factor-based sigma.
    sig=[0]*(n+1)
    for d in range(1,n+1):
        for k in range(d,n+1,d): sig[k]+=d
    return sig

def fiber(m, return_stats=False):
    """Entire s^-1(m), for m != 1; at 1 raise for infinite prime fiber.
    Split sources into an exponent-one prime extension and squarefull source.
    """
    if m<0: return ([],{}) if return_stats else []
    if m==0: return ([1],{}) if return_stats else [1]
    if m==1: raise ValueError('s^-1(1) is the infinite set of all primes')
    spf=spf_sieve(m)
    out=set(); core_hits=[]; full_hits=[]
    for a in range(2,m):
        sa=sigma_factorization(factors(a,spf)); denom=sa-a
        diff=m-sa
        if diff>0 and diff%denom==0:
            p=diff//denom
            if p>=2 and p<=m and spf[p]==p and a%p:
                out.add(a*p); core_hits.append([a,p,a*p])
    B=(m-1)**2
    squarefull_candidates=0
    for v in range(1,icbrt(B)+1):
        fv=factors(v,spf)
        if any(e>1 for e in fv.values()): continue
        for u in range(1,isqrt(B//(v**3))+1):
            n=u*u*v**3
            if n<4: continue
            squarefull_candidates += 1
            fs={p:2*e for p,e in factors(u,spf).items()}
            for p in fv: fs[p]=fs.get(p,0)+3
            if sigma_factorization(fs)-n==m:
                out.add(n); full_hits.append([u,v,n])
    stats={'target':m,'bound':B,'core_candidates':max(0,m-2),
           'squarefull_candidates':squarefull_candidates,
           'core_hits':core_hits,'squarefull_hits':full_hits}
    return (sorted(out),stats) if return_stats else sorted(out)

def even_prefix(X):
    """All even target fibers <=X using Pomerance--Yang's parity reduction.
    Sources even: <=2X (safe at small X). Odd squares: prime roots <=X-1,
    composite roots <=floor(cuberoot(X^2)). Return fibers, with no conjecture.
    """
    sig=sigma_sieve(max(2,2*X)); spf=spf_sieve(max(2,X))
    out={m:[] for m in range(2,X+1,2)}
    for n in range(2,2*X+1,2):
        m=sig[n]-n
        if m in out: out[m].append(n)
    for p in range(3,X,2):
        if spf[p]==p: out[p+1].append(p*p)
    root_bound=icbrt(X*X)
    for a in range(9,root_bound+1,2):
        if spf[a]==a: continue
        fs={p:2*e for p,e in factors(a,spf).items()}
        n=a*a; m=sigma_factorization(fs)-n
        if m in out: out[m].append(n)
    for m in out: out[m]=sorted(set(out[m]))
    return out

def goldbach_witness(m,spf):
    """Only a witness search, never an untouchability test."""
    if m==1: return 2
    if m==3: return 4
    if m==7: return 8
    for p in range(3,(m-1)//2,2):
        q=m-1-p
        if spf[p]==p and spf[q]==q: return p*q
    return None

def two_prime_core(m,a):
    """All guarded n=a*p*q with p<q primes coprime to a.
    Hyperbola identity (s(a)p+sigma(a))(s(a)q+sigma(a))=s(a)m+a*sigma(a).
    """
    if a<2: raise ValueError('a=1 degenerates to Goldbach; require a>=2')
    t=sigma_factorization(factors(a)); s=t-a; D=s*m+a*t
    out=[]
    for d in range(1,isqrt(D)+1):
        if D%d: continue
        e=D//d
        if (d-t)%s or (e-t)%s: continue
        p=(d-t)//s; q=(e-t)//s
        if p<2 or q<=p or gcd(a,p*q)!=1: continue
        if factors(p)!={p:1} or factors(q)!={q:1}: continue
        out.append(a*p*q)
    return sorted(set(out))
