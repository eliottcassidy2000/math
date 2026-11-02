#!/usr/bin/env python3
"""Analytic finite-field cycle counts and exact padding matrices; no tournament search."""
from math import gcd,lcm,factorial
import json

def check(ok,label):
    if not ok:raise RuntimeError(label)

def factors(n):
    fs=[];p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:n//=p;a+=1
            fs.append((p,a))
        p+=1
    if n>1:fs.append((n,1))
    return fs

def divisors(n):return [d for d in range(1,n+1) if n%d==0]

def mobius(n):
    fs=factors(n)
    return 0 if any(a>1 for p,a in fs) else (-1)**len(fs)

def fixed(p,d):return 1 if d%2 else 2*p**(gcd(d,p-1)//2)-1

def exact_periods(p):
    return {ell:sum(mobius(ell//d)*fixed(p,d) for d in divisors(ell)) for ell in divisors(p-1)}

def order(x,p):
    r=1;y=x%p
    while y!=1:y=y*x%p;r+=1
    return r

def support_periods(p):
    d=(p-1)//2;orders=[lcm(2,order(x,p)) for x in range(1,d+1)];out={ell:0 for ell in divisors(p-1)};out[1]=1
    for mask in range(1,1<<d):
        period=1
        for i,a in enumerate(orders):
            if mask>>i&1:period=lcm(period,a)
        out[period]+=2*(p-1)**mask.bit_count()
    return out

def clock(M):
    if M==1:return 1,1
    fs=factors(M);A=max(a for p,a in fs);L=1
    for p,a in fs:
        ell=(p-1)*p**(a-1) if p!=2 else 1 if a<=2 else 2**max(1,a-3)
        L=lcm(L,ell)
    return A,L

def precisions(M):
    c=1;d=1;out=[]
    while d%M:
        out.append(M//gcd(M,d));c+=1;d*=c*(c-1)
    return out

def pad(v,mods):return tuple(((i+1)*v[i]+(v[i-1] if i else 0))%m for i,m in enumerate(mods))

def matrix_checks():
    ms=list(range(1,129))+[275,1331,3025];basis_checks=0;clock_checks=0;weighted_checks=0
    for M in ms:
        A,L=clock(M);mods=precisions(M);k=len(mods);weights=[factorial(c)*factorial(c-1)%M for c in range(1,k+1)]
        check(all(mods[i-1]%mods[i]==0 for i in range(1,k)),('precision',M))
        for j in range(k):
            v=tuple(int(i==j)%mods[i] for i in range(k));base=None
            for r in range(A+L):
                if r==A-1:base=v
                if r==A-1+L:check(v==base,('matrix law',M,j));basis_checks+=1;break
                oldq=tuple(w*z%M for w,z in zip(weights,v))
                nv=pad(v,mods);newq=tuple(w*z%M for w,z in zip(weights,nv))
                expected=tuple(((i+1)*oldq[i]+((i+1)*i*oldq[i-1] if i else 0))%M for i in range(k))
                check(newq==expected,('weighted matrix',M,j,r));weighted_checks+=1;v=nv
        # Exact onset/period witnessed by K1; padding exponent starts at zero.
        v=tuple(int(i==0)%mods[i] for i in range(k));seen={};r=0
        while v not in seen:seen[v]=r;v=pad(v,mods);r+=1
        check((seen[v],r-seen[v])==(A-1,L),('sharp clock',M));clock_checks+=1
    return {'moduli':len(ms),'formal_basis_checks':basis_checks,'weighted_matrix_checks':weighted_checks,'sharp_seed_clocks':clock_checks}

def main():
    primes=[n for n in range(3,98,2) if factors(n)==[(n,1)]]
    support_primes=[3,5,7,11,13,17,19,23,29,31];root_checks=0;tables=[]
    for p in primes:
        E=exact_periods(p)
        check(sum(E.values())==2*p**((p-1)//2)-1,('total',p))
        check(all(z>=0 and z%ell==0 for ell,z in E.items()),('cycles',p))
        check(all(z==0 for ell,z in E.items() if ell>1 and ell%2),('odd period',p))
        for d in range(1,2*(p-1)+1):
            roots=sum(pow(x,d,p)==1 for x in range(1,p));check(roots==gcd(d,p-1),('root count',p,d));root_checks+=1
        if p in support_primes:check(E==support_periods(p),('support inversion',p))
        if p in [3,5,7,11]:tables.append({'prime':p,'actual_complete':True,'periods':[{'period':ell,'states':z,'cycles':z//ell} for ell,z in E.items() if z]})
    check(fixed(11,1)==1 and fixed(11,2)==21,'odd fixed/zero control')
    indicator=[0]*5;indicator[1]=1
    check(all(x*indicator[x]%5==indicator[x] for x in range(5)),'full-space hostile fixed')
    check(not any(all(indicator[-x%5]==e*indicator[x]%5 for x in range(5)) for e in [-1,1]),'full-space hostile excluded')
    print(json.dumps({'status':'PASS','formal_prime_checks':len(primes),'support_crosscheck_primes':support_primes,'root_of_unity_checks':root_checks,'actual_cycle_tables':tables,'padding_matrix':matrix_checks(),'limits':'No new actual-state realization asserted outside primes3,5,7,11.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
