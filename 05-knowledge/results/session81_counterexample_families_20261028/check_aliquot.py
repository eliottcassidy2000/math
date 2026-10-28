#!/usr/bin/env python3
"""Exact integer checks; no assert, no external dependencies."""
from fractions import Fraction
from itertools import product
from math import isqrt
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out
def sigma(n):
    z=1
    for p,a in factors(n).items():z*=(p**(a+1)-1)//(p-1)
    return z
def s(n):return sigma(n)-n
def divisor_s(n):return sum(d+(n//d if d*d!=n else 0) for d in range(1,isqrt(n)+1) if n%d==0)-n
def valuation(n,p):
    a=0
    while n%p==0:n//=p;a+=1
    return a

def pp_inverse(x,y,p,b):
    q=p**b;h=(q-1)//(p-1)
    if (y-q*x)%h or (y-x)%h:return None
    c=(y-q*x)//h;d=(y-x)//h
    if c<=0 or d<=0:return None
    a=valuation(d,p);num=p*c-d//p**a
    if num%(p-1):return None
    m=num//(p-1)
    if m<=0 or m%p==0 or sigma(m)!=c:return None
    return p**a*m

def inverse_pair(x,y,k):
    """Complete fiber for (s(n),s(k*n)), k>=2; profile enumeration."""
    if k<2:raise ValueError('k >= 2 required')
    if x<0:return []
    if x==0:return [1] if y==s(k) else []
    if x==1:
        cand=set(factors(k));num=y-sigma(k);den=s(k)
        if num%den==0 and num//den>1:cand.add(num//den)
        return sorted(n for n in cand if s(n)==1 and s(k*n)==y)
    N=(x-1)**2;pf=list(factors(k).items());ranges=[]
    for p,b in pf:
        A=0;power=1
        while power*p<=N:A+=1;power*=p
        ranges.append(range(A+1))
    out=set()
    for aa in product(*ranges):
        R=Fraction(1);support=1
        for (p,b),a in zip(pf,aa):
            support*=p**a
            R*=Fraction(p**(a+b+1)-1,p**(a+1)-1)
        if support>N:continue
        need(R>k,'positive profile slope')
        n=(y-R*x)/(R-k)
        if n.denominator!=1 or not 1<=n<=N:continue
        n=int(n)
        if any(valuation(n,p)!=a for (p,b),a in zip(pf,aa)):continue
        if s(n)==x and s(k*n)==y:out.add(n)
    return sorted(out)

def main():
    for n in range(1,10001):need(s(n)==divisor_s(n),'independent divisor sum')
    pp=0
    for p,b in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(5,2),(7,2)]:
        for n in range(1,5001):
            need(pp_inverse(s(n),s(p**b*n),p,b)==n,'prime-power recovery');pp+=1
    need(pp_inverse(1,16,2,2)==3,'composite-four repaired')
    need(pp_inverse(1,4,2,1) is None,'image guard hostile')
    need(s(22655)==s(24275)==5857,'collision x')
    need(s(40*22655)==s(40*24275)==1303480,'collision y')
    complete=inverse_pair(5857,1303480,40)
    need(complete==[22655,24275],'complete collision pair fiber')
    # Exhaustive arbitrary-pair image test: x>=2 gives N<=841;
    # x=1 candidates include all prime sources because y<=200 implies n<=200.
    checked=0
    for k in [4,6,12,40]:
        image={}
        for n in range(1,842):image.setdefault((s(n),s(k*n)),[]).append(n)
        for x in range(0,31):
            for y in range(0,201):
                need(inverse_pair(x,y,k)==image.get((x,y),[]),'complete inverse image');checked+=1
    for k in range(2,201):
        if all(k%p for p in [3,5,7]):need(s(k*35)-s(k*27)==8*s(k),'coprime fiber transport')
    print(json.dumps({'status':'PASS','independent_sigma_inputs':10000,'prime_power_reconstructions':pp,'arbitrary_pair_checks':checked,'multiplier40_complete_pair_fiber':complete},indent=2))
if __name__=='__main__':main()
