#!/usr/bin/env python3
"""Exact bounded replay of continuation spectra. Standard library only."""
import json, math
from itertools import product
from fractions import Fraction

def require(x,msg):
    if not x: raise RuntimeError(msg)
def divisors(n): return [d for d in range(1,math.isqrt(n)+1) if n%d==0 for d in ({d,n//d})]
def phi(n): return sum(math.gcd(j,n)==1 for j in range(1,n+1))
def order2(n):
    if n==1:return 1
    x=2%n; t=1
    while x!=1:x=x*2%n;t+=1
    return t
def carry(w,b=3):
    C=0;H=0
    for a in w:C=b*C+2**H;H+=a
    return C,H

def spectrum(w,b=3):
    C,H=carry(w,b); k=len(w); N=b*C+2**H; A=b**(k+1)
    rows=[]
    for d in sorted(divisors(N)):
        t=order2(d)
        hits=[v for v in range(1,t+1) if (2**(H+v)-A)%d==0]
        require(len(hits)<=1,'cyclic phase uniqueness')
        if hits:rows.append({'divisor':d,'period':t,'phase':hits[0],'weight':phi(d)})
    return N,H,A,rows

def run():
    gates=0; anchors={}
    for w in [(1,1,5),(1,3,3)]:
        N,H,A,rows=spectrum(w); L=order2(N)
        hist={}; values=[]
        for v in range(1,3*L+1):
            D=2**(H+v)-A; g=math.gcd(N,D)
            predicted=sum(r['weight'] for r in rows if (v-r['phase'])%r['period']==0)
            require(g==predicted,'divisor expansion');gates+=1
            if v<=L:hist[g]=hist.get(g,0)+1
            values.append(D//g)
        for v in range(L):
            require(values[v+2*L]-(2**L+1)*values[v+L]+2**L*values[v]==0,'denominator recurrence');gates+=1
        anchors[str(w)]={'N':N,'period':L,'rows':rows,'gcd_histogram':hist,'mean_gcd':str(Fraction(sum(g*n for g,n in hist.items()),L)), 'first_denominators':values[:12]}
    # General direct expansion, prime powers included; modest explicit universe.
    tested=0; prime_power_examples=[]
    for b in (3,5,7):
      for k in range(1,4):
       for w in product(range(1,5),repeat=k):
        C,H=carry(w,b);N=b*C+2**H; A=b**(k+1)
        ds=sorted(divisors(N)); phis={d:phi(d) for d in ds}
        if any(N%(p*p)==0 for p in range(3,math.isqrt(N)+1,2)) and len(prime_power_examples)<8:
            prime_power_examples.append([b,list(w),N])
        for v in range(1,41):
            D=2**(H+v)-A
            require(math.gcd(N,D)==sum(phis[d] for d in ds if D%d==0),'general totient identity');gates+=1
            x=Fraction(N,D)
            require(x.denominator==abs(D)//math.gcd(N,D),'reduced denominator');gates+=1
        tested+=1
    require(pow(2,9,13)==81*pow(128,-1,13)%13,'phase13')
    require(pow(2,7,17)==81*pow(128,-1,17)%17,'phase17')
    require((9-7)%math.gcd(12,8)!=0,'incompatible CRT phases')
    require(anchors['(1, 1, 5)']['gcd_histogram']=={1:27,5:8,185:1},'anchor185 hist')
    require(anchors['(1, 3, 3)']['gcd_histogram']=={1:19,13:2,17:3},'anchor221 hist')
    # Rational GF equality tested coefficientwise, including a finite all-clock recurrence.
    return {'status':'FINITE-EXACT','gates':gates,'word_cases':tested,'append_range':[1,40], 'multipliers':[3,5,7],'letters':[1,4],'word_lengths':[1,3], 'prime_power_examples':prime_power_examples,'anchors':anchors,'CRT_221':'9 mod 12 and 7 mod 8 incompatible modulo 4','scope':'Algebraic theorem in companion proof; finite checks are not Collatz convergence.'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
