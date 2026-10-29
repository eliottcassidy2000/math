#!/usr/bin/env python3
"""Independent stdlib audit. All gates remain enabled with python -O."""
from fractions import Fraction
from itertools import product
from math import gcd, lcm
from functools import lru_cache
import json

checks=0
def check(ok, tag):
    global checks
    checks+=1
    if not ok: raise RuntimeError(tag)

def coef(w):
    # Direct affine composition, rather than the numerator recurrence.
    slope, shift=Fraction(1), Fraction(0)
    for v in w:
        slope, shift=3*slope/2**v,(3*shift+1)/2**v
    x=shift/(1-slope)
    H=sum(w); D=2**H-3**len(w)
    return int(shift*2**H), D, x

@lru_cache(None)
def divisors(n):
    lo=[d for d in range(1, int(n**.5)+1) if n%d==0]
    return sorted(set(lo+[n//d for d in lo]))
@lru_cache(None)
def phi(n):
    ans=n; t=n; p=2
    while p*p<=t:
        if t%p==0:
            ans-=ans//p
            while t%p==0:t//=p
        p+=1
    if t>1: ans-=ans//t
    return ans
@lru_cache(None)
def powers(d):
    if d==1:return {0:0}
    out={}; t=1
    while t not in out:
        out[t]=len(out);t=t*2%d
    return out

def channels(w):
    C,D,x=coef(w);H=sum(w); k=len(w);N=3*C+2**H
    ans=[]
    for d in divisors(N):
        tab=powers(d);r=(pow(3,k+1,d)*pow(pow(2,H,d),-1,d))%d if d>1 else 0
        if r in tab: ans.append((d,len(tab),tab[r]))
    return N,ans

summary={'universe':'all nonempty words length <=4, letters in 1..4, sum<=12; all append letters 1..80','words':0}
prime_power=None
for k in range(1,5):
 for w in product(range(1,5),repeat=k):
    if sum(w)>12:continue
    summary['words']+=1
    C,D,x=coef(w);H=sum(w)
    check(C%2==1 and D%2==1,'odd carrier')
    phase=x
    for i,v in enumerate(w):
        check(phase.denominator==x.denominator,'orbit denominator invariant')
        phase=(3*phase+1)/2**v
        check(phase.numerator%2==1 and phase.denominator%2==1,'exact valuation')
    check(phase==x,'closure')
    N,chs=channels(w)
    check(gcd(N,6)==1,'N unit conditions')
    T=len(powers(N))
    for v in range(1,81):
        delta=2**(H+v)-3**(k+1);g=gcd(N,delta)
        alt=sum(phi(d) for d,t,r in chs if v%t==r)
        check(g==alt,'totient channels')
        check(g==gcd(N,2**(H+v+T)-3**(k+1)),'period')
        cc,dd,xx=coef(w+(v,))
        check((cc,dd,xx)==(N,delta,Fraction(N,delta)),'append affine closure')
        check(xx.denominator==abs(delta)//g,'signed reduced denominator')
        # Coefficient extraction from z^first/(1-z^period), v>=1.
        gf=sum(phi(d) for d,t,r in chs if v>=(r or t) and (v-(r or t))%t==0)
        check(gf==g,'GF coefficients')
        if prime_power is None:
            for d in divisors(g):
                if any(d==p*p for p in range(2,int(d**.5)+1)):
                    prime_power={'word':w,'append':v,'N':N,'delta':delta,'gcd':g,'square_divisor':d};break
    # CRT equivalence tested pairwise at a complete joint period.
    for (d,t,r),(e,u,s) in product(chs,repeat=2):
        compatible=(r-s)%gcd(t,u)==0
        # Modular inverse implementation independently tests actual existence.
        actual=any((r+j*t)%u==s for j in range(u//gcd(t,u)))
        check(compatible==actual,'CRT compatibility')

for w,expected in [((1,1,5),[(1,1,0),(5,4,1),(37,36,25),(185,36,25)]),((1,3,3),[(1,1,0),(13,12,9),(17,8,7)])]:
 N,chs=channels(w);check(chs==expected,'anchor channels')
 C,D,x=coef(w)
 check(x.denominator==101,'anchor denominator')
 check(((3*x+1)/2).denominator==101,'true anchor step')
summary['prime_power_control']=prime_power
summary['anchor_115_channels']=channels((1,1,5))
summary['anchor_133_channels']=channels((1,3,3))
summary['checks']=checks
print(json.dumps(summary,indent=2,sort_keys=True))
