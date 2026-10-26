#!/usr/bin/env python3
from fractions import Fraction
from pathlib import Path
import json
from aliquot_fibers import factors, sigma_factorization, aliquot, spf_sieve, goldbach_witness

def phi(n):
    z=n
    for p in factors(n): z=z//p*(p-1)
    return z

def h(d): return max(Fraction(0),Fraction(1,2*d)-Fraction(1,aliquot(2*d)))

def direct_g(M):
    return sum((Fraction(phi(M//d),M//d)*h(d) for d in range(1,M+1) if M%d==0),Fraction(0))

cases=0; strict=0
for M in range(1,101):
    for q in [3,5,7,11,13]:
        if M%q==0: continue
        gain=direct_g(M*q)-direct_g(M)
        rhs=sum((Fraction(phi(M//d),M//d)*(h(q*d)-h(d)/q) for d in range(1,M+1) if M%d==0),Fraction(0))
        if gain!=rhs or gain<0: raise RuntimeError('prime-addition identity/positivity failed')
        if M%6==0 and gain<=0: raise RuntimeError('strict gain failed')
        for d in range(1,M+1):
            if M%d: continue
            S=aliquot(2*d); T=S+2*d
            if S>2*d and h(q*d)-h(d)/q != Fraction(T,q*S*(q*S+T)):
                raise RuntimeError('active-term exact gain failed')
        strict+=gain>0; cases+=1

spf=spf_sieve(19562)
boundary_primes=[p for p in range(5,9782) if spf[p]==p]
for p in boundary_primes:
    n=goldbach_witness(2*p+1,spf)
    if n is None or aliquot(n)!=2*p+1: raise RuntimeError('Booker distinctness boundary failure')
out={'prime_addition_identity':{'M':[1,100],'new_primes':[3,5,7,11,13],
      'cases':cases,'strict_gains':strict,'arithmetic':'exact rational, separate direct divisor enumeration'},
 'Booker_distinctness_boundary':{'prime_p_range':[5,9781],'checked_primes':len(boundary_primes),
      'failures':0,'scope':'Only small supplemental check; large4e18 Goldbach computation cited, not rerun'}}
p=Path(__file__).resolve().parents[1]/'results'/'density_monotonicity.json'
p.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
