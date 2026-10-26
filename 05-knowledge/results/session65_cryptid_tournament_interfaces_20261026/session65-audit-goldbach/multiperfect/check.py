#!/usr/bin/env python3
"""Independent divisor-sieve and exact family audit; no audited code execution.
Usage: python3 -B check.py [path/to/session65-root/multiperfect]
"""
from fractions import Fraction
from math import comb, isqrt
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent.parent/'session65-root/multiperfect'

def need(ok,why):
    if not ok:
        raise ValueError(why)

def valuation(n):
    need(n>0,'positive valuation argument')
    e=0
    while n%2==0:
        n//=2;e+=1
    return e

def factors(n):
    out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            power=0
            while n%p==0:
                n//=p;power+=1
            out.append((p,power))
        p=3 if p==2 else p+2
    if n>1:
        out.append((n,1))
    return out

def divisors_from_factors(n):
    values=[1]
    for p,e in factors(n):
        values=[d*p**a for d in values for a in range(e+1)]
    need(len(values)==len(set(values)) and all(n%d==0 for d in values),'distinct actual divisors')
    return values

def paired_divisor_sum(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            total+=d
            if d!=n//d:
                total+=n//d
    return total

def main():
    sums=(SOURCE/'SHA256SUMS').read_text().splitlines()
    for line in sums:
        h,name=line.split('  ',1)
        p=SOURCE/name
        need(p.resolve().is_relative_to(SOURCE.resolve()) and not p.is_symlink(),'safe input')
        need(hashlib.sha256(p.read_bytes()).hexdigest()==h,'sealed input hash: '+name)
    need((SOURCE/'normal.json').read_bytes()==(SOURCE/'optimized.json').read_bytes(),'source normal/-O equality')
    N=100000
    sigma=[0]*(N+1)
    for d in range(1,N+1):
        for multiple in range(d,N+1,d):
            sigma[multiple]+=d
    for n in range(1,10001):
        need(sigma[n]==paired_divisor_sum(n),'independent sieve vs paired divisors')
    violations=[]
    for n in range(1,N+1,2):
        fac=factors(n)
        odd=[(p,e) for p,e in fac if e%2]
        cost=sum(valuation(p+1)+valuation(e+1)-1 for p,e in odd)
        need(cost==valuation(sigma[n]),'odd valuation ledger')
        need((sigma[n]%2==1)==(isqrt(n)**2==n),'odd sigma parity iff square')
        need(len(odd)<=cost,'every odd exponent costs at least one')
        if cost==1:
            need(len(odd)==1 and odd[0][0]%4==odd[0][1]%4==1,'Euler budget-one form')
        if sigma[n]%n==0 and n>1:
            violations.append([n,sigma[n]//n])
    # Independently enumerate all divisors of both M and D, not just substitute identity.
    family=[]
    for r in range(2,25):
        M=2**r-1
        D=comb(2**r,2)
        sm=sum(divisors_from_factors(M))
        sd=sum(divisors_from_factors(D))
        need(D==2**(r-1)*M and sd==M*sm,'arc-family exact identity')
        ratio=Fraction(sd,D)
        need(ratio==Fraction(sm,2**(r-1)),'abundancy transport')
        need((ratio.denominator==1)==(valuation(sm)>=r-1),'necessary AND sufficient in family')
        prime=factors(M)==[(M,1)]
        need((ratio==2)==prime,'perfect iff prime apex')
        if ratio.denominator==1:
            need(ratio>=2,'multiperfect multiplier at least two')
        family.append({'r':r,'D':D,'M':M,'sigma_D':sd,'sigma_M':sm,'abundancy':str(ratio),'v2_sigma_M':valuation(sm),'multiperfect':ratio.denominator==1})
    baseline=json.loads((SOURCE/'normal.json').read_text())['family']
    for own,old in zip(family,baseline):
        need((own['r'],own['D'],own['M'],own['sigma_M'],own['abundancy'],own['v2_sigma_M'],own['multiperfect'])==(old['exponent'],old['arc_dimension'],old['apex'],old['sigma_apex'],old['abundancy'],old['binary_valuation'],old['multiperfect']),'all baseline family rows independently matched')
    need(len(family)==len(baseline)==23,'family census extent')
    need(sigma[120]==360 and sigma[21]==32 and 32%21!=0,'scope hostiles')
    # Euler form is not sufficient even when the entire binary budget matches k=2.
    need(sigma[5]==6 and valuation(sigma[5])==1 and 5%4==1 and sigma[5]!=2*5,'Euler-form nonsufficiency')
    # A square need not be multiperfect for any odd k>=3.
    need(sigma[9]==13 and isqrt(9)**2==9 and 13%9!=0,'odd-k square nonsufficiency')
    perfect=[n for n in range(2,N+1) if sigma[n]==2*n]
    multiperfect=[[n,sigma[n]//n] for n in range(2,N+1) if sigma[n]%n==0]
    need(perfect==[6,28,496,8128],'small perfect census control')
    need(multiperfect==[[6,2],[28,2],[120,3],[496,2],[672,3],[8128,2],[30240,4],[32760,4]],'small multiperfect census control')
    need(not violations,'bounded odd census only')
    print(json.dumps({'status':'PASS FINITE-EXACT; no odd nonexistence inference','sealed_hashes':len(sums),'sieve_limit':N,'paired_divisor_checks':10000,'odd_ledger_checks':N//2,'family':family,'multiperfect_census_through_100000':multiperfect,'odd_witnesses_through_99999':violations,'hostiles':{'120':{'sigma':360,'k':3},'21':{'sigma':32,'does_not_divide_sigma':True},'5':{'sigma':6,'Euler_form_not_sufficient':True},'9':{'sigma':13,'square_not_sufficient':True}}},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
