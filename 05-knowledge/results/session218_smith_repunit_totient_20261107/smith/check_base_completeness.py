#!/usr/bin/env python3
"""Independent exact finite verification of the moving-base classification."""
import json
LIMIT=100000
spf=list(range(LIMIT+1))
for p in range(2,int(LIMIT**0.5)+1):
    if spf[p]==p:
        for n in range(p*p,LIMIT+1,p):
            if spf[n]==n:spf[n]=p

def factors(n):
    out=[]
    while n>1:out.append(spf[n]);n//=spf[n]
    return out

def digitsum(n,b):
    total=0
    while n:n,r=divmod(n,b);total+=r
    return total

def smith(n,b):
    ff=factors(n)
    return len(ff)>=2 and digitsum(n,b)==sum(digitsum(p,b) for p in ff)

def divisors(n):
    assert n>0
    out=set()
    d=1
    while d*d<=n:
        if n%d==0:out.update((d,n//d))
        d+=1
    return sorted(out)

exceptions={6,8,9,12}
construction_counts={'four':0,'twice_prime':0,'factor_sum_base':0}
for n in range(4,LIMIT+1):
    ff=factors(n)
    if len(ff)<2 or n in exceptions:continue
    F=sum(ff)
    if n==4:b=5;tag='four'
    elif n%2==0 and spf[n//2]==n//2:b=n//2-1;tag='twice_prime'
    else:
        b=n-F+1;tag='factor_sum_base'
        assert n>12 and 2*F<n+2 and 2*b>n and max(ff)<b
        assert divmod(n,b)==(1,F-1)
    assert smith(n,b),(n,b,ff)
    construction_counts[tag]+=1

failed=[];exception_details=[]
for n in range(6,5001):
    ff=factors(n)
    if len(ff)<2:continue
    D=n-sum(ff)
    bases=[d+1 for d in divisors(D)]
    hits=[b for b in bases if smith(n,b)]
    if not hits:failed.append(n)
    if n in exceptions:
        exception_details.append({'n':n,'F':sum(ff),'possible_bases':bases,'digit_and_factor_sums':[[digitsum(n,b),sum(digitsum(p,b) for p in ff)] for b in bases]})
assert failed==sorted(exceptions)
for n in range(6,301):
    ff=factors(n)
    if len(ff)<2:continue
    direct=[b for b in range(2,n+2) if smith(n,b)]
    filtered=[d+1 for d in divisors(n-sum(ff)) if smith(n,d+1)]
    assert direct==filtered
assert all(smith(4,b) for b in range(5,101))
assert not any(smith(4,b) for b in range(2,5))
print(json.dumps({'status':'PASS','construction_limit':LIMIT,'construction_counts':construction_counts,'exhaustive_divisor_sieve_limit':5000,'independent_all_bases_limit':300,'exceptions':exception_details,'claims':'Finite controls only; all-integer result has an elementary proof.'},indent=2))
