#!/usr/bin/env python3
"""Independent exact checks; no assertions, floating arithmetic, or external modules."""
from fractions import Fraction
from itertools import product
from math import factorial, gcd
import json

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def cycles(p):
    unseen=set(range(len(p))); result=[]
    while unseen:
        i=min(unseen); orbit=[]
        while i in unseen:
            unseen.remove(i); orbit.append(i); i=p[i]
        result.append(orbit)
    return result

def group_order(gens):
    identity=tuple(range(len(gens[0]))); seen={identity}; stack=[identity]
    while stack:
        q=stack.pop()
        for g in gens:
            r=compose(g,q)
            if r not in seen: seen.add(r); stack.append(r)
    return len(seen)

primes=(2,3,5,7,11,13); exponents=(5,7,3,4,1,2)
N=1
for p,e in zip(primes,exponents): N*=p**e
check(N==39046339332000,'input integer')
divisors=[]; heights={}; brute={}
for powers in product(*(range(e+1) for e in exponents)):
    b=1
    for p,e in zip(primes,powers): b*=p**e
    if b==1: continue
    h=min(e//f for e,f in zip(exponents,powers) if f)
    k=0; power=b
    while N%power==0: k+=1; power*=b
    check(h==k,('height disagreement',b,h,k))
    divisors.append(b); heights[b]=h; brute[h]=brute.get(h,0)+1
check(len(divisors)==5759,'divisor count')
check(brute=={1:5616,2:120,3:16,4:4,5:2,7:1},'height distribution')
cycle=[4,2,5,3,7]
check([heights[b] for b in cycle]==cycle[1:]+cycle[:1],'literal cycle')
check(heights[13]==heights[4]==2,'fold compatibility')
cyclebasin=0
for b in divisors:
    seen=set(); x=b
    while x!=1 and x not in seen:
        seen.add(x); x=heights[x]
    check((x!=1)==(heights[b]>=2),('basin classification',b))
    cyclebasin+=x!=1
check(cyclebasin==143,'basin count')
# Exact divisor tail counts.
tails={}
for k in range(1,9):
    t=1
    for e in exponents: t*=e//k+1
    tails[k]=t-1
    check(t-1==sum(h>=k for h in heights.values()),('tail',k))
# Mod9 and its only nontrivial quotient mod3 both fail the prime map.
check(11%9==2%9 and 1%9!=5%9,'mod9 hostile')
check(2%3==5%3 and 5%3!=3%3,'mod3 hostile')
# Two-block monodromy permutations. Joining transposition gives one full cycle.
monodromy=[]
for a in range(2,8):
    for b in range(2,8):
        n=a+b; sigma=list(range(n)); tau=list(range(n))
        for i in range(a): sigma[i]=(i+1)%a
        for i in range(a,n): sigma[i]=a+(i-a+1)%b
        tau[0],tau[a]=a,0
        sigma=tuple(sigma); tau=tuple(tau); c=compose(sigma,tau)
        check(sorted(map(len,cycles(sigma)))==sorted([a,b]),'sigma passport')
        check(len(cycles(c))==1,'joining full cycle')
        check(len(cycles(sigma))+len(cycles(tau))+len(cycles(c))==n+2,'genus zero')
        d=gcd(a,b); predicted=factorial(n//d)**d*d
        # Independent BFS finite group for degree <=8 only.
        actual=group_order([sigma,tau]) if n<=8 else None
        if actual is not None: check(actual==predicted,('group order',a,b,actual,predicted))
        # Polynomial critical value by exact substitution.
        norm=Fraction(n**n,a**a*b**b); q=Fraction(a,n)
        check(norm*q**a*(1-q)**b==1,'critical normalization')
        monodromy.append({'a':a,'b':b,'d':d,'predicted_group_order':predicted,'enumerated_order':actual})
# The specific cover is hyperbolic as a triangle orbifold despite sphere domain.
orbifold=Fraction(1,10)+Fraction(1,2)+Fraction(1,7)-1
check(orbifold==Fraction(-9,35),'orbifold characteristic')
check(7*orbifold==2-Fraction(1,2)-Fraction(4,5)-5*Fraction(1,2),'orbifold pullback')
print(json.dumps({'status':'PASS','N':N,'nonunit_divisors':len(divisors),'height_distribution':brute,'tail_counts':tails,'cycle':cycle,'cycle_basin':cyclebasin,'terminal_basin':len(divisors)-cyclebasin,'belyi_checks':monodromy,'specific_triangle_orbifold_characteristic':str(orbifold)},sort_keys=True,indent=2))
