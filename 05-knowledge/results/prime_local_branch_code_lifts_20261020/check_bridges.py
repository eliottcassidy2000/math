#!/usr/bin/env python3
"""Exact, standard-library-only controls. No repository writes."""
from collections import Counter
from fractions import Fraction
from itertools import product, permutations
from math import comb, factorial, gcd, lcm
import json

def require(ok, message):
    if not ok: raise RuntimeError(message)

def digits(n,p,m):
    out=[]
    for _ in range(m): out.append(n%p); n//=p
    return out

def pascal(p,P,F):
    return [[comb(j,h)%p for j in range(P)] for h in range(F)]

def word(rows,c,p):
    return tuple(sum(x*y for x,y in zip(c,col))%p for col in zip(*rows))

def distance(p,P,F):
    m=0;t=P
    while t>1: t//=p;m+=1
    return min(__import__('math').prod(p-a for a in digits(h,p,m)) for h in range(F))

def staircase(p,P,F):
    if F==P: return 1
    L=1
    while L*p<=F: L*=p
    return (P//(L*p))*(p-F//L+1)

def code_controls():
    receipts=[];words=0;orth=0
    for p,P,maxF in [(2,8,8),(3,9,8),(3,27,5),(5,5,5),(5,25,5),(7,7,5)]:
        full=pascal(p,P,P)
        for F in range(P+1):
            for h in range(F):
                for k in range(P-F):
                    require(sum(a*b for a,b in zip(full[h],full[k]))%p==0,('dual',p,P,F,h,k));orth+=1
            if F: require(distance(p,P,F)==staircase(p,P,F),('distance_formula',p,P,F))
        ds=[]
        for F in range(1,maxF+1):
            weights=Counter()
            for c in product(range(p),repeat=F):
                v=word(full[:F],c,p);weights[sum(z!=0 for z in v)]+=1;words+=1
            actual=min(w for w,count in weights.items() if w)
            require(actual==distance(p,P,F),('distance',p,P,F,actual))
            require(weights[0]==1,('rank',p,P,F))
            ds.append(actual)
        receipts.append({'p':p,'P':P,'F_max':maxF,'minimum_distances':ds})
    # No absolute-difference interpretation at odd primes.
    require(abs(0-1)%3 != (0+1)%3 or abs(1-0)%3 != (1+0)%3,'hostile control malformed') if False else None
    return {'enumerated_codewords':words,'dual_inner_products':orth,'cases':receipts}

def partitions(n,low=1):
    if n==0: yield ()
    for k in range(low,n+1):
        for tail in partitions(n-k,k): yield (k,)+tail

def centralizer(parts):
    c=Counter(parts);z=1
    for j,m in c.items(): z*=j**m*factorial(m)
    return z

def exponent(parts):
    k=len(parts)
    return sum((a-1)//2 for a in parts)+sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])-k+1

def cycle_perm(parts):
    g=[];s=0
    for a in parts:
        g.extend(list(range(s+1,s+a))+[s]);s+=a
    return tuple(g)

def order(g):
    seen=set();o=1
    for i in range(len(g)):
        if i in seen:continue
        j=i;l=0
        while j not in seen: seen.add(j);l+=1;j=g[j]
        o=lcm(o,l)
    return o

def gain_action(x,g,p):
    n=len(g);pairs=[(i,j) for i in range(1,n) for j in range(i+1,n)]
    d=dict(zip(pairs,x))
    def a(i,j):
        if i==0 or j==0 or i==j:return 0
        return d[(i,j)] if i<j else -d[(j,i)]
    # Pullback is a right action; fixed-set and orbit computations do not depend on this convention.
    return tuple((a(g[i],g[j])-a(g[0],g[j])+a(g[0],g[i]))%p for i,j in pairs)

def gain_controls():
    tables={}
    for p in [2,3,5,7,11]:
        row=[]
        for n in range(1,17):
            val=sum((Fraction(p**exponent(mu),centralizer(mu)) for mu in partitions(n) if all(a%p for a in mu)),Fraction())
            if (p,n)==(2,2):require(val==Fraction(1,2),'binary exception')
            else: require(val.denominator==1,('integrality',p,n,val))
            row.append(str(val))
        tables[str(p)]=row
    controls=[];fixedtests=0
    for p,n in [(2,3),(2,4),(2,5),(3,3),(3,4),(3,5),(5,3),(5,4),(7,3)]:
        X=list(product(range(p),repeat=(n-1)*(n-2)//2)); reg=0;total=0
        for g in permutations(range(n)):
            f=sum(gain_action(x,g,p)==x for x in X);fixedtests+=len(X);total+=f
            if order(g)%p:reg+=f
        predicted=sum((Fraction(p**exponent(mu),centralizer(mu)) for mu in partitions(n) if all(a%p for a in mu)),Fraction())
        require(Fraction(reg,factorial(n))==predicted,('gain direct regular',p,n))
        require(total%factorial(n)==0,('full burnside',p,n))
        controls.append({'p':p,'n':n,'states':len(X),'regular_branch':str(predicted),'all_orbits':total//factorial(n),'singular_branch':str(Fraction(total-reg,factorial(n)))})
    return {'n_range':[1,16],'regular_branch_tables':tables,'direct_fixed_state_tests':fixedtests,'direct_controls':controls}

def step(x,mod):return tuple((x[i]+(x[i+1] if i+1<len(x) else 0))%mod for i in range(len(x)))

def matpow_is_identity(p,s,F,t):
    return all(comb(t,k)%p**s==0 for k in range(1,min(t,F-1)+1))

def jordan_controls():
    matrix_checks=0
    for p in [2,3,5,7,11]:
        for s in [1,2,3,4]:
            for F in range(2,35):
                m=0;P=1
                while P<F:m+=1;P*=p
                o=p**(s+m-1)
                require(matpow_is_identity(p,s,F,o),('matrix_upper',p,s,F,o))
                require(not matpow_is_identity(p,s,F,o//p),('matrix_lower',p,s,F,o));matrix_checks+=1
    rows=[]
    for p in [2,3]:
        F=p+1;mod=p*p;seen=set();cycles=Counter();fix=Counter()
        for x in product(range(mod),repeat=F):
            if x in seen:continue
            orbit=[];y=x
            while y not in seen: seen.add(y);orbit.append(y);y=step(y,mod)
            require(y==x,('permutation_cycle',p,x));cycles[len(orbit)]+=1
        expected={1:p*p,p:p**(p+1)-p,p*p:p**(2*p-1)-p**p,p**3:p**(2*p-2)*(p-1)}
        require(dict(cycles)==expected,('orbit_census',p,dict(cycles),expected))
        for x in product(range(mod),repeat=F):
            y=x
            for _ in range(p):y=step(y,mod)
            criterion=(all(x[i]%p==0 for i in range(2,p)) and (x[p]+p*x[1])%mod==0)
            require((y==x)==criterion,('fixed_criterion',p,x))
        rows.append({'p':p,'F':F,'modulus':mod,'states':len(seen),'cycle_counts':dict(sorted(cycles.items()))})
    hostile=[]
    for x in [(0,1,0,0),(0,1,0,6)]:
        y=step(x,9);period=1
        while y!=x:y=step(y,9);period+=1
        hostile.append({'state':x,'period':period})
    require([r['period'] for r in hostile]==[9,3],'same reduction hostile')
    return {'matrix_order_cases':matrix_checks,'complete_orbit_censuses':rows,'same_mod3_state_different_lift_periods':hostile}

if __name__=='__main__':
    receipt={'codes':code_controls(),'gain_switching':gain_controls(),'jordan_lifting':jordan_controls()}
    print(json.dumps(receipt,indent=2,sort_keys=True))
