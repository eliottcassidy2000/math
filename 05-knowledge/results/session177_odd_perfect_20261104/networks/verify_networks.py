#!/usr/bin/env python3
"""Exact stdlib-only replay. This is not a universal odd-perfect search."""
from math import gcd, prod, isqrt
from itertools import combinations
from fractions import Fraction
import json
from pathlib import Path

def require(condition):
    if not condition:
        raise AssertionError("Exact certificate check failed")

def factor(n):
    out={}; p=2
    while p*p<=n:
        while n%p==0: out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def sigma_factor(f):return prod(sum(p**j for j in range(e+1)) for p,e in f.items())
def sigma(n):return sigma_factor(factor(n))
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def sc(f):
    edges={p:set(factor(sigma_factor({p:e}))) & set(f) for p,e in f.items()}
    for p in f:
        reached={p};todo=[p]
        while todo:
            for q in edges[todo.pop()]-reached:reached.add(q);todo.append(q)
        if reached!=set(f):return False
    return True

def branch(q):
    # Hypothesis: q||N is the Euler component; every newly forced prime
    # has exact exponent 2. Exiting this grammar gives a fourth-power clause.
    f={q:1}; ordinary=set(); steps=[]
    while True:
        A=prod(p**e for p,e in f.items()); s=sigma_factor(f)
        R=s//gcd(s,2*A); rf=factor(R)
        row={'factors':dict(f),'A':A,'sigma_A':s,'residual':R}
        steps.append(row)
        if s>2*A:row['stop']='abundancy exceeds two';break
        if gcd(R,A)>1:row['stop']='fixed exponent capacity exceeded';break
        high=[p for p,e in rf.items() if e>2]
        if high:ordinary.update(high);row['stop']='new prime forced to exponent at least four';break
        require(rf and 2 not in rf and (q not in rf))
        ordinary.update(rf)
        for p in rf:f[p]=2
        require(len(steps) < 20)
    return {'Euler_prime':q,'Euler_exponent':1,'fourth_power_disjunction':sorted(ordinary),'steps':steps}

def main():
    # Independent divisor-enumeration arithmetic oracle.
    for n in range(1,501):require(sigma(n) == sum(divisors(n)))
    cut_checks=0
    for n in [6,28,496,8128]:
        f=factor(n);require(sigma(n) == 2 * n and sc(f))
        keys=list(f)
        for mask in range(1,2**len(keys)-1):
            sf={p:f[p] for i,p in enumerate(keys) if mask>>i&1}
            A=prod(p**e for p,e in sf.items());B=n//A
            R=sigma(A)//gcd(sigma(A),2*A)
            require(B % R == 0 and gcd(A, R) == 1)
            outgoing=prod(r**v for p,e in f.items() if p in sf
                          for r,v in factor(sigma_factor({p:e})).items() if r not in sf)
            incoming=prod(r**v for p,e in f.items() if p not in sf
                          for r,v in factor(sigma_factor({p:e})).items() if r in sf)
            # Even-perfect analogue: the extra incoming 2 is at vertex 2.
            epsilon=int(2 in sf)
            require(Fraction(outgoing, incoming) == Fraction(sigma(A), 2 ** epsilon * A))
            require(gcd(outgoing, incoming) == 1)
            require(Fraction(sigma(A), A) * Fraction(sigma(R), R) <= 2)
            cut_checks+=1
    require(sc(factor(120)) and sigma(120) == 3 * 120)  # SC is not sufficient.
    require(factor(sigma(17)) == {2: 1, 3: 2})
    require(factor(sigma(3 ** 2)) == {13: 1})
    require(factor(sigma(13 ** 2)) == {3: 1, 61: 1})
    require(factor(sigma(61 ** 2)) == {3: 1, 13: 1, 97: 1})
    require(sigma(13 ** 4) % 3 != 0) # changed exponent breaks the earlier edge.
    require(factor(sigma(3 ** 4)) == {11: 2}) # do not retain the old 3^2 -> 13 edge.
    qs=[q for q in range(5,98,4) if factor(q)=={q:1}]
    clauses=[branch(q) for q in qs]
    require(next((c for c in clauses if c['Euler_prime'] == 17))['fourth_power_disjunction'] == [3, 13])
    require(next((c for c in clauses if c['Euler_prime'] == 5))['fourth_power_disjunction'] == [3, 13, 61])
    out={'status':'FINITE-EXACT; not an all-exponent or all-prime exclusion',
         'independent_sigma_checks':500,'perfect_unitary_cut_checks':cut_checks,
         'strong_connectivity_hostile':120,'clauses':clauses}
    path=Path(__file__).with_name('network_certificates.json');path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'sigma_checks':500,'cut_checks':cut_checks,'clauses':len(clauses),'all_assertions_passed':True}))
if __name__=='__main__':main()
