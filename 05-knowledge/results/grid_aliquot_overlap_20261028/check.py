"""Exact controls for divisor shape, additive overlap, and the two grids."""
from itertools import combinations_with_replacement
from math import isqrt, prod
import json
import sys

sys.stdout.reconfigure(encoding='utf-8',newline='\n')

def need(ok,message):
    if not ok: raise RuntimeError(message)

def factor(n):
    out={}; p=2
    while p*p<=n:
        while n%p==0: out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out

def sigma(n):
    return prod((p**(a+1)-1)//(p-1) for p,a in factor(n).items())

def proper(n):
    return [d for d in range(1,n) if n%d==0]

def mul(x,y,m=None):
    a,b=x; c,d=y
    out=(a*c-b*d,a*d+b*c-b*d)
    return tuple(t%m for t in out) if m else out

def fourth(z,m=None):
    zz=mul(z,z,m); return mul(zz,zz,m)

def main():
    N=1_000_000
    # Direct additive sieve and a separate least-prime multiplicative evaluator.
    s=[0]*(N+1); divisor_count=[0]*(N+1)
    for d in range(1,N//2+1):
        for n in range(2*d,N+1,d): s[n]+=d; divisor_count[n]+=1
    lp=[0]*(N+1)
    for p in range(2,N+1):
        if lp[p]: continue
        for n in range(p,N+1,p):
            if not lp[n]: lp[n]=p
    sig=[1]*(N+1); tau=[1]*(N+1); omega=[0]*(N+1); power=[1]*(N+1); exponent=[0]*(N+1)
    for n in range(2,N+1):
        p=lp[n]; rest=n//p
        power[n]=power[rest]*p if rest%p==0 else p
        exponent[n]=exponent[rest]+1 if rest%p==0 else 1
        base=n//power[n]
        sig[n]=sig[base]*(p*power[n]-1)//(p-1)
        tau[n]=tau[base]*(exponent[n]+1)
        omega[n]=omega[base]+1
    need(all(s[n]==sig[n]-n and divisor_count[n]==tau[n]-1 for n in range(1,N+1)),
         'independent additive and multiplicative evaluations')
    excess_two=[n for n in range(1,N+1) if s[n]==n+2]
    two_prime=[n for n in excess_two if omega[n]==2]
    formula=[]
    for a in range(2,30):
        q=2**(a+1)-3; n=2**a*q
        if n>N: break
        if factor(q)=={q:1}: formula.append(n)
    need(two_prime==formula,'complete two-prime excess-two formula in census')
    need(all(next(p for p in factor(n) if p!=2)%8==5 for n in two_prime),
         'near-perfect odd prime parameters lie in the excluded Eisenstein residue')
    need([n for n in range(1,N+1) if tau[n]==6 and s[n]>=n]==[12,18,20,28],
         'all nondeficient five-proper-divisor sources')
    need([n for n in excess_two if tau[n]==6]==[20],'unique five-divisor overlap')
    need(all(s[n]<=n for n in range(2,N+1) if omega[n]==2 and tau[n]==4),
         'distinct-semiprime boundary')
    need(all(isqrt(sig[n])**2!=sig[n] for n in two_prime),'two-prime excess-two family has no square-sigma edge')
    rows=[]
    for n in [4,9,25,8,49,16,6,10,15,14,21,35,22,12,18,20,28,24,30,100,650,1001,343,57]:
        ds=proper(n)
        need(sum(ds)==s[n] and len(ds)==tau[n]-1,'literal divisor list')
        rows.append({'n':n,'factorization':factor(n),'proper_divisors':ds,'s':s[n],
                     'nontrivial_sum':s[n]-1,'excess':s[n]-n,'tau':tau[n],
                     'near_perfect':s[n]>n and n%(s[n]-n)==0})

    reps=sorted({sum(x**4 for x in xs) for xs in combinations_with_replacement(range(3),4)
                 if sum(x**4 for x in xs)<=50})
    budget=[n for n in range(51) if n//16+n%16<=4]
    need(reps==budget and len(reps)==14 and len([n for n in reps if n])==13,'fourth-power window')
    missing=[n for n in range(1,51) if n not in reps]
    need(len(missing)==37 and 37%16==5,'thirteen positives and thirty-seven omissions')
    need([n for n in range(1,81) if n//16+n%16<=4][-1]==64,'finite head boundary')
    need(81//16+81%16==6 and 3**4==81,'hostile outside two-denomination range')

    quartics={fourth((a,b),8) for a in range(8) for b in range(8)}
    need(quartics=={(0,0),(1,0),(0,1),(7,7)},'quartic residues are zero and cubic units')
    reachable={(0,0)}; minima={0:0}; rational={0:[0]}
    for k in range(1,6):
        reachable={((a+c)%8,(b+d)%8) for a,b in reachable for c,d in quartics}
        rational[k]=sorted(a for a,b in reachable if b==0)
        for a in rational[k]: minima.setdefault(a,k)
    need(rational[4]==[0,1,2,3,4,6,7] and rational[5]==list(range(8)),'Eisenstein residue obstruction and boundary')
    need([minima[r] for r in range(8)]==[0,1,2,3,4,5,4,2],'exact modular term costs')
    need(all(minima[next(p for p in factor(n) if p!=2)%8]==5 for n in two_prime),
         'the odd prime factor needs at least five Eisenstein fourth powers')
    need(fourth((1,2))==(9,0) and fourth((2,0))==(16,0),'first grid supplies nine')
    need(tuple(2*x+2*y for x,y in zip(fourth((1,2)),fourth((2,0))))==(50,0),
         '50 has four fourth powers already on the first rectangular grid')
    grid=[]
    for a in range(-5,6):
        for b in range(-5,6):
            u,v=2*a-b,b
            need((u-v)%2==0 and (u*u+3*v*v)==4*(a*a-a*b+b*b),'half-shift norm identity')
            need((u+v)//2==a,'coordinate recovery')
    for n in [13,37]:
        witness=next((x,y) for x in range(isqrt(n)+1) for y in range(isqrt(n//3)+1) if x*x+3*y*y==n)
        grid.append({'norm':n,'rectangular_first_grid':witness,'fourth_power_residue':n%8})

    orbits=[]
    for seed in [1001,1000001]:
        trace=[seed]
        while trace[-1]:
            need(len(trace)<=30,'unexpected orbit length')
            trace.append(sigma(trace[-1])-trace[-1])
        orbits.append(trace)
    bound=342**2
    fiber=[n for n in range(1,bound+1) if s[n]==343]
    triples=[n for n in fiber if omega[n]==3 and tau[n]==8]
    need(len(fiber)==19 and triples==[1001],'complete inherited anchor fiber')
    need(sigma(343)==400 and sigma(1001)==1344 and 1001+1001000==1001**2,'two square-sum predicates')
    need(sigma(1001000)-1001000!=1001,'external square edge is not an aliquot edge in either direction')
    need(orbits[0]==[1001,343,57,23,1,0] and orbits[1][-4:]==[21,11,1,0],'short orbit controls')
    family5=[{'n':n,'s':sigma(n)-n} for n in [5,10,25,50,100]]
    fibers={str(t):[n for n in range(1,(t-1)**2+1) if s[n]==t] for t in [11,12,13,37]}
    need(fibers=={'11':[21],'12':[121],'13':[27,35],'37':[155,203,299,323]},'complete bounded inverse fibers')
    print(json.dumps({'status':'FINITE-EXACT controls of displayed proofs and inherited claims',
        'census_bound':N,'excess_two_all':excess_two,'excess_two_two_primes':two_prime,
        'nondeficient_tau6':[12,18,20,28],'divisor_rows':rows,'fourth_power_values_0_50':reps,
        'fourth_power_omissions_1_50':missing,'Eisenstein_fourth_residues_mod8':sorted(quartics),
        'rational_residues_by_term_budget':rational,'minimum_terms_mod8':minima,
        'Eisenstein_50_witness':{'bases':[[1,2],[1,2],[2,0],[2,0]],'fourth_powers':[[9,0],[9,0],[16,0],[16,0]]},
        'norm_witnesses':grid,'doubling_squaring_family':family5,'inverse_fibers':fibers,
        'short_aliquot_orbits':orbits,'complete_343_fiber':fiber,'three_prime_sources_of_343':triples,
        's_1001000':sigma(1001000)-1001000},indent=2,sort_keys=True))

if __name__=='__main__': main()
