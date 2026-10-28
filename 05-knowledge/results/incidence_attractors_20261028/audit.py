"""Independent additive-divisor and full-context audit; imports no producer code."""
from itertools import combinations
from math import gcd,isqrt,prod
from pathlib import Path
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    bank=(2,3,4,5,7,11,13);N=60060
    sums=[0]*(N+1);divs=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):sums[n]+=d;divs[n]+=1
    def square(n):return isqrt(n)**2==n
    subsets=[tuple(bank[i] for i in range(7) if mask>>i&1) for mask in range(128)]
    legal=[s for s in subsets if all(gcd(a,b)==1 for a,b in combinations(s,2))]
    zero=[s for s in subsets if square(prod(sums[b] for b in s))]
    actual=sorted(prod(s) for s in legal if square(sums[prod(s)]))
    # Full finite residual table, rather than iterative automaton refinement.
    groups={}
    for s in legal:
        n=prod(s);response=[]
        for t in subsets:
            good=all(gcd(a,b)==1 for a,b in combinations(t,2)) and gcd(n,prod(t))==1
            response.append(0 if not good else 1+int(square(sums[n*prod(t)])))
        groups.setdefault(tuple(response),[]).append(n)
    merged=sorted(sorted(g) for g in groups.values() if len(g)>1)
    # Each merged pair has the same outputs under every appended guard,
    # and its successful successors stay in a common residual class.
    owners={n:key for key,ns in groups.items() for n in ns}
    for ns in groups.values():
        for x,y in combinations(ns,2):
            for b in bank:
                need((gcd(x,b)==1)==(gcd(y,b)==1),'guard congruence')
                if gcd(x,b)==1:need(owners[x*b]==owners[y*b],'successor congruence')
    a=json.loads((Path(__file__).parent/'normal.json').read_text(encoding='utf-8'))['arithmetic']
    need(actual==a['square_sources'],'square support list')
    need(len(groups)==a['minimal_guarded_states']==86,'full-context class count')
    need(merged==a['merged_pairs'],'merged class equality')
    small_primes=[p for p in range(2,14) if divs[p]==2]
    def fus(n):
        if n==1:return True
        r=sum(n%p==0 for p in small_primes)
        sqfree=all(n%(p*p)!=0 for p in small_primes)
        return divs[n]-2==(2**r-1-int(sqfree))+(r-int(divs[n]==2))
    additional={}
    for title,observation in [('with_FUS',fus),('with_tau',lambda n:divs[n])]:
        sigs=set()
        for s in legal:
            n=prod(s);response=[]
            for t in subsets:
                good=all(gcd(a,b)==1 for a,b in combinations(t,2)) and gcd(n,prod(t))==1
                response.append(None if not good else (square(sums[n*prod(t)]),observation(n*prod(t))))
            sigs.add(tuple(response))
        additional[title]=len(sigs)
    need(additional==a['additional_observer_class_counts'],'independent observer hierarchy')
    # Direct divisor counts and prime counts identify exponent pattern (2,1,1):
    # tau=12, omega=3. Count primes by trial primality, not a factor recurrence.
    primes=[p for p in range(2,365) if divs[p]==2]
    candidates=[n for n in range(2,365) if divs[n]==12 and sum(n%p==0 for p in primes)==3 and square(sums[n])]
    need(candidates==[364],'complete minimum bound')
    context=[m for m in range(1,206) if gcd(m,52)==1 and square(sums[4*m])!=square(sums[13*m])]
    need(context==[205],'common context bound')
    fano_primes=(2,5,7,19,29,59,89);weights={}
    # Independent direct divisor-pair enumeration, without a parity matrix.
    for mask in range(128):
        ps=[p for i,p in enumerate(fano_primes) if mask>>i&1];n=prod(ps)
        total=sum(d+n//d if d*d!=n else d for d in range(1,isqrt(n)+1) if n%d==0)
        if square(total):weights[len(ps)]=weights.get(len(ps),0)+1
    need(weights=={0:1,3:7,4:7,7:1},'independent Hamming census')
    # Norm-one shell audit by brute coordinates and explicit unique descent.
    roots={(1,1,1),(2,1,2),(1,2,2)};counts={r:0 for r in roots}
    def parent(x):
        a,b,c=x;r=a+2*b-2*c;s=2*a+b-2*c;h=3*c-2*a-2*b
        options=[(r,-s,h),(r,s,h),(-r,s,h)]
        good=[y for y in options if min(y)>0 and y[2]<c]
        need(len(good)==1,'unique positive parent');return good[0]
    for c in range(1,151):
        for a0 in range(1,c+1):
            b2=c*c+1-a0*a0;b=isqrt(b2)
            if b>0 and b*b==b2:
                x=(a0,b,c)
                while x not in roots:x=parent(x)
                counts[x]+=1
    return {'sigma_method':'add every divisor into its multiples','sigma_universe':[1,N],
            'subset_universe':128,'legal_states':len(legal),'formal_zero_states':len(zero),
            'additional_observer_class_counts':additional,
            'full_residual_table_entries':len(legal)*len(subsets),'minimal_classes':len(groups),
            'merged_pairs':merged,'square_sources':actual,'least_p2qr_square_sigma':364,
            'least_common_context':205,'norm_one_shell_c_bound':150,
            'norm_one_shell_root_counts':{str(k):v for k,v in sorted(counts.items())},
            'Fano_direct_divisor_subset_count':128,'Fano_square_weight_counts':weights}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
