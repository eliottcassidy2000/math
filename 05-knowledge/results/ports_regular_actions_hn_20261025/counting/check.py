"""A smaller fixed-target response graph for unrestricted tournament counts."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import factorial, gcd, isqrt, prod
import sys


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pow2(e):
    return Fraction(1 << e) if e >= 0 else Fraction(1, 1 << (-e))


def count_target(n, gauge=True):
    if type(n) is not int or n < 0:
        raise ValueError('target must be a nonnegative integer')
    if n == 0:
        return Fraction(1), {'states': 1, 'transitions': 0, 'rows': []}
    top=n if n%2 else n-1
    states={(0,(0,)*((top+1)//2)): Fraction(1)}
    rows=[]
    for ell in range(top,0,-2):
        new=defaultdict(Fraction)
        transitions=0
        for (used,response), weight in states.items():
            remaining=n-used
            choices=[remaining] if ell==1 else range(remaining//ell+1)
            for k in choices:
                rest=remaining-ell*k
                selected=response[(ell-1)//2] if ell<=remaining else 0
                exponent=k*(ell-1)//2+ell*k*(k-1)//2
                exponent += k*(rest-selected) if gauge else k*selected
                transfer=pow2(exponent)/((ell**k)*factorial(k))
                future=range(1,min(ell-2,rest)+1,2)
                next_response=tuple(response[(d-1)//2]+k*((d-gcd(ell,d)) if gauge else gcd(ell,d)) for d in future)
                if gauge:
                    need(all(0<=x<=n for x in next_response), 'defect packing bound')
                new[(n-rest,next_response)]+=weight*transfer
                transitions+=1
        rows.append({'length':ell,'input_states':len(states),'output_states':len(new),'transitions':transitions})
        states=new
    need(set(states)=={(n,())}, 'all terminal paths fill target')
    value=states[(n,())]
    need(value.denominator==1, 'unrestricted count is integral')
    return value, {'states':sum(r['input_states'] for r in rows),
                   'transitions':sum(r['transitions'] for r in rows),
                   'peak_states':max(r['input_states'] for r in rows),'rows':rows}


def partitions(n, low=1):
    if n==0:
        yield ()
    for first in range(low,n+1,2):
        for tail in partitions(n-first,first):
            yield (first,)+tail


def partition_count(n):
    total=Fraction(0)
    for p in partitions(n):
        exponent=sum((a-1)//2 for a in p)+sum(gcd(a,b) for a,b in combinations(p,2))
        denominator=prod(a**k*factorial(k) for a,k in Counter(p).items())
        total+=Fraction(1<<exponent,denominator)
    return total


def literal_edge_orbits(n):
    edges=list(combinations(range(n),2))
    total=0
    for perm in permutations(range(n)):
        unseen=set(edges); valid=True; cycles=0
        while unseen:
            edge=min(unseen); cur=edge; parity=0
            while True:
                unseen.remove(cur)
                a,b=perm[cur[0]],perm[cur[1]]
                parity^=a>b
                cur=tuple(sorted((a,b)))
                if cur==edge:
                    break
            if parity:
                valid=False;break
            cycles+=1
        if valid:
            total+=1<<cycles
    need(total%factorial(n)==0,'literal Burnside divisibility')
    return total//factorial(n)


def response(p, d):
    return sum(gcd(a,d) for a in p)


def defect(p, d):
    return d*len(p)-response(p,d)


def prefix_weight(p):
    e=sum((a-1)//2 for a in p)+sum(gcd(a,b) for a,b in combinations(p,2))
    z=prod(a**k*factorial(k) for a,k in Counter(p).items())
    return Fraction(1<<e,z)


def controls():
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    for n in range(28):
        a,sa=count_target(n,True); b,sb=count_target(n,False)
        need(a==b==partition_count(n), ('independent partition agreement',n))
        need(sa['states']<=sb['states'],('state quotient',n))
        if n<len(known):
            need(a==known[n],('known prefix',n))
        if n<=7:
            need(a==literal_edge_orbits(n),('literal unordered-edge action',n))
    # A strict, realizable same-degree merger: processed lengths >3.
    a=(27,);b=(9,9,9)
    need(sum(a)==sum(b)==27,'strict witness degree')
    need([response(a,d) for d in (1,3)]!=[response(b,d) for d in (1,3)],'old distinguishes')
    need([defect(a,d) for d in (1,3)]==[defect(b,d) for d in (1,3)]==[0,0],'new merges')
    c=(23,);d=(5,9,9)
    need(sum(c)==sum(d)==23,'smaller strict witness degree')
    need([response(c,j) for j in (1,3)]==[1,1] and [response(d,j) for j in (1,3)]==[3,7], 'smaller old profiles')
    need([defect(c,j) for j in (1,3)]==[defect(d,j) for j in (1,3)]==[0,2],'smaller projective merger')
    # The fixed ACTUAL future degree is essential. Core substitution excess breaks it.
    need(pow2(1*(1-4)) != pow2(3*(1-4)), 'one actual point plus three excess points hostile')
    wheel_a=(30013,30047);wheel_b=(30029,30031)
    need(sum(wheel_a)==sum(wheel_b)==60060,'wheel degree sidecar')
    need(all(response(wheel_a,j)==response(wheel_b,j)==2 for j in range(1,59,2)), 'wheel horizon through57')
    need(response(wheel_a,59)==2 and response(wheel_b,59)==60, 'wheel first new response at59')
    need(30031==59*509 and all(all(p%q for q in range(2,isqrt(p)+1)) for p in (30013,30047,30029,59,509)), 'exact wheel factor controls')
    # Finite complete normalized-future controls, including non-equivalent histories.
    comparison_count=0
    for remaining in range(1,13):
        for ell in (1,3,5):
            future=[p for p in partitions(remaining) if not p or max(p)<=ell]
            alphabet=list(range(1,min(ell,remaining)+1,2))
            past=[p for mass in range(0,22) for p in partitions(mass,ell+2)]
            # Fix past degree, as required for DP state comparison.
            by_degree=defaultdict(list)
            for p in past: by_degree[sum(p)].append(p)
            for histories in by_degree.values():
                normalized={p:tuple(sum(defect(p,d) for d in q) for q in future) for p in histories}
                for p,q in combinations(histories,2):
                    eq=all(defect(p,d)==defect(q,d) for d in alphabet)
                    need(eq==(normalized[p]==normalized[q]),'coarsest normalized response')
                    comparison_count+=1
    for bad in (-1,True,1.5,'3'):
        try: count_target(bad)
        except ValueError: pass
        else: raise RuntimeError('bad input accepted')
    return comparison_count


def main():
    comparisons=controls()
    benchmarks=[]
    for n in (30,40,60,80,100):
        value,new=count_target(n,True)
        oldvalue,old=count_target(n,False)
        need(value==oldvalue,('larger response agreement',n))
        need(new['states']<old['states'],('strict graph shrink',n))
        benchmarks.append({'target':n,'old_states':old['states'],'new_states':new['states'],
            'old_transitions':old['transitions'],'new_transitions':new['transitions'],
            'old_peak':old['peak_states'],'new_peak':new['peak_states'],
            'exact_count':str(value.numerator)})
    print(json.dumps({'status':'PASS','normalized_future_comparisons':comparisons,
        'unmerged_partition_range':[0,27],'literal_edge_action_range':[0,7],
        'strict_merger':{'prefixes':[[23],[5,9,9]],'used_degree':23,'future_lengths':[1,3]},
        'zero_defect_merger':[[27],[9,9,9]],
        'wheel_horizon':{'histories':[[30013,30047],[30029,30031]],'common_degree':60060,'first_differing_odd_response':59},
        'graph_comparisons':benchmarks,
        'scope':'fixed-target ordinary A000568 only; no pair-free substitution or no-fixed mode'},indent=2))


if __name__=='__main__':
    main()
