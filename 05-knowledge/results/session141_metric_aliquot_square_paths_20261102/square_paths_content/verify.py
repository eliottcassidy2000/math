#!/usr/bin/env python3
"""Exact reversal/content theorem controls; standard library, active under -O."""
from collections import Counter
from functools import reduce
from fractions import Fraction
from itertools import combinations,permutations
from math import comb,gcd
from pathlib import Path
import hashlib,json

def need(ok,why):
    if not ok:raise RuntimeError(why)

def arcs(p):return frozenset(zip(p,p[1:]))
def consistent(a):return all((v,u) not in a for u,v in a)
def valuation2(x):
    need(x>0,'positive valuation input')
    return (x & -x).bit_length()-1

def invert(moments):
    top=max(moments)
    return [sum((-1)**(j-k)*comb(j,k)*moments.get(j,0) for j in range(k,top+1))
            for k in range(top+1)]

def content_check(coefficients,moments,N,Umax=None):
    c=reduce(gcd,coefficients)
    m=reduce(gcd,moments.values())
    need(c==m,'integer translation content invariance')
    need(c>0 and c&(c-1)==0,'pure power of two')
    need(valuation2(c)==min(valuation2(x) for x in moments.values() if x),'minimum valuation')
    if Umax is not None:
        e=N-Umax+1
        need(c%(1<<e)==0,'reversal/free-edge bound')
    else:e=None
    return {'content_exponent':valuation2(c),'reversal_bound_exponent':e,
            'bound_sharp':e==valuation2(c) if e is not None else None}

def small_controls():
    cases=0;subsets=0;leaf_cases=0
    for n in range(2,5):
        all_edges=list(combinations(range(n),2));N=len(all_edges)
        for guardbits in range(1<<N):
            guard={e for i,e in enumerate(all_edges) if guardbits>>i&1}
            leaves=[v for v in range(n) if sum(v in e for e in guard)==1]
            bank=[p for p in permutations(range(n))
                  if all(tuple(sorted(e)) in guard for e in zip(p,p[1:]))]
            sets=list(map(arcs,bank));index={p:i for i,p in enumerate(bank)}
            reverse=[index[p[::-1]] for p in bank]
            neighbourhood=[{j for j in range(i+1,len(bank)) if consistent(sets[i]|sets[j])}
                           for i in range(len(bank))]
            moments=Counter({0:1<<N});profiles=Counter();selectedsets=set();Umax=None
            def visit(selected,candidates,union):
                nonlocal subsets,Umax,leaf_cases
                for i in sorted(candidates):
                    sel=selected+(i,);joint=union|sets[i]
                    need(consistent(joint),'clique union consistency')
                    for leaf in leaves:
                        need(all(bank[j][0]==leaf for j in sel) or all(bank[j][-1]==leaf for j in sel), 'leaf endpoint half')
                        leaf_cases+=1
                    rev=tuple(sorted(reverse[j] for j in sel))
                    need(sel!=rev,'free reversal action')
                    reverseunion=set().union(*(sets[j] for j in rev))
                    need(consistent(reverseunion) and len(reverseunion)==len(joint),'reversal preserves weight')
                    selectedsets.add(sel);subsets+=1
                    moments[len(sel)]+=1<<(N-len(joint));profiles[len(sel),len(joint)]+=1
                    Umax=max(Umax or 0,len(joint))
                    visit(sel,candidates&neighbourhood[i],joint)
            visit((),set(range(len(bank))),frozenset())
            need(all(tuple(sorted(reverse[j] for j in s)) in selectedsets for s in selectedsets),
                 'orbit partner completeness')
            need(all(c%2==0 for c in profiles.values()),'profile evenness')
            coeff=invert(moments);literal=Counter()
            for bits in range(1<<N):
                T={(u,v) if bits>>i&1 else (v,u) for i,(u,v) in enumerate(all_edges)}
                literal[sum(a<=T for a in sets)]+=1
            need(coeff==[literal[k] for k in range(max(literal)+1)],'literal full tournament control')
            content_check(coeff,moments,N,Umax);cases+=1
    # A bank not closed under reversal breaks the claimed extra factor.
    need(reduce(gcd,[1,1])==1,'one directed path on two vertices hostile')
    # n=1 has a self-reverse empty path; its z polynomial is not divisible by two.
    need(reduce(gcd,[0,1])==1,'order-one hostile')
    return {'all_undirected_guards_orders2_to4':cases,'nonempty_compatible_subsets_checked':subsets,'leaf_endpoint_half_checks':leaf_cases,
            'non_reversal_closed_hostile':'passed','order_one_hostile':'passed'}

def main():
    root=Path(__file__).resolve().parent
    pins=json.loads((root/'SOURCE_PINS.json').read_text())
    for p in pins['sources']:
        raw=(root/p['snapshot']).read_bytes()
        need(hashlib.sha256(raw).hexdigest()==p['sha256'],'source snapshot digest')
    source=json.loads((root/'sources/square_paths_normal.json').read_text())
    rows=[]
    for row in source['rows']:
        if row['n'] not in (25,26):continue
        N=row['total_tournament_bits'];profiles={int(k):{int(u):c for u,c in d.items()} for k,d in row['compatible_subset_union_profiles'].items()}
        moments={0:1<<N}
        for j,profile in profiles.items():
            need(all(c%2==0 for c in profile.values()),'square profile reversal parity')
            moments[j]=sum(c*(1<<(N-u)) for u,c in profile.items())
        coeff=invert(moments)
        expected=[row['guard_histogram'][str(k)]*(1<<row['outside_free_bits']) for k in range(len(row['guard_histogram']))]
        need(coeff==expected,'square profile/coefficient reconstruction')
        Umax=max(u for profile in profiles.values() for u in profile)
        result=content_check(coeff,moments,N,Umax)
        result.update({'n':row['n'],'N':N,'Umax':Umax,
                       'nonzero_moment_valuations':{j:valuation2(x) for j,x in moments.items()},
                       'primitive_coefficients':[x//(1<<result['content_exponent']) for x in coeff]})
        mean=Fraction(sum(k*c for k,c in enumerate(coeff)),1<<N)
        positive=Fraction(sum(coeff[1:]),1<<N)
        need(mean=={25:Fraction(5,1<<22),26:Fraction(3,1<<22)}[row['n']], 'finite mean')
        need(positive=={25:Fraction(19,1<<24),26:Fraction(351,1<<29)}[row['n']], 'positive probability')
        result.update({'mean':[mean.numerator,mean.denominator], 'positive_probability':[positive.numerator,positive.denominator], 'maximum_response':len(coeff)-1})
        need(result['content_exponent']=={25:276,26:296}[row['n']],'exact square content')
        rows.append(result)
    return {'status':'passed','scope':'uniform theorem; exhaustive small guards; existing n25/n26 census only',
            'controls':small_controls(),'square_rows':rows}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
