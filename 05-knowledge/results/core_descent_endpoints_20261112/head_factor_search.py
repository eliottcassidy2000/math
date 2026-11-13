"""Whole head identities with every attained content factor and the new birth bank."""
from context import *
from fractions import Fraction
from math import gcd,isqrt
from collections import defaultdict,Counter
import time

def divisors(n):
    small=[d for d in range(1,isqrt(n)+1) if n%d==0]
    return sorted(set(small+[n//d for d in small]))

def run():
    start=time.time();bank=load('hamilton_clocks_9.json.gz')['records'];factors={};stats=Counter()
    for r in bank:
        h=value(r['exponential_coefficients'],1)
        if h not in factors:factors[h]=r['rows']
    queries=defaultdict(list)
    for r in load('even_valuations.json.gz')['records']:
        for d in range(9):
            for shift in range(13):
                t=[Fraction(c*(i+1)**shift,2**d) for i,c in enumerate(r['exponential_coefficients'])];t[0]+=1
                if any(c.denominator!=1 for c in t):continue
                t=[int(c) for c in t];G=abs(gcd(*t))
                for factor in divisors(G):
                    if factor not in factors:continue
                    profile=trim([c//factor for c in t]);queries[profile].append({'birth_id':r['birth_id'],'source_port':r['port'],'division_exponent':d,'shift':shift,'factor':factor,'factor_rows':factors[factor]})
    hits=[];seen_profiles=set();firsts={value(p,1) for p in queries}
    for r in bank:
        p=tuple(r['exponential_coefficients'])
        if p in queries and p not in seen_profiles:
            seen_profiles.add(p)
            hits.extend({**q,'repair':r,'source_bank':'complete pointed through9'} for q in queries[p])
    del bank
    stats['old_matching_profiles']=len(seen_profiles);stats['queries']=sum(map(len,queries.values()))
    stats['distinct_target_profiles']=len(queries);seen=set()
    names=['birth_extensions_11.json.gz','birth_extensions_12.json.gz','birth_extensions_13.json.gz','birth_extensions_13_synchrony.json.gz','birth_extensions_14_rank5.json.gz','birth_sources_synchrony_12.json.gz']
    for name in names:
        for r in load(name)['records']:
            rows=r['rows'];key=tuple(rows)
            if key in seen:continue
            seen.add(key);stats['new_source_presentations']+=1
            h=r.get('H')
            if h is None:h=int(tables(np.array(rows,dtype=np.int64))[0])
            if h not in firsts:continue
            stats['H_filtered_new_sources']+=1
            for p in range(len(rows)):
                e=trim(expcoeff([int(c) for c in skeleton(np.array(rows,dtype=np.int64),p)]));stats['new_full_pointed_profiles']+=1
                if e in queries and e not in seen_profiles:
                    seen_profiles.add(e)
                    hits.extend({**q,'repair':{'rows':rows,'port':p,'exponential_coefficients':list(e)},'source_bank':name} for q in queries[e])
    oldclosed={(r['source']['birth_id'],r['source']['port']) for r in load('repaired_rays.json.gz')['families']}
    fresh=sorted({(h['birth_id'],h['source_port']) for h in hits}-oldclosed)
    out={'scope':'Every478 old even axis; divisions0..8; shifts0..12; every divisor of the integer profile content having a physical factor witness through9. Whole pointed H kernels from every ordinary tournament through9 plus all listed retained birth layers. New-source duplicate row presentations are suppressed; this is not an unmarked size census. An identity is not by itself full ray closure.','statistics':dict(stats),'source_receipts':names,'hits':hits,'matched_old_axes':sorted({(h['birth_id'],h['source_port']) for h in hits}),'matched_axes_outside_previous_five':fresh}
    save('head_factor_search.json.gz',out)
    print(json.dumps({k:v for k,v in out.items() if k!='hits'}),flush=True)
    print('seconds',round(time.time()-start,2),'hits',len(hits),flush=True)
if __name__=='__main__':run()
