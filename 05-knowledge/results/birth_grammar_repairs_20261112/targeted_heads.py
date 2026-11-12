"""Larger-core search for two early hostiles and all eight pure3 axes."""
from common import *
from hamilton_clocks import skeleton,trim
from math import gcd
from collections import defaultdict,Counter
from fractions import Fraction
import time

def run():
    chosen={3,4,14,42,28,35,84,90,92,97,98,206};queries=defaultdict(list)
    for r in read('even_valuations.json.gz')['records']:
        if r['birth_id'] not in chosen:continue
        j=r['eventual_valuations_by_parity'][0]
        for d in range(j):
            for shift in range(7):
                e=[Fraction(c*(b+1)**shift,2**d) for b,c in enumerate(r['exponential_coefficients'])];e[0]+=1
                if any(c.denominator!=1 for c in e):continue
                e=[int(c) for c in e];G=abs(gcd(*e))
                for factor in sorted({1,G}):
                    profile=trim([c//factor for c in e]);queries[value(profile,1)].append({'birth_id':r['birth_id'],'source_port':r['port'],'division_exponent':d,'shift':shift,'factor':factor,'profile':profile})
    sources=read('unit_ten.json.gz',True)['records']
    for n in [11,12]:sources+=read(f'birth_extensions_{n}.json.gz')['records']
    hits=[];stats=Counter();seen=set();start=time.time()
    for r in sources:
        if r['H'] not in queries or tuple(r['rows']) in seen:continue
        seen.add(tuple(r['rows']));stats['filtered_sources']+=1
        for p in range(len(r['rows'])):
            e=trim(expcoeff([int(c) for c in skeleton(np.array(r['rows'],dtype=np.int64),p)]));stats['full_pointed_profiles']+=1
            for q in queries[r['H']]:
                if e==q['profile']:hits.append({**q,'rows':r['rows'],'port':p})
        if stats['filtered_sources']%1000==0:print(json.dumps({'stats':dict(stats),'hits':len(hits),'seconds':round(time.time()-start,1)}),flush=True)
    out={'scope':'Twelve declared original even axes (two early hostiles, two positive-control families, all eight pure3 axes); all eventual unpaid divisions, shifts0..6, factor1 or full coefficient content. Source graphs: complete strong order10 graphs admitting a unit selection, plus retained birth extension layers11,12. This filtered source universe is not all tournaments at these sizes.','selected_birth_ids':sorted(chosen),'statistics':dict(stats),'hits':hits}
    save('targeted_heads.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='hits'}),flush=True);print('hits',len(hits),flush=True)
if __name__=='__main__':run()
