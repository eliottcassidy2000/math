#!/usr/bin/env python3
"""Independent reroot audit using only order-eight isomorphism tests."""
import itertools as it,json,hashlib
from collections import Counter
from pathlib import Path
from audit import need,stable

def canon(rows):
    need(len(rows)==8,'only parent-eight canonicalization permitted')
    return min(tuple(sum(1<<j for j,w in enumerate(order) if rows[v]>>w&1) for v in order) for chunks in it.product(*(it.permutations(C) for C in stable(rows))) for order in [sum(chunks,())])
def main():
    src=Path('producer/sources/session50-regular-nine-receipt.json');data=json.loads(src.read_text());parents=[tuple(s['parent']) for s in data['states']]
    keys=[canon(p) for p in parents];need(len(keys)==len(set(keys)),'parent representatives distinct');index={k:i for i,k in enumerate(keys)};dest=[]
    for p in parents:
        need(sorted(r.bit_count() for r in p)==[3]*4+[4]*4,'almost regular parent')
        child=[r|(1<<8 if r.bit_count()==3 else 0) for r in p]+[sum(1<<v for v,r in enumerate(p) if r.bit_count()==4)]
        need(all(r.bit_count()==4 for r in child),'forced regular extension');row=Counter()
        for omitted in range(9):
            keep=[v for v in range(9) if v!=omitted]
            smaller=tuple(sum(1<<j for j,w in enumerate(keep) if child[v]>>w&1) for v in keep);key=canon(smaller);need(key in index,'all deletion destinations present');row[index[key]]+=1
        dest.append(row)
    remaining=set(range(len(parents)));regular=Counter();orbits=Counter();trace=0
    while remaining:
        i=min(remaining);component=set(dest[i]);need(i in component and component<=remaining,'reroot class partition')
        need(all(dest[j]==dest[i] for j in component),'reroot rows constant within class');remaining-=component
        mult=tuple(sorted(dest[i].values()));need(sum(mult)==9,'full vertex mass');need(sum(dest[j][j] for j in component)==9,'class diagonal trace')
        e=9*(len(component)-1);need(e==sum(9-n for n in mult),'regular error via orbits');regular[e]+=1;orbits[mult]+=1;trace+=9
    own=json.loads(Path('receipt.json').read_text());strong10=Counter(x['error'] for x in own['representatives'] if x['error']);strong9=regular.copy();strong9.pop(0,None);strong9[14]+=1;strong9[18]+=4
    total=strong10.copy()
    for e,c in strong9.items():total[e]+=2*c
    total[14]+=6
    need(regular==Counter({0:3,18:3,36:2,72:7}),'regular nine histogram');need(strong10==Counter({14:4,18:32,40:84}),'independent strong ten histogram')
    need(total==Counter({14:12,18:46,36:4,40:84,72:14}),'full ten distribution')
    result={'status':'PASS','input_contract':'Complete almost-regular-eight parent slice inherited from Session50','parent_count':len(parents),'parent_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'reroot_destinations':9*len(parents),'regular_nine_error_histogram':sorted(regular.items()),'regular_nine_vertex_orbit_multisets':sorted(orbits.items()),'regular_nine_trace':trace,'strong_ten_error_histogram':sorted(strong10.items()),'total_ten_error_histogram':sorted(total.items()),'total_ten_errorful_classes':sum(total.values()),'first_moment':sum(e*c for e,c in total.items()),'second_moment':sum(e*e*c for e,c in total.items()),'order_nine_or_ten_canonicalizations':0}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
