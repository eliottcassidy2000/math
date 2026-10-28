"""Parent-only canonical-root generation without materializing root-switch rows.

Uses the preceding session's exact parent isomorphism implementation. The new
step is early rejection at a smaller pointed root and compact pointed indices.
Canonical augmentation itself is classical; no algorithmic priority claim.
"""
import argparse
from collections import Counter,defaultdict
from functools import lru_cache
import gzip,hashlib,importlib.util,json
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'relation_tower_frobenius_ports_20261028'
spec=importlib.util.spec_from_file_location('parent_oracle',PRIOR/'check_root_switch.py')
old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)

def generate(parents,N,cache_size,degree_first=False):
    started=perf_counter(); n=N-1
    groups=[old.isomorphisms(a,a,True) for a in parents]
    mask_ids=[]; reps=[]; offsets=[]; total=0
    buckets=defaultdict(list)
    for k,(a,group) in enumerate(zip(parents,groups)):
        cm=[min(old.move_mask(s,g) for g in group) for s in range(1<<n)]
        representatives=sorted(set(cm)); rank={s:i for i,s in enumerate(representatives)}
        mask_ids.append(tuple(rank[s] for s in cm)); reps.append(representatives)
        offsets.append(total); total+=len(representatives)
        buckets[old.colours(a)[1]].append(k)
    lookups=0
    @lru_cache(maxsize=cache_size)
    def classify(a):
        nonlocal lookups
        lookups+=1
        for k in buckets[old.colours(a)[1]]:
            maps=old.isomorphisms(a,parents[k])
            if maps:
                p=maps[0]
                old.need(all((a[i]>>j&1)==(parents[k][p[i]]>>p[j]&1) for i in range(n) for j in range(n)),'parent transport valid')
                return k,p
        raise RuntimeError('missing parent')
    children=[]; auts=Counter(); tested=0; rejected=0; accepted_trace=0
    processed=0; rejection_positions=Counter(); degree_rejections=0; degree_skips=0
    for k,ss in enumerate(reps):
        for rank,s in enumerate(ss):
            idx=offsets[k]+rank
            child=old.extend(parents[k],s); diagonal=1
            if degree_first:
                degrees=[row.bit_count() for row in child]
                if min(degrees)<s.bit_count():
                    degree_rejections+=1; rejected+=1; processed+=1
                    continue
            for v in range(n):
                if degree_first and degrees[v]>s.bit_count():
                    degree_skips+=1
                    continue
                tested+=1
                a,mask=old.delete(child,v); parent,p=classify(a)
                key=offsets[parent]+mask_ids[parent][old.move_mask(mask,p)]
                if key<idx:
                    rejected+=1; rejection_positions[v+1]+=1
                    break
                diagonal+=(key==idx)
            else:
                children.append(child)
                aut=diagonal*sum(old.move_mask(s,g)==s for g in groups[k])
                auts[aut]+=1; accepted_trace+=diagonal
                old.need(aut%2==1,'odd automorphism order')
                if N<=7: old.need(aut==len(old.isomorphisms(child,child,True)),'independent child automorphism audit')
            processed+=1
            if N==9 and processed%250000==0:
                print(json.dumps({'order':N,'processed':processed,'total':total,'seconds':round(perf_counter()-started,2)}),flush=True)
    old.need(processed==total and rejected+len(children)==total,'partition of candidates')
    old.need(len(children)==old.burnside(N),'independent ordinary count')
    row={'order':N,'classes':len(children),'pointed_states':total,'reroots_examined':tested,'old_all_rows_reroots':N*total,'own_root_queries_saved':total,'rejected_candidates':rejected,'degree_rejections':degree_rejections,'higher_degree_roots_skipped':degree_skips,'reject_position_histogram':dict(sorted(rejection_positions.items())),'automorphism_order_histogram':dict(sorted(auts.items())),'parent_cache_capacity':cache_size,'parent_cache':classify.cache_info()._asdict(),'exact_parent_isomorphism_lookups':lookups,'seconds':round(perf_counter()-started,3)}
    print(json.dumps({k:v for k,v in row.items() if k not in ('automorphism_order_histogram','reject_position_histogram')}),flush=True)
    return children,row

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--max-order',type=int,default=9); ap.add_argument('--cache-size',type=int,default=65536); ap.add_argument('--degree-first',action='store_true'); args=ap.parse_args()
    if not 1<=args.max_order<=9: raise ValueError('bounded execution 1..9')
    parents=[()]; rows=[]
    for N in range(1,args.max_order+1):
        parents,row=generate(parents,N,args.cache_size,args.degree_first); rows.append(row)
    raw=(json.dumps({'order':args.max_order,'rows':parents},separators=(',',':'))+'\n').encode()
    sha=hashlib.sha256(raw).hexdigest()
    prior=json.loads((PRIOR/'root_switch.json').read_text(encoding='utf-8'))
    for a,b in zip(rows,prior['results']):
        old.need(a['classes']==b['classes'] and a['pointed_states']==b['pointed_states'],'prior count agreement')
        old.need({str(k):v for k,v in a['automorphism_order_histogram'].items()}==b['automorphism_order_histogram'],'all automorphism orders agree')
    if args.max_order==9:
        reference=gzip.decompress((PRIOR/'tournament_bank.json.gz').read_bytes())
        if not args.degree_first:
            old.need(raw==reference,'entire bank byte identity')
        else:
            # Different chosen roots change representatives. Match every new
            # representative to the prior bank by an exact bijective class map.
            original=[tuple(a) for a in json.loads(reference)['rows']]
            buckets=defaultdict(list)
            for i,a in enumerate(original): buckets[old.colours(a)[1]].append(i)
            seen=set()
            for a in parents:
                found=[i for i in buckets[old.colours(a)[1]] if old.isomorphisms(a,original[i])]
                old.need(len(found)==1 and found[0] not in seen,'exact reference-bank class bijection')
                seen.add(found[0])
            old.need(len(seen)==len(original),'reference coverage')
    hostile=old.hostile()
    report={'status':'FINITE-EXACT','generation':'least pointed-root certificate; early reject at first smaller root','degree_first':args.degree_first,'results':rows,'bank_sha256':sha,'byte_identical_to_prior_bank':args.max_order==9 and not args.degree_first,'exact_class_bijection_to_prior_bank':args.max_order==9,'hostile_81_27':hostile,'source_oracle':str((PRIOR/'check_root_switch.py').relative_to(HERE.parents[2])),'checks':old.CHECKS}
    (HERE/('root_degree.json' if args.degree_first else 'root_minimum.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'bank_sha256':sha,'bank_classes':len(parents),'status':'PASS'}))

if __name__=='__main__':main()
