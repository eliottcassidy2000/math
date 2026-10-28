#!/usr/bin/env python3
"""Optional second-stage saving: certified subgroup orbits within each first row."""
from collections import Counter
import argparse
import json
from time import perf_counter
from baseline_sources import generation as old
import check_blockfill as base


def closure(generators, n):
    identity=tuple(range(n)); group={identity}; stack=[identity]
    while stack:
        p=stack.pop()
        for g in generators:
            q=base.compose(g,p)
            if q not in group:
                group.add(q); stack.append(q)
    return tuple(sorted(group))


def adaptive_row(engine, child, source, parent_group):
    N=len(child); n=N-1
    # This source is already a canonical parent with the least mask in its orbit.
    stabilizer=tuple(p+(n,) for p in parent_group if old.move_mask(source[1],p)==source[1])
    old.need(stabilizer, 'missing identity stabilizer')
    records=[None]*N; records[n]=(source,tuple(sorted(stabilizer)))
    group=stabilizer; queries=0; propagated=0; enlargements=0
    while True:
        # Any certified subgroup transporter gives ALL pointed transports at v,
        # because the already known record at u is complete.
        for u in range(N):
            if records[u] is None:
                continue
            key,maps=records[u]
            for g in group:
                v=g[u]
                if records[v] is None:
                    records[v]=(key,tuple(sorted(base.compose(g,p) for p in maps)))
                    propagated+=1
        missing=[v for v in range(N) if records[v] is None]
        if not missing:
            break
        v=missing[0]
        record=base.aligned_root(engine,child,v); queries+=1
        key,maps=record
        generators=list(group)
        inv=base.inverse(maps[0])
        generators.extend(base.compose(p,inv) for p in maps) # stabilizer at v
        matches=[u for u,r in enumerate(records) if r is not None and r[0]==key]
        if matches:
            u=matches[0]
            generators.append(base.compose(maps[0],base.inverse(records[u][1][0])))
        records[v]=record
        next_group=closure(generators,N)
        enlargements+=len(next_group)>len(group)
        group=next_group
    records=tuple(records)
    all_automorphisms=tuple(sorted(p for key,maps in records if key==source for p in maps))
    old.need(group==all_automorphisms,'adaptive subgroup failed to recover full group')
    row=tuple(sorted(Counter(key for key,maps in records).items()))
    return row,records,{'classified_roots':queries,'inferred_roots':propagated,
                        'known_original_root':1,'subgroup_enlargements':enlargements}


def build(parents,n):
    start=perf_counter()
    engine,states,parent_data,counts=base.prepare(parents,n)
    parent_groups={code:group for code,raw,incoming,group in parent_data}
    rows,roots,groups,representatives={},{},{},{}
    measures=Counter(); N=n+1
    for source,child in sorted(states.items()):
        if source in rows: continue
        row,records,metrics=adaptive_row(engine,child,source,parent_groups[source[0]])
        measures.update(metrics)
        base.check_row(row,source,states,N)
        representatives[row]=child
        for target,w in row:
            old.need(target not in rows,'overlapping blocks')
            target_records=records if target==source else base.fill_from_anchor(records,target,N)
            rows[target]=row; roots[target]=target_records
            groups[target]=tuple(sorted(p for key,maps in target_records if key==target for p in maps))
    old.need(set(rows)==set(states),'uncovered states')
    old.need(measures['classified_roots']+measures['inferred_roots']+measures['known_original_root']==N*len(representatives),
             'root accounting')
    counts.update(dict(measures)); counts.update({'order':N,'ordinary_classes':len(representatives),
        'pointed_states':len(states),'all_canonicalizer_requests':engine.requests,
        'total_searches':engine.searches,'total_orders_tested':engine.orders_tested,
        'seconds':perf_counter()-start})
    return {'states':states,'Q':rows,'root_transports':roots,'automorphisms':groups,
            'parent_data':parent_data,'bank':tuple(representatives.values()),'counts':counts}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--max-order',type=int,default=8); args=ap.parse_args()
    old.need(1<=args.max_order<=8,'finite scope 1..8')
    parents=[()]; results=[]
    for n in range(args.max_order):
        fixed=base.build(parents,n,'blockfill'); adaptive=build(parents,n)
        hashes=base.compare(fixed,adaptive)
        results.append({'order':n+1,'fixed_first_row':fixed['counts'],
                        'adaptive':adaptive['counts'],'equal_complete_output_sha256':hashes})
        parents=list(adaptive['bank'])
    return {'status':'FINITE-EXACT','source_pin':base.PIN,'results':results,
            'timing_scope':'One illustrative pass, not a repeated benchmark; same raw mask preparation and all outputs. Subgroup closure overhead may outweigh saved classifications.'}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
