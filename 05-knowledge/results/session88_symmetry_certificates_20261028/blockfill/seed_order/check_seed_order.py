#!/usr/bin/env python3
"""Bounded follow-up: choose each block's seed by decreasing known stabilizer size.

Does not modify the published block-fill/adaptive algorithms. No timing claim.
"""
from collections import Counter
from pathlib import Path
import argparse
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from baseline_sources import generation as old
import check_blockfill as base
import check_adaptive as adaptive


def stabilizer_order(source,parent_group):
    return sum(old.move_mask(source[1],p)==source[1] for p in parent_group)


def build(parents,n):
    engine,states,parent_data,counts=base.prepare(parents,n)
    parent_groups={code:group for code,raw,incoming,group in parent_data}
    stabilizers={key:stabilizer_order(key,parent_groups[key[0]]) for key in states}
    traversal=sorted(states,key=lambda key:(-stabilizers[key],key))
    rows,roots,groups,representatives={},{},{},{}
    measures=Counter(); anchors=[]; N=n+1
    for source in traversal:
        if source in rows: continue
        child=states[source]
        row,records,metrics=adaptive.adaptive_row(engine,child,source,parent_groups[source[0]])
        measures.update(metrics)
        base.check_row(row,source,states,N)
        weights=dict(row)
        old.need(weights[source]==min(weights.values()),'largest stabilizer failed to select smallest orbit')
        full_group_size=sum(len(maps) for key,maps in records if key==source)
        old.need(all(weights[key]*stabilizers[key]==full_group_size for key in weights),
                 'orbit-stabilizer seed accounting')
        canonical_seed=row[0][0]
        # Keep exactly the old output bank and ordering despite a different seed.
        representatives[canonical_seed]=states[canonical_seed]
        anchors.append({'block_minimum':canonical_seed,'source':source,
            'original_source_orbit':weights[canonical_seed],'chosen_source_orbit':weights[source],
            'original_stabilizer_order':stabilizers[canonical_seed],
            'chosen_stabilizer_order':stabilizers[source],**metrics})
        for target,w in row:
            old.need(target not in rows,'overlapping blocks')
            target_records=records if target==source else base.fill_from_anchor(records,target,N)
            rows[target]=row; roots[target]=target_records
            groups[target]=tuple(sorted(p for key,maps in target_records if key==target for p in maps))
    old.need(set(rows)==set(states),'uncovered states')
    old.need(measures['classified_roots']+measures['inferred_roots']+measures['known_original_root']==N*len(representatives),
             'root accounting')
    # Costs depend on the selected rooted presentation. Compare changed blocks
    # directly, outside all production counters; unchanged seeds are identical.
    audit=old.ParentCanonicalizer(max_order=n)
    differences=Counter(); delta_total=0; old_queries=0; examples={}
    for anchor in anchors:
        k=anchor['block_minimum']; source=anchor['source']
        new=anchor['classified_roots']
        if source==k:
            previous=new
        else:
            unused_row,unused_records,metrics=adaptive.adaptive_row(audit,states[k],k,parent_groups[k[0]])
            previous=metrics['classified_roots']
        delta=new-previous; delta_total+=delta; old_queries+=previous
        differences[delta]+=1
        if delta and str(delta) not in examples:
            examples[str(delta)]={**anchor,'previous_classified_roots':previous,'query_change':delta}
    counts.pop('preparation_seconds')
    counts.update(dict(measures)); counts.update({'order':N,'ordinary_classes':len(representatives),
        'pointed_states':len(states),'all_canonicalizer_requests':engine.requests,
        'total_searches':engine.searches,'total_orders_tested':engine.orders_tested,
        'changed_seed_blocks':sum(a['source']!=a['block_minimum'] for a in anchors),
        'strictly_smaller_seed_orbit_blocks':sum(a['chosen_source_orbit']<a['original_source_orbit'] for a in anchors),
        'per_block_query_change_histogram':dict(sorted(differences.items())),
        'total_query_change':delta_total,'independently_reconstructed_original_queries':old_queries,
        'query_change_examples':examples})
    return {'states':states,'Q':rows,'root_transports':roots,'automorphisms':groups,
        'parent_data':parent_data,'bank':tuple(representatives[k] for k in sorted(representatives)),
        'counts':counts}


def adverse_seed_explanation():
    explanation=[]
    for source in ((131074,103),(155648,2)):
        child=old.extend(old.matrix(7,source[0]),source[1])
        row,records=base.complete_row(old.ParentCanonicalizer(),child)
        parent_group=old.ParentCanonicalizer().maps(old.matrix(7,source[0]))[1]
        H=tuple(p+(7,) for p in parent_group if old.move_mask(source[1],p)==source[1])
        first_key,maps=records[0]
        generators=list(H)
        inv=base.inverse(maps[0])
        generators.extend(base.compose(p,inv) for p in maps)
        if first_key==source:
            generators.append(base.compose(maps[0],base.inverse(H[0])))
        G=adaptive.closure(generators,8)
        old.need(len(G)==15,'adverse seed full-group recovery after first query')
        explanation.append({'source':source,'literal_child_code':old.encode(child,tuple(range(8))),
            'row':row,'initial_stabilizer_order':len(H),
            'first_classified_vertex':0,'first_classified_type':first_key,
            'group_order_after_first_classification':len(G),
            'distinct_root_types_known_after_first_classification':len({source,first_key})})
    old.need([x['distinct_root_types_known_after_first_classification'] for x in explanation]==[2,1],
             'adverse seed first-query explanation changed')
    return explanation


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-order',type=int,default=8);args=ap.parse_args()
    old.need(1<=args.max_order<=8,'finite scope 1..8')
    parents=[()];results=[]
    for n in range(args.max_order):
        original=adaptive.build(parents,n); reordered=build(parents,n)
        hashes=base.compare(original,reordered)
        previous=original['counts']['classified_roots']; current=reordered['counts']['classified_roots']
        old.need(reordered['counts']['independently_reconstructed_original_queries']==previous,
                 'per-block original cost reconstruction mismatch')
        old.need(current-previous==reordered['counts']['total_query_change'],'query delta accounting')
        results.append({'order':n+1,'original_adaptive_queries':previous,
                        'stabilizer_first':reordered['counts'],'equal_complete_output_sha256':hashes})
        parents=list(original['bank'])
    old.need([r['stabilizer_first']['classified_roots'] for r in results]==[0,1,3,8,37,242,2586,47377][:args.max_order],
             'post-construction query-count regression mismatch')
    return {'status':'FINITE-EXACT','source_pin':base.PIN,'results':results,
            'adverse_seed_explanation':adverse_seed_explanation(),
            'scope':'Same complete Q, named states, raw bank ordering, all root transports, full automorphisms and parent maps. Exact query-count comparison only; no runtime benchmark.',
            'selection_theorem':'Within a block |H_j|=|Aut(U)|/m_j, so decreasing known |H_j| makes the first seed a smallest vertex orbit. This does not imply minimum adaptive row cost.'}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
