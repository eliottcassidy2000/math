#!/usr/bin/env python3
"""Compare independent and frozen producer facts; test deliberately bad receipts."""
from copy import deepcopy
import argparse
import json
from pathlib import Path

H = Path(__file__).resolve().parent

def require(ok, why):
    if not ok:
        raise ValueError(why)


def compare(a,p):
    require(a['status'] == p['status'] == 'PASS','status')
    require(a['total_labelled_tournaments'] == p['total_labelled_tournaments'] == 1099,'universe')
    require(len(a['exhaustive_rows']) == len(p['exhaustive_rows']) == 5,'row universe')
    for ar,pr in zip(a['exhaustive_rows'],p['exhaustive_rows']):
        for key,value in pr.items():
            require(ar.get(key) == value,'census disagreement: '+key)
    for key in ('ordinary_collision_order_4','merged_collision_order_5','merged_collision_order_8'):
        for field,value in a[key].items():
            require(p[key][field] == value,'witness disagreement: '+key+'.'+field)
    for key in ('all_small_recenter_checks_Q_orders_0_through_5','all_vs_distinguished_cycle_control','no_extreme_class_count_and_regular_lower_bound_controls'):
        require(a[key] == p[key],'control disagreement: '+key)
    return True


def hostile(a,p):
    cases=[]
    def reject(label, mutate):
        q=deepcopy(p)
        mutate(q)
        try:
            compare(a,q)
        except (ValueError,KeyError,TypeError):
            cases.append(label)
        else:
            raise ValueError('hostile result accepted: '+label)
    reject('inflate target class count',lambda q:q['exhaustive_rows'][4].__setitem__('target_unmarked_classes',12))
    reject('erase merged loss',lambda q:q['exhaustive_rows'][4].__setitem__('target_unmarked_converse_merged_classes',10))
    reject('drop a labelled order',lambda q:q['exhaustive_rows'].pop())
    reject('misstate full automorphism histogram',lambda q:q['exhaustive_rows'][1].__setitem__('target_automorphism_order_histogram',{'1':2}))
    reject('corrupt target adjacency',lambda q:q['merged_collision_order_5']['D_1_plus_Q'].__setitem__(0,'00000000000'))
    reject('corrupt explicit recenter map',lambda q:q['merged_collision_order_5']['explicit_D_isomorphism_old_to_new'].__setitem__(0,0))
    reject('confuse ordinary and merged witness',lambda q:q['ordinary_collision_order_4'].__setitem__('source_isomorphic_to_other_converse',False))
    reject('inflate path-cycle bijection to all cycles',lambda q:q['all_vs_distinguished_cycle_control'].__setitem__('distinguished_cycles',2))
    reject('inflate regular lower bound',lambda q:q['no_extreme_class_count_and_regular_lower_bound_controls'][-1].__setitem__('ordinary_lower_bound',12))
    return cases


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    a=json.loads((H/'normal.json').read_text())
    p=json.loads((H/'producer/normal.json').read_text())
    compare(a,p)
    checked=[]
    for rel, reference in [('optimized.json',a),('producer/optimized.json',p),('producer-replay-normal.json',p),('producer-replay-optimized.json',p)]:
        require(json.loads((H/rel).read_text()) == reference,'replay differs: '+rel)
        checked.append(rel)
    result={'status':'PASS','census_rows_compared':5,'complete_collision_witnesses_compared':3,
            'replay_receipts_compared':checked,'hostile_receipts_rejected':hostile(a,p)}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    Path(args.output).write_text(text)
    print(text,end='')

if __name__ == '__main__':
    main()
