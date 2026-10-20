#!/usr/bin/env python3
"""Replay each row using the uncompressed vertex-pair action engine.
Shares normalizer generation and graph-isomorphism/global-sign primitives.
Independent of coordinate-action construction, action deduplication, and cube order.
"""
import fixed_cube_engine as m
uncompressed=m.relative_records
original_classify=m.classify
import extension_census as e
import json,time,sys

def main():
    expected=json.load(open(__file__.replace('pair_action_replay.py','extension_counts.json')))
    bykey={(r['n'],r['class']):r for r in expected['rows']}
    m.relative_records=uncompressed
    # Replaying the original classifier additionally guards the WL bucket filter.
    m.classify=original_classify
    fields=['order','normalizer_order','pair_orbits','fixed_objects_each','relative_tournament_orbits','relative_even_graph_orbits','tau_T','tau_E','difference','tournament_multiplicity_histogram','even_graph_multiplicity_histogram','escape_graph_multiplicity_histogram','tournament_fusion_defect','even_graph_fusion_defect','sign_escape_relative_orbits','sign_escape_global_types']
    receipts=[]
    for n,name,gens in e.cases():
        start=time.monotonic();r=m.census_case(n,name,gens);want=bykey[n,name]
        for key in fields:
            # JSON normalizes integer histogram keys to strings.
            e.need(json.loads(json.dumps(r[key]))==want[key],f'pair-action replay mismatch: {n} {name} {key}')
        receipts.append({'n':n,'class':name,**{key:r[key] for key in fields}})
        print(n,name,'PASS',round(time.monotonic()-start,2),'seconds',file=sys.stderr,flush=True)
    print(json.dumps({'status':'PASS','scope':'all 18 classified mixed-prime rows at n=11,12','method':'uncompressed pair actions and original unbucketed isomorphism classification','shared_primitives':'normalizer candidates, pair-orbit construction, full isomorphism search, twin-quotient global sign','checks':m.CHECKS,'rows':receipts},sort_keys=True,indent=2))
if __name__=='__main__':main()
