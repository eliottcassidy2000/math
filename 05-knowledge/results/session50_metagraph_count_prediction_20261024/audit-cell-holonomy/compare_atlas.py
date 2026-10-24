#!/usr/bin/env python3
"""Authenticate frozen producer parent records against independent exact atlases."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import independent as audit
from independent_predictor import canonical_parent,attach

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--receipt',required=True); args=parser.parse_args()
    source=Path(__file__).resolve().parent/'source/predictor/parents-through-seven.json'
    raw=source.read_bytes(); payload=json.loads(raw)
    audit.require(hashlib.sha256(raw).hexdigest()=='df16881db1060e0423ed73658bdcb5923266f839681b74ac5b1bddd5a232f3ad','frozen atlas source hash')
    expected={}; six=None
    for n in range(7):
        _,classes=audit.full_classes(n)
        expected[n]={canonical_parent(a) for a,g in classes}
        if n==6: six=classes
    expected[7]={canonical_parent(attach(a,mask)) for a,g in six for mask in range(64)}
    summary=[]
    for level in payload['levels']:
        n=level['order']; observed=set(); groups=Counter()
        for record in level['parents']:
            rows=record['rows']; a=[[(rows[v]>>w)&1 for w in range(n)] for v in range(n)]
            audit.require(all(a[v][v]==0 for v in range(n)) and all(a[v][w]+a[w][v]==1 for v in range(n) for w in range(v)),'input actual tournament')
            key=canonical_parent(a)
            audit.require(key not in observed,'input classes pairwise nonisomorphic')
            observed.add(key)
            literal,_=audit.degree_group_and_key(a)
            supplied={tuple(g) for g in record['actions']}
            audit.require(set(literal)==supplied,'input actions are full literal automorphism group')
            groups[len(literal)]+=1
        audit.require(observed==expected[n],'complete independent parent atlas match')
        summary.append({'order':n,'classes':len(observed),'group_order_histogram':dict(sorted(groups.items()))})
    report={'status':'PASS','atlas_sha256':hashlib.sha256(raw).hexdigest(),'source':'Complete independently generated atlases; literal degree-preserving arc tests authenticate every supplied group.','levels':summary,'checks':dict(sorted(audit.CHECKS.items()))}
    Path(args.receipt).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','complete_atlas_levels':len(summary),'records':sum(x['classes'] for x in summary),'full_actions_authenticated':True},sort_keys=True))
if __name__=='__main__': main()
