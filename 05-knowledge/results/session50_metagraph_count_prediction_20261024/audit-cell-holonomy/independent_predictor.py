#!/usr/bin/env python3
"""Independent atlas-through-seven construction and child-refinement-only sum.

Imports only this audit's elementary matrix utilities, never producer code.
The supplied mathematical zero-error theorem at strong order eight is explicit.
"""
import argparse
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import json
import sys
sys.dont_write_bytecode=True
import independent as audit

CANON_ORDERS=Counter()
def canonical_parent(a):
    n=len(a)
    audit.require(n<=7,'no child-order canonicalization')
    CANON_ORDERS[n]+=1
    degrees=[sum(row) for row in a]
    cells=[tuple(v for v in range(n) if degrees[v]==d) for d in sorted(set(degrees))]
    return min(audit.code(a,tuple(v for cell in ps for v in cell)) for ps in product(*(permutations(c) for c in cells)))

def attach(a,mask):
    n=len(a); b=[row[:]+[(mask>>v)&1] for v,row in enumerate(a)]
    b.append([1-((mask>>v)&1) for v in range(n)]+[0]); return b

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--receipt',required=True); args=parser.parse_args()
    _,six=audit.full_classes(6)
    keys=set()
    for a,g in six:
        for mask in range(64): keys.add(canonical_parent(attach(a,mask)))
    parents=[audit.decode(7,key) for key in sorted(keys)]
    audit.require(len(parents)==456,'complete independently generated parent atlas')
    pointed=0; upper=0; cells=Counter(); parent_group_orders=Counter(); strong7=0
    input_mass=0
    for a in parents:
        g,_=audit.degree_group_and_key(a)
        parent_group_orders[len(g)]+=1
        input_mass+=5040//len(g)
        if len(audit.strong_components(a))==1: strong7+=audit.excess(a,g)
        unseen=set(range(128))
        while unseen:
            mask=min(unseen)
            orbit={sum(1<<p[v] for v in range(7) if (mask>>v)&1) for p in g}
            audit.require(orbit<=unseen,'independent parent subset orbit partition')
            unseen.difference_update(orbit)
            b=attach(a,mask)
            size=next(len(c) for c in audit.stable(b) if 7 in c)
            pointed+=1; upper+=size; cells[size]+=1
    audit.require(input_mass==2**21,'parent labelled mass')
    audit.require(strong7==14,'strong-seven discrepancy from actual parent groups')
    correction=2*strong7
    numerator=upper-correction
    audit.require(numerator%8==0,'corrected upper trace integral')
    report={'schema':'session50-independent-refinement-prediction-v1','parent_source':'Complete full-permutation class partition of labelled order six, then all 64 attachments per class and degree-cell canonicalization only at order seven.',
            'parent_classes':len(parents),'parent_automorphism_order_histogram':dict(sorted(parent_group_orders.items())),
            'parent_labelled_mass':input_mass,'pointed_states':pointed,'root_refinement_size_histogram':dict(sorted(cells.items())),
            'root_cell_upper_sum':upper,'strong_error_seven':strong7,'strong_error_eight':{'value':0,'source':'Odd equitable-cell size proof in audited theorem, not a target census.'},
            'SCC_error_eight':correction,'corrected_trace':numerator,'predicted_order_eight_classes':numerator//8,
            'canonicalized_order_histogram':dict(sorted(CANON_ORDERS.items())),'deletion_or_reroot_isomorphism_tests':0,
            'checks':dict(sorted(audit.CHECKS.items()))}
    audit.require((pointed,upper,numerator//8)==(54256,55068,6880),'postcomputation comparison only')
    Path(args.receipt).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','parent_classes':len(parents),'root_upper_sum':upper,'E8':correction,'prediction':numerator//8,'max_canonicalized_order':max(CANON_ORDERS)},sort_keys=True))
if __name__=='__main__': main()
