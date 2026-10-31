#!/usr/bin/env python3
"""Hostile and positive controls for the one-run spectrum bound."""
from itertools import product
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'square_sum_spectrum_overnight_20261030'))
from audit_feedback import gaps,require,subset_dp,tournament


def main():
    contexts=[()]+[(a,) for a in range(4)]
    cyclic=acyclic=checks=independent=0
    for p,q,a in product(contexts,contexts,range(4)):
        one=p+(a,)+q;rows=tournament(one);b=2+len(p)
        triangle=any((rows[b]>>x&1) and (rows[x]>>y&1) and (rows[y]>>b&1)
                     for x in range(len(rows)) for y in range(len(rows)))
        if triangle:cyclic+=1
        else:acyclic+=1
        for m in range(7):
            w=p+(a,)*m+q;h=gaps(w)
            if triangle and m:require(h>=2**m+1,'cyclic module exponential lower bound')
            if not triangle:require(h==gaps(p+q),'acyclic module erasure')
            if m<=3:require(subset_dp(tournament(w))==h,'independent literal count');independent+=1
            checks+=1
    print(json.dumps({'status':'passed','fixed_context_length_at_most':1,'run_lengths':list(range(7)),
                      'templates':100,'cyclic_templates':cyclic,'acyclic_templates':acyclic,
                      'counts_checked':checks,'independent_subset_dp_counts':independent,
                      'scope':'Finite controls; the all-run bound is proved by an explicit odd-cycle injection.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
