#!/usr/bin/env python3
"""C7 strong fixed fibres below225; completeness uses cited Busch f(11)=225."""
from collections import Counter
import json
from verify import fixed_masks, rows, path_count, scc_count, require
out={}
for n in range(7,11):
    parts=(1,)*(n-7)+(7,)
    vals=Counter(); strong=Counter(); below=Counter(); total=0
    for mask in fixed_masks(parts):
        total+=1; r=rows(n,mask); h=path_count(r)
        vals[h]+=1
        if scc_count(r)==1:
            strong[h]+=1
            if h<225: below[h]+=1
    out[n]={'cycle_type':parts,'fixed_labelled_count':total,'strong_count':sum(strong.values()),'strong_minimum':min(strong) if strong else None,'strong_H_below_225':dict(sorted(below.items()))}
require(out[7]['strong_H_below_225']=={175:6,189:2},'inherited circulants')
require(all(not out[n]['strong_H_below_225'] for n in range(8,11)),'new bounded fibre check')
print(json.dumps({'status':'PASS','table':out,'scope':'All strong tournaments with a specified C7 action and H<225, using the cited all-order Busch floor to exclude n>=11.'},indent=2,sort_keys=True))
