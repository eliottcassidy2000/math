#!/usr/bin/env python3
"""Independent audit of the p=29 hostile and adjacent new finite minima."""
from pathlib import Path
from collections import defaultdict
from math import comb
import argparse,json
from verify_period_ranks import matrix_mod,rank_forward
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--base',type=Path);args=parser.parse_args()
HERE=Path(__file__).resolve().parent
candidates=[args.base]if args.base is not None else[HERE.parent,HERE.parent/'session6-periods']
BASE=next((q for q in candidates if(q/'next_thresholds.json').is_file()),None)
if BASE is None:raise FileNotFoundError('Pass --base with next_thresholds.json directory')
data={row['p']:row for row in json.loads((BASE/'next_thresholds.json').read_text())}
result={'status':'PASS','cases':[]}
for p,delta in [(29,42),(23,58),(31,74)]:
    row=data[p];f={tuple(map(int,k.split(','))):v%p for k,v in row['witness'].items()}
    if row['degree']!=delta or not f:raise RuntimeError(('witness metadata',p))
    if any(min(t)<=0 or sum(t)!=delta or v==0 for t,v in f.items()):raise RuntimeError(('interior homogeneous',p))
    if any(f.get((b,c,a),0)!=v for(a,b,c),v in f.items()):raise RuntimeError(('cyclicity',p))
    norm=defaultdict(int)
    for(a,b,c),v in f.items():
        norm[a,b,c]+=v
        for j in range(b+1):norm[j,a+b-j,c]+=v*(-1)**(a+b-j)*comb(b,j)
        for j in range(a+1):norm[j+b,a-j,c]+=v*(-1)**(b+j)*comb(a,j)
    if any(v%p for v in norm.values()):raise RuntimeError(('literal norm',p))
    degrees=[]
    for n in range(delta+1):
        A=matrix_mod(n,'interior',p);rank=rank_forward(A,p);null=A.shape[1]-rank
        if null!=int(n==delta):raise RuntimeError(('minimum',p,n,null))
        degrees.append({'n':n,'columns':A.shape[1],'rank':rank,'nullity':null})
    result['cases'].append({'p':p,'minimum':delta,'witness_terms':len(f),'degrees':degrees})
    print('PASS complete minimum and literal witness',p,delta,'terms',len(f),flush=True)
    with open('period_next_minima_receipt.json','w')as out:json.dump(result,out,indent=2)
