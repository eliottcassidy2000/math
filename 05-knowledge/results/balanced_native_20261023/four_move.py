#!/usr/bin/env python3
"""All original-pool four-slot trades and a rational linear obstruction, if present."""
from itertools import combinations
from collections import defaultdict
from fractions import Fraction as F
from math import lcm,gcd
from functools import reduce
from pathlib import Path
import json,subprocess
def main():
    root=Path(__file__).resolve().parents[3]
    path='0f9dec9c3e1d:05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'
    inv=json.loads(subprocess.check_output(['git','show',path],cwd=root));options={r['source']:r for r in inv['records']}
    sources=sorted(options);labels=[(sum(b<<j for j,b in enumerate(options[s]['source_label'][:4])),options[s]['source_label'][4]) for s in sources]
    if len(set(labels))!=28:raise RuntimeError('two-slot trade omitted: labels repeat')
    bank=defaultdict(list)
    for i,j in combinations(range(28),2):bank[labels[i][0]^labels[j][0],(labels[i][1]+labels[j][1])%128].append((i,j))
    trades=[]
    for bucket in bank.values():
        for a,b in combinations(bucket,2):
            if set(a)&set(b):raise RuntimeError('distinct labels should forbid overlap in pair-sum collision')
            row=[0]*28
            for i in a:row[i]=1
            for i in b:row[i]=-1
            trades.append(row)
    A=[[F(x) for x in row] for row in trades];rank=0;pivots=[]
    for c in range(28):
        p=next((i for i in range(rank,len(A)) if A[i][c]),None)
        if p is None:continue
        A[rank],A[p]=A[p],A[rank];d=A[rank][c];A[rank]=[x/d for x in A[rank]]
        for i in range(len(A)):
            if i!=rank and A[i][c]:
                d=A[i][c];A[i]=[x-d*y for x,y in zip(A[i],A[rank])]
        pivots.append(c);rank+=1
    plus={176,1394,1457,1477,1532,1923,2087,2176,3369};minus={1425,1549,1841,1921,2185,2553,2835,2905,3000}
    delta=[int(s in plus)-int(s in minus) for s in sources];witnesses=[]
    for c in range(28):
        if c in pivots:continue
        v=[F(0)]*28;v[c]=1
        for i,p in enumerate(pivots):v[p]=-A[i][c]
        scale=lcm(*(x.denominator for x in v));v=[int(x*scale) for x in v];div=reduce(gcd,v);v=[x//div for x in v]
        if any(sum(x*y for x,y in zip(v,row)) for row in trades):raise RuntimeError('invariant failed actual trade')
        change=sum(x*y for x,y in zip(v,delta))
        if change:witnesses.append({'coefficients':v,'endpoint_change':change})
    witnesses.sort(key=lambda w:(max(map(abs,w['coefficients'])),sum(map(abs,w['coefficients'])),w['coefficients']))
    # Found by an optional integer optimizer; accepted ONLY by exact integer checks here.
    small=[1,0,3,0,2,-1,1,0,3,-3,-2,0,0,3,0,0,-2,0,-1,-2,0,1,0,-2,-3,3,-1,-3]
    if any(sum(a*b for a,b in zip(small,row)) for row in trades):raise RuntimeError('small invariant violates a complete four-slot trade')
    endpoint_values=[sum(v for s,v in zip(sources,small) if s in E) for E in (minus,plus)]
    if endpoint_values!=[-3,-1]:raise RuntimeError('small invariant does not separate marked endpoints')
    inherited=json.loads(subprocess.check_output(['git','show','a793424e20:05-knowledge/results/native_packet_decomposition_20261023/native/four_move_receipt.json'],cwd=root))
    inherited_values=[sum(v for s,v in zip(inherited['sources'],inherited['coefficients']) if s in E) for E in (minus,plus)]
    if inherited_values!=[-6,-6]:raise RuntimeError('inherited single invariant comparison')
    return {'status':'FINITE-EXACT complete four-slot trade universe','sources':sources,'trades':trades,'rank':rank,
      'invariant_dimension':28-rank,'separating_invariant':witnesses[0] if witnesses else None,
      'small_exact_invariant':small,'small_endpoint_values':endpoint_values,'inherited_invariant_values':inherited_values,
      'scope':'all paired-charge-safe four-slot trades in the complete28-source binary supply, no one-active restriction'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
