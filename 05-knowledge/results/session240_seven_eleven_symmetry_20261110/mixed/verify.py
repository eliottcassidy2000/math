#!/usr/bin/env python3
"""Exact mixed 7+11+1 kernel enumeration from replayed finite core inputs."""
import json
from pathlib import Path
from collections import Counter

def require(q,msg):
    if not q: raise RuntimeError(msg)

base=Path(__file__).resolve().parent.parent
raw=json.loads((base/'gaps/probe.json').read_text())
data=raw.get('primes',raw)
r7=data['7']['rows'];r11=data['11']['rows']
hist=Counter(); best=[]
for a in r7:
    for b in r11:
        h=3*a['H']*b['H']+a['F2']*b['H']+a['H']*b['F2']+a['F2']*b['F2']
        require(h%77==0,'action divisibility')
        hist[h]+=2
        best.append((h,a['mask'],b['mask']))
m=min(hist)
require(m==3049061631,'minimum')
require(hist[m]==120,'minimizer count')
require(sum(hist.values())==512,'universe')
require(min(r['F2'] for r in r11)==1571526,'eleven F2 minimum')
print(json.dumps({'status':'FINITE-EXACT inherited kernel; independently audit vertex DP',
 'universe':'all orientations fixed by specified (7)(11)(1), filtered strong',
 'strong_orientations':sum(hist.values()),'minimum_order':19,'minimum_H':m,
 'minimum_H_over_77':m//77,'minimizing_orientations':hist[m],
 'histogram':dict(sorted(hist.items())),
 'minimizing_pairs':[[a,b] for h,a,b in best if h==m]},indent=2,sort_keys=True))
