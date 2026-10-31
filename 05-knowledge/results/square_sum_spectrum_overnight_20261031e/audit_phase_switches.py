#!/usr/bin/env python3
"""Independent permutation consumer of the finite parity-preserving switch bank."""
import json
from math import isqrt
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'session124_mixed_radix_square_rays_20261031'))
from core import Constructor
bank=json.loads((HERE/'phase_switches.json').read_text())
first=bank['first'];last=bank['last']
if first!=112 or last!=2811 or sorted(map(int,bank['moves']))!=list(range(first,last+1)):
    raise RuntimeError('complete declared base interval')
c=Constructor();paths=0;vertices=0;types={}
for n in range(first,last+1):
    pair=c.pair(n);changed=[]
    for i,(path,move) in enumerate(zip(pair,bank['moves'][str(n)])):
        length=n+i
        if move[0]=='reverse':
            _,a,b=move
            if not(2<=a<b<length-1 and (b-a)%2==0):raise RuntimeError('odd reversal domain')
            order=list(range(a))+list(range(b,a-1,-1))+list(range(b+1,length))
        else:
            if move[0]!='exchange':raise RuntimeError('move type')
            _,a,b,z=move
            if not(1<=a<b<z<length-1 and (b-a)%2==0 and (z-b)%2==0):raise RuntimeError('even exchange domain')
            order=list(range(a+1))+list(range(b+1,z+1))+list(range(a+1,b+1))+list(range(z+1,length))
        if sorted(order)!=list(range(length)) or any(j%2!=v%2 for j,v in enumerate(order)):
            raise RuntimeError('individual positional parity')
        out=[path[j] for j in order]
        if out==path or out[:2]!=[1,3] or out[-1]!=(8 if length%2==0 else 10):raise RuntimeError('nontrivial typed path')
        if any(isqrt(u+v)**2!=u+v for u,v in zip(out,out[1:])):raise RuntimeError('new edge')
        changed.append(out);paths+=1;vertices+=length;types[move[0]]=types.get(move[0],0)+1
    colors={v:j%2 for j,v in enumerate(changed[0])}
    if any(colors[v]!=j%2 for j,v in enumerate(changed[1]) if v<=n):raise RuntimeError('weak pair phase')
print(json.dumps({'status':'passed','first':first,'last':last,'paths':paths,'vertices_checked':vertices,'move_types':types,
                  'scope':'All2700 base pairs; each changed member is distinct and preserves every source-label position parity.'},indent=2,sort_keys=True))
