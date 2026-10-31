#!/usr/bin/env python3
"""Explicit parity-preserving local alternatives in the fixed weak-B base bank."""
import json
from math import isqrt
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'session124_mixed_radix_square_rays_20261031'))
from core import Constructor


def square(x):return isqrt(x)**2==x


def find(path):
    n=len(path)
    for i in range(2,n-3):
        for j in range(i+2,n-1,2):
            if square(path[i-1]+path[j]) and square(path[i]+path[j+1]):
                return ['reverse',i,j]
    for a in range(1,n-4):
        for b in range(a+2,n-2,2):
            if not square(path[a]+path[b+1]):continue
            for c in range(b+2,n-1,2):
                if square(path[c]+path[a+1]) and square(path[b]+path[c+1]):
                    return ['exchange',a,b,c]
    return None


def apply(path,move):
    if move[0]=='reverse':
        _,i,j=move
        return path[:i]+path[i:j+1][::-1]+path[j+1:]
    _,a,b,c=move
    return path[:a+1]+path[b+1:c+1]+path[a+1:b+1]+path[c+1:]


def main():
    constructor=Constructor();records={};failures=[];checks=0
    for n in range(42,1062):
        pair=constructor.bases.base(n);moves=[find(p) for p in pair]
        if None in moves:
            failures.append({'n':n,'moves':moves})
            continue
        changed=[apply(p,m) for p,m in zip(pair,moves)]
        constructor.bases.check_pair(changed)
        for p,q in zip(pair,changed):
            colors={v:j%2 for j,v in enumerate(p)}
            if p==q or any(colors[v]!=j%2 for j,v in enumerate(q)):raise RuntimeError('switch changes individual phase or nothing')
            checks+=1
        records[str(n)]=moves
    last=max(r['n'] for r in failures)
    bank={k:v for k,v in records.items() if int(k)>last}
    out={'status':'passed','first_uniform_base_index':last+1,'finite_bank_end':1061,
         'moves':bank,'failed_local_move_tests':failures,'individually_checked_switches':checks,
         'scope':'Failures concern only these two move types in the supplied base paths, not all weak pairs.'}
    (HERE/'leaf_switches.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('moves','failed_local_move_tests')},sort_keys=True))


if __name__=='__main__':main()
