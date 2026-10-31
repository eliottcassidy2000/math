#!/usr/bin/env python3
"""Finite switches closing a uniform SAME-PHASE two-choice supply."""
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031d'))
from leaf_switches import Constructor,find,apply

def main():
    c=Constructor();bank={};failures=[];types={};vertices=0
    for n in range(112,2812):
        pair=c.pair(n);moves=[find(p) for p in pair]
        if None in moves:
            failures.append([n,moves]);continue
        other=[apply(p,m) for p,m in zip(pair,moves)]
        c.bases.check_pair(other)
        for p,q,m in zip(pair,other,moves):
            colors={v:j%2 for j,v in enumerate(p)}
            if p==q or any(colors[v]!=j%2 for j,v in enumerate(q)):
                raise RuntimeError('individual phase or distinctness')
            types[m[0]]=types.get(m[0],0)+1;vertices+=len(q)
        bank[str(n)]=moves
    out={'status':'passed' if not failures else 'failed','first':112,'last':2811,
         'moves':bank,'failures':failures,'move_types':types,'vertices_checked':vertices,
         'source':'Constructor.pair(n, empty word), including its A refinements above1061'}
    (HERE/'phase_switches.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='moves'},sort_keys=True,indent=2))
    if failures:raise RuntimeError('incomplete same-phase base interval')

if __name__=='__main__':main()
