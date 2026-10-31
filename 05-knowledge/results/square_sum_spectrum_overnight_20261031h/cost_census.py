#!/usr/bin/env python3
"""Exact positive-deficit costs, and the arbitrary-slope hostile reroute."""
import json,sys
from pathlib import Path
from bounded_sources import search,word
HERE=Path(__file__).resolve().parent
def main():
    Q=int(sys.argv[1]) if len(sys.argv)>1 else 256
    name=f'strong_profiles{Q}.json'
    bank=json.loads((HERE/name).read_text())['profiles']
    prev,t=search(bank,Q,Q//2-1,2*Q)
    positive=[];targets={}
    for state in sorted(prev):
        F,q,g,h=state
        if q and F+q+g<Q and F+q+g not in targets:
            targets[F+q+g]={'state':state,'word':word(prev,state)}
        if q and not q&(q-1) and q>F+g:
            positive.append({'state':state,'c':q-F-g,'word':word(prev,state)})
    print(json.dumps({'states':len(prev),'transitions':t,'profiles':len(bank),'Q':Q,
          'positive_sources':positive,'arbitrary_slope_one_zero_targets':targets,
          'missing_odd_targets':[n for n in range(3,Q,2) if n not in targets]},indent=2,sort_keys=True))
if __name__=='__main__':main()
