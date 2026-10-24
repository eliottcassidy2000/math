#!/usr/bin/env python3
"""Two independent dynamic programs for ordinary and anti-palindromic path fibres."""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import json
import independent as G
HERE=Path(__file__).resolve().parent

def run():
    report=[]
    for rec in json.loads((HERE/'normal.json').read_text())['order_ten']['representatives']:
        if not rec['anti_automorphism_order']:continue
        a=tuple(rec['rows']);n=len(a);q=G.maps(a,G.reverse(a))[0]
        G.need(rec['automorphism_order']==rec['anti_automorphism_order']==1,'rigid unique anti-map for fibre controls')
        pairs=[(v,q[v]) for v in range(n) if v<q[v]];pair_of={v:i for i,p in enumerate(pairs) for v in p}
        @lru_cache(None)
        def h(mask,last):
            if not mask:return 1
            return sum(h(mask^(1<<v),v) for v in range(n) if (mask>>v)&1 and G.arc(a,last,v))
        @lru_cache(None)
        def blue(mask,last):
            if not mask:return G.arc(a,last,q[last])
            total=0
            for i,p in enumerate(pairs):
                if not (mask>>i)&1:continue
                for v in p:
                    if G.arc(a,last,v):total+=blue(mask^(1<<i),v)
            return total
        H=sum(h(((1<<n)-1)^(1<<v),v) for v in range(n))
        B=sum(blue(((1<<len(pairs))-1)^(1<<pair_of[v]),v) for v in range(n))
        G.need(0<B<H and B%2==H%2==1,'each rigid SC error node has both anti-palindromic and other path fibres')
        report.append({'rows':a,'hamiltonian_paths':H,'blue_ports':B,'black_ports':H-B})
    G.need(len(report)==24,'all self-converse positive-error fibre tests')
    return {'status':'PASS','scope':'24 rigid self-converse class path fibres, not line pairing or metagraph adjacency','blue_port_sum':sum(r['blue_ports'] for r in report),'black_port_sum':sum(r['black_ports'] for r in report),'blue_histogram':sorted(Counter(r['blue_ports'] for r in report).items()),'records':report}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
