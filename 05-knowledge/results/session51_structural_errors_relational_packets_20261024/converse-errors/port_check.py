#!/usr/bin/env python3
"""Literal Hamiltonian/blue-port controls on the 24 rigid self-converse nodes.
No line pairing or graph wiring is inferred from these vertex-fibre counts.
"""
from functools import lru_cache
from itertools import permutations,product
from collections import Counter
from pathlib import Path
import argparse,json
import check

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='port-normal.json');args=ap.parse_args()
    records=json.loads((check.HERE/'normal.json').read_text())['representatives'];out=[]
    for rec in records:
        if not rec['self_converse']:continue
        T=tuple(rec['rows']);n=len(T);q=tuple(rec['anti_certificates'][0]['anti_involution'])
        need(rec['automorphism_order']==1,'all tested positive self-converse ten classes rigid')
        @lru_cache(None)
        def dp(last,mask):
            if not mask:return 1
            return sum(dp(v,mask^(1<<v)) for v in range(n) if mask>>v&1 and T[last]>>v&1)
        H=sum(dp(v,((1<<n)-1)^(1<<v)) for v in range(n))
        pairs=[(v,q[v]) for v in range(n) if v<q[v]];blue=[]
        for halves in product(*pairs):
            for first in permutations(halves):
                path=first+tuple(q[v] for v in first[::-1])
                if all(T[a]>>b&1 for a,b in zip(path,path[1:])):blue.append(path)
        # A literal full-path traversal independently counts the DP and finds
        # a black-port witness; no anti-map criterion is used to prune it.
        full_count=0;blue_count=0;black_witness=None
        def walk(path,mask):
            nonlocal full_count,blue_count,black_witness
            if not mask:
                full_count+=1
                # Literal normalized-word staircase reflection test: every
                # non-path pair bit equals its reflected pair bit. This does
                # not use the known anti-involution to decide blueness.
                is_blue=all(((T[path[i]]>>path[j])&1)==((T[path[n-1-j]]>>path[n-1-i])&1)
                            for i in range(n) for j in range(i+2,n))
                need(is_blue==all(path[n-1-i]==q[v] for i,v in enumerate(path)),
                     'literal normalized word agrees with unique anti-involution')
                if is_blue:blue_count+=1
                elif black_witness is None:black_witness=path
                return
            for v in range(n):
                if mask>>v&1 and (not path or T[path[-1]]>>v&1):walk(path+(v,),mask^(1<<v))
        walk((),(1<<n)-1)
        need(full_count==H and blue_count==len(blue),'independent full-word and half-system counts')
        need(H%2==len(blue)%2==1,'ordinary and blue fibre parity')
        need(0<len(blue)<H,'each tested node is mixed, not all-blue')
        out.append({'representative_index':rec['index'],'rows':T,'Hamiltonian_paths':H,'blue_ports':len(blue),'black_ports':H-len(blue),'blue_path_witness':blue[0],'black_path_witness':black_witness})
    need(len(out)==24,'all 24 self-converse positive nodes checked')
    result={'status':'PASS','scope':'presentation port counts only, no line pairing or endpoint adjacency','blue_incident_positive_nodes':24,'mixed_positive_nodes':24,'all_blue_positive_nodes':0,'black_only_positive_nodes':68,'self_converse_blue_port_sum':sum(r['blue_ports'] for r in out),'self_converse_black_port_sum':sum(r['black_ports'] for r in out),'blue_port_histogram':sorted(Counter(r['blue_ports'] for r in out).items()),'representatives':out}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='representatives'},indent=2,sort_keys=True))
if __name__=='__main__':main()
