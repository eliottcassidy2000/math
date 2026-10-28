#!/usr/bin/env python3
"""Direct arc tests and sparse Q-square tests, separately from classification."""
from collections import Counter
import json
import check_blockfill as b


def main():
    parents=[()]; out=[]
    for n in range(8):
        r=b.build(parents,n,'blockfill'); N=n+1; maps_checked=0; entries=0
        for target,child in r['states'].items():
            row=r['Q'][target]; square=Counter()
            for mid,w in row:
                for dst,v in r['Q'][mid]: square[dst]+=w*v
            b.need(dict(square)=={dst:N*w for dst,w in row},'direct sparse Q square fails')
            for vertex,(source,maps) in enumerate(r['root_transports'][target]):
                origin=r['states'][source]
                b.need(maps,'empty transporter fiber')
                for p in maps:
                    b.need(sorted(p)==list(range(N)) and p[-1]==vertex,'bad port bijection')
                    # Deliberately do not call inherited encode/move_mask.
                    for u in range(N):
                        for v in range(N):
                            b.need(((origin[u]>>v)&1)==((child[p[u]]>>p[v])&1),
                                   'transporter fails direct adjacency test')
                    maps_checked+=1
            entries+=len(row)
        out.append({'order':N,'pointed_states':len(r['states']),
                    'nonzero_Q_entries':entries,'full_arc_preserving_maps_checked':maps_checked,
                    'sparse_Q_square_checked_for_every_state':True})
        parents=list(r['bank'])
    return {'status':'FINITE-EXACT','checks':out}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
