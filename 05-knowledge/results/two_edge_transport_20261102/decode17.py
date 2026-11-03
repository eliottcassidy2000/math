"""Construct a count-17 tournament from any six-vertex orientation, in <=3 flips.

The radius-three guarantee is FINITE-EXACT for the complete 15-bit cube.
This is not an all-order or arbitrary-target inverse.
"""
from collections import deque
from itertools import combinations
import argparse
import json
from scout import rows,h_count

def decode(mask,max_error=None):
    if not 0<=mask<1<<15:raise ValueError('source must encode six vertices in fifteen bits')
    q=deque([mask]);previous={mask:None};counts={mask:h_count(rows(6,mask))}
    if max_error is not None and (counts[mask]-17)**2>max_error:raise ValueError('source exceeds requested barrier')
    while q:
        x=q.popleft()
        if counts[x]==17:
            path=[]
            while x is not None:path.append(x);x=previous[x]
            path.reverse();edges=list(combinations(range(6),2))
            return {'source':mask,'masks':path,'counts':[counts[x] for x in path],
                    'edge_flips':[edges[(a^b).bit_length()-1] for a,b in zip(path,path[1:])],
                    'target_out_rows':rows(6,path[-1])}
        for e in range(15):
            y=x^(1<<e)
            if y not in previous:
                value=h_count(rows(6,y))
                if max_error is None or (value-17)**2<=max_error:
                    previous[y]=x;counts[y]=value;q.append(y)
    raise RuntimeError('no target in the requested source component')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mask',type=int);parser.add_argument('--max-error',type=int);args=parser.parse_args()
    print(json.dumps(decode(args.mask,args.max_error),indent=2))
