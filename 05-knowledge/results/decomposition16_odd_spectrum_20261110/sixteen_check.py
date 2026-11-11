#!/usr/bin/env python3
"""Independent all64-matrix audit of the incoming five-bin candidate."""
from collections import Counter
from itertools import permutations
from hashlib import sha256
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'05-knowledge/results/session250_sixteen_20261110/geometry/DIGRAPH_CANDIDATE_SIXTEEN.md'

def require(x,msg):
    if not x:raise RuntimeError(msg)

def main():
    pairs=[(i,j) for i in range(3) for j in range(3) if i!=j]
    P=list(permutations(range(3)))
    classes={}
    for mask in range(64):
        edges={p for bit,p in enumerate(pairs) if mask>>bit&1}
        key=min(tuple(int((p[i],p[j]) in edges) for i,j in pairs) for p in P)
        H=sum(all((p[i],p[i+1]) in edges for i in range(2)) for p in P)
        tournament=all(int((i,j) in edges)+int((j,i) in edges)==1
                       for i in range(3) for j in range(i+1,3))
        value=(len(edges),H,tournament)
        require(key not in classes or classes[key]==value,'isomorphism invariance')
        classes[key]=value
    fine=Counter(e for e,H,t in classes.values())
    group=lambda e:0 if e==0 else 1 if e<3 else 2 if e==3 else 3 if e<6 else 4
    coarse=Counter(group(e) for e,H,t in classes.values())
    middle=Counter(H for e,H,t in classes.values() if e==3)
    require([fine[e] for e in range(7)]==[1,1,4,4,4,1,1],'fine grading')
    require([coarse[k] for k in range(5)]==[1,5,4,5,1],'coarse grading')
    require(middle=={1:3,3:1},'middle H loss')
    require(sum(t for e,H,t in classes.values())==2,'tournament middle subset')
    require(sum(e==3 and not t for e,H,t in classes.values())==2,'nontournament middle')
    # A single transposition swaps the six arc positions in three pairs;
    # a three-cycle produces two cycles of length3. Direct polynomial coefficients.
    from math import comb
    burnside=[(comb(6,e)+3*(comb(3,e//2) if e%2==0 else 0)
               +2*(comb(2,e//3) if e%3==0 else 0))//6 for e in range(7)]
    require(burnside==[fine[e] for e in range(7)],'independent Burnside grading')
    print(json.dumps({'status':'PASS','universe':'all64 loopless binary adjacency matrices of order3',
        'ordinary_isomorphism_classes':len(classes),'fine_arc_histogram':[fine[e] for e in range(7)],
        'five_bin_histogram':[coarse[k] for k in range(5)],'middle_H_histogram':dict(middle),
        'tournament_classes':2,'nontournament_middle_classes':2,
        'historical_status':'incoming candidate reconstruction, original referent still unconfirmed',
        'source_sha256':{str(SOURCE.relative_to(ROOT)):sha256(SOURCE.read_bytes()).hexdigest()},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':main()
