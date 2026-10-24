#!/usr/bin/env python3
"""Independent matrix-orbit and fixed-position derivation for the two C5 cells."""
from collections import Counter
from itertools import permutations,product
import json,argparse
from pathlib import Path

def need(ok,label):
    if not ok:raise RuntimeError(label)
def translate(M,a,b):return tuple(tuple(M[(i+a)%5][(j+b)%5] for j in range(5)) for i in range(5))
def twist(M):return tuple(tuple(M[-j%5][-i%5] for j in range(5)) for i in range(5))
def canonical(M):return min(translate(M,a,b) for a,b in product(range(5),repeat=2))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='matrix-normal.json');args=ap.parse_args()
    P=list(permutations(range(5)));universes={1:{tuple(tuple(int(j==p[i]) for j in range(5)) for i in range(5)) for p in P},2:set()}
    for p in P:
        for q in P:
            if any(p[i]==q[i] for i in range(5)):continue
            universes[2].add(tuple(tuple(int(j==p[i] or j==q[i]) for j in range(5)) for i in range(5)))
    report={}
    for d,U in universes.items():
        need(all(all(sum(M[i][j] for i in range(5))==d for j in range(5)) for M in U),'column regularity')
        orbits={}
        for M in U:orbits.setdefault(canonical(M),set()).add(M)
        symmetric=sum(all(M[i][j]==M[j][i] for i in range(5) for j in range(5)) for M in U)
        fixed=Counter();twisted_counts={}
        for M,orb in orbits.items():
            need(len(orb) in (5,25),'nonuniform orbit size five or twenty-five')
            is_sc=canonical(twist(M))==M;fixed[(len(orb),is_sc)]+=1
        for a,b in product(range(5),repeat=2):
            twisted_counts[(a,b)]=sum(translate(twist(M),a,b)==M for M in U)
        need(all(v==symmetric for key,v in twisted_counts.items() if key[0]==key[1]),'diagonal twisted maps fix symmetric count')
        non_diag={v for key,v in twisted_counts.items() if key[0]!=key[1]}
        need(non_diag==({1} if d==1 else {2}),'off-diagonal twisted fixed counts')
        fixed_orbits=sum(c for (size,is_sc),c in fixed.items() if is_sc)
        need(sum(twisted_counts.values())==25*fixed_orbits,'twisted Burnside full averaging')
        report[d]={'labelled_matrices':len(U),'symmetric_matrices':symmetric,'ordinary_classes':len(orbits),'self_converse_classes':fixed_orbits,'merged_classes':(len(orbits)+fixed_orbits)//2,'class_histogram_orbit_size_selfconverse':[[list(k),v] for k,v in sorted(fixed.items())],'twisted_fixed_counts':{str(k):v for k,v in twisted_counts.items()}}
    need([report[d]['self_converse_classes'] for d in (1,2)]==[6,24],'complete self-converse counts')
    result={'status':'PASS','sectors':report};Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
