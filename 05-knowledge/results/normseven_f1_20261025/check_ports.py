#!/usr/bin/env python3
"""Exact joint-fibre transport; counting marginals and adjacency are separate."""
import hashlib
import itertools as it
import json
from collections import Counter
from pathlib import Path


def shear(x):
    return x[0],x[1],x[0]^x[2]


def marginal(table, positions):
    return {key:sum(1 for x in table if tuple(x[i] for i in positions)==key)
            for key in it.product(range(2),repeat=len(positions))}


def main():
    cube=list(it.product(range(2),repeat=3));checks=0
    for x in cube:
        if shear(shear(x))!=x or sum(x)%2!=(shear(x)[1]^shear(x)[2]):
            raise RuntimeError('shear inverse or character map')
        checks+=2
    parity=[[x for x in cube if sum(x)%2==e] for e in range(2)]
    for size in range(3):
        for subset in it.combinations(range(3),size):
            if marginal(parity[0],subset)!=marginal(parity[1],subset):
                raise RuntimeError('old proper marginals differ')
            checks+=1
    new=[[shear(x) for x in table] for table in parity]
    if marginal(new[0],(1,2))==marginal(new[1],(1,2)):
        raise RuntimeError('new pair query failed to distinguish')
    checks+=1
    if sum(shear((1,0,0)))!=2:
        raise RuntimeError('cube adjacency hostile')
    checks+=1
    for epsilon in range(2):
        decoded=[(y[0],y[1],y[0]^y[2]) for y in new[epsilon]]
        if Counter(decoded)!=Counter(parity[epsilon]):
            raise RuntimeError('old ports not restored')
        checks+=1
    print(json.dumps({'status':'FINITE-EXACT','checks':checks,'states':8,
                      'old_parity_relations':parity,'new_relations':new,
                      'new_pair_counts':[{''.join(map(str,k)):v for k,v in marginal(t,(1,2)).items()} for t in new],
                      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
