#!/usr/bin/env python3
"""Finite controls for the elementary affine-gauge proof; no assert checks."""
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def collatz(x):
    return 3*x+1 if x % 2 else x//2

def main():
    branches = 0
    for delta in (0, 2, 4):
        a = delta+2
        for x in range(1, 4097):
            label = a*x+1
            expected = 3*label+delta if x % 2 else (label+1)//2
            require(a*collatz(x)+1 == expected, ('branch', delta, x))
            branches += 1
        candidates = [(aa, b) for aa in range(1, 21) for b in range(-5, 6)
                      if 2*(aa+b) == 2*aa+b+1
                      and 4*aa+b == 3*(aa+b)+delta]
        require(candidates == [(a, 1)], ('coefficient equations', delta, candidates))
    walls = [(0, 6, 13, 7), (0, 20, 41, 21),
             (2, 10, 41, 21), (4, 2, 13, 7)]
    for delta, x, before, after in walls:
        require((delta+2)*x+1 == before, 'source wall')
        require((delta+2)*collatz(x)+1 == after, 'destination wall')
    forbidden = [(0,3,7), (0,10,21), (2,5,21), (4,1,7)]
    for delta,x,value in forbidden:
        require((delta+2)*x+1 == value, 'forbidden preimage')
    source = '01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md'
    print(json.dumps({
        'status':'PASS; FINITE-EXACT arithmetic controls',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256':{source:sha256((ROOT/source).read_bytes()).hexdigest()},
        'branch_checks':branches, 'walls':walls, 'forbidden_preimages':forbidden,
        'scope':'Affine gauges only; no convergence or spectrum completeness claim'
    }, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
