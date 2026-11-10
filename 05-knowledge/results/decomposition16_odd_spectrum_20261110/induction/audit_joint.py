#!/usr/bin/env python3
"""Independent DP and three-root marked-class controls for JOINT_CORES.md."""
from collections import Counter, defaultdict
from functools import reduce
from hashlib import sha256
from itertools import permutations
from math import comb, gcd
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_one_root', HERE/'check.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def require(x, reason):
    if not x:
        raise RuntimeError(reason)

def fixed_key(rows, keep, marks):
    rest = tuple(v for v in keep if v not in marks)
    return min(m.code_at(rows, tuple(marks)+q) for q in permutations(rest))

def main():
    actual = defaultdict(Counter); kernels = defaultdict(set)
    comparisons = 0; labelled = 0
    for n in range(3,6):
        marked_classes = {}
        for code in range(1 << (n*(n-1)//2)):
            rows=m.rows_from_code(n,code); V=tuple(range(n)); roots=(0,1,2)
            cube=tuple(m.hdp(m.induced(rows,tuple(v for v in V if v>=3 or not(mask>>v&1)))) for mask in range(8))
            parts=m.sccs(rows); touched=[c for c in parts if set(c)&set(roots)]
            outside=[c for c in parts if not(set(c)&set(roots))]
            factor=1
            for c in outside: factor*=m.hdp(m.induced(rows,c))
            d=reduce(gcd,cube)
            require(d%factor==0,'three-root outside gcd')
            if d==1:
                require(all(len(c)==1 for c in outside),'three-root primitive kernel')
                keep=tuple(sorted(v for c in touched for v in c)); k=len(touched)
                require(5**(len(keep)-k)<=cube[0]**3,'three-root integer core bound')
                key=fixed_key(rows,V,roots)
                val=(cube,fixed_key(rows,keep,roots),len(keep),k)
                require(key not in marked_classes or marked_classes[key]==val,'marked canonical consistency')
                marked_classes[key]=val
                for mask in range(8):
                    vertices=tuple(v for v in keep if v>=3 or not(mask>>v&1))
                    require(m.hdp(m.induced(rows,vertices))==cube[mask],'entire cube preserved')
                    comparisons+=1
            labelled+=1
        for cube,key,b,k in marked_classes.values():
            actual[cube][n]+=1
            kernels[cube].add((key,b,k))
    for cube, observed in actual.items():
        predicted=Counter()
        for key,b,k in kernels[cube]:
            for n in range(b,6): predicted[n]+=comb(n-b+k,k)
        require(predicted==observed,('three-root rational series',cube))
    hostile=(16,9,3,5,14)
    cube=tuple(m.hdp(m.induced(hostile,tuple(v for v in range(5) if v>=2 or not(mask>>v&1)))) for mask in range(4))
    require(cube==(15,3,5,1) and len(m.sccs(hostile))==1 and cube[0]*cube[3]==cube[1]*cube[2], 'strong zero determinant witness')
    files=[HERE/'../JOINT_CORES.md',HERE/'../joint_check.py',HERE/'check.py']
    print(json.dumps({'status':'PASS','independence':'subset DP and Floyd reachability versus reviewed literal permutations and DFS reachability',
        'universe':'all labelled tournaments n3..5 with roots0,1,2; marked quotient fixes all three names',
        'labelled_three_root_graphs':labelled,'deletion_cube_entry_checks':comparisons,
        'three_root_marked_classes':sum(sum(p.values()) for p in actual.values()),'three_root_rational_series':len(actual),
        'strong_two_root_zero_determinant_hostile':{'rows':hostile,'cube':cube},
        'reviewed_sha256':{str(p.resolve().relative_to(HERE.parent)):sha256(p.read_bytes()).hexdigest() for p in files},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':main()
