#!/usr/bin/env python3
"""Differential control of root residual-antichain solver against subset enumeration.
Loads only the callable module; never runs its full n7 main or writes producer files.
"""
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'audit/anchors.py'
spec=importlib.util.spec_from_file_location('root_anchor_audited', SOURCE)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
checks=0

def union(xs):
    result=0
    for x in xs: result |= x
    return result

def reference(cols,target,budget):
    return any(union(cc)&target==target
               for k in range(min(budget,len(cols))+1)
               for cc in combinations(cols,k))

def check(cols,target,budget):
    global checks
    got=module.can_cover(cols,target,budget)[0]
    expected=reference(cols,target,budget)
    if got != expected:
        raise RuntimeError((cols,target,budget,got,expected))
    checks+=1

for size in range(6):
    for cols in combinations(range(1,16),size):
        target=union(cols)
        for budget in range(5): check(cols,target,budget)
extra=[((),0), ((),1), ((0,),0), ((0,),1), ((1,),2), ((7,5),2),
       ((1,1,0),1), ((1,2,4,8,15),15), ((5,10,3,12),15),
       ((12,3,10,5),15), ((1<<17,(1<<17)|4,4), (1<<17)|4)]
for cols,target in extra:
    for budget in range(5):check(cols,target,budget)
print(json.dumps({'status':'PASS: differential finite controls',
                  'checks':checks,'families_four_bit':4944,'budgets_per_family':5,
                  'extra_families':len(extra),
                  'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
                  'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},
                 sort_keys=True,indent=2))
