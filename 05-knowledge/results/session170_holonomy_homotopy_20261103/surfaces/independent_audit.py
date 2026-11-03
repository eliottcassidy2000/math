import importlib.util
from functools import lru_cache
from collections import deque
from pathlib import Path
spec=importlib.util.spec_from_file_location('subject',str(Path(__file__).with_name('check.py')))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
@lru_cache(None)
def direct(n,a,x):
    if n==0:return 0
    h=1<<(n-1); u,A=a%h,a//h; y,t=x%h,x//h
    if not A and not t:return direct(n-1,u,y)
    if not A and t:return direct(n-1,y,u)
    if A and not t:return direct(n-1,u,y)^int(y!=0)
    return 1^int(y!=0)^direct(n-1,y,u)
for n in range(3,9):
    if any(m.twist(n,a,x)!=direct(n,a,x) for a in range(1<<n) for x in range(1<<n)):
        raise RuntimeError('Cayley-Dickson recursion mismatch')
print('Independent direct four-case Cayley-Dickson recursion agrees at every basis product for n=3..8.')
def closure(gs):
    id=tuple(range(len(gs[0]))); seen={id}; q=deque([id])
    while q:
        p=q.popleft()
        for g in gs:
            z=tuple(g[p[x]] for x in range(len(p)))
            if z not in seen:seen.add(z);q.append(z)
    return seen
for n in (3,4):
    actual=closure([m.left(n,a) for a in range(1,1<<n)])
    paired=closure([g for pair in m.surface_generators(n) for g in pair])
    if actual!=paired:raise RuntimeError('Different permutation subgroups')
    print('n=',n,'actual all-basis left group equals paired subgroup, order=',len(actual))
