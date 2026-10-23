#!/usr/bin/env python3
"""Finite probe only: one auxiliary bit for heterogeneous cube-face menus."""
from collections import deque
from itertools import combinations_with_replacement, product
import json
from native_cube import assignments

def faces(r):
    return [{x for x in range(1<<r) if all(b==2 or (x>>j)&1==b for j,b in enumerate(p))}
            for p in product(range(3),repeat=r)]

def main():
    receipts=[]
    for r,maxs in ((1,2),(2,4),(3,3)):
        menus=faces(r);count=0;visited=0;base=0;maxdist=0
        for s in range(1,maxs+1):
            for indices in combinations_with_replacement(range(len(menus)),s):
                E=[menus[i] for i in indices];ends=assignments(E)
                if not ends:continue
                count+=1;base+=len(ends);start=ends[0];q=deque([start]);d={start:0};mask=(1<<r)-1
                while q:
                    x=q.popleft();occupied=set(x)
                    for i in range(s):
                        for j in range(r+1):
                            z=x[i]^(1<<j)
                            if z in occupied or (z&mask) not in E[i]:continue
                            y=x[:i]+(z,)+x[i+1:]
                            if y not in d:d[y]=d[x]+1;q.append(y)
                missing=[x for x in ends if x not in d]
                if missing:raise RuntimeError(json.dumps(dict(r=r,E=[sorted(e) for e in E],start=start,unreachable=missing[0])))
                visited+=len(d);maxdist=max(maxdist,max(d[x] for x in ends))
        receipts.append(dict(r=r,max_slots=maxs,feasible_menu_multisets=count,
                             boundary_assignments=base,lifted_states_visited=visited,
                             max_distance_from_selected_start=maxdist))
    return dict(status='FINITE-EXACT positive probe; general one-bit conjecture OPEN',
                quotient='slot permutations only; actual labels recover by permutation',cases=receipts)

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
