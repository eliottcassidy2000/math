#!/usr/bin/env python3
"""Independent certificate checks: cut/rejoin lower bound and length-one fibers.
Imports only the full-universe enumerators; uses different exchange logic.
"""
from graceful_square_exchange_20261019 import square_paths, graceful_trees, check
from itertools import combinations, permutations
from collections import defaultdict
from math import isqrt


def square(x,y): return isqrt(x+y)**2 == x+y

def rejoin_exists(parts):
    def visit(last,remaining):
        if not remaining:return True
        for i in remaining:
            for a,b in {(parts[i][0],parts[i][-1]),(parts[i][-1],parts[i][0])}:
                if last is None or square(last,a):
                    if visit(b,remaining-{i}):return True
        return False
    return visit(None,set(range(len(parts))))

paths=square_paths(28); check(len(paths)==52,'source universe')
checks=0
for p in paths:
    for k in range(3):
        for cuts in combinations(range(1,28),k):
            borders=(0,)+cuts+(28,)
            parts=[p[a:b] for a,b in zip(borders,borders[1:])]+[(29,)]
            check(not rejoin_exists(parts),'unexpected <=2-cut square extension')
            checks+=1
print('Independent cut-and-rejoin rejection:',checks,'source/cut sets; all 52 Q28 paths; 0,1,2 cuts.')
print('Every orientation and ordering of the remaining path segments and singleton29 is searched.')
for n in range(3,10):
    codes=graceful_trees(n); fibers=defaultdict(list)
    for t in codes:fibers[t[1:]].append(t[0])
    check(all(len(v)>=2 and len(v)%2==0 for v in fibers.values()),'length-one even fiber')
    print('Length-one fibers:',n,'trees',len(codes),'fibers',len(fibers),'all positive even cliques')
# Four K4 stars are exactly the rainbow spanning trees of its proper1factorization.
color={(0,1):0,(2,3):0,(0,2):1,(1,3):1,(0,3):2,(1,2):2}
rainbow=[]
for ee in combinations(color,3):
    if len({color[e] for e in ee})<3:continue
    reached={0}
    while True:
        rr=reached|{v for e in ee if set(e)&reached for v in e}
        if rr==reached:break
        reached=rr
    if len(reached)==4:rainbow.append(frozenset(ee))
check(len(rainbow)==4,'K4 count')
check(all(len(a-b)>1 for a,b in combinations(rainbow,2)),'K4 isolation')
print('General-rainbow hostile: proper1factorized K4 has exactly4 isolated rainbow spanning trees.')
print('ALL AUDIT CHECKS PASSED.')
