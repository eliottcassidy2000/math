#!/usr/bin/env python3
"""Anonymous [7] rows forget actual pointed keys, groups and joint ports."""
from collections import Counter
from itertools import permutations
import json
from baseline_sources import generation as old
from check_blockfill import complete_row, inverse, compose


def cyclic(steps):
    return tuple(sum(1 << ((v+s)%7) for s in steps) for v in range(7))


def hamilton_paths(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            for w in range(n):
                if (rows[v]>>w)&1 and not ((mask>>w)&1):
                    dp[mask|(1<<w)][w]+=dp[mask][v]
    return sum(dp[-1])


def pair_orbits(group,n):
    unseen={(i,j) for i in range(n) for j in range(n) if i!=j}; sizes=[]
    while unseen:
        i,j=min(unseen)
        orbit={(p[i],p[j]) for p in group}
        old.need(orbit <= unseen,'pair orbits overlap')
        unseen-=orbit; sizes.append(len(orbit))
    return sorted(sizes)


def main():
    out=[]
    for steps in ((1,2,3),(1,2,4)):
        rows=cyclic(steps); engine=old.ParentCanonicalizer()
        row, records=complete_row(engine,rows)
        code, maps=old.ParentCanonicalizer().maps(rows)
        base=inverse(maps[0]); group=tuple(sorted(compose(p,base) for p in maps))
        direct=tuple(p for p in permutations(range(7)) if old.encode(rows,p)==old.encode(rows,tuple(range(7))))
        old.need(group==direct,'independent permutation aut mismatch')
        old.need(len(row)==1 and row[0][1]==7,'not anonymous [7]')
        out.append({'steps':steps, 'canonical_code':code,'named_row':row,'anonymous_row':[7],
                    'H':hamilton_paths(rows),'automorphism_order':len(group),
                    'ordered_distinct_pair_orbit_sizes':pair_orbits(group,7)})
    old.need(out[0]['named_row']!=out[1]['named_row'],'named dictionaries unexpectedly equal')
    old.need(out[0]['canonical_code']!=out[1]['canonical_code'],'hostiles are isomorphic')
    return {'status':'FINITE-EXACT','same_anonymous_row_different_named_states':out,
            'interpretation':'A complete row identifies its block only in the common named pointed-state dictionary. Anonymous orbit-weight profiles do not recover the child, group action, or pair ports.'}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
