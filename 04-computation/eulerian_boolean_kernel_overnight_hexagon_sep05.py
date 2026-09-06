#!/usr/bin/env python3
"""Bounded native-kernel probe using complete generator orbits, not samples.

Universe: every labelled Eulerian graph through n=7, quotienting by
adjacent vertex transpositions generating S_n. Every triangle neighbor
is retained; the actual Boolean matrix is distinguished from multiplicity.
"""
from collections import deque
import hashlib
import json
import sympy as sp
import even_graph_triangle_quotient_spectrum_thm4078 as old
import eulerian_uniform_orbit_sampler_overnight_hexagon_sep05 as fixed


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def generator_orbits(n):
    edges,index = old.edge_system(n)
    states = set(old.cycle_space(n,index))
    maps = []
    for v in range(n-1):
        perm = list(range(n))
        perm[v],perm[v+1] = perm[v+1],perm[v]
        maps.append(tuple(index[tuple(sorted((perm[i],perm[j])))] for i,j in edges))
    orbit = {}
    reps,sizes = [],[]
    for start in sorted(states):
        if start in orbit:
            continue
        label = len(reps)
        reps.append(start)
        orbit[start] = label
        queue = deque([start])
        size = 0
        while queue:
            state = queue.popleft()
            size += 1
            for mapping in maps:
                target = old.permute_mask(state,mapping)
                require(target in states, "permutation preserves cycle space")
                if target not in orbit:
                    orbit[target] = label
                    queue.append(target)
                else:
                    require(orbit[target]==label, "generator orbit partition")
        sizes.append(size)
    require(len(reps)==fixed.OrbitSampler(n).orbit_count, "independent Burnside orbit count")
    return reps,orbit,sizes,index


def main():
    rows=[]
    for n in range(3,8):
        reps,orbit,sizes,index = generator_orbits(n)
        if n<=5:
            r2,o2,s2,_ = old.quotient_data(n)
            require((reps,orbit,sizes)==(r2,o2,s2), "literal full-permutation control")
        triangles = old.simple_cycles(n,3,index)
        M = old.weighted_operator(reps,orbit,triangles)
        B = old.boolean_support(M)
        even = [i for i,r in enumerate(reps) if r.bit_count()%2==0]
        odd = [i for i,r in enumerate(reps) if r.bit_count()%2]
        for state,label in orbit.items():
            row = [0]*len(reps)
            for triangle in triangles:
                row[orbit[state ^ triangle]] += 1
            require(row==M[label], "all-labelled-state row and representative row agree")
        require(all(sizes[i]*M[i][j]==sizes[j]*M[j][i] for i in range(len(reps)) for j in range(len(reps))), "reversible multiplicity control")
        block=sp.Matrix([[B[i][j] for j in odd] for i in even])
        rank=block.rank()
        nullity=len(reps)-2*rank
        require(nullity>=len(even)-len(odd), "chiral index lower bound")
        payload={"n":n,"labels":len(orbit),"classes":len(reps),"even":len(even),
                 "odd":len(odd),"block_rank":rank,"index":len(even)-len(odd),
                 "nullity":nullity,"boolean_edges":sum(map(sum,B))//2}
        if n>=5:
            c4=old.edge_mask(index,[(0,1),(1,2),(2,3),(3,0)])
            c5=old.edge_mask(index,[(0,1),(1,2),(2,3),(3,4),(4,0)])
            bow=old.edge_mask(index,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)])
            h7=old.edge_mask(index,[(i,j) for i in range(5) for j in range(i+1,5)
                                    if (i,j) not in {(0,1),(0,2),(1,2)}])
            gauge=[[M[orbit[x]][orbit[y]] for y in (c5,h7)] for x in (c4,bow)]
            require(gauge==[[4*(n-4),2*(n-4)],[4,4]], "all-height diagonal-gauge hostile")
            require(all(B[orbit[x]][orbit[y]]==1 for x in (c4,bow) for y in (c5,h7)), "Boolean hostile block")
            payload["gauge_hostile_block"]=gauge
        if len(even)==len(odd):
            payload["block_determinant"]=str(block.det())
        print("COMPLETE NATIVE KERNEL",json.dumps(payload,sort_keys=True),flush=True)
        rows.append(payload)
        if rank<len(odd):
            print("EXTRA ODD-SIDE KERNEL",[[str(x) for x in v] for v in block.nullspace()],flush=True)
    print("SEMANTIC SHA256",hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest())
    print("FINITE-EXACT n3..7; no all-order nullity or gap claim")


if __name__=="__main__":
    main()
