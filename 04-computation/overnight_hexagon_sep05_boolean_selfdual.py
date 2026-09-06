#!/usr/bin/env python3
"""All-n trace dictionary: exact bounded primal/dual and complement controls.

Only n=3,...,6 cycle spaces are materialized. Dual classes are built by
adjacent-transposition regauging; complement-fixed status is independently
tested on literal odd-triangle sets under all vertex permutations.
"""
from collections import Counter
from itertools import combinations, permutations
from math import factorial
import hashlib
import json
import sympy as sp

import even_graph_triangle_quotient_spectrum_thm4078 as primal
import even_graph_triangle_quotient_spectrum_thm4078_independent_audit as dual
import self_complementary_graphs_burnside_s566 as old_complement
import eulerian_uniform_orbit_sampler_overnight_hexagon_sep05 as fixed
from overnight_hexagon_sep05_boolean_peck import parity_trivial


def require(ok, data):
    if not ok:
        raise RuntimeError(repr(data))


def dual_orbits(n, edges, index):
    size = 1 << len(edges)
    maps = []
    for j in range(n-1):
        p = list(range(n))
        p[j],p[j+1] = p[j+1],p[j]
        maps.append([dual.regauge_after_permutation(x,n,tuple(p),edges,index) for x in range(size)])
    orbit = {}
    classes = []
    for seed in range(size):
        if seed in orbit:
            continue
        label = len(classes)
        orbit[seed] = label
        stack, members = [seed], []
        while stack:
            x = stack.pop()
            members.append(x)
            for transform in maps:
                y = transform[x]
                if y not in orbit:
                    orbit[y] = label
                    stack.append(y)
        classes.append(sorted(members))
    return classes,orbit


def triangle_set(mask,n,gedges):
    signs = {edge for j,edge in enumerate(gedges) if mask >> j & 1}
    return frozenset(t for t in combinations(range(n),3)
                     if sum(tuple(sorted(e)) in signs for e in combinations(t,2))%2)


def main():
    gates = 0
    rows = []
    for n in range(3,7):
        gedges,gindex = dual.gauge_edges(n)
        classes,where = dual_orbits(n,gedges,gindex)
        all_negative = (1 << len(gedges))-1
        involution = [where[cl[0]^all_negative] for cl in classes]
        fixed_orbits = [i for i,j in enumerate(involution) if i==j]
        all_triangles = frozenset(combinations(range(n),3))
        literal_fixed = []
        for i,cl in enumerate(classes):
            T = triangle_set(cl[0],n,gedges)
            target = all_triangles-T
            is_fixed = any(frozenset(tuple(sorted(p[v] for v in t)) for t in T)==target
                           for p in permutations(range(n)))
            require(is_fixed==(i in fixed_orbits), ('literal two-graph complement',n,i))
            if is_fixed:
                literal_fixed.append(i)
            gates += 1
        edges,ei = primal.edge_system(n)
        reps,orbit,sizes,_ = primal.quotient_data(n)
        parity = [(-1)**rep.bit_count() for rep in reps]
        difference = sum(parity)
        require(len(classes)==len(reps) and len(fixed_orbits)==difference, ('invariant trace identity',n))
        coordinates = [sum(((rep >> ei[e]) & 1) << j for j,e in enumerate(gedges)) for rep in reps]
        P = sp.Matrix([[sum((-1)**((x&h).bit_count()) for h in cl) for cl in classes] for x in coordinates])
        require(P.det()!=0, ('full orbit Fourier change of basis',n))
        for i in range(len(reps)):
            for j in range(len(classes)):
                require(parity[i]*P[i,j]==P[i,involution[j]], ('multiplication translates dual',n,i,j))
                gates += 1
        triangles = primal.simple_cycles(n,3,ei)
        M = sp.Matrix(primal.weighted_operator(reps,orbit,triangles))
        B = sp.Matrix(primal.boolean_support(M.tolist()))
        failing = 0
        for j in fixed_orbits:
            require(M*P[:,j]==sp.zeros(len(reps),1), ('weighted fixed-orbit zero mode',n,j))
            failing += int(B*P[:,j]!=sp.zeros(len(reps),1))
            gates += 1
        # Literal Eulerian members of every individual labelled cut class.
        even_member_counts = []
        for cl in classes:
            H = sum((1 << ei[e]) for j,e in enumerate(gedges) if cl[0] >> j & 1)
            count = 0
            for subset in range(1 << (n-1)):
                switched = H
                for i,j in edges:
                    inside_i = int(i>0 and (subset >> (i-1)) & 1)
                    inside_j = int(j>0 and (subset >> (j-1)) & 1)
                    if inside_i != inside_j:
                        switched ^= 1 << ei[i,j]
                degrees = [sum((switched >> ei[min(i,j),max(i,j)])&1 for j in range(n) if j!=i)%2 for i in range(n)]
                count += int(not any(degrees))
            require(count==1 if n%2 else count in (0,2**(n-2)), ('odd/even representative boundary',n,count))
            even_member_counts.append(count)
            gates += 1
        if n%2:
            complement_mask = (1 << len(edges))-1
            eulerian_self = sum(orbit[rep^complement_mask]==i for i,rep in enumerate(reps))
            require(eulerian_self==difference, ('odd-order self-complementary Eulerian count',n))
            gates += 1
        if n==4:
            j = fixed_orbits[0]
            require(tuple(P[:,j])==(6,0,-2) and tuple(B*P[:,j])==(0,4,0), 'minimal Boolean transfer hostile')
            require(even_member_counts[j]==0, 'fixed dual orbit need not have Eulerian representative')
            gates += 2
        row = {'n':n,'orbits':len(reps),'parity_index':difference,
               'self_complementary_two_graph_orbits':len(fixed_orbits),
               'fixed_modes_failing_Boolean':failing,
               'Eulerian_members_per_cut_class':dict(sorted(Counter(even_member_counts).items()))}
        rows.append(row)
        print('EXACT PRIMAL/DUAL',row)
        gates += 2
    for k in range(1,6):
        n = 4*k+1
        numerator = sum(fixed.class_size(n,p)*2**fixed.closed_dimension(p)
                        for p in fixed.partitions(n) if parity_trivial(p))
        difference = numerator//factorial(n)
        ordinary = old_complement.burnside_undirected(4*k)[1]
        require(difference==ordinary, ('recovered ordinary self-complementary count',k,difference,ordinary))
        print('RECOVERED SELF-COMPLEMENTARY n',4*k,'COUNT',ordinary,'= PARITY INDEX n',n)
        gates += 1
    require(old_complement.brute_force_undirected(5)[1]==2, 'ordinary SC5 differs from Eulerian SC5=1')
    gates += 1
    print('GATES',gates,'SEMANTIC SHA256',hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest())
    print('PROVED all-n invariant trace dictionary; self-complementary zero modes are WEIGHTED only')
    print('OPEN constructive native Boolean kernel dictionary and exact all-n nullity')


if __name__=='__main__':
    main()
