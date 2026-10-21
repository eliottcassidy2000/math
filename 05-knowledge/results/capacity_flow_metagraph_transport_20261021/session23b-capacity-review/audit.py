#!/usr/bin/env python3
"""Independent capacity audit: path-tree elimination and direct orientation census.
Does not import producer code and never touches Git.
"""
from itertools import combinations, product, permutations
from collections import Counter
from fractions import Fraction
import json


def check(value, message):
    if not value:
        raise RuntimeError(message)


def pairs(n):
    return list(combinations(range(n), 2))


def cyclic(n, steps):
    return [(u, v) if (v-u) % n in steps else (v, u) for u, v in pairs(n)]


def score(n, arcs):
    result = [0] * n
    for u, v in arcs:
        result[u] += 1
    return tuple(result)


def divergence(n, arcs, values):
    result = [0] * n
    for (u, v), x in zip(arcs, values):
        result[u] -= x
        result[v] += x
    return tuple(result)


def tree_flows(n, arcs, M, demand=None, zero_path=False):
    """Assign nonpath edges, solve each path edge from its prefix cut.
    Unlike producer, the eliminated tree is a path, not a star.
    The function allows either orientation of any path-tree edge.
    """
    demand = (0,) * n if demand is None else demand
    check(sum(demand) == 0, 'demand sum')
    path_ids = [i for i, (u, v) in enumerate(arcs) if abs(u-v) == 1]
    path_ids.sort(key=lambda i: min(arcs[i]))
    check(len(path_ids) == n-1, 'path spanning tree')
    chord_ids = [i for i in range(len(arcs)) if i not in path_ids]
    cut = [[(-1 if u <= j < v else 1 if v <= j < u else 0)
            for i in chord_ids for u, v in [arcs[i]]] for j in range(n-1)]
    bp = [sum(demand[:j+1]) for j in range(n-1)]
    signs = [1 if arcs[i][0] < arcs[i][1] else -1 for i in path_ids]
    for chord_values in product(range(M+1), repeat=len(chord_ids)):
        path_values = []
        for row, b, sg in zip(cut, bp, signs):
            x = sg * (sum(a*y for a, y in zip(row, chord_values)) - b)
            if x < 0 or x > M or (zero_path and x != 0):
                break
            path_values.append(x)
        else:
            f = [0] * len(arcs)
            for i, v in zip(chord_ids, chord_values):
                f[i] = v
            for i, v in zip(path_ids, path_values):
                f[i] = v
            yield tuple(f)


def normalized_regular(n):
    base = pairs(n)
    free = [i for i, (u, v) in enumerate(base) if v != u+1]
    result = []
    for bits in product((0, 1), repeat=len(free)):
        arcs = base[:]
        for i, bit in zip(free, bits):
            if bit:
                arcs[i] = arcs[i][::-1]
        if score(n, arcs) == ((n-1)//2,) * n:
            result.append(arcs)
    return result


def hamiltonian_paths(n, arcs):
    arc_set = set(arcs)
    return sum(all((p[i], p[i+1]) in arc_set for i in range(n-1))
               for p in permutations(range(n)))


def automorphisms(n, arcs):
    arc_set = set(arcs)
    return sum(all((p[u], p[v]) in arc_set for u, v in arcs)
               for p in permutations(range(n)))


def finite_diff(values):
    diffs = []
    while values:
        diffs.append(values[0])
        values = [b-a for a, b in zip(values, values[1:])]
    return diffs


def main():
    result = {}
    # Independent full n5 counts and exact h numerator from differences.
    A5 = cyclic(5, {1, 2})
    full5 = [sum(1 for f in tree_flows(5, A5, M)) for M in range(11)]
    expected = [(M+1)**2*(3*(M+1)**4+5*(M+1)**2+4)//12 for M in range(11)]
    check(full5 == expected, 'full n5 polynomial')
    check(finite_diff(full5)[7:] == [0]*4, 'degree-six differences')
    from math import comb
    h = [sum((-1)**j*comb(7, j)*full5[k-j] for j in range(min(7, k)+1)) for k in range(11)]
    check(h == [1,17,72,72,17,1,0,0,0,0,0], 'full n5 numerator')
    result['full_n5'] = {'M0_through_10': full5, 'h_star': h[:6], 'route': 'Hamiltonian-path-tree elimination'}

    # Direct path-normalized tournament census, independent of flow census.
    regular = {n: normalized_regular(n) for n in (3, 5, 7)}
    check([len(regular[n]) for n in (3,5,7)] == [1,3,91], 'fixed-path regular census')
    result['regular_fixed_path_census'] = {str(n): len(regular[n]) for n in regular}

    # Independent restricted n5 count and complete set-level transport.
    restricted_rows = []
    restricted_base = {}
    for M in range(7):
        source = set(tree_flows(5, A5, M, zero_path=True))
        check(len(source) == comb(M+2, 2), 'restricted n5 simplex count')
        for target in regular[5]:
            D = {i for i in range(10) if A5[i] != target[i]}
            image = {tuple(M-x if i in D else x for i, x in enumerate(f)) for f in source}
            truth = set(tree_flows(5, target, M, zero_path=True))
            check(image == truth, 'regular restricted affine transport')
            check(all(abs(A5[i][0]-A5[i][1]) != 1 for i in D), 'fixed path is not reversed')
        restricted_rows.append(len(source))
        restricted_base[M] = source
    result['restricted_n5'] = {'M0_through_6': restricted_rows, 'formula': 'binomial(M+2,2)', 'three_regular_bases': True}

    # Fixed-path complement has changed demand, with all path coordinates zero.
    path_transport = []
    for n in (3,5,7):
        A = cyclic(n, set(range(1, (n+1)//2)))
        D = {i for i,(u,v) in enumerate(A) if abs(u-v) != 1}
        target = [arc[::-1] if i in D else arc for i, arc in enumerate(A)]
        for M in ((1,2,3) if n<7 else (1,)):
            src = set(tree_flows(n,A,M,zero_path=True))
            b = (-M,)+(0,)*(n-2)+(M,)
            tgt = set(tree_flows(n,target,M,b,zero_path=True))
            image = {tuple(M-x if i in D else x for i,x in enumerate(f)) for f in src}
            check(image == tgt, 'zero-path source-to-sink transport')
            path_transport.append({'n':n, 'M':M, 'count':len(src)})
    result['zero_path_nonpath_complement'] = path_transport

    # Strong deletion: inspect every cut of every path-normalized regular n5/n7 T.
    cut_checks=0; min_out={}
    for n in (5,7):
        minimum=10**9
        for A in regular[n]:
            for mask in range(1,(1<<n)-1):
                S={v for v in range(n) if (mask>>v)&1}; k=len(S)
                out=sum(u in S and v not in S for u,v in A)
                removed=sum(i in S and i+1 not in S for i in range(n-1))
                check(out == k*(n-k)//2, 'regular cut law')
                check(removed <= min(k,n-k), 'path crossing bound')
                check(out-removed > 0, 'T minus P is strong')
                minimum=min(minimum,out-removed); cut_checks+=1
        min_out[str(n)]=minimum
    result['strong_deletion']={'cuts_checked':cut_checks,'minimum_outgoing_remaining':min_out}

    # n7 map across all 91 regular path presentations; no code quotienting.
    A7=cyclic(7,{1,2,3}); F7=list(tree_flows(7,A7,1,zero_path=True))
    check(len(F7)==91,'restricted n7 flow census')
    images=0
    for target in regular[7]:
        D={i for i in range(21) if A7[i]!=target[i]}
        for f in F7:
            g=tuple(1-x if i in D else x for i,x in enumerate(f))
            check(divergence(7,target,g)==(0,)*7,'n7 transported circulation')
            check(all(g[i]==0 for i,(u,v) in enumerate(target) if abs(u-v)==1),'n7 fixed path coordinates')
            images+=1
    result['n7_regular_basis_transport']={'base_flows':len(F7),'bases':91,'images_checked':images}

    # Restricted codegree3: B_D 1 = -p, so interior3 = 1 + binary p-demand flows.
    interior3={}
    for n,A in ((5,A5),(7,A7)):
        p=(-1,)+(0,)*(n-2)+(1,)
        G=list(tree_flows(n,A,1,p,zero_path=True))
        check(G, 'restricted codegree at most three')
        interior3[str(n)]=len(G)
        for g in G:
            f=tuple(0 if abs(u-v)==1 else x+1 for (u,v),x in zip(A,g))
            check(divergence(n,A,f)==(0,)*n,'interior-three translation')
            check(all(0<x<3 for (u,v),x in zip(A,f) if abs(u-v)!=1),'strict interior bounds')
    result['restricted_interior_at_M3']=interior3

    # Unequal capacities, arbitrary demands, signed-current and energy invariants.
    transport_checks=0; capacities=(1,2,3); edge_list=pairs(3)
    for signs in product((0,1),repeat=3):
        A=[e[::-1] if s else e for e,s in zip(edge_list,signs)]
        for flips in product((0,1),repeat=3):
            Ap=[e[::-1] if s else e for e,s in zip(A,flips)]
            shift=divergence(3,A,[c*s for c,s in zip(capacities,flips)])
            for f in product(*(range(c+1) for c in capacities)):
                g=tuple(c-x if s else x for c,x,s in zip(capacities,f,flips))
                b=divergence(3,A,f); bp=divergence(3,Ap,g)
                check(bp==tuple(x-y for x,y in zip(b,shift)),'unequal-capacity demand transport')
                j=tuple(2*x-c for x,c in zip(f,capacities)); jp=tuple(2*x-c for x,c in zip(g,capacities))
                check(divergence(3,A,j)==divergence(3,Ap,jp),'centered-current divergence')
                check(sum((i+1)*x*x for i,x in enumerate(j))==sum((i+1)*x*x for i,x in enumerate(jp)),'fixed-edge quadratic energy')
                transport_checks+=1
    result['unequal_capacities_and_centered_currents']={'cases':transport_checks,'capacities':capacities}

    # H and total-use information genuinely differs within the regular fibre.
    info=[]
    for steps in ({1,2,3},{1,2,4}):
        A=cyclic(7,steps); F=list(tree_flows(7,A,1)); histogram=Counter(map(sum,F))
        info.append({'steps':sorted(steps),'H':hamiltonian_paths(7,A),'automorphisms':automorphisms(7,A),'full_M1':len(F),'q4':histogram[4]})
    check([x['H'] for x in info]==[175,189],'H hostile')
    check([x['q4'] for x in info]==[28,21],'total-use hostile')
    check([x['full_M1'] for x in info]==[2640,2640],'full n7 score fibre')
    result['information_loss']=info
    result['status']='PASS'
    return result

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
