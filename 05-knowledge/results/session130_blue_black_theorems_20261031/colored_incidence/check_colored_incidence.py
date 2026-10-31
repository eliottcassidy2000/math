#!/usr/bin/env python3
"""Exact colored-line/deletion transport. Standard library; no input files.
Explicit failure gates survive python -O. The n<=5 atlas is rebuilt from all
fixed-path masks and independently checked against all labelled tournaments.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations, permutations
from math import gcd
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def tiles(n):
    return tuple((x, y) for y in range(1, n-1) for x in range(n, y+1, -1))


@lru_cache(None)
def reflection_indices(n):
    ts = tiles(n)
    index = {xy: k for k, xy in enumerate(ts)}
    return tuple(index[n-y+1, n-x+1] for x, y in ts)


def reflect(n, t):
    return sum(((t >> k) & 1) << j for k, j in enumerate(reflection_indices(n)))


def full(n):
    return (1 << len(tiles(n))) - 1


def line(n, t):
    return min(t, t ^ full(n))


def inner(n, t):
    """Restrict to vertices n-1,...,2 and renumber them n-2,...,1."""
    index = {xy: k for k, xy in enumerate(tiles(n-2))}
    return sum(((t >> k) & 1) << index[x-1, y-1]
               for k, (x, y) in enumerate(tiles(n)) if x < n and y > 1)


def central_of_lower(n, s):
    """Delete the first path vertex from an order-(n-1) mask."""
    # Lower vertices n-1,...,1 represent original vertices n,...,2.
    index = {xy: k for k, xy in enumerate(tiles(n-2))}
    return sum(((s >> k) & 1) << index[x, y]
               for k, (x, y) in enumerate(tiles(n-1)) if x < n-1)


def adjacency(n, t):
    a = [[0]*n for _ in range(n)]
    for i in range(n-1):
        a[i][i+1] = 1
    for k, (x,y) in enumerate(tiles(n)):
        i,j = n-x,n-y
        a[i][j] = 1-((t >> k)&1)
        a[j][i] = (t >> k)&1
    return tuple(tuple(row) for row in a)


@lru_cache(None)
def perms(n):
    return tuple(permutations(range(n)))


@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n),2))


def canon(a):
    n = len(a)
    return min(sum(a[p[i]][p[j]] << (len(pairs(n))-1-k)
                   for k, (i,j) in enumerate(pairs(n))) for p in perms(n))


def transpose(a):
    return tuple(zip(*a))


def atlas(n, full_labelled=True):
    arrays = [adjacency(n,t) for t in range(full(n)+1)]
    ordinary = [canon(a) for a in arrays]
    merged = [min(canon(a),canon(transpose(a))) for a in arrays]
    labels = {c:i for i,c in enumerate(sorted(set(merged)))}
    owner = [labels[c] for c in merged]
    fibers = [[t for t in range(full(n)+1) if owner[t] == u]
              for u in range(len(labels))]
    all_labelled_classes = set()
    for b in range((1 << len(pairs(n))) if full_labelled else 0):
        a = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs(n)):
            a[i][j] = (b >> k)&1
            a[j][i] = 1-a[i][j]
        all_labelled_classes.add(canon(a))
    if full_labelled:
        check(set(ordinary) == all_labelled_classes, 'Fixed-path atlas incomplete')
    check(len(fibers) == {3:2,4:3,5:10,6:34}[n], 'Merged class count')
    rows = []
    for u, f in enumerate(fibers):
        a = arrays[f[0]]
        h = sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in perms(n))
        aut = sum(all(a[p[i]][p[j]] == a[i][j] for i,j in pairs(n)) for p in perms(n))
        sheet_count = len({ordinary[t] for t in f})
        b = sum(reflect(n,t) == t for t in f)
        check(len(f)*aut == sheet_count*h, 'Hamiltonian orbit-stabilizer law')
        check(h % 2 == 1 and aut % 2 == 1, 'Odd path/aut invariant')
        check((b % 2 == 1) == (sheet_count == 1), 'Self-converse blue parity')
        check((b != 0) == (sheet_count == 1), 'Blue support is self-converse')
        rows.append(dict(node=u,H=h,aut=aut,ordinary_sheets=sheet_count,
                         fiber_size=len(f),blue_size=b,masks=f))
    return owner,rows


def component_count(vertices, edges):
    parent = list(range(vertices))
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for u,v in edges:
        parent[root(u)] = root(v)
    return len({root(v) for v in range(vertices)})


def main():
    lift_audits = []
    for n in range(4,8):
        total,blue,black = Counter(),Counter(),Counter()
        for t in range((full(n)+1)//2):
            r = reflect(n,t)
            check(r != (t ^ full(n)), 'Forbidden anti-fixed mask exists')
            child = line(n-1,t >> (n-2))
            total[child] += 1
            (blue if r == t else black)[child] += 1
            check(((t ^ full(n)) >> (n-2)) == ((t >> (n-2)) ^ full(n-1)), 'Deletion not C-equivariant')
            check(inner(n,t) == central_of_lower(n,t >> (n-2)), 'Central restriction mismatch')
        eligible = []
        for s in range((full(n-1)+1)//2):
            c = central_of_lower(n,s)
            good = reflect(n-2,c) == c
            if good: eligible.append(s)
            check(total[s] == 2**(n-2), 'Nonuniform line lift count')
            check(blue[s] == (2 if good else 0), 'Blue lift theorem')
            check(black[s] == 2**(n-2)-blue[s], 'Black lift theorem')
        lift_audits.append(dict(n=n,lower_lines=len(total),eligible_lower_lines=len(eligible),
                               upper_blue_lines=sum(blue.values()),upper_black_lines=sum(black.values()),
                               total_lifts_per_line=2**(n-2),blue_lifts_per_eligible_line=2,
                               blue_lifts_per_ineligible_line=0))
    own3,rows3 = atlas(3)
    own4,rows4 = atlas(4)
    own5,rows5 = atlas(5)
    check(len({(r['H'],r['blue_size']) for r in rows5}) == 10, 'n5 H/blue classifier')
    own6,rows6 = atlas(6, full_labelled=False)
    collision_groups = defaultdict(list)
    for r in rows6:
        collision_groups[r['H'],r['blue_size']].append(r['node'])
    collisions = [dict(H=h,blue_size=b,nodes=us) for (h,b),us in sorted(collision_groups.items()) if len(us)>1]
    check(any(c['nodes']==[11,16] and c['H']==17 and c['blue_size']==5 for c in collisions), 'n6 H/blue hostile pair')
    colored_edges = Counter()
    literal_edges,blue_edges,halfblack_edges,black_edges = [],[],[],[]
    for t in range(32):
        uv = tuple(sorted((own5[t],own5[t^63])))
        literal_edges.append(uv)
        if reflect(5,t) == t:
            color='blue';blue_edges.append(uv)
        else:
            color='black';black_edges.append(uv)
            if t < line(5,reflect(5,t)):
                halfblack_edges.append(uv)
        colored_edges[uv,color] += 1
    check(len(blue_edges)==8 and len(black_edges)==24 and len(halfblack_edges)==12,'Colored line counts')
    for (uv,c),count in colored_edges.items():
        check(c!='black' or count%2==0,'Black parallel-pair parity')
    profiles=[]
    for row in rows5:
        u=row['node'];raw=[0]*3;b=[0]*3;k=Counter()
        for t in row['masks']:
            v=own4[t >> 3];raw[v]+=1
            r=reflect(5,t)
            if r==t:b[v]+=1
            elif t<r:k[tuple(sorted((v,own4[r >> 3])))]+=1
        degree=[0]*3
        for (v,w),count in k.items():
            degree[v]+=count;degree[w]+=count
        check(raw==[b[i]+degree[i] for i in range(3)],'Colored deletion decomposition')
        g=reduce(gcd,raw)
        primitive=[x//g for x in raw]
        code=7**raw[0]*17**raw[1]*55**raw[2]
        blue_code=7**b[0]*17**b[1]*55**b[2]
        factors=1
        for (v,w),count in k.items():factors*=((7,17,55)[v]*(7,17,55)[w])**count
        check(code==blue_code*factors,'Arithmetic color factorization')
        profiles.append(dict(node=u,raw=raw,primitive=primitive,blue=b,black_degree=degree,
                             black_parent_edges=[dict(parents=vw,multiplicity=count) for vw,count in sorted(k.items())],
                             raw_code_7_17_55=code,blue_code_7_17_55=blue_code,gcd=g))
    check(len({tuple(p['primitive']) for p in profiles})==10,'Primitive fingerprint completeness')
    check([sum(p['blue'][j] for p in profiles) for j in range(3)]==[2,4,10],'Blue column mass')
    # A literal edge observable, not a node-code cardinality match.
    ratio_counts={'blue':Counter(),'black':Counter()}
    lower_ratios=[]
    for s in range(4):
        lower_ratios.append(str(Fraction((7,17,55)[own4[s^7]],(7,17,55)[own4[s]])))
    for t in range(32):
        s=t >> 3
        ratio=Fraction((7,17,55)[own4[s^7]],(7,17,55)[own4[s]])
        ratio_counts['blue' if reflect(5,t)==t else 'black'][str(ratio)]+=1
    check(lower_ratios==['55/7','1','55/17','17/55'],'Selected ratio dictionary changed')
    check(all(ratio_counts['blue'][q]==2 and ratio_counts['black'][q]==6 for q in lower_ratios),'Color-independent lower ratio law')
    # d does not descend through the R quotient, even if the target is R-quotiented.
    orbit=lambda n,t:min(t,t^full(n),reflect(n,t),reflect(n,t)^full(n))
    check(reflect(5,2)==8,'Witness gauge changed')
    check(orbit(4,2>>3)!=orbit(4,8>>3),'Unexpected line-quotient descent')
    ranks={}
    for name,edges in [('full_literal',literal_edges),('blue',blue_edges),('black_literal',black_edges),('black_R_quotient',halfblack_edges),('full_R_quotient',blue_edges+halfblack_edges)]:
        c=component_count(10,edges)
        ranks[name]=dict(vertices=10,edges=len(edges),components=c,cycle_rank=len(edges)-10+c)
    check([ranks[k]['cycle_rank'] for k in ('full_literal','blue','black_literal','black_R_quotient','full_R_quotient')]==[23,1,18,6,11],'Cycle ranks')
    print(json.dumps(dict(status='FINITE-EXACT checks; symbolic proofs in REPORT.md',
        lift_audits=lift_audits,n5_nodes=rows5,n5_profiles=profiles,
        n5_colored_edges=[dict(nodes=uv,color=c,multiplicity=count) for (uv,c),count in sorted(colored_edges.items())],
        n5_graph_ranks=ranks,lower_line_ratios=lower_ratios,
        n5_H_blue_classifier=[dict(node=r['node'],H=r['H'],blue_size=r['blue_size']) for r in rows5],
        n6_H_blue_collisions=collisions,
        n6_H_blue_hostile_witnesses=[dict(node=r['node'],representative_mask=r['masks'][0],H=r['H'],blue_size=r['blue_size'],aut=r['aut'],fiber_size=r['fiber_size'],raw_deletion_profile=[sum(own5[t >> 4] == v for t in r['masks']) for v in range(10)]) for r in rows6 if r['node'] in [11,16,27,28]],
        n5_ratio_counts={c:dict(v) for c,v in ratio_counts.items()},
        no_line_R_quotient_descent=dict(upper_masks=[2,8],lower_masks=[0,1],lower_line_R_orbits=[0,1]),
        controls=dict(independent_full_labelled_n3_n5=True,optimized_safe_failure_gates=True,no_external_data=True)),indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
