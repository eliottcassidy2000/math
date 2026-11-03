#!/usr/bin/env python3
"""Exact controls for Session 170 topology compiler. Standard library only."""
from itertools import combinations, permutations
from collections import defaultdict
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)


def rank2(columns):
    pivots = {}
    for column in columns:
        while column:
            p = column.bit_length() - 1
            if p not in pivots:
                pivots[p] = column
                break
            column ^= pivots[p]
    return len(pivots)


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub-1) & mask


def cube_complex(k, facets):
    """Cells are (free-coordinate mask, 1-bits among fixed coordinates)."""
    faces = set()
    for facet in facets:
        faces.update(subsets(facet))
    cells = [[] for _ in range(max(f.bit_count() for f in faces)+1)]
    full = (1 << k)-1
    for free in sorted(faces):
        for fixed in subsets(full ^ free):
            cells[free.bit_count()].append((free, fixed))
    indices = [{cell:i for i,cell in enumerate(layer)} for layer in cells]
    boundaries = [None]
    columns = [[]]
    for d in range(1, len(cells)):
        bd = []
        for free, fixed in cells[d]:
            one = {}
            for j, coord in enumerate(i for i in range(k) if free >> i & 1):
                sign = (-1)**j
                for bit, factor in [(0,-sign),(1,sign)]:
                    face = (free ^ (1 << coord), fixed | (bit << coord))
                    one[indices[d-1][face]] = factor
            bd.append(one)
        boundaries.append(bd)
        columns.append([sum(1 << i for i in x) for x in bd])
    # Verify d^2=0 over the integers, not only mod 2.
    for d in range(2, len(cells)):
        for boundary in boundaries[d]:
            twice = defaultdict(int)
            for i, coeff in boundary.items():
                for j, coeff2 in boundaries[d-1][i].items():
                    twice[j] += coeff*coeff2
            require(all(x == 0 for x in twice.values()), 'boundary squared')
    ranks = [0] + [rank2(cols) for cols in columns[1:]] + [0]
    betti = [len(cells[d])-ranks[d]-ranks[d+1] for d in range(len(cells))]
    return cells, boundaries, betti


def surface_test(k):
    pairs = [(i,(i+1)%k) for i in range(k)]
    cells, bd, betti = cube_complex(k, [(1 << i)|(1 << j) for i,j in pairs])
    incidences = [0]*len(cells[1])
    orient = [0]*len(cells[1])
    for (free,fixed),boundary in zip(cells[2],bd[2]):
        pair = next((i,j) for i,j in pairs if free == ((1 << i)|(1 << j)))
        # Positive cyclic coordinate order; convert to increasing coordinate order.
        coefficient = (-1)**fixed.bit_count() * (1 if pair[0] < pair[1] else -1)
        for edge,sign in boundary.items():
            incidences[edge] += 1
            orient[edge] += coefficient*sign
    require(all(x == 2 for x in incidences), 'surface edge incidence')
    require(all(x == 0 for x in orient), 'surface integral orientation cycle')
    # At every corner, links have exactly these coordinate vertices and cyclic edges.
    require(all(sum(i in edge for edge in pairs) == 2 for i in range(k)), 'vertex links')
    genus = 1 + (k-4)*(2**(k-3))
    require(betti == [1,2*genus,1], 'surface Betti formula')
    counts = list(map(len,cells))
    require(counts == [2**k,k*2**(k-1),k*2**(k-2)], 'surface cell count')
    return {'coordinates':k,'cells':counts,'genus':genus,'betti_mod2':betti,
            'tournament_order':6*k,'hamilton_paths':15**k}


def torus_test(d):
    k = 2*d
    facets = []
    for bits in range(1 << d):
        facets.append(sum(1 << (2*i+((bits >> i)&1)) for i in range(d)))
    cells,_,betti = cube_complex(k,facets)
    from math import comb
    require(betti == [comb(d,j) for j in range(d+1)], 'torus Betti formula')
    require(list(map(len,cells)) == [4**d*comb(d,j) for j in range(d+1)], 'torus cells')
    return {'dimension':d,'coordinates':k,'cells':list(map(len,cells)),
            'betti_mod2':betti,'tournament_order':12*d,'hamilton_paths':15**(2*d)}


def sphere_test(k):
    full = (1 << k)-1
    cells,_,betti = cube_complex(k,[full ^ (1 << i) for i in range(k)])
    expected = [1] + [0]*(k-2) + [1]
    require(betti == expected, 'cube boundary sphere')
    return {'coordinates':k,'dimension':k-1,'cells':list(map(len,cells)),'betti_mod2':betti}


def tournament_orbits(n):
    pairs = list(combinations(range(n),2))
    pindex = {pair:i for i,pair in enumerate(pairs)}
    actions = []
    for perm in permutations(range(n)):
        mapping = []
        for i,j in pairs:
            x,y = perm[i],perm[j]
            mapping.append((pindex[tuple(sorted((x,y)))],int(x > y)))
        actions.append(mapping)
    size = 1 << len(pairs)
    canonical = [-1]*size
    representatives = []
    for mask in range(size):
        if canonical[mask] >= 0:
            continue
        representatives.append(mask)
        for mapping in actions:
            moved = 0
            for src,(dst,flip) in enumerate(mapping):
                moved |= (((mask >> src)&1)^flip) << dst
            require(canonical[moved] in (-1,mask),'overlapping tournament orbit')
            canonical[moved] = mask
    return pairs,canonical,representatives


def corner_test(n):
    pairs,canon,reps = tournament_orbits(n)
    chords = [pairs.index((i,j)) for i,j in pairs if j >= i+2]
    m = len(chords)
    backbone = sum(1 << pairs.index((i,i+1)) for i in range(n-1))
    # S151 flag bit 1 reverses the forward chord; tournament mask bit 1 is forward.
    flags = [canon[backbone | sum((1 ^ ((t >> i)&1)) << chord for i,chord in enumerate(chords))]
             for t in range(1 << m)]
    classes = sorted(set(flags))
    require(classes == reps,'Hamilton marked flags cover every class')
    index = {x:i for i,x in enumerate(classes)}
    flags = [index[x] for x in flags]
    full = (1 << m)-1
    lines = [(t,t^full) for t in range(1 << m) if t < (t^full)]
    parent = list(range(len(classes)))
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for t,u in lines:
        parent[root(flags[t])] = root(flags[u])
    component_roots = {root(i) for i in range(len(classes))}
    components = len(component_roots)
    all_arcs = (1 << len(pairs))-1
    converse = [index[canon[code ^ all_arcs]] for code in classes]
    sc_count = sum(converse[i] == i for i in range(len(classes)))
    chord_pairs = [pairs[i] for i in chords]
    reflection = [chord_pairs.index((n-1-j,n-1-i)) for i,j in chord_pairs]
    reflected = [sum(((t >> i)&1) << reflection[i] for i in range(m))
                 for t in range(1 << m)]
    require(all(flags[reflected[t]] == converse[flags[t]] for t in range(1 << m)),
            'reflection agrees with converse class')
    fixed_flags = [t for t in range(1 << m) if reflected[t] == t]
    d = (n-1)**2 // 4
    require(len(fixed_flags) == 2**d, 'reflection fixed cube dimension')
    require(set(flags[t] for t in fixed_flags) == {i for i in range(len(classes)) if converse[i] == i},
            'fixed flags cover precisely SC classes')
    fixed_components = sum(root(converse[r]) == r for r in component_roots)
    blue_lines = sum(reflected[t] == t for t,u in lines)
    require(blue_lines == 2**(d-1), 'blue fixed line count')
    require(not any(reflected[t] == u for t,u in lines), 'no line inversion')
    k,e = len(classes),len(lines)
    graph_rank = e-k+components
    thick_rank = 2*e-k
    extra_rank = e-components
    require(graph_rank+extra_rank == thick_rank,'free-factor rank formula')
    # Exact chain complex of flag hub graph plus line graph plus triangular cells.
    # vertices: k classes and hub; first 2e edges are flags, remaining e edges lines.
    bd1 = [(1 << k) ^ (1 << c) for c in flags]
    bd1 += [(1 << flags[t]) ^ (1 << flags[u]) for t,u in lines]
    bd2 = [(1 << t) ^ (1 << u) ^ (1 << (2*e+i)) for i,(t,u) in enumerate(lines)]
    require(rank2(bd1) == k,'connected thickening')
    require(rank2(bd2) == e,'independent triangular boundaries')
    require(len(bd1)-rank2(bd1)-rank2(bd2) == thick_rank,'thickening H1')
    cubical = None
    if n <= 5:
        cells,boundaries,_ = cube_complex(m,[full])
        # Identify only cube corners. All positive-dimensional cells remain distinct.
        lookup = {cell:i for i,cell in enumerate(cells[0])}
        corner_lookup = {i:flags[fixed] for i,(free,fixed) in enumerate(cells[0])}
        quot1 = []
        for boundary in boundaries[1]:
            col = 0
            for i in boundary:
                col ^= 1 << corner_lookup[i]
            quot1.append(col)
        ranks = [0,rank2(quot1)]
        ranks += [rank2([sum(1 << i for i in x) for x in boundaries[d]])
                  for d in range(2,len(cells))]
        ranks += [0]
        ncells = [k]+list(map(len,cells[1:]))
        betti = [ncells[d]-ranks[d]-ranks[d+1] for d in range(len(cells))]
        require(betti == [1,thick_rank]+[0]*(len(cells)-2),'corner cube quotient homology')
        cubical = {'cells':ncells,'betti_mod2':betti}
    fixed_rank = 2**d-sc_count
    gamma_trace = blue_lines-sc_count+fixed_components
    extra_trace = blue_lines-fixed_components
    require(gamma_trace+extra_trace == fixed_rank, 'relative reflection character')
    require((thick_rank-fixed_rank)%2 == 0, 'reflection eigenspace integrality')
    reflection_record = {'fixed_cube_dimension':d,'sc_classes':sc_count,
        'fixed_flags':len(fixed_flags),'fixed_locus_rank':fixed_rank,
        'total_eigenspace_dimensions':[(thick_rank+fixed_rank)//2,(thick_rank-fixed_rank)//2],
        'graph_eigenspace_dimensions':[(graph_rank+gamma_trace)//2,(graph_rank-gamma_trace)//2],
        'extra_trace':extra_trace}
    return {'n':n,'chords':m,'flags':2*e,'classes':k,'lines':e,'components':components,
            'metagraph_cycle_rank':graph_rank,'corner_quotient_rank':thick_rank,
            'extra_free_rank':extra_rank,'direct_cubical_check':cubical,
            'reflection':reflection_record}


def seed_profile(mask):
    n=6
    pairs=list(combinations(range(n),2))
    arcs=[[False]*n for _ in range(n)]
    for z,(i,j) in enumerate(pairs):
        arcs[i][j]=bool(mask >> z & 1)
        arcs[j][i]=not arcs[i][j]
    reach=[row[:] for row in arcs]
    for h in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][h] and reach[h][j]
    require(all(reach[i][j] for i in range(n) for j in range(n)),'seed strong')
    h=[0]*(1 << n)
    for size in range(1,n+1):
        for vertices in combinations(range(n),size):
            h[sum(1 << x for x in vertices)] = sum(
                all(arcs[p[i]][p[i+1]] for i in range(size-1)) for p in permutations(vertices))
    profiles=[[0]*(n+1) for _ in range(1 << n)]
    profiles[0][0]=1
    for available in range(1,1 << n):
        first=available & -available
        for part in subsets(available):
            if not part & first:
                continue
            for j in range(1,n+1):
                profiles[available][j] += h[part]*profiles[available ^ part][j-1]
    return profiles[-1]


def main():
    p = seed_profile(18)
    require(p == seed_profile(22),'equal seed profiles')
    require(p == [0,15,93,138,73,15,1],'expected seed profile')
    falling=[1]
    shared=[0]*7
    for j in range(7):
        for degree,coefficient in enumerate(falling):
            shared[degree] += p[j]*coefficient
        if j < 6:
            next_falling=[0]*(len(falling)+1)
            for degree,coefficient in enumerate(falling):
                next_falling[degree] -= j*coefficient
                next_falling[degree+1] += coefficient
            falling=next_falling
    require(shared == [0,0,6,0,8,0,1], 'shared slot polynomial')
    _,canon,_=tournament_orbits(6)
    require(canon[18] != canon[22],'nonisomorphic seed pair')
    require((18 ^ 22).bit_count() == 1,'one flip seed pair')
    out={'seed':{'masks':[18,22],'profile_coefficients':p,'shared_slot_coefficients':shared},
         'surfaces':[surface_test(k) for k in range(3,10)],
         'tori':[torus_test(d) for d in range(1,5)],
         'spheres':[sphere_test(k) for k in range(2,8)],
         'corner_quotients':[corner_test(n) for n in range(3,7)]}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
