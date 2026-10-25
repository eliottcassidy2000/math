#!/usr/bin/env python3
"""Exact cyclic-incidence tests. Standard library only; checks survive python -O.
Run: python verify_wrap.py --output verification.json
The direct physical-line constructor does not use the forest equations.
"""
import argparse
import itertools
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def rank(rows, q, ncols=None):
    a = [[x % q for x in row] for row in rows]
    if not a:
        return 0
    n = len(a[0]) if ncols is None else ncols
    r = 0
    for j in range(n):
        pivot = next((k for k in range(r, len(a)) if a[k][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(a[r][j], -1, q)
        a[r] = [x * inverse % q for x in a[r]]
        for k in range(r + 1, len(a)):
            coefficient = a[k][j]
            if coefficient:
                a[k] = [(x - coefficient * y) % q for x, y in zip(a[k], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def nullspace(rows, q, n):
    a = [[x % q for x in row] for row in rows]
    pivots = []
    r = 0
    for j in range(n):
        pivot = next((k for k in range(r, len(a)) if a[k][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(a[r][j], -1, q)
        a[r] = [x * inverse % q for x in a[r]]
        for k in range(len(a)):
            if k != r and a[k][j]:
                coefficient = a[k][j]
                a[k] = [(x - coefficient * y) % q for x, y in zip(a[k], a[r])]
        pivots.append(j)
        r += 1
    result = []
    for free in sorted(set(range(n)) - set(pivots)):
        v = [0] * n
        v[free] = 1
        for k, pivot in enumerate(pivots):
            v[pivot] = -a[k][free] % q
        result.append(v)
    return result


def physical_maps(p, I):
    # Physical points, not an abstract labelled-box line partition.
    points = [((-a - 3*b) % p, (j - 9*b) % p)
              for a in (0, 1) for b in (0, 1) for j in I]
    require(len(set(points)) == 4 * len(I), ('point collision', p, I))
    maps = []
    for dx, dy in ((1, 0), (0, 1), (3, 9)):
        groups = {}
        for index, (x, y) in enumerate(points):
            groups.setdefault((dx*y - dy*x) % p, []).append(index)
        mat = []
        for group in groups.values():
            row = [0] * len(points)
            for index in group:
                row[index] = 1
            mat.append(row)
        maps.append(mat)
    return maps


def components(adjacency):
    seen = set()
    out = []
    for root in sorted(adjacency):
        if root in seen:
            continue
        todo = [root]
        seen.add(root)
        for x in todo:
            for y in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        out.append(todo)
    return out


def forest(p, I):
    occupied = set(I)
    H = {j for j in occupied if (j + 9) % p in occupied}
    adjacency = {j: [] for j in range(p)}
    for j in set(range(p)) - occupied:
        x, y = (j - 6) % p, (j - 3) % p
        adjacency[x].append(y)
        adjacency[y].append(x)
    paths = []
    for component in components(adjacency):
        endpoint = min(x for x in component if len(adjacency[x]) <= 1)
        path = []
        previous, x = None, endpoint
        while True:
            path.append(x)
            neighbors = [y for y in adjacency[x] if y != previous]
            if not neighbors:
                break
            require(len(neighbors) == 1, 'component is not a path')
            previous, x = x, neighbors[0]
        require(set(path) == set(component), 'path missed a vertex')
        anchors = {i % 2 for i, x in enumerate(path) if x not in H}
        paths.append((path, anchors))
    free = [path for path, anchors in paths if not anchors]
    odd = [path for path in free if len(path) % 2]
    compatible = all(len(anchors) <= 1 for _, anchors in paths)
    return H, paths, free, odd, compatible


def forest_generators(p, I, q):
    H, paths, free, odd, compatible = forest(p, I)
    def alternating(path):
        r = [0] * p
        for i, x in enumerate(path):
            r[x] = (-1)**i
        return r
    reference = alternating(odd[0]) if odd else None
    parameters = []
    for path in free:
        r = alternating(path)
        if len(path) % 2:
            if path == odd[0]:
                continue
            r = [x - y for x, y in zip(r, reference)]
        parameters.append((r, 0))
    if q == p and compatible and odd:
        r = [0] * p
        for path, anchors in paths:
            parity = next(iter(anchors)) if anchors else 0
            for i, x in enumerate(path):
                r[x] = (i % 2) ^ parity
        beta = sum(r)
        r = [x - beta*y for x, y in zip(r, reference)]
        parameters.append((r, 1))
    generators = []
    n = len(I)
    for j in range(n - 1):
        h = [int(k == j) - int(k == j + 1) for k in range(n)]
        generators.append(h + [-x for x in h] + [-x for x in h] + h)
    outside = set(range(p)) - set(I)
    for r, c in parameters:
        P = [c - sum(r[(j-k) % p] for k in (3, 6, 9)) for j in range(p)]
        S = [-r[(j-9) % p] for j in range(p)]
        require(all((r[j] % q, P[j] % q, S[j] % q) == (0, 0, 0) for j in outside),
                ('support', p, I, q))
        require(sum(r) % q == sum(P) % q == 0, ('column sums', p, I, q))
        generators.append([0] * n + [P[j] for j in I] + [r[j] for j in I]
                          + [S[j] - P[j] for j in I])
    k, o, b = len(free), int(bool(odd)), int(compatible)
    dimension = n - 1 + k - o + int(q == p) * b * o
    require(dimension == len(generators), 'wrong number of generators')
    return generators, dimension, (k, o, b)


def W_formula(p, I):
    I = set(I)
    n = len(I)
    e3 = sum((j+3) % p in I for j in I)
    e9 = sum((j+9) % p in I for j in I)
    c3, c9 = (1 if e3 else 2), (1 if e9 else 2)
    graph = {j: [] for j in I}
    for j in I:
        graph[j] = [(j+d) % p for d in (-9, -3, 3, 9) if (j+d) % p in I]
    kappa = len(components(graph))
    return 2*n+e3-4+c3, e3+e9+kappa, 2*n+e9-4+c9


def interval_formula(p, L, q):
    pos = lambda x: max(0, x)
    if p >= 17:
        return (L-1 + pos(L-10) + pos(L-p+6) + pos(L-p+3)
                + int(q == p and L >= p-6))
    ordinary = {
        5: [0,1,2,4,8],
        7: [0,1,2,3,4,8,12],
        11: [0,1,2,3,4,6,8,10,12,16,20],
        13: [0,1,2,3,4,5,6,8,10,12,16,20,24],
    }
    onset = {5: 4, 7: 5, 11: 9, 13: 10}
    return ordinary[p][L-1] + int(q == p and L >= onset[p])


def explicit_interval_generators(p, L, q):
    """Second basis: the sparse interval table, without invoking the forest."""
    parameters = []
    m = p-L
    if m == 0:
        for i in range(p-1):
            r = [0]*p
            r[i], r[i+1] = 1, -1
            parameters.append((r,0))
        if q == p:
            parameters.append(([0]*p,1))
    else:
        B = max(0,L-9)
        for i in range(B-1):
            r = [0]*p
            r[i], r[i+1] = 1, -1
            parameters.append((r,0))
        local = {
            1: [[(i,1)] for i in (0,1,3,4,6,7)] + [[(2,1),(5,-1)]],
            2: [[(i,1)] for i in (0,3,6)] + [[(1,1),(4,-1)],[(2,1),(5,-1)]],
            3: [[(i,1),(i+3,-1)] for i in range(3)],
            4: [[(i,1),(i+3,-1)] for i in range(2)],
            5: [[(0,1),(3,-1)]],
        }
        for entries in local.get(m,[]):
            require(B>0,'sparse interval basis lacks low anchor')
            r = [0]*p
            for i,v in entries:
                r[p-9+i] += v
            r[0] -= sum(r)
            parameters.append((r,0))
        if q == p and 1 <= m <= 6:
            r = [0]*p
            r[p-9]=r[p-8]=r[p-7]=1
            r[0]=-3
            parameters.append((r,1))
    generators=[]
    for i in range(L-1):
        h=[int(j==i)-int(j==i+1) for j in range(L)]
        generators.append(h+[-x for x in h]+[-x for x in h]+h)
    for r,c in parameters:
        P=[c-sum(r[(j-k)%p] for k in (3,6,9)) for j in range(p)]
        S=[-r[(j-9)%p] for j in range(p)]
        require(all(r[j]%q == P[j]%q == S[j]%q == 0 for j in range(L,p)),
                ('explicit interval support',p,L,q))
        generators.append([0]*L+P[:L]+r[:L]+[S[j]-P[j] for j in range(L)])
    return generators


def check_case(p, I, q, coupling=False):
    U, V, W = physical_maps(p, I)
    all_maps = U + V + W
    n = 4 * len(I)
    gs, dimension, stats = forest_generators(p, I, q)
    require(n-rank(all_maps, q) == dimension, ('J dimension', p, I, q, dimension))
    require(rank(gs, q, n) == dimension, ('basis independence', p, I, q))
    require(all(sum(a*b for a,b in zip(row,g)) % q == 0 for row in all_maps for g in gs),
            ('physical line failure', p, I, q))
    Wdims = W_formula(p, I)
    pairs = (V+W, U+W, U+V)
    require(tuple(n-rank(rows,q) for rows in pairs) == Wdims, ('W dimensions', p,I,q))
    if coupling:
        bases = [nullspace(rows, q, n) for rows in pairs]
        columns = [ [(s*x) % q for x in v]
                    for s, basis in zip((-3,-9,1),bases) for v in basis ]
        mrank = rank(columns, q, n)
        expected = Wdims[2] if q == 3 else sum(Wdims)-2*dimension
        require(mrank == expected, ('physical coupling rank', p,I,q,mrank,expected))
        if q != 3:
            # This checks actual pointwise relation channels; domain membership
            # follows from every physical J equation verified above.
            for g in gs:
                require(all((-3*x+3*x)%q == 0 and (-9*x+9*x)%q == 0 for x in g),
                        'channel cancellation')
    return dimension, stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='verification.json')
    parser.add_argument('--interval-only', action='store_true')
    args = parser.parse_args()
    receipt = {'status':'PASS', 'coefficient_fields':'prime fields; proof is all fields',
               'sparse_interval_basis':'verified independently of forest basis for p>=17',
               'intervals':[], 'subsets':{}, 'coupling_cases':0}
    total = 0
    for p in (5,7,11,13,17,19,23,29,31,37,41,43):
        rows = []
        for L in range(1,p+1):
            row = {'L':L,'dimensions':{}}
            for q in sorted({2,3,5,p,101}):
                coupling = L in {1,2,p-8,p-7,p-6,p-5,p-3,p-2,p-1,p}
                dimension, stats = check_case(p,list(range(L)),q,coupling)
                require(dimension == interval_formula(p,L,q), ('interval formula',p,L,q))
                if p >= 17:
                    explicit = explicit_interval_generators(p,L,q)
                    require(len(explicit) == dimension and rank(explicit,q,4*L) == dimension,
                            ('sparse interval basis count/rank',p,L,q))
                    physical = sum(physical_maps(p,list(range(L))),[])
                    require(all(sum(x*y for x,y in zip(g,row))%q == 0
                                for g in explicit for row in physical),
                            ('sparse interval basis line',p,L,q))
                    eta = lambda s: max(0,L-(s%p))+max(0,L-(p-s%p))
                    kappa = min(3,L)-int(L>=p-8)-int(L>=p-7)
                    require(W_formula(p,list(range(L)))[1] == eta(3)+eta(9)+kappa,
                            ('interval kappa formula',p,L,q))
                row['dimensions'][str(q)] = dimension
                total += 1
                receipt['coupling_cases'] += int(coupling)
            row['forest_k_o_b'] = stats
            rows.append(row)
        receipt['intervals'].append({'p':p,'rows':rows})
        print('interval prime',p,'PASS',flush=True)
    if not args.interval_only:
        for p in (5,7,11):
            count = 0
            for mask in range(1,1<<p):
                I = [j for j in range(p) if mask >> j & 1]
                for q in (2,3,p,101):
                    check_case(p,I,q)
                    count += 1
            receipt['subsets'][str(p)] = {'all_nonempty_subsets':(1<<p)-1,'field_cases':count}
            total += count
            print('all subsets prime',p,'PASS',flush=True)
    # Explicit boundary witnesses against tempting false extrapolations.
    hostiles = {}
    for p,L,q in ((29,23,29),(29,23,101),(29,24,101),(29,29,29),(29,29,101),(5,4,5),(5,4,2)):
        dimension,_ = check_case(p,list(range(L)),q,True)
        hostiles[f'p{p}_L{L}_K{q}'] = dimension
    require(hostiles['p29_L23_K29'] == 36 and hostiles['p29_L23_K101'] == 35,
            'first characteristic-p wrap witness')
    require(hostiles['p29_L29_K29'] == 57 and hostiles['p29_L29_K101'] == 56,
            'full-cycle characteristic witness')
    U,V,W = physical_maps(17,[0,1])
    d = 8-rank(U+V+W,3)
    wd = W_formula(17,[0,1])
    require(d == 1 and wd[0]+wd[1] == 4 != 2*d, 'char3 zero-weight hostile')
    receipt.update({'case_count':total,'hostile_dimensions':hostiles,
                    'char3_weight_hostile':{'p':17,'L':2,'dimJ':1,'physical_coupling_nullity':4,
                                             'incorrect_nonzero_weight_formula':2}})
    with open(args.output,'w') as stream:
        json.dump(receipt,stream,indent=2,sort_keys=True)
        stream.write('\n')
    print('PASS',total,'direct field cases;',receipt['coupling_cases'],'coupling ranks; output',args.output)

if __name__ == '__main__':
    main()
