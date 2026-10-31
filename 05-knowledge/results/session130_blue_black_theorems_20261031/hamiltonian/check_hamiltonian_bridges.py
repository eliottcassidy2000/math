#!/usr/bin/env python3
"""Stdlib exact checker. No repository, network, or assert dependency."""
import argparse
import itertools as it
import json
from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, mask):
    a = [[0] * n for _ in range(n)]
    for k, (u, v) in enumerate(it.combinations(range(n), 2)):
        a[u][v] = (mask >> k) & 1
        a[v][u] = 1 - a[u][v]
    return a


def encode(a, order=None):
    order = list(range(len(a))) if order is None else order
    return sum(a[order[u]][order[v]] << k for k, (u, v) in
               enumerate(it.combinations(range(len(a)), 2)))


def canonical(a):
    return min(encode(a, p) for p in it.permutations(range(len(a))))


def ham_dp(a):
    n = len(a)
    if n == 0:
        return 1
    d = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for s in range(1, 1 << n):
        for u in range(n):
            if d[s][u]:
                for v in range(n):
                    if not ((s >> v) & 1) and a[u][v]:
                        d[s | (1 << v)][v] += d[s][u]
    return sum(d[-1])


def paths(a):
    return [p for p in it.permutations(range(len(a)))
            if all(a[u][v] for u, v in zip(p, p[1:]))]


def normal_cycle(p):
    k = p.index(min(p))
    return tuple(p[k:] + p[:k])


def odd_cycles(a):
    """Each directed cycle once, oriented; rotational equivalence only."""
    n = len(a)
    result = []
    for k in range(3, n + 1, 2):
        for ss in it.combinations(range(n), k):
            for tail in it.permutations(ss[1:]):
                p = (ss[0],) + tail
                if all(a[p[j]][p[(j + 1) % k]] for j in range(k)):
                    result.append(p)
    return result


def cycle_data(a, rho=None):
    cc = odd_cycles(a)
    supports = [sum(1 << v for v in c) for c in cc]
    alpha2 = sum(not (u & v) for k, u in enumerate(supports) for v in supports[k+1:])
    result = {'cycles': len(cc), 'c3': sum(len(c) == 3 for c in cc),
              'c5': sum(len(c) == 5 for c in cc), 'alpha2': alpha2}
    if rho is not None:
        lookup = {c: k for k, c in enumerate(cc)}
        q = intersecting = disjoint = 0
        for k, c in enumerate(cc):
            image = normal_cycle(tuple(rho[v] for v in reversed(c)))
            j = lookup[image]
            if j == k:
                q += 1
            elif j > k:
                if supports[k] & supports[j]:
                    intersecting += 1
                else:
                    disjoint += 1
        result.update(q=q, mirror_intersecting=intersecting, mirror_disjoint=disjoint)
        fixed_vertices = [v for v in range(len(a)) if rho[v] == v]
        if fixed_vertices:
            center = fixed_vertices[0]
            center_cycles = sum(center in c for c in cc)
            require((center_cycles-q) % 2 == 0, 'center orbit parity')
            result['center_nonfixed_pairs'] = (center_cycles-q)//2
        require(alpha2 % 2 == disjoint % 2, 'fixed two-packings mismatch')
    return result


def tile_matrix(n, mask):
    a = [[0] * n for _ in range(n)]
    tiles = [(u, v) for u in range(n) for v in range(u+2, n)]
    for u in range(n-1):
        a[u][u+1] = 1
    for k, (u, v) in enumerate(tiles):
        a[u][v] = (mask >> k) & 1
    for u in range(n):
        for v in range(u+1, n):
            a[v][u] = 1 - a[u][v]
    return a


def blue_masks(n):
    tiles = [(u, v) for u in range(n) for v in range(u+2, n)]
    indices = {p:k for k,p in enumerate(tiles)}
    groups = []
    for k, (u, v) in enumerate(tiles):
        j = indices[(n-1-v, n-1-u)]
        if k <= j:
            groups.append((1 << k) | (1 << j))
    for bits in range(1 << len(groups)):
        yield sum(g for k, g in enumerate(groups) if (bits >> k) & 1)


def blue_check(n, masks):
    distribution = Counter()
    count = 0
    witnesses = {}
    for mask in masks:
        a = tile_matrix(n, mask)
        rho = tuple(reversed(range(n)))
        h = ham_dp(a)
        data = cycle_data(a, rho)
        rhs = 1 + 2 * data['q'] + 4 * data['mirror_intersecting']
        require((h-rhs) % 8 == 0, ('mod8', n, mask, h, data))
        require(h % 8 == (1 + 2*data['cycles'] + 4*data['alpha2']) % 8,
                ('OCF truncation', n, mask))
        if n % 2 == 0:
            require(data['q'] == 0 and h % 4 == 1, ('even stratum', n, mask))
        if n <= 5:
            require(h == 1 + 2*data['cycles'], ('exact low order', n, mask))
        if n % 2:
            z = n//2
            vertices = [v for v in range(n) if v != z]
            base = [[a[u][v] for v in vertices] for u in vertices]
            hb = ham_dp(base)
            require((h-hb-2*data['q']-4*data['center_nonfixed_pairs'])%8 == 0,
                    ('center deletion mod8',n,mask))
            b = pair_cover_weights(base)
            ss = {k for k,v in enumerate(vertices) if a[z][v]}
            crossing = sum(b[u][v] for u in ss for v in range(n-1)
                           if v not in ss and base[u][v])
            require(h-hb == 2*crossing, ('canonical cut',n,mask))
            mirror = sum(b[u][n-2-u] for u in ss if base[u][n-2-u])
            require((mirror-data['q'])%2==0, ('cycle-cut parity bridge',n,mask))
        distribution[h % 8] += 1
        witnesses.setdefault(str(h % 8), {'tile_mask': mask, 'H': h, **data})
        count += 1
    return {'n': n, 'tested': count, 'H_mod8': dict(sorted(distribution.items())),
            'first_residue_witnesses': witnesses}


def pair_cover_weights(a):
    """Independent enumeration of ordered two-path covers, allowing the cut arc."""
    n=len(a)
    w=[[0]*n for _ in range(n)]
    for p in it.permutations(range(n)):
        bad=[k for k in range(n-1) if not a[p[k]][p[k+1]]]
        if len(bad)>1:
            continue
        for k in (bad if bad else range(n-1)):
            w[p[k]][p[k+1]]+=1
    b=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1,n):
            require((w[u][v]+w[v][u])%2==0,'pair integrality')
            b[u][v]=b[v][u]=(w[u][v]+w[v][u])//2
    return b


def n5_table():
    n = 5
    tiles = [(u, v) for u in range(n) for v in range(u+2, n)]
    index = {e:k for k,e in enumerate(tiles)}
    def reflect(mask):
        return sum(((mask >> k)&1) << index[(n-1-v,n-1-u)]
                   for k,(u,v) in enumerate(tiles))
    cl = {}
    row = {}
    for mask in range(64):
        a = tile_matrix(n,mask)
        c = canonical(a)
        cl[mask] = c
        if c not in row:
            aa = tournament(n,c)
            pp = paths(aa)
            require(len(pp)==ham_dp(aa), ('path-dp',c))
            aut = sum(encode(aa,p)==c for p in it.permutations(range(n)))
            converse = canonical([[aa[v][u] for v in range(n)] for u in range(n)])
            row[c] = {'canonical': c, 'score': sorted(map(sum,aa)), 'H':len(pp),
                      'aut':aut, 'converse':converse, 'fibre':0, 'blue':0,
                      **cycle_data(aa)}
        row[c]['fibre'] += 1
        if reflect(mask)==mask:
            row[c]['blue'] += 1
            data = cycle_data(a, tuple(reversed(range(n))))
            if 'q' in row[c]:
                require(row[c]['q']==data['q'], ('fixed cycle class invariant',c))
                require(row[c]['mirror_intersecting']==data['mirror_intersecting'],
                        ('mirror intersection class invariant',c))
            row[c].update({k:data[k] for k in ('q','mirror_intersecting','mirror_disjoint')})
    merged = lambda c: min(c,row[c]['converse'])
    black = Counter()
    blue = Counter()
    for mask in range(64):
        comp = 63-mask
        if mask > comp:
            continue
        endpoints=tuple(sorted((merged(cl[mask]),merged(cl[comp]))))
        (blue if reflect(mask)==mask else black)[endpoints] += 1
    require(all(v%2==0 for v in black.values()), 'black twins')
    fundamental = {p:v//2 for p,v in black.items()}
    for c,r in row.items():
        require(r['H']==r['aut']*r['fibre'], ('fibre identity',c))
        if r['converse']==c:
            d=(r['fibre']-r['blue'])//2
            r['fundamental_black_degree']=d
            require((2*r['aut']*d+r['aut']*r['blue']-1-2*r['q']-
                     4*r['mirror_intersecting'])%8==0, ('normalized mod8',c))
        else:
            r['fundamental_black_degree']=r['fibre']
        if merged(c)==c:
            measured=sum(v*((p[0]==c)+(p[1]==c)) for p,v in fundamental.items())
            require(measured==r['fundamental_black_degree'], ('degree',c,measured,r))
    require(len(row)==12 and sum(r['blue'] for r in row.values())==16,'n5 controls')
    return {'classes': list(row.values()),
            'merged_nodes': sorted(set(map(merged,row))),
            'blue_lines': [{'endpoints':p,'multiplicity':v} for p,v in sorted(blue.items())],
            'fundamental_black_edges': [{'endpoints':p,'multiplicity':v}
                                        for p,v in sorted(fundamental.items())]}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-blue-n',type=int,default=7)
    args=parser.parse_args()
    out={'status':'PASS','conventions':{'vertices':'0..n-1','base_path':'0->1->...->n-1',
         'tile_bits':'lex pairs (u,v), v>=u+2; bit1 means u->v',
         'canonical':'minimum integer upper-triangle code over all label permutations',
         'cycles':'directed simple odd cycles length>=3, modulo rotations, not reversal'},
         'blue_exhaustive':[blue_check(n,blue_masks(n)) for n in range(3,args.max_blue_n+1)],
         'n5':n5_table()}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
