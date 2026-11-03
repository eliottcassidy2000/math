#!/usr/bin/env python3
"""Exact finite controls for S170 holonomy tomography. Standard library only."""
import itertools as it
import json
import math
from collections import Counter
from functools import lru_cache

def require(ok, label):
    if not ok:
        raise RuntimeError(label)

@lru_cache(None)
def sign(n, a, x):
    if n == 0:
        return 0
    half = 1 << (n - 1)
    u, r = a % half, a // half
    y, t = x % half, x // half
    c = 0 if not u else 1 + bool(y) + bool(y ^ u)
    return (sign(n - 1, u, y) + r * bool(y) + t * (c + r)) % 2

def left(n, a):
    size = 1 << n
    return tuple((x ^ a) + size * (s ^ sign(n, a, x))
                 for s in range(2) for x in range(size))

def compose(a, b):
    return tuple(a[x] for x in b)

def inverse(a):
    out = [0] * len(a)
    for x, y in enumerate(a):
        out[y] = x
    return tuple(out)

def commutator(a, b):
    return compose(compose(compose(a, b), inverse(a)), inverse(b))

def anf(values):
    out = list(values)
    for i in range(len(out).bit_length() - 1):
        for s in range(len(out)):
            if s >> i & 1:
                out[s] ^= out[s ^ (1 << i)]
    return out

def rank(rows):
    pivots = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)

def mask(values):
    return sum(int(b) << i for i, b in enumerate(values))

def monomial(d, a):
    return mask([(y & a) == a for y in range(1 << d)])

def rm_basis(d, k):
    return [monomial(d, a) for a in range(1 << d) if a.bit_count() <= k]

def diagonal_table(n, perm):
    size, M = 1 << n, 1 << (n - 1)
    h = [perm[y] // size for y in range(M)]
    require(all(perm[s * size + x] == (s ^ h[x % M]) * size + x
                for s in range(2) for x in range(size)), 'not a repeated diagonal')
    return h

def kernel_basis(d):
    M = 1 << d
    return [((1 << 0) | (1 << y)) * ((1 << M) + 1) for y in range(1, M)] + [((1 << M) - 1) << M]

def evaluation_rank(basis, points):
    return rank([sum(((b >> x) & 1) << i for i, b in enumerate(basis)) for x in points])

def expected_nullity(d, points):
    M = 1 << d
    owners = {}
    for x in points:
        owners.setdefault(x % M, set()).add(x // M)
    q = len(owners)
    doubled = any(len(ts) == 2 for ts in owners.values())
    if q < M:
        return M - q - int(doubled)
    if doubled:
        return 0
    return 1 - (sum(next(iter(ts)) for ts in owners.values()) % 2)

def seed_edges(bits):
    return tuple(bool(bits >> i & 1) for i, _ in enumerate(it.combinations(range(5), 2)))

def edge(matrix, a, b):
    if a > b:
        return not edge(matrix, b, a)
    return matrix[list(it.combinations(range(5), 2)).index((a, b))]

def seed_profile(bits):
    matrix = seed_edges(bits)
    covers = set()
    for p in it.permutations(range(5)):
        for cuts in range(16):
            paths, start, valid = [], 0, True
            for j in range(4):
                if cuts >> j & 1:
                    paths.append(p[start:j + 1])
                    start = j + 1
                elif not edge(matrix, p[j], p[j + 1]):
                    valid = False
                    break
            if valid:
                paths.append(p[start:])
                covers.add(tuple(sorted(paths)))
    profile = Counter(map(len, covers))
    require(profile == {1:9,2:33,3:31,4:10,5:1}, 'seed path profile')
    for start in range(5):
        seen = {start}
        for _ in range(5):
            seen |= {b for a in seen for b in range(5) if a != b and edge(matrix, a, b)}
        require(len(seen) == 5, 'seed is not strong')
    for u in range(8):
        value = sum(profile[j] * math.prod(range(u-j+1,u+1)) for j in profile if j <= u)
        require(value == u**5+6*u**3+2*u, 'slot polynomial')
    return profile

def tournament(word):
    N = 5 * len(word)
    out = [0] * N
    for a, b in it.combinations(range(N), 2):
        forward = a // 5 < b // 5 or edge(seed_edges(8 if not word[a//5] else 10), a%5, b%5)
        tail, head = (a,b) if forward else (b,a)
        out[tail] |= 1 << head
    return tuple(out)

def physical_trial(start):
    n, x, s = 4, start, 0
    control = [0,0,0]
    state = control + [s] + [(x >> i) & 1 for i in range(n)]
    initial = tuple(state)
    prev = tournament(state)
    arcs = []
    itinerary = [(0,4,1),(1,1,1),(2,10,1),(1,1,-1),(2,10,-1),(0,4,-1)]
    for owner,a,direction in itinerary:
        require(control[owner] == (direction < 0), 'illegal binary traversal')
        control[owner] ^= 1
        s ^= sign(n,a,x) ^ int(direction < 0)
        x ^= a
        target = control + [s] + [(x >> i) & 1 for i in range(n)]
        for j in range(len(state)):
            if state[j] != target[j]:
                state[j] ^= 1
                nxt = tournament(state)
                changed = [(u,v) for u,v in it.combinations(range(40),2)
                           if bool(prev[u] >> v & 1) != bool(nxt[u] >> v & 1)]
                require(changed == [(5*j,5*j+2)], 'not the claimed literal arc')
                arcs.extend(changed)
                prev = nxt
        require(state == target, 'register update failure')
    require(control == [0,0,0] and x == start, 'index/controller did not close')
    require(s == int(start in range(4,8)), 'shifted plane sign')
    return {'start_index':start,'final_sign':s,'literal_arc_reversals':len(arcs),
            'vertices':40,'blocks':8,'hamilton_count':9**8,'initial_word':initial,'final_word':state}

def main():
    report = {'source':'S170 holonomy code tests','levels':[]}
    for n in range(3,8):
        d, M = n-1, 1 << (n-1)
        L = [left(n,a) for a in range(1<<n)]
        comms = 0
        for a in range(1,1<<n):
            for b in range(a+1,1<<n):
                h = diagonal_table(n,commutator(L[a],L[b]))
                require(sum(h)%2 == 0, 'commutator parity')
                require(all(sum(h[y] for y in range(M) if y>>j&1)%2 == 0 for j in range(d)), 'commutator first moment')
                require(all(not c or a.bit_count() <= d-2 for a,c in enumerate(anf(h))), 'commutator degree')
                comms += 1
        # Every two-flat is the sum of three square commutator tables.
        flat_controls = 0
        Z = commutator(L[1], L[2])
        for u in range(1,M):
            for v in range(u+1,M):
                C = compose(compose(commutator(L[u], L[M+v]), commutator(L[u], L[v])), Z)
                h = diagonal_table(n,C)
                require(h == [int(y in (0,u,v,u^v)) for y in range(M)], 'affine plane circuit')
                flat_controls += 1
        layers = []
        W = commutator(L[1],L[M+2])
        for r in range(2,d+1):
            if r > 2:
                W = commutator(L[1 << (r-1)],W)
            h = diagonal_table(n,W)
            require(h == [int(y < (1<<r)) for y in range(M)], 'minimal-support word')
            basis = rm_basis(d,d-r)
            K = sum(math.comb(d,j) for j in range(r,d+1))
            info = [y for y in range(M) if y.bit_count() <= d-r]
            require(len(basis) == K and len(info) == K and rank(basis) == K, 'RM dimension')
            require(evaluation_rank(basis,info) == K, 'information set')
            # I^r is generated by translates of coordinate r-flats.
            flats = []
            for directions in it.combinations(range(d),r):
                subspace = [sum(((z>>j)&1)<<i for j,i in enumerate(directions)) for z in range(1<<r)]
                for offset in range(M):
                    flats.append(sum(1 << (offset^y) for y in subspace))
            require(rank(flats) == K and rank(flats+basis) == K, 'augmentation/RM spaces')
            minimum = None
            if d <= 4:
                words = [0]
                for b in basis:
                    words += [v^b for v in words]
                minimum = min(v.bit_count() for v in words if v)
                require(minimum == 1<<r, 'minimum code distance')
            layers.append({'r':r,'probe_count':K,'distance':1<<r,'moved_signed_sheets':1<<(r+2),
                           'enumerated_minimum_distance':minimum})
        kb = kernel_basis(d)
        bank = list(range(M-1))+[M]
        require(evaluation_rank(kb,bank) == M, 'full group optimal bank')
        require(evaluation_rank(kb,range(M)) == M-1, 'lower half hostile')
        subsets = range(1<<(2*M)) if n==3 else range(1024)
        checked = 0
        for v in subsets:
            word = v if n==3 else ((v*0x9e3779b97f4a7c15)^((v+1)*0xbf58476d1ce4e5b9)) & ((1<<(2*M))-1)
            points = [x for x in range(2*M) if word>>x&1]
            require(M-evaluation_rank(kb,points) == expected_nullity(d,points), 'probe classification')
            checked += 1
        if n <= 4:
            kernel = [0]
            for b in kb:
                kernel += [v^b for v in kernel]
            responses=set()
            for a in range(2*M):
                for g in kernel:
                    response=(a,tuple(sign(n,a,x)^((g>>(x^a))&1) for x in bank))
                    require(response not in responses, 'full-group response collision')
                    responses.add(response)
            require(len(responses) == 1 << (M+n), 'group response cardinality')
        report['levels'].append({'n':n,'signed_sheets':4*M,'generator_pairs':comms,
                                 'full_group_probes':M,'probe_sets_checked':checked,'flat_circuits_checked':flat_controls,'layers':layers})
    for bits in [8,10]:
        seed_profile(bits)
    require(not any(all(edge(seed_edges(8),a,b)==edge(seed_edges(10),p[a],p[b])
                        for a,b in it.combinations(range(5),2)) for p in it.permutations(range(5))), 'seeds isomorphic')
    report['physical_trials']=[physical_trial(0),physical_trial(4)]
    report['status']='PASS'
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
