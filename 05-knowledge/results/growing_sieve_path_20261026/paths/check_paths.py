#!/usr/bin/env python3
"""Exact complete labelled tournaments n<=6; no assertions disabled by -O."""
import hashlib
import itertools
import json
from collections import Counter

GATES = 0

def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)

def tournament(n, mask):
    a = [0] * n
    for k, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if mask >> k & 1:
            a[u] |= 1 << v
        else:
            a[v] |= 1 << u
    return a

def reach(a, vertices=None):
    n = len(a)
    bits = (1 << n) - 1 if vertices is None else vertices
    r = [(a[u] | 1 << u) & bits if bits >> u & 1 else 0 for u in range(n)]
    for k in range(n):
        for u in range(n):
            if r[u] >> k & 1:
                r[u] |= r[k]
    return r

def components(a, vertices=None):
    n = len(a)
    bits = (1 << n) - 1 if vertices is None else vertices
    r = reach(a, bits)
    todo = bits
    cc = []
    while todo:
        u = (todo & -todo).bit_length() - 1
        c = sum(1 << v for v in range(n) if todo >> v & 1 and r[u] >> v & 1 and r[v] >> u & 1)
        cc.append(c)
        todo &= ~c
    cc.sort(key=lambda c: -r[(c & -c).bit_length() - 1].bit_count())
    return cc, r

def independent_paths(a):
    n = len(a)
    paths = []
    def walk(path, left):
        if not left:
            paths.append(tuple(path))
            return
        choices = left if not path else a[path[-1]] & left
        while choices:
            bit = choices & -choices
            v = bit.bit_length() - 1
            walk(path + [v], left ^ bit)
            choices ^= bit
    walk([], (1 << n) - 1)
    return paths

def reconstruct(u):
    """Promise input U(T). Only u is used by this reconstruction."""
    n = len(u)
    cc, r = components(u)
    owner = {v: k for k, c in enumerate(cc) for v in range(n) if c >> v & 1}
    cutreach = {k: [(z, reach(u, c & ~(1 << z))) for z in range(n) if c >> z & 1]
                for k, c in enumerate(cc) if c.bit_count() >= 3}
    a = u.copy()
    certificates = []
    for x, y in itertools.combinations(range(n), 2):
        if u[x] >> y & 1 or u[y] >> x & 1:
            continue
        if owner[x] != owner[y]:
            if r[x] >> y & 1:
                a[x] |= 1 << y
                certificates.append([x, y, 'outer'])
            elif r[y] >> x & 1:
                a[y] |= 1 << x
                certificates.append([y, x, 'outer'])
            else:
                raise RuntimeError('unordered outer SCCs')
            continue
        orientations = []
        for z, rz in cutreach[owner[x]]:
            if z in (x, y):
                continue
            xy, yx = bool(rz[x] >> y & 1), bool(rz[y] >> x & 1)
            if xy != yx:
                orientations.append((x, y, z) if xy else (y, x, z))
        if not orientations or len({v[:2] for v in orientations}) != 1:
            raise RuntimeError('no consistent one-vertex cut orientation')
        v, w, z = orientations[0]
        a[v] |= 1 << w
        certificates.append([v, w, z])
    return a, certificates

def subcube_languages(n):
    pairs = list(itertools.combinations(range(n), 2))
    indices = {p: k for k, p in enumerate(pairs)}
    full = (1 << len(pairs)) - 1
    unions, languages = [0] * (full + 1), [0] * (full + 1)
    perms = list(itertools.permutations(range(n)))
    for pi, p in enumerate(perms):
        fixed = ones = used = 0
        for x, y in zip(p, p[1:]):
            bit = 1 << indices[min(x, y), max(x, y)]
            fixed |= bit
            if x < y:
                ones |= bit
            used |= 1 << (x * n + y)
        free = full ^ fixed
        s = free
        while True:
            m = ones | s
            unions[m] |= used
            languages[m] |= 1 << pi
            if not s:
                break
            s = (s - 1) & free
    return perms, unions, languages

def run():
    tables, witnesses = [], {}
    for n in range(1, 7):
        perms, unions, languages = subcube_languages(n)
        index = {p: k for k, p in enumerate(perms)}
        require(len(set(unions)) == len(unions), ('union injectivity', n))
        require(len(set(languages)) == len(languages), ('language injectivity', n))
        counts = Counter(tournaments=len(unions), languages=len(set(languages)), unions=len(set(unions)))
        for mask, union in enumerate(unions):
            a = tournament(n, mask)
            u = [(union >> (v * n)) & ((1 << n) - 1) for v in range(n)]
            paths = independent_paths(a)
            direct_union = [0] * n
            direct_language = 0
            for p in paths:
                direct_language |= 1 << index[p]
                for x, y in zip(p, p[1:]):
                    direct_union[x] |= 1 << y
            require(direct_union == u, ('DFS union', n, mask))
            require(direct_language == languages[mask], ('DFS language', n, mask))
            counts['paths'] += len(paths)
            restored, cert = reconstruct(u)
            require(restored == a, ('cut reconstruction', n, mask))
            counts['missing_arcs'] += len(cert)
            ca, _ = components(a)
            cu, _ = components(u)
            require(ca == cu, ('outer SCCs', n, mask))
            if len(ca) == 1 and n >= 3:
                counts['strong'] += 1
            for c in ca:
                if c.bit_count() < 3:
                    continue
                for z in range(n):
                    if c >> z & 1:
                        keep = c & ~(1 << z)
                        require(reach(a, keep) == reach(u, keep), ('deletion lemma', n, mask, c, z))
                        counts['deletion_checks'] += 1
                for x in range(n):
                    if not c >> x & 1:
                        continue
                    for y in range(n):
                        if not (c >> y & 1 and a[x] >> y & 1 and not u[x] >> y & 1):
                            continue
                        good = []
                        for z in range(n):
                            if z in (x, y) or not c >> z & 1:
                                continue
                            blocks, _ = components(a, c & ~(1 << z))
                            if len(blocks) >= 4 and blocks[0] >> x & 1 and blocks[-1] >> y & 1:
                                second, penult = blocks[1], blocks[-2]
                                if all(a[v] >> z & 1 for v in range(n) if second >> v & 1) and a[z] & penult == penult:
                                    good.append((z, blocks))
                        require(bool(good), ('classical 4.12 certificate', n, mask, x, y))
                        counts['internal_missing_arcs'] += 1
                        if 'first_strong_missing' not in witnesses:
                            witnesses['first_strong_missing'] = dict(n=n, mask=mask, rows=a, used_rows=u,
                                arc=[x, y], separator=good[0][0], blocks=good[0][1], paths=paths)
        tables.append(dict(n=n, **dict(sorted(counts.items()))))
    # Hostile controls distinguishing full language, one path, reachability, and carrier.
    trans = [6, 4, 0]
    chain = [2, 4, 0]
    cycle = [2, 4, 1]
    require(independent_paths(trans) == independent_paths(chain), 'arbitrary digraph noninjectivity')
    require(independent_paths(trans)[0] in independent_paths(cycle), 'one path insufficient')
    require(reach(trans, 5) != reach(chain, 5), 'strong hypothesis needed')
    require(reconstruct(chain)[0] == trans, 'transitive3 hidden long arc')
    missing = witnesses['first_strong_missing']
    require(missing['n'] == 5 and missing['mask'] == 10, 'minimal strong hostile')
    require(all(r == 31 for r in reach(missing['used_rows'])), 'ordinary reachability loses direction')
    require(any(p.index(1) < p.index(4) for p in missing['paths']), 'unused does not mean fixed relative path order')
    result = dict(status='FINITE-EXACT PASS', universe='all labelled tournaments 1<=n<=6; no quotient or sampling',
        gates=GATES, tables=tables, witnesses=witnesses,
        encoding='lexicographic unordered pairs, low bit first, 1 means smaller label -> larger label',
        algorithms=['path-constraint subcube filling', 'independent vertex DFS', 'SCC and vertex-deletion cut closure'])
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    run()
