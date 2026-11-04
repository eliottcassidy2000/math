#!/usr/bin/env python3
"""Independent first-two-layer birth audit and small one-exit C2 controls.

This imports no research producer. Python integers and recursive endpoint
deletion rebuild every cut on the actual first two layers. Literal path
permutations verify the one-exit joint identity via actual C2 outputs.
"""
from functools import lru_cache
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path

GATES = 0
SOURCE = (1790,1668,1672,1778,1670,1942,1718,1280,1631,128,512)


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def birth(rows, incoming):
    n = len(rows)
    return tuple(row | ((1 << n) if incoming >> i & 1 else 0)
                 for i, row in enumerate(rows)) + (((1 << n)-1)^incoming,)


def literal(rows):
    n = len(rows)
    H = [0] * (1 << n)
    E = [[0] * n for _ in H]
    H[0] = 1
    for length in range(1, n+1):
        for path in permutations(range(n), length):
            if all(rows[a] >> b & 1 for a, b in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                H[mask] += 1
                E[mask][path[-1]] += 1
    return H, E


def F(tab, root):
    h, e = tab
    full = len(h)-1
    return sum(e[a][root]*h[full^a] for a in range(full+1))


def opposite(rows):
    full = (1 << len(rows))-1
    return tuple(full ^ (1 << i) ^ row for i, row in enumerate(rows))


def allcuts(rows):
    n, full = len(rows), (1 << len(rows))-1
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1)
                for v in range(n)]
    @lru_cache(None)
    def end(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(end(rest, u) for u in range(n) if rest & incoming[v] >> 0 & (1 << u))
    @lru_cache(None)
    def start(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(start(rest, u) for u in range(n) if rest & rows[v] & (1 << u))
    ee = [[end(mask, v) if mask >> v & 1 else 0 for v in range(n)]
          for mask in range(full+1)]
    ss = [[start(mask, v) if mask >> v & 1 else 0 for v in range(n)]
          for mask in range(full+1)]
    K = [[sum(ee[a][i]*ss[full^a][j] for a in range(full+1))
          for j in range(n)] for i in range(n)]
    h = sum(ee[full])
    # Different evaluation path: update incoming/outgoing cut sums directly.
    values = []
    for cut in range(full+1):
        inside = [i for i in range(n) if cut >> i & 1]
        outside = [i for i in range(n) if not (cut >> i & 1)]
        value = sum(ee[full][i] for i in inside)+sum(ss[full][j] for j in outside)
        value += sum(K[i][j] for i in inside for j in outside)
        need(value >= h, 'literal-cut monotonicity')
        values.append(value)
    need(values[0] == values[-1] == h, 'universal cut boundary')
    return values


def main():
    cases = 0
    for n in range(2, 5):
        pairs = list(combinations(range(n), 2))
        for code in range(1 << len(pairs)):
            rows = [0]*n
            for k, (i,j) in enumerate(pairs):
                if code >> k & 1: rows[i] |= 1 << j
                else: rows[j] |= 1 << i
            original = literal(rows)
            for a in range(n):
                if rows[a].bit_count() != 1: continue
                b = rows[a].bit_length()-1
                dual = literal(opposite(rows))
                sigma = F(dual,a)-original[0][-1]
                U = ((1 << n)-1) ^ (1 << a)
                need(sigma == original[1][U][b], 'one-exit start deficit')
                first = birth(rows, ((1 << n)-1) ^ (1 << a))
                child = birth(first, (1 << n)-1)
                observed = literal(child)
                for r in range(n):
                    if r == a: continue
                    J = sum(original[1][m][r]*original[1][U^m][b]
                            for m in range(U+1) if m & U == m)
                    need(F(observed,r)-3*F(original,r) == 2*J,
                         'joint identity through actual C2 graph')
                    need(J == 0 if r == b else J >= sigma, 'wrong-exit joint floor')
                    cases += 1
    parents = [SOURCE]
    levels = []
    for depth in (1, 2):
        retained, raw_hits, cuts = [], [0,0,0], 0
        for parent in parents:
            values = allcuts(parent)
            cuts += len(values)
            for cut, h in enumerate(values):
                if h <= 4291:
                    retained.append(birth(parent, cut))
                    if h in (4287,4289,4291):raw_hits[(h-4287)//2] += 1
        need(cuts == (2048 if depth == 1 else 151552), 'complete independent cut universe')
        need(len(retained) == (37 if depth == 1 else 164), 'independent retained universe')
        need(raw_hits[1] == 0, 'independent no-plus-two boundary')
        levels.append({'depth':depth,'cuts':cuts,'retained':len(retained),'raw_hits_0_2_4':raw_hits})
        parents = retained
    out = {'status':'PASS','gates':GATES,'one_exit_literal_cases':cases,
           'independent_induced_levels':levels,
           'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    raw=(json.dumps(out,sort_keys=True,indent=2)+'\n').encode()
    Path(__file__).with_suffix('.json').write_bytes(raw)
    print(json.dumps({'gates':GATES,'cases':cases,'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
