"""Exact R4 boundary detours. Explicit checks survive python -O."""
from itertools import combinations, permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R4 = (2, 4, 9, 3)
CASES = [
    ((2,4,9,3), [('A',1,2,3,0),(2,'A',1,3,0),(2,'A',0,1,3)]),
    ((2,12,9,1), [('A',2,0,1,3),(1,3,'A',2,0),(1,2,3,'A',0)]),
    ((6,12,8,1), [('A',1,2,3,0),(3,'A',1,2,0),(0,3,'A',1,2),
                   (0,2,3,'A',1),(2,3,'A',0,1)]),
    ((6,4,8,3), [('A',1,2,3,0),(2,'A',1,3,0),(2,'A',0,1,3)]),
    ((10,4,9,2), [('A',1,2,0,3),(2,0,3,'A',1),(2,3,'A',0,1)]),
]
R5 = tuple(sum(1 << v for v in range(5) if (v-u) % 5 in (1,2))
           for u in range(5))
CASES_R5 = [
    ((6,12,24,17,3), [('A',3,0,1,2,4),(1,2,4,'A',3,0),(1,2,3,4,'A',0)]),
    ((18,12,9,17,6), [('A',0,1,3,4,2),(4,'A',0,3,2,1),(2,4,'A',0,3,1)]),
    ((10,20,9,18,5), [('A',0,1,4,2,3),(4,3,'A',0,2,1),(3,'A',0,2,4,1)]),
]
EXPECTED = [[(1,-5),(-1,6)], [(1,-4),(-1,4)],
            [(1,-5),(-1,4),(1,-4),(-1,4)], [(1,-5),(-1,5)],
            [(1,-4),(-1,4)]]


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def induced(a, p):
    return tuple(sum(1 << j for j,v in enumerate(p) if a[u] >> v & 1) for u in p)


def transpose(a):
    n = len(a)
    return tuple(sum(1 << j for j in range(n) if a[n-1-j] >> (n-1-i) & 1)
                 for i in range(n))


def join(parts):
    n = sum(map(len, parts))
    out, offset = [], 0
    for a in parts:
        later = sum(1 << j for j in range(offset + len(a), n))
        out.extend((r << offset) | later for r in a)
        offset += len(a)
    return tuple(out)


def is_path(a, p):
    return sorted(p) == list(range(len(a))) and all(a[u] >> v & 1 for u,v in zip(p,p[1:]))


def blue(a, p):
    read = induced(a, p)
    return read == transpose(read)


def edge_set(p):
    return set(tuple(sorted(e)) for e in zip(p,p[1:]))


def complement(a, p):
    need(is_path(a,p), 'K input is not a Hamiltonian path')
    edges = edge_set(p)
    out = list(a)
    for u,v in combinations(range(len(a)),2):
        if (u,v) not in edges:
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    return tuple(out)


def expand(word, size):
    return tuple(v for token in word for v in
                 (range(size) if token == 'A' else [size+token]))


def blue_tilings(n):
    free = [(u,v) for u,v in combinations(range(n),2) if v != u+1]
    seen, orbits = set(), []
    for pair in free:
        if pair in seen:
            continue
        u,v = pair
        orbit = {pair, (n-1-v,n-1-u)}
        seen.update(orbit)
        orbits.append(orbit)
    for mask in range(1 << len(orbits)):
        out = [0] * n
        for u in range(n-1):
            out[u] |= 1 << (u+1)
        for j,orbit in enumerate(orbits):
            for u,v in orbit:
                a,b = (u,v) if mask >> j & 1 else (v,u)
                out[a] |= 1 << b
        yield tuple(out)


def symbolic_checks():
    hp = [p for p in permutations(range(4)) if is_path(R4,p)]
    need(hp == [(0,1,2,3),(1,2,3,0),(2,3,0,1),(3,0,1,2),(3,1,2,0)],
         'R4 path list changed')
    need([induced(R4,p) for p in hp] == [b for b,_ in CASES], 'R4 presentation list')
    need(sum(induced(R4,p) == R4 for p in permutations(range(4))) == 1,
         'R4 has a nontrivial automorphism')
    result = []
    for index,(b,words) in enumerate(CASES):
        # Two distinct formal A ports suffice for all path and edge-XOR identities.
        initial = join([(2,0),b])
        current = initial
        parity = set()
        for word in words:
            p = expand(word,2)
            parity ^= edge_set(p)
            current = complement(current,p)
        need(parity == edge_set(tuple(range(6))), 'formal boundary edge XOR failed')
        need(current == complement(initial,tuple(range(6))), 'formal target differs')
        # Affine scores alpha*n+beta; A denotes its first vertex.
        scores = {i:(0,b[i].bit_count()) for i in range(4)}
        scores['A'] = (1,-1-b[3].bit_count())
        differences = []
        for word in words:
            first,last = word[0],word[-1]
            need(last != 'A', 'affine score proxy omits the last A vertex')
            differences.append((scores[first][0]+scores[last][0]-1,
                                scores[first][1]+scores[last][1]+1))
            scores = {v:(1-a,-1-c+int(v==first)-int(v==last))
                      for v,(a,c) in scores.items()}
        need(differences[1:] == EXPECTED[index], 'symbolic endpoint formula differs')
        result.append(dict(B=b, words=words, endpoint_differences=differences))
    return result


def r5_symbolic_checks():
    hp = [p for p in permutations(range(5)) if is_path(R5,p)]
    need(len(hp) == 15, 'R5 Hamiltonian path count')
    need(set(induced(R5,p) for p in hp) == {b for b,_ in CASES_R5},
         'R5 path orbit presentations incomplete')
    need(all(sum(induced(R5,p) == b for p in hp) == 5 for b,_ in CASES_R5),
         'R5 path orbit size differs')
    result = []
    for b,words in CASES_R5:
        changed_right = tuple(words[0][1:])
        need(induced(b,changed_right) != b, 'first right path remains in its old orbit')
        initial = join([(2,0),b])
        current = initial
        parity = set()
        for word in words:
            p = expand(word,2)
            parity ^= edge_set(p)
            current = complement(current,p)
        need(parity == edge_set(tuple(range(7))), 'R5 formal boundary XOR failed')
        need(current == complement(initial,tuple(range(7))), 'R5 formal target differs')
        scores = {i:(0,2) for i in range(5)}
        scores['A'] = (1,-3)
        differences = []
        for word in words:
            first,last = word[0],word[-1]
            need(last != 'A', 'R5 score proxy omits last A vertex')
            differences.append((scores[first][0]+scores[last][0]-1,
                                scores[first][1]+scores[last][1]+1))
            scores = {v:(1-a,-1-c+int(v==first)-int(v==last))
                      for v,(a,c) in scores.items()}
        need(differences == [(0,0),(1,-5),(-1,5)], 'R5 endpoint formulas differ')
        result.append(dict(B=b, words=words, endpoint_differences=differences))
    return result


def full_controls(cases, family):
    records = []
    for middle_n in range(8):
        middles = certificates = steps = 0
        for z in blue_tilings(middle_n):
            middles += 1
            need(blue(z,tuple(range(middle_n))), 'middle generator not blue')
            for b,words in cases:
                source = join([transpose(b),z,b])
                n = len(source)
                standard = tuple(range(n))
                need(blue(source,standard), 'source path not blue')
                target = complement(source,standard)
                actual = source
                for word in words:
                    p = expand(word,n-len(b))
                    need(is_path(actual,p), 'lifted macro has absent arc')
                    need(not blue(actual,p), 'lifted macro has blue step')
                    actual = complement(actual,p)
                    steps += 1
                need(actual == target, 'literal labelled target failed')
                certificates += 1
        records.append(dict(family=family, middle_order=middle_n, blue_middle_presentations=middles,
                            all_outer_presentation_certificates=certificates,
                            black_steps=steps))
    return records


def main():
    report = dict(formal_R4_boundary_identities=symbolic_checks(),
                  formal_R5_boundary_identities=r5_symbolic_checks(),
                  complete_blue_middle_controls=full_controls(CASES,'R4') +
                                                full_controls(CASES_R5,'R5'))
    result = json.dumps(report,indent=2) + '\n'
    (HERE / 'moves.json').write_text(result,encoding='utf-8')
    print(result,end='')


if __name__ == '__main__':
    main()
