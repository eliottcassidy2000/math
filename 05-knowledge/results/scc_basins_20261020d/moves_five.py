"""Uniform three-black detours for all strong order-five outer SCCs.

The search is finite on a two-port prefix and five right vertices. Every
accepted score certificate is symbolic for all n >= 10. Explicit checks
remain active under python -O.
"""
from collections import Counter
from itertools import combinations, permutations
import json
from pathlib import Path
from moves import (need, induced, transpose, join, is_path, blue, edge_set,
                   complement, expand, blue_tilings)


def strong(a):
    n = len(a)
    full = (1 << n) - 1
    for v in range(n):
        seen = 1 << v
        while True:
            new = seen
            for u in range(n):
                if seen >> u & 1:
                    new |= a[u]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def tilings(n):
    free = [(u,v) for u,v in combinations(range(n),2) if v != u+1]
    for mask in range(1 << len(free)):
        out = [1 << (u+1) if u+1 < n else 0 for u in range(n)]
        for i,(u,v) in enumerate(free):
            if mask >> i & 1:
                out[u] |= 1 << v
            else:
                out[v] |= 1 << u
        yield tuple(out)


def score_step(scores, word):
    first,last = word[0],word[-1]
    return {v:(1-a,-1-c+int(v == first)-int(v == last))
            for v,(a,c) in scores.items()}


def score_difference(scores, word):
    first,last = word[0],word[-1]
    return (scores[first][0]+scores[last][0]-1,
            scores[first][1]+scores[last][1]+1)


def no_integer_zero(difference, lower):
    a,c = difference
    if a == 0:
        return c != 0
    return c % a != 0 or -c // a < lower


def compile_certificates():
    h = 5
    words = [w for w in permutations(['A']+list(range(h))) if w[-1] != 'A']
    paths = {w:expand(w,2) for w in words}
    edges = {w:frozenset(edge_set(paths[w])) for w in words}
    byedge = {e:w for w,e in edges.items()}
    standard = tuple(range(h+2))
    e0 = frozenset(edge_set(standard))
    certificates = []
    all_right = list(tilings(h))
    need(len(all_right) == 64, 'fixed-path universe is not all 64 presentations')
    for b in all_right:
        if not strong(b):
            continue
        initial = join([(2,0),b])
        scores = {i:(0,b[i].bit_count()) for i in range(h)}
        scores['A'] = (1,-1-b[-1].bit_count())
        found = None
        for p in permutations(range(h)):
            if not is_path(b,p) or induced(b,p) == b:
                continue
            w1 = ('A',)+p
            state1 = complement(initial,paths[w1])
            scores1 = score_step(scores,w1)
            for w2 in words:
                w3 = byedge.get(e0 ^ edges[w1] ^ edges[w2])
                if w3 is None or not is_path(state1,paths[w2]):
                    continue
                d2 = score_difference(scores1,w2)
                if not no_integer_zero(d2,10):
                    continue
                scores2 = score_step(scores1,w2)
                d3 = score_difference(scores2,w3)
                if not no_integer_zero(d3,10):
                    continue
                state2 = complement(state1,paths[w2])
                if not is_path(state2,paths[w3]):
                    continue
                need(complement(state2,paths[w3]) == complement(initial,standard),
                     'exact seven-vertex target mismatch')
                need(edges[w1] ^ edges[w2] ^ edges[w3] == e0,
                     'boundary edge-XOR mismatch')
                found = dict(B=b, words=(w1,w2,w3),
                             endpoint_differences=(d2,d3))
                break
            if found is not None:
                break
        need(found is not None, 'no uniform certificate for a strong presentation')
        certificates.append(found)
    need(len(certificates) == 50, 'strong five-vertex presentation count changed')
    return certificates


def literal_controls(certificates):
    records = []
    for middle_n in range(6):
        count = steps = 0
        for z in blue_tilings(middle_n):
            for cert in certificates:
                b,words = cert['B'],cert['words']
                original = join([transpose(b),z,b])
                standard = tuple(range(len(original)))
                need(blue(original,standard), 'original presentation is not blue')
                state = original
                for word in words:
                    p = expand(word,len(original)-5)
                    need(is_path(state,p), 'literal lifted path failed')
                    need(not blue(state,p), 'literal lifted move is blue')
                    state = complement(state,p)
                    steps += 1
                need(state == complement(original,standard), 'literal target failed')
                count += 1
        records.append(dict(middle_order=middle_n,certificates=count,black_steps=steps))
    return records


def main():
    certs = compile_certificates()
    differences = Counter(tuple(d) for c in certs for d in c['endpoint_differences'])
    report = dict(universe='all 64 five-vertex tournaments with path 01234',
                  strong_presentations=50, uniform_three_move_certificates=certs,
                  endpoint_difference_histogram={str(k):v for k,v in sorted(differences.items())},
                  literal_controls=literal_controls(certs))
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()
                      if k != 'uniform_three_move_certificates'},indent=2))


if __name__ == '__main__':
    main()
