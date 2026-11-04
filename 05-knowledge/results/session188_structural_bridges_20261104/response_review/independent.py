#!/usr/bin/env python3
"""Direct Hamilton-word audit, independent of the cover normal form."""
from itertools import combinations, permutations
from collections import defaultdict, Counter
from pathlib import Path
import json


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


PAIRS = list(combinations(range(5), 2))
EDGE_BIT = {p:i for i,p in enumerate(PAIRS)}
PAIR_MASKS = [(1 << i) | (1 << j) for i,j in PAIRS]
BASIS = sorted(a | (c << 5) for a in PAIR_MASKS for c in PAIR_MASKS)
INDEX = {s:i for i,s in enumerate(BASIS)}


def arc(code, u, v):
    if u < v:
        return (code >> EDGE_BIT[u,v]) & 1
    return 1-((code >> EDGE_BIT[v,u]) & 1)


def rows(code):
    return [sum(1 << v for v in range(5) if v != u and arc(code,u,v)) for u in range(5)]


def is_strong(code):
    reach = [[u == v or bool(arc(code,u,v)) for v in range(5)] for u in range(5)]
    for middle in range(5):
        for u in range(5):
            for v in range(5):
                reach[u][v] |= reach[u][middle] and reach[middle][v]
    return all(all(row) for row in reach)


def quartics_from_words():
    vectors = [[0]*100 for _ in range(1024)]
    skeletons = 0
    for word in permutations(range(7)):
        ix, iy = word.index(5), word.index(6)
        if ix in (0,6) or iy in (0,6) or abs(ix-iy) == 1:
            continue
        a = (1 << word[ix-1]) | (1 << word[ix+1])
        c = (1 << word[iy-1]) | (1 << word[iy+1])
        mon = a | (c << 5)
        req_mask = req_value = 0
        for u,v in zip(word,word[1:]):
            if u < 5 and v < 5:
                bit = 1 << EDGE_BIT[min(u,v),max(u,v)]
                req_mask |= bit
                if u < v:
                    req_value |= bit
        check(req_mask.bit_count() == 2, 'literal quartic word must test two old arcs')
        for code in range(1024):
            if code & req_mask == req_value:
                vectors[code][INDEX[mon]] += 1
        skeletons += 1
    return vectors, skeletons


def full_query(code, amask, cmask):
    def edge(u,v):
        if u < 5 and v < 5:
            return arc(code,u,v)
        if u >= 5 and v >= 5:
            return u == 5 and v == 6
        old, port = (u,v) if u < 5 else (v,u)
        bit = ((amask if port == 5 else cmask) >> old) & 1
        return bit if u < 5 else 1-bit
    return sum(all(edge(u,v) for u,v in zip(w,w[1:])) for w in permutations(range(7)))


def full_polynomial(code):
    result = defaultdict(int)
    for word in permutations(range(7)):
        factors = []
        for u,v in zip(word,word[1:]):
            if u < 5 and v < 5:
                if not arc(code,u,v):
                    break
            elif u >= 5 and v >= 5:
                if (u,v) != (5,6):
                    break
            else:
                old, port = (u,v) if u < 5 else (v,u)
                mask = 1 << (old + 5*(port-5))
                factors.append({mask:1} if u < 5 else {0:1, mask:-1})
        else:
            terms = {0:1}
            for factor in factors:
                expanded = defaultdict(int)
                for a,x in terms.items():
                    for b,y in factor.items():
                        expanded[a | b] += x*y
                terms = expanded
            for s,v in terms.items():
                result[s] += v
    return {s:v for s,v in result.items() if v}


def main():
    vectors,skeletons = quartics_from_words()
    groups = defaultdict(list)
    for code,vector in enumerate(vectors):
        groups[tuple(vector)].append(code)
    histogram = dict(Counter(map(len,groups.values())))
    check(len(groups) == 502 and histogram == {2:492,4:10}, 'exact class census')
    exceptions = []
    strong_count = 0
    for vector,codes in groups.items():
        if len(codes) == 2:
            check(codes[0] ^ codes[1] == 1023, 'nonconverse ordinary pair')
        else:
            shapes = []
            triangles = set()
            for code in codes:
                adj = rows(code)
                source = [i for i,r in enumerate(adj) if r.bit_count() == 4]
                sink = [i for i,r in enumerate(adj) if r == 0]
                check(len(source) == len(sink) == 1, 'exception exterior shape')
                tri = tuple(i for i in range(5) if i not in source+sink)
                check(all(sum(arc(code,u,v) for v in tri if v != u) == 1 for u in tri), 'middle not cyclic')
                triangles.add(tri)
                shapes.append((source[0],sink[0],arc(code,tri[0],tri[1])))
            check(len(triangles) == 1, 'named triangle not fixed')
            tri = next(iter(triangles))
            outside = [i for i in range(5) if i not in tri]
            check(set(shapes) == {(outside[0],outside[1],0),(outside[0],outside[1],1),
                                  (outside[1],outside[0],0),(outside[1],outside[0],1)},
                  'four choices not all present')
            exceptions.append(dict(triangle=tri,codes=codes,shapes=shapes))
        for code in codes:
            if is_strong(code):
                check(len(codes) == 2, 'strong exception')
                strong_count += 1
    check(strong_count == 544, 'strong count')
    check({tuple(x['triangle']) for x in exceptions} == set(combinations(range(5),3)), 'not all ten triangle sets')
    # Compare every coordinate of every code against published exact records.
    published = json.loads((Path(__file__).resolve().parent.parent/'response'/'classification.json').read_text())
    check(published['quartic_basis'] == [[s&31,s>>5] for s in BASIS], 'basis ordering differs')
    checked_codes = set()
    for rec in published['records']:
        for code,adj in zip(rec['codes'],rec['rows']):
            check(vectors[code] == rec['quartic_vector'], ('full exact vector mismatch',code))
            check(rows(code) == adj, ('adjacency mismatch',code))
            checked_codes.add(code)
    check(checked_codes == set(range(1024)), 'published record coverage')
    T,U = 5,131
    p,q = full_polynomial(T),full_polynomial(U)
    difference = {s:p.get(s,0)-q.get(s,0) for s in set(p)|set(q) if p.get(s,0)!=q.get(s,0)}
    counts = dict(Counter(s.bit_count() for s in difference))
    check(counts == {2:6,3:36}, 'difference degree counts')
    query = full_query(T,1,4),full_query(U,1,4)
    check(query == (31,29), 'literal query separation')
    old_counts = [sum(all(arc(code,u,v) for u,v in zip(w,w[1:])) for w in permutations(range(5))) for code in (T,U)]
    check(old_counts == [3,3], 'old Hamilton counts')
    check(rows(T) == [10,0,3,6,15] and rows(U) == [6,0,10,3,15], 'named example rows')
    check(vectors[T] == vectors[U], 'example quartics differ')
    return dict(status='PASS', literal_quartic_word_skeletons=skeletons,
                full_named_coordinates_compared=102400, quartic_classes=502,
                class_histogram=histogram,strong_tournaments=strong_count,
                exceptional_classes=sorted(exceptions,key=lambda x:x['triangle']),
                example=dict(T=rows(T),U=rows(U),old_hamilton_counts=old_counts,
                             literal_query_counts=query,difference_degree_counts=counts,
                             complete_difference=[[s&31,s>>5,v] for s,v in sorted(difference.items())]))


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
