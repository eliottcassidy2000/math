#!/usr/bin/env python3
"""Exact, source-marked Gray / component-chart / tournament controls.

Standard library only. Every check remains enabled under python -O.
No unit-distance realization or universal number-theoretic claim is inferred.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0
PERM = (0, 1, 4, 7, 2, 3, 6, 5, 8, 9)
GRAY = ((0, 0), (0, 1), (1, 1), (1, 0))


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def span(rows):
    words = {0}
    for row in rows:
        words |= {word ^ row for word in words}
    return tuple(sorted(words))


def gray_word(a, b, t):
    values = tuple((a + j*b) % 4 for j in range(4)) + (2*t,)
    bits = sum((GRAY[v] for v in values), ())
    return sum(bit << PERM[i] for i, bit in enumerate(bits))


def encode(a, b, u, v, w):
    return (170*a) ^ (204*b) ^ (15*u) ^ (240*v) ^ (768*w)


def decode(word):
    bit = lambda i: (word >> i) & 1
    return bit(0)^bit(1), bit(0)^bit(2), bit(0), bit(4), bit(8)


def maximum_pack(code):
    """DP on attainable unions, all of which are themselves codewords.

    A disjoint union is XOR, so at most |C| states are needed, not 2^length.
    Integer order is topological because adjoining a nonempty disjoint mask
    strictly increases the integer. Runtime O(|C|^2) apart from bit costs.
    """
    best = {0: ()}
    for used in code:
        if used not in best:
            continue
        for word in code[1:]:
            if used & word == 0:
                merged = used | word
                candidate = best[used] + (word,)
                if len(candidate) > len(best.get(merged, ())):
                    best[merged] = candidate
    return max(best.values(), key=len)


def cosets(code, subcode):
    remaining = set(code)
    result = []
    while remaining:
        root = min(remaining)
        cell = tuple(sorted(root ^ word for word in subcode))
        need(set(cell) <= remaining, "disjoint actual coset")
        result.append((root, cell))
        remaining.difference_update(cell)
    return result


def all_codes(n):
    """Every distinct binary linear code on n named coordinates, incl zero."""
    frontier = { (0,) }
    seen = set(frontier)
    while frontier:
        new = set()
        for code in frontier:
            for v in range(1, 1 << n):
                if v not in code:
                    candidate = tuple(sorted(set(code) | {v ^ w for w in code}))
                    if candidate not in seen:
                        new.add(candidate)
        seen.update(new)
        frontier = new
    return sorted(seen, key=lambda c: (len(c), c))


def literal_min_cover(code):
    """Independent small-case set cover over ALL disjoint-generator charts."""
    nonzero = code[1:]
    spaces = set()

    def packs(start, used, rows):
        spaces.add(span(rows))
        for i in range(start, len(nonzero)):
            if used & nonzero[i] == 0:
                packs(i+1, used | nonzero[i], rows + (nonzero[i],))

    packs(0, 0, ())
    positions = {word: i for i, word in enumerate(code)}
    cells = set()
    for subspace in spaces:
        for root in code:
            cells.add(sum(1 << positions[root ^ word] for word in subspace))
    # Non-maximal charts cannot improve a covering by positive-cost charts.
    cells = [s for s in cells if not any(s != t and s & t == s for t in cells)]
    by_bit = [[s for s in cells if s & (1 << i)] for i in range(len(code))]

    @lru_cache(None)
    def cover(left):
        if not left:
            return 0
        i = (left & -left).bit_length() - 1
        return 1 + min(cover(left & ~s) for s in by_bit[i])

    return cover((1 << len(code))-1)


def rm1(m, pairs=0):
    size = 1 << m
    rows = [(1 << size)-1]
    rows += [sum(((x >> j) & 1) << x for x in range(size)) for j in range(m)]
    rows += [3 << (size + 2*j) for j in range(pairs)]
    return span(rows)


def boolean_atoms(code, n):
    """Coordinate signature classes are atoms of the generated Boolean ring."""
    groups = {}
    for i in range(n):
        signature = tuple((word >> i) & 1 for word in code)
        if any(signature):
            groups[signature] = groups.get(signature, 0) | (1 << i)
    return tuple(sorted(groups.values()))


def cycle_code(n):
    edges = tuple(combinations(range(n), 2))
    code = []
    for mask in range(1 << len(edges)):
        parity = 0
        for i, (a, b) in enumerate(edges):
            if mask & (1 << i):
                parity ^= (1 << a) | (1 << b)
        if not parity:
            code.append(mask)
    return edges, tuple(code)


def graph_chart(n, root, generators):
    """A binary inequality graph, with pins; no geometric embedding asserted."""
    edges = []
    vertex_count = n
    used = 0
    for row in generators:
        need(row and not used & row, "nonempty disjoint generator")
        used |= row
        support = [i for i in range(n) if row & (1 << i)]
        hub = support[0]
        for v in support[1:]:
            if ((root >> hub) ^ (root >> v)) & 1:
                edges.append(tuple(sorted((hub, v))))
            else:
                hidden = vertex_count
                vertex_count += 1
                edges.extend([(hub, hidden), (v, hidden)])
    pins = tuple((i, (root >> i) & 1) for i in range(n) if not used & (1 << i))
    return vertex_count, tuple(edges), pins


def brute_graph_support(n, graph):
    """Enumerate every vertex assignment, independently of parity elimination."""
    total, edges, pins = graph
    support = Counter()
    for word in range(1 << total):
        if any(((word >> i) & 1) != bit for i, bit in pins):
            continue
        if any(((word >> a) ^ (word >> b)) & 1 == 0 for a, b in edges):
            continue
        support[word & ((1 << n)-1)] += 1
    return support


def tournament(a, b):
    edges = ((0, 1), (1, 2), (2, 3), (0, 3),
             (2, 0) if a else (0, 2), (3, 1) if b else (1, 3))
    degree = tuple(sorted(sum(x == v for x, y in edges) for v in range(4)))
    paths = sum(all((p[i], p[i+1]) in edges for i in range(3))
                for p in permutations(range(4)))
    return degree, paths


def majority(x, y, z):
    return (x & y) | (x & z) | (y & z)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true",
                        help="compare with the committed JSON, without writing")
    args = parser.parse_args()
    code = span((15, 240, 768, 170, 204))
    base = span((15, 240, 768))
    need(len(code) == 32 and len(base) == 8, "inherited dimensions")
    need(set(code) == {gray_word(a, b, t) for a in range(4)
                      for b in range(4) for t in range(2)}, "literal Z4 image")
    receipts = []
    for a, b, t in product(range(4), range(4), range(2)):
        a0, a1, b0, b1 = a % 2, a // 2, b % 2, b // 2
        u, v = a1, a1 ^ b1 ^ (a0*b0)
        word = gray_word(a, b, t)
        need(word == encode(a0, b0, u, v, t), "source-marked Gray formula")
        need(decode(word) == (a0, b0, u, v, t), "bit inverse")
        need((a0+2*u, b0+2*(v ^ u ^ (a0*b0)), t) == (a, b, t),
             "quaternary inverse with carry")
        need(decode(word ^ 170) == (a0^1, b0, u, v, t), "named ac flip")
        need(decode(word ^ 204) == (a0, b0^1, u, v, t), "named bd flip")
        receipts.append([a, b, t, word])
    for x, y in product(range(4), repeat=2):
        expected = tuple(GRAY[x][i] ^ GRAY[y][i] ^ ((x % 2)*(y % 2))
                         for i in range(2))
        need(GRAY[(x+y) % 4] == expected, "Gray addition carry")
        lee = min((x-y) % 4, (y-x) % 4)
        need(sum(u != v for u, v in zip(GRAY[x], GRAY[y])) == lee, "Lee isometry")
    for a, b, t in product(range(4), range(4), range(2)):
        for aa, bb, tt in product(range(4), range(4), range(2)):
            first = [(a+j*b) % 4 for j in range(4)] + [2*t]
            second = [(aa+j*bb) % 4 for j in range(4)] + [2*tt]
            correction = sum(((first[j] % 2)*(second[j] % 2)) << PERM[2*j+eta]
                             for j in range(5) for eta in range(2))
            need(gray_word((a+aa) % 4, (b+bb) % 4, t^tt) ==
                 gray_word(a,b,t) ^ gray_word(aa,bb,tt) ^ correction,
                 "full Z4 code addition with carry")
    need(gray_word(1,1,0) == 150 and encode(1,1,0,0,0) == 102,
         "missing carry source hostile")
    # Omitting the carry changes the decoder, but NOT aggregate charts.
    for a0, b0 in product(range(2), repeat=2):
        actual = {gray_word(a0+2*u, b0+2*v, w)
                  for u,v,w in product(range(2), repeat=3)}
        no_carry = {encode(a0,b0,u,u^v,w) for u,v,w in product(range(2), repeat=3)}
        need(actual == no_carry, "aggregate counts conceal missing source carry")
    need(tuple(w for w in code if w.bit_count() == 2) == (768,), "unique marked pair")
    rotated = ((768 << 1) & 1023) | (768 >> 9)
    need(rotated == 513 and rotated not in code, "mod10 shift hostile")

    # Actual component graphs: all assignments, not just proposed certificates.
    charts = []
    for root in (0, 170, 204, 102):
        graph = graph_chart(10, root, (15,240,768))
        support = brute_graph_support(10, graph)
        expected = {root ^ w for w in base}
        need(set(support) == expected, "exact graph support")
        need(set(support.values()) == {1}, "unique hidden extension")
        charts.append({"root": root, "vertices": graph[0], "edges": len(graph[1]),
                       "support": sorted(support)})
    odd_graph = graph_chart(10, 0, (15,240,768))
    odd_graph = (odd_graph[0], odd_graph[1]+((0,1),), odd_graph[2])
    need(not brute_graph_support(10, odd_graph), "odd cycle guard")
    pinned = graph_chart(3, 4, (3,))
    need(set(brute_graph_support(3,pinned)) == {4,7}, "fixed-coordinate pins retained")
    need(set(brute_graph_support(3,(pinned[0], pinned[1], ()))) == {0,3,4,7},
         "dropping pin doubles support")

    # Universal finite test universe, independently checked by literal set cover.
    small_counts, cover_counts, histogram = {}, {}, Counter()
    for n in range(7):
        universe = all_codes(n)
        small_counts[str(n)] = len(universe)
        for c in universe:
            packing = maximum_pack(c)
            nu, k = len(packing), len(c).bit_length()-1
            for x,y in combinations(packing,2):
                need(not x & y, "packing disjoint")
            need(all(x in c and x for x in packing), "packing in code")
            cells = cosets(c, span(packing))
            predicted = 1 << (k-nu)
            need(len(cells) == predicted, "attained chart bound")
            need(sum(len(cell) for _,cell in cells) == len(c), "full cover")
            closed = all((x & y) in c for x in c for y in c)
            need(closed == (predicted == 1), "one chart iff Boolean ring")
            atoms = boolean_atoms(c,n)
            generated = span(atoms)
            need(set(c) <= set(generated), "Boolean completion contains code")
            if closed:
                need(generated == c, "Boolean atoms reconstruct closed code")
            if (1 << n)-1 in c:
                need(sum(packing) == (1 << n)-1, "all-ones maximal pack covers")
            if n <= 4:
                need(literal_min_cover(c) == predicted, "independent exact set cover")
                cover_counts[str(n)] = cover_counts.get(str(n),0)+1
            histogram[f"{n},{k},{nu}"] += 1
    need(small_counts == {"0":1,"1":2,"2":5,"3":16,"4":67,"5":374,"6":2825},
         "all-code universe sizes")

    rm_rows = []
    for m in range(1,7):
        for q in range(4):
            c = rm1(m,q)
            packing = maximum_pack(c)
            need(len(c) == 1 << (m+1+q), "RM dimension")
            need(len(packing) == q+2, "RM packing")
            rm_rows.append({"m":m,"pairs":q,"length":(1<<m)+2*q,
                            "dimension":m+1+q,"nu":q+2,"charts":1<<(m-1)})
    # Both recursive operations, at actual word-set level.
    for m in range(1,6):
        c, size = rm1(m), 1 << m
        doubled = tuple(sorted({u | ((u ^ (v*((1<<size)-1))) << size)
                                for u in c for v in range(2)}))
        need(doubled == rm1(m+1), "Plotkin one-new-selector doubling")
        diagonal = tuple(sorted(w | (w << size) for w in c))
        need(len(maximum_pack(diagonal)) == len(maximum_pack(c)), "duplication preserves nu")
    for a in all_codes(3):
        for b in all_codes(3):
            joined = tuple(sorted(x | (y<<3) for x in a for y in b))
            need(len(maximum_pack(joined)) == len(maximum_pack(a))+len(maximum_pack(b)),
                 "direct sum packing additivity")

    graph_rows = []
    for n in (3,4,5):
        edges,c = cycle_code(n)
        packing = maximum_pack(c)
        k = len(edges)-n+1
        graph_rows.append({"graph":f"K{n}","edges":len(edges),"dimension":k,
                           "nu":len(packing),"charts":1<<(k-len(packing)),
                           "packed_edge_sets":[[list(e) for i,e in enumerate(edges)
                                              if row & (1<<i)] for row in packing]})
        need(len(c) == 1<<k, "cycle dimension")
        need(len(packing) == {3:1,4:1,5:3}[n], "graphic packing")
    need(all(w.bit_count() != 2 for w in cycle_code(5)[1]), "K5 weight-two hostile")

    # Exact old blue involution, including the displaced affine coset.
    swap = (1,0,3,2,5,4,7,6,8,9)
    sigma = lambda w: sum(((w>>i)&1) << swap[i] for i in range(10))
    fixed = tuple(w for w in code if sigma(w) == w)
    need(fixed == span((15,240,768,204)), "blue intersection")
    for w in code:
        need(w ^ sigma(w) == 255*decode(w)[0], "retained blue defect bit")
    need(set(code) == set(fixed) | {170^w for w in fixed}, "two-stage fold")
    collapsed = {sum(((w>>i)&1)<<j for j,i in enumerate((0,2,4,6))) for w in fixed}
    need(collapsed == {w for w in range(16) if w.bit_count()%2 == 0}, "even-parity fold")
    need(len(maximum_pack(fixed)) == 3 and len(fixed) == 16, "two charts in blue part")
    # Blue B=I2^4+F2^2 forces its unique self-dual subcode to be I2^5.
    blue = span((3,12,48,192,256,512))
    blue_dual = span((3,12,48,192))
    self_dual_extensions = set()
    for w in blue:
        c = span(blue_dual+(w,))
        if len(c) == 32 and all((x&y).bit_count()%2 == 0 for x in c for y in c):
            self_dual_extensions.add(c)
    need(self_dual_extensions == {span((3,12,48,192,768))}, "unique blue self-dual")
    need(set(code) != set(next(iter(self_dual_extensions))), "Hamming escapes blue fixed locus")
    atoms = boolean_atoms(code,10)
    boolean_completion = span(atoms)
    need(len(atoms) == 9 and len(boolean_completion) == 512, "Boolean completion dimension")
    need((170 & 204) == 136 and 136 not in code, "quadratic interaction hostile")
    need(set(boolean_completion) == {w for w in range(1024)
                                     if ((w >> 8) ^ (w >> 9)) & 1 == 0},
         "completion preserves only distinguished repetition pair")
    # Incoming majority obstruction: every binary constraint is preserved by
    # coordinatewise majority, hence every 2-CNF network and its projection is.
    for relation_mask in range(16):
        relation = tuple(w for w in range(4) if relation_mask & (1 << w))
        for x,y,z in product(relation,repeat=3):
            need(majority(x,y,z) in relation, "every binary relation majority closed")
    need(majority(0,170,204) == 136 and 136 not in code, "arbitrary two-list obstruction")
    # A fixed three-colour list permits an NAE primitive; forgetting that list
    # admits every binary port tuple in ordinary six-colouring.
    nae_rows = []
    for ports in product(range(2),repeat=3):
        listed = sum(all(x != y for x,y in zip(inner,ports))
                     for inner in permutations(range(3)))
        unrestricted = sum(all(x != y for x,y in zip(inner,ports))
                           for inner in permutations(range(6),3))
        need(bool(listed) == (len(set(ports)) > 1), "three-list triangle NAE")
        need(unrestricted > 0, "unrestricted six colours repair NAE")
        nae_rows.append({"ports":ports,"listed_extensions":listed,
                         "six_colour_extensions":unrestricted})

    types = {(0,0):"T",(1,0):"+",(0,1):"-",(1,1):"S"}
    tournament_rows = []
    for a,b in product(range(2),repeat=2):
        degrees,h = tournament(a,b)
        need(h == {(0,0):1,(1,0):3,(0,1):3,(1,1):5}[a,b], "marked tournament H")
        tournament_rows.append({"a":a,"b":b,"type":types[a,b],
                                "degrees":degrees,"H":h})
    need(tournament(1,0)[1] == tournament(0,1)[1], "H loses two selector states")
    need(1 not in code, "ten-edge tournament relabelling 0-to-1 hostile")
    payload = {"status":"FINITE-EXACT; abstract charts, not unit-distance gadgets",
               "permutation":PERM,"code":code,"source_receipts":receipts,
               "component_charts":charts,"all_binary_codes":small_counts,
               "independent_set_covers":cover_counts,
               "packing_histogram":dict(sorted(histogram.items())),
               "rm_family":rm_rows,"graph_codes":graph_rows,
               "tournament_controller":tournament_rows,
               "carry_hostile":{"source":[1,1,0],"correct":150,"forgotten_carry":102},
               "mod10_hostile":{"word":768,"rotation":rotated},
               "boolean_completion":{"atoms":atoms,"dimension":9,"words":512,
                                     "added_words":480,"missing_product":136},
               "majority":{"inputs":[0,170,204],"output":136,
                           "three_list_nae":nae_rows},
               "blue_fixed_intersection":fixed,"checks":CHECKS}
    encoded = json.dumps(payload,sort_keys=True,separators=(",",":"))
    payload["content_sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    target = ROOT/"05-knowledge/results/gray_selector_recursion_20261018e.json"
    if args.verify:
        need(json.loads(target.read_text(encoding="utf-8")) ==
             json.loads(json.dumps(payload)), "frozen artifact replay")
    else:
        target.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"checks":CHECKS,"sha256":payload["content_sha256"],
                      "all_codes":sum(small_counts.values()),
                      "literal_cover_checks":sum(cover_counts.values()),
                      "path":str(target),"mode":"verify" if args.verify else "write"}))


if __name__ == "__main__":
    main()
