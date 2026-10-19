#!/usr/bin/env python3
"""All optimal H8+I2 component atlases, with the retained primal/dual map.

Standard library; all checks active under -O. Incoming SL2 port permutations
are read as data and independently checked on each affected codeword.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def span(rows):
    out = {0}
    for r in rows:
        out |= {r ^ w for w in out}
    return tuple(sorted(out))


def matchings(items):
    if not items:
        yield ()
        return
    for b in items[1:]:
        for tail in matchings(tuple(x for x in items[1:] if x != b)):
            yield ((items[0],b),) + tail


def linear_word(slope):
    return (170*(slope&1)) ^ (204*((slope>>1)&1)) ^ (240*((slope>>2)&1))


def slope(word):
    bit = lambda i: (word >> i) & 1
    return (bit(0)^bit(1)) + 2*(bit(0)^bit(2)) + 4*(bit(0)^bit(4))


def relabel(word, perm):
    return sum(((word >> i) & 1) << perm[i] for i in range(len(perm)))


def key_partition(cells):
    return tuple(sorted(tuple(sorted(cell)) for cell in cells))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify",action="store_true")
    args = ap.parse_args()
    H = span((170,204,240,255))
    C = span((170,204,240,255,768))
    E = span((255,768))
    fibre = {a:tuple(sorted(linear_word(a)^e for e in E)) for a in range(8)}
    need(set().union(*map(set,fibre.values())) == set(C), "slope quotient exhausted")
    cells = {pair:tuple(sorted(fibre[pair[0]]+fibre[pair[1]]))
             for pair in combinations(range(8),2)}
    need(len(set(cells.values())) == 28, "28 distinct pair preimages")

    # Independent maximum-chart enumeration from actual disjoint codewords.
    directions = set()
    for x,y,z in combinations(C[1:],3):
        if not (x&y or x&z or y&z):
            need(sorted((x.bit_count(),y.bit_count(),z.bit_count())) == [2,4,4],
                 "only possible maximum packing weights")
            directions.add(span((x,y,z)))
    need(len(directions) == 7, "seven maximum direction spaces")
    maximum_charts = {tuple(sorted(r^x for x in D)) for D in directions for r in C}
    need(maximum_charts == set(cells.values()), "all and only 28 maximum charts")
    need(all(set(E) <= set(D) for D in directions), "common phase subspace")

    match_list = list(matchings(tuple(range(8))))
    need(len(match_list) == 105, "all perfect matchings of eight named points")
    records, nodes, edges = [], set(), set()
    hist = Counter()
    partitions = set()
    for matching in match_list:
        atlas = key_partition(cells[p] for p in matching)
        partitions.add(atlas)
        need(set().union(*map(set,atlas)) == set(C), "full boundary support")
        need(sum(map(len,atlas)) == len(C), "four disjoint maximum charts")
        decoded = tuple(sorted(tuple(sorted({slope(w) for w in chart})) for chart in atlas))
        need(decoded == matching, "matching recovered without forgetting source")
        differences = tuple(sorted(a^b for a,b in matching))
        distinct = len(set(differences))
        need(distinct in (1,2,4), "one two or four directions")
        hist[distinct] += 1

        # The SAME marked matching, now on physical coordinate addresses.
        pair_rows = tuple((1<<a)|(1<<b) for a,b in matching)
        even_pairs = span(pair_rows[0]^x for x in pair_rows[1:])
        need(len(even_pairs) == 8, "even-pair code dimension three")
        orthogonal = [v for v in range(256)
                      if all((v&e).bit_count()%2 == 0 for e in even_pairs)]
        extensions = {span(even_pairs+(v,)) for v in orthogonal if v not in even_pairs}
        need(len(extensions) == 3, "three self-dual extensions")
        hamming = tuple(sorted(c for c in extensions
                               if min(w.bit_count() for w in c if w) == 4))
        need(len(hamming) == 2, "two Hamming extensions per matching")
        need((H in hamming) == (distinct == 1), "parallel iff edge incident with H")
        need(hamming not in edges, "matching edge injective")
        edges.add(hamming); nodes.update(hamming)
        records.append({"matching":matching,"direction_multiset":differences,
                        "atlas":atlas,"hamming_edge":hamming})
    need(len(partitions) == 105, "all optimal partitions distinct")
    need(hist == {1:7,2:42,4:56}, "sharp 7+42+56 decomposition")
    need(len(nodes) == 30 and len(edges) == 105, "old Hamming incidence carrier")
    node_list = sorted(nodes)
    degrees = Counter(v for edge in edges for v in edge)
    need(set(degrees.values()) == {7}, "all thirty Hamming vertices degree seven")

    # The primal/dual identification has a sharp affine-translation boundary.
    hostile = ((0,1),(2,3),(4,6),(5,7))
    moved = tuple(sorted(tuple(sorted((a^4,b^4))) for a,b in hostile))
    need(moved == ((0,2),(1,3),(4,5),(6,7)) and moved != hostile,
         "physical translation moves primal matching")
    for b in range(8):
        perm = tuple(x^b for x in range(8))+(8,9)
        for cell in maximum_charts:
            need(tuple(sorted(relabel(w,perm) for w in cell)) == cell,
                 "physical coordinate translations fix each maximum chart")

    # The seven parallel atlases: SL2(3) acts through A4 with orbits 4+3.
    source = ROOT/"05-knowledge/results/marked_geometries_selectors_20261018o.json"
    incoming = json.loads(source.read_text(encoding="utf-8"))
    groups = [g for g in incoming["witt"]["GL2_3"] if g["determinant"] == 1]
    need(len(groups) == 24, "all determinant-one port actions")
    actions, action_records = set(), []
    for g in groups:
        perm = tuple(g["port_permutation"])+(8,9)
        need({relabel(w,perm) for w in C} == set(C), "actual code permutation")
        action = tuple(slope(relabel(linear_word(a),perm)) for a in range(8))
        need(sorted(action) == list(range(8)) and action[0] == 0, "quotient automorphism")
        for a in range(8):
            for b in range(8):
                need(action[a^b] == action[a]^action[b], "linear quotient action")
        actions.add(action)
        action_records.append({"matrix":g["matrix"],"ports":perm,"slopes":action})
    need(len(actions) == 12, "central order-two kernel")
    outer = tuple(sorted({g[4] for g in actions}))
    inner = tuple(x for x in range(1,8) if x not in outer)
    plane = (0,)+inner
    need(len(outer) == 4 and len(inner) == 3, "four plus three orbit sizes")
    need(all(x^y in plane for x in plane for y in plane), "inner directions form plane")
    need({outer[0]^w for w in plane} == set(outer), "outer directions affine coset")
    k4 = {w:tuple((a,b) for a,b in combinations(outer,2) if a^b == w) for w in inner}
    for w,pairs in k4.items():
        need(len(pairs) == 2 and {x for p in pairs for x in p} == set(outer), "K4 perfect matching")
    for g in actions:
        for w in inner:
            transported = tuple(sorted(tuple(sorted((g[a],g[b]))) for a,b in k4[w]))
            need(transported == k4[g[w]], "equivariant K4 matching dictionary")
    need(len({tuple(g[a] for a in outer) for g in actions}) == 12, "faithful A4 on four atlases")
    minus = next(g for g in action_records if g["matrix"] == [2,0,0,2])
    need(minus["slopes"] == tuple(range(8)), "minus-I disappears on quotient")
    need(relabel(15,minus["ports"]) == 240, "minus-I still changes component phase")

    # Incoming four H8 atlases lift to exactly our four parallel atlases.
    def parallel_atlas(d):
        return key_partition(cells[(a,a^d)] for a in range(8) if a < (a^d))
    incoming_atlases = set()
    for atlas in incoming["witt"]["four_atlases"]:
        incoming_atlases.add(key_partition(
            tuple(sorted([w for w in cell]+[w^768 for w in cell])) for cell in atlas))
    need(incoming_atlases == {parallel_atlas(d) for d in outer},
         "same actual four partitions, not just cardinality")

    # Keep the phase of the primal/dual operators. Their commutator is an
    # actual complement, and products can have order four.
    def act(word,a,b,c):
        correction = ((slope(word)&b).bit_count() % 2) ^ c
        return word ^ linear_word(a) ^ (255*correction)

    params = [(a,b,c) for a in range(8) for b in range(8) for c in range(2)]
    actions128 = {g:tuple(act(w,*g) for w in C) for g in params}
    need(len(set(actions128.values())) == 128, "faithful binary phase group")
    order_hist = Counter()
    for a,b,c in params:
        order_hist[1 if (a,b,c)==(0,0,0) else
                   (4 if (a&b).bit_count()%2 else 2)] += 1
        for w in C:
            need(act(act(w,a,b,c),a,b,c) == w ^ (255*((a&b).bit_count()%2)),
                 "square is exact mod-four phase")
    for a,b,c in params:
        for aa,bb,cc in params:
            combined = (a^aa,b^bb,c^cc^((b&aa).bit_count()%2))
            for w in C:
                need(act(act(w,aa,bb,cc),a,b,c) == act(w,*combined),
                     "ordered central cocycle composition")
    need(order_hist == {1:1,2:71,4:56}, "phase-group orders")
    need(act(act(0,1,0,0),0,1,0) == 170^255 and
         act(act(0,0,1,0),1,0,0) == 170,
         "smallest noncommuting source hostile")

    payload = {"status":"FINITE-EXACT; optimal abstract component atlases",
               "common_phase_subspace":E,"slope_fibres":fibre,
               "maximum_charts":[{"pair":p,"words":w} for p,w in cells.items()],
               "matching_atlas_records":records,
               "direction_counts":dict(sorted(hist.items())),
               "hamming_vertices":node_list,
               "hamming_edges":[[node_list.index(c) for c in edge] for edge in sorted(edges)],
               "incoming_source":source.name,
               "incoming_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
               "SL2_actions":action_records,"outer_orbit":outer,"inner_orbit":inner,
               "K4_matchings":k4,
               "translation_hostile":{"original":hostile,"moved":moved},
               "central_phase_hostile":{"source":15,"image":240},
               "binary_phase_group":{"order":128,"element_orders":dict(order_hist),
                                     "law":"(a,b,c)(a',b',c')=(a+a',b+b',c+c'+b.dot(a'))",
                                     "noncommuting_source":0,"TS_image":85,"ST_image":170},
               "checks":CHECKS}
    text = json.dumps(payload,sort_keys=True,separators=(",",":"))
    payload["content_sha256"] = hashlib.sha256(text.encode()).hexdigest()
    target = ROOT/"05-knowledge/results/hamming_matching_atlases_20261018e.json"
    if args.verify:
        need(json.loads(target.read_text(encoding="utf-8")) ==
             json.loads(json.dumps(payload)), "artifact replay")
    else:
        target.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"checks":CHECKS,"sha256":payload["content_sha256"],
                      "charts":28,"atlases":105,"parallel":7,"directions":dict(hist)}))


if __name__ == "__main__":
    main()
