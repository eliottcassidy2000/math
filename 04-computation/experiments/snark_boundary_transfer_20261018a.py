"""Exact marked snark/vertex-six interfaces; no geometric realization inference.

Standard library only. Run python -B [-O] this_file.py for deterministic JSON.
The inherited four-pole atlas is checked by colouring and by matchings here.
"""
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product
import json


def require(ok, message):
    if not ok:
        raise ValueError(message)


def nat(x, minimum=0):
    require(type(x) is int and x >= minimum, "exact integer out of range")
    return x


@dataclass(frozen=True)
class Graph:
    n: int
    edges: tuple


def audit_graph(g):
    require(type(g) is Graph, "Graph required")
    nat(g.n)
    require(type(g.edges) is tuple, "canonical tuple of edges required")
    for e in g.edges:
        require(type(e) is tuple and len(e) == 2, "edge shape")
        a, b = e
        nat(a); nat(b)
        require(a < b < g.n, "simple ordered edge required")
    require(tuple(sorted(set(g.edges))) == g.edges, "canonical edge bank required")
    return g


def graph(n, edges):
    nat(n)
    raw = tuple(edges)
    for e in raw:
        require(type(e) is tuple and len(e) == 2, "edge shape")
        nat(e[0]); nat(e[1])
        require(e[0] != e[1], "loop")
    return audit_graph(Graph(n, tuple(sorted({tuple(sorted(e)) for e in raw}))))


def adjacency(g):
    audit_graph(g)
    result = [set() for _ in range(g.n)]
    for a, b in g.edges:
        result[a].add(b); result[b].add(a)
    return tuple(frozenset(x) for x in result)


def cubic(g):
    require(all(len(x) == 3 for x in adjacency(g)), "cubic graph required")
    return g


def line_graph(g, deleted_vertices=(), active=None):
    audit_graph(g)
    require(type(deleted_vertices) is tuple, "vertex tuple required")
    for v in deleted_vertices:
        nat(v); require(v < g.n, "vertex outside graph")
    ids = tuple(range(len(g.edges))) if active is None else active
    require(type(ids) is tuple, "edge index tuple required")
    for i in ids:
        nat(i); require(i < len(g.edges), "edge index outside graph")
    require(ids == tuple(sorted(set(ids))), "canonical edge indices")
    live = set(range(g.n)) - set(deleted_vertices)
    return graph(len(ids), ((i, j) for i, j in combinations(range(len(ids)), 2)
                           if set(g.edges[ids[i]]) & set(g.edges[ids[j]]) & live))


def audit_colouring(g, k, colours):
    audit_graph(g); nat(k, 1)
    require(type(colours) is tuple and len(colours) == g.n, "colour vector shape")
    for c in colours:
        nat(c); require(c < k, "colour outside palette")
    require(all(colours[a] != colours[b] for a, b in g.edges), "improper colouring")
    return colours


def colourings(g, k):
    """Complete finite DSATUR enumeration, not a capped search."""
    adj = adjacency(g); nat(k, 1)
    values = [-1] * g.n
    def visit(left):
        if not left:
            yield tuple(values)
            return
        v = max(left, key=lambda x: (len({values[y] for y in adj[x] if values[y] >= 0}),
                                      len(adj[x]), -x))
        used = {values[y] for y in adj[v]}
        for c in range(k):
            if c not in used:
                values[v] = c
                yield from visit(left - {v})
        values[v] = -1
    yield from visit(set(range(g.n)))


def greedy_extend(g, k, prescribed):
    adj = adjacency(g); nat(k, 1)
    require(type(prescribed) is tuple, "prescribed tuple required")
    require(k > max(map(len, adj), default=-1), "degree-plus-one guard")
    colours = [-1] * g.n
    for item in prescribed:
        require(type(item) is tuple and len(item) == 2, "marked assignment shape")
        v, c = item; nat(v); nat(c)
        require(v < g.n and c < k and colours[v] == -1, "invalid marked assignment")
        colours[v] = c
    require(all(colours[a] < 0 or colours[b] < 0 or colours[a] != colours[b]
                for a, b in g.edges), "incompatible marked assignment")
    for v in range(g.n):
        if colours[v] < 0:
            used = {colours[w] for w in adj[v]}
            colours[v] = next(c for c in range(k) if c not in used)
    return audit_colouring(g, k, tuple(colours))


@dataclass(frozen=True)
class Pole:
    source: Graph
    removed_edge: int
    active: tuple
    ports: tuple
    constraints: Graph


def make_pole(g, edge_index):
    cubic(g); nat(edge_index)
    require(edge_index < len(g.edges), "missing source edge")
    u, v = g.edges[edge_index]
    active = tuple(i for i in range(len(g.edges)) if i != edge_index)
    ports = tuple(i for w in (u, v) for i, e in enumerate(g.edges)
                  if w in e and i != edge_index)
    require(len(ports) == 4 and len(set(ports)) == 4, "four distinct marked ports")
    return Pole(g, edge_index, active, ports, line_graph(g, (u, v), active))


def audit_pole(p):
    require(type(p) is Pole, "Pole required")
    cubic(p.source); nat(p.removed_edge)
    for values in (p.active, p.ports):
        require(type(values) is tuple, "marked indices tuple required")
        for v in values:
            nat(v)
    audit_graph(p.constraints)
    require(p == make_pole(p.source, p.removed_edge), "source/port reconstruction mismatch")
    return p


def boundary_word(p, labels):
    audit_pole(p); audit_colouring(p.constraints, 3, labels)
    return tuple(labels[p.active.index(i)] for i in p.ports)


def profile(p):
    audit_pole(p)
    counts = Counter(); witnesses = {}
    positions = tuple(p.active.index(i) for i in p.ports)
    for labels in colourings(p.constraints, 3):
        word = tuple(labels[i] for i in positions)
        counts[word] += 1
        witnesses.setdefault(word, labels)
    return counts, witnesses


def matching_count(p, word):
    """Independent boundary count by three disjoint perfect matchings."""
    audit_pole(p)
    require(type(word) is tuple and len(word) == 4, "boundary word shape")
    for c in word:
        nat(c); require(c < 3, "boundary colour")
    g = p.source; u, v = g.edges[p.removed_edge]
    live = frozenset(range(g.n)) - {u, v}
    internal = frozenset(e for e in g.edges if set(e) <= live)
    adj = adjacency(g)
    boundary_live = tuple(next(x for x in g.edges[i] if x in live) for i in p.ports)
    by_colour = [tuple(x for x, c in zip(boundary_live, word) if c == t) for t in range(3)]
    if any(len(set(x)) < len(x) for x in by_colour):
        return 0
    @lru_cache(None)
    def match(remaining):
        if not remaining:
            return (frozenset(),)
        if len(remaining) % 2:
            return ()
        a = min(remaining)
        return tuple(rest | {tuple(sorted((a, b)))} for b in sorted(adj[a] & remaining)
                     for rest in match(remaining - {a, b}))
    banks = [match(live - set(xs)) for xs in by_colour]
    return sum(1 for a in banks[0] for b in banks[1] if not a & b for c in banks[2]
               if not (a | b) & c and a | b | c == internal)


def defects(g, colours):
    cubic(g); audit_colouring(line_graph(g), 6, colours)
    result = [0] * g.n
    for e, c in zip(g.edges, colours):
        for v in e:
            result[v] ^= 1 + c // 2
    return tuple(result)


def shade_fibre(g, pairs):
    """Exact source-specific lifting obstruction and count after forgetting shades."""
    cubic(g)
    require(type(pairs) is tuple and len(pairs) == len(g.edges), "pair vector shape")
    for c in pairs:
        nat(c); require(c < 3, "pair label outside palette")
    constraints = []
    for v in range(g.n):
        incident = tuple(i for i,e in enumerate(g.edges) if v in e)
        if len({pairs[i] for i in incident}) == 1:
            return {"status":"monochromatic_star", "vertex":v, "count":0}
        constraints.extend((a,b) for a,b in combinations(incident,2) if pairs[a] == pairs[b])
    cg = graph(len(g.edges),constraints)
    adj = adjacency(cg)
    require(max(map(len,adj),default=0) <= 2, "cubic shade degree invariant")
    phase = [-1]*cg.n; parents = {}; components = []
    for seed in range(cg.n):
        if phase[seed] >= 0:
            continue
        phase[seed] = 0; parents[seed] = None
        todo = [seed]; component = []
        for a in todo:
            component.append(a)
            for b in sorted(adj[a]):
                if phase[b] < 0:
                    phase[b] = phase[a]^1; parents[b] = a; todo.append(b)
                elif phase[b] == phase[a]:
                    left = [a]
                    while parents[left[-1]] is not None:
                        left.append(parents[left[-1]])
                    right = [b]
                    while right[-1] not in left:
                        right.append(parents[right[-1]])
                    cycle = tuple(left[:left.index(right[-1])+1]+right[-2::-1])
                    require(len(cycle)%2 == 1 and len(set(cycle)) == len(cycle), "odd cycle reader")
                    require(all(tuple(sorted(e)) in cg.edges
                                for e in zip(cycle,cycle[1:]+cycle[:1])), "actual shade cycle")
                    return {"status":"odd_shade_cycle", "constraint_graph":cg,
                            "cycle":cycle, "count":0}
        components.append(tuple(sorted(component)))
    require(len(cg.edges) == sum(d != 0 for d in (
        __import__('functools').reduce(int.__xor__,(pairs[i]+1 for i,e in enumerate(g.edges)
                                                  if v in e),0) for v in range(g.n))),
            "one constraint per defective star")
    return {"status":"liftable", "constraint_graph":cg, "phase":tuple(phase),
            "components":tuple(components), "count":2**len(components)}


def shade_lift(g, pairs, component_bits):
    fibre = shade_fibre(g,pairs)
    require(fibre["status"] == "liftable", "pair assignment has no six-colour lift")
    require(type(component_bits) is tuple and len(component_bits) == len(fibre["components"]),
            "one bit per marked shade component")
    colours = [2*c+s for c,s in zip(pairs,fibre["phase"])]
    for vertices, bit in zip(fibre["components"],component_bits):
        nat(bit); require(bit < 2, "component shade bit")
        for v in vertices:
            colours[v] ^= bit
    return audit_colouring(line_graph(g),6,tuple(colours))


def lift_two_defects(p, labels, central_pair, shades):
    """Source-bound full-colour lift; no presumed global snark theorem needed."""
    word = boundary_word(p, labels); nat(central_pair)
    require(central_pair < 3, "pair label")
    a, aa, b, bb = word
    require(a == aa and b == bb and central_pair not in (a, b), "two-defect pair guard")
    require(type(shades) is tuple and len(shades) == len(p.source.edges), "shade vector")
    for s in shades:
        nat(s); require(s < 2, "binary shade")
    require(shades[p.ports[0]] != shades[p.ports[1]] and
            shades[p.ports[2]] != shades[p.ports[3]], "marked shade XOR constraints")
    pair = [central_pair] * len(shades)
    for i, c in zip(p.active, labels):
        pair[i] = c
    result = tuple(2 * c + s for c, s in zip(pair, shades))
    actual = defects(p.source, result)
    endpoints = set(p.source.edges[p.removed_edge])
    require(all(d == (central_pair + 1 if i in endpoints else 0)
                for i, d in enumerate(actual)), "full defect reconstruction failed")
    return result


def palette_graph(p):
    audit_pole(p)
    return graph(p.constraints.n + 3, (*combinations(range(3), 2),
                 *((a + 3, b + 3) for a, b in p.constraints.edges),
                 *((a, b + 3) for a in range(3) for b in range(p.constraints.n))))


def palette_lift(p, labels, palette):
    boundary_word(p, labels)
    require(type(palette) is tuple and len(palette) == 3, "palette triple")
    for c in palette:
        nat(c); require(c < 6, "palette colour")
    require(len(set(palette)) == 3, "palette must be injective")
    remaining = sorted(set(range(6)) - set(palette))
    result = palette + tuple(remaining[c] for c in labels)
    return audit_colouring(palette_graph(p), 6, result)


def wire_relation(length):
    """length is number of suppressed zero-defect degree-two vertices."""
    nat(length)
    return tuple(tuple(int(b == (a ^ (length % 2))) for b in range(6)) for a in range(6))


def compose_relation(a, b):
    for matrix in (a, b):
        require(type(matrix) is tuple and len(matrix) == 6, "six-port matrix")
        for row in matrix:
            require(type(row) is tuple and len(row) == 6, "six-port row")
            for x in row:
                nat(x)
    return tuple(tuple(sum(a[i][j] * b[j][k] for j in range(6)) for k in range(6))
                 for i in range(6))


PETERSEN_EDGES = ((0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,8),(4,9),
                  (5,7),(5,8),(6,8),(6,9),(7,9))


def audit_subdivision(g, paths):
    audit_graph(g)
    require(type(paths) is tuple and len(paths) == 15, "one path per marked Petersen edge")
    branch = {}; interiors = set(); edges = set(g.edges)
    for (a, b), path in zip(PETERSEN_EDGES, paths):
        require(type(path) is tuple and len(path) >= 2, "nonempty path")
        for x in path:
            nat(x); require(x < g.n, "path source vertex")
        require(len(set(path)) == len(path), "path is not simple")
        require(all(tuple(sorted(e)) in edges for e in zip(path, path[1:])), "missing path edge")
        for v, x in ((a, path[0]), (b, path[-1])):
            require(v not in branch or branch[v] == x, "inconsistent branch mark")
            branch[v] = x
        require(not interiors & set(path[1:-1]), "paths have overlapping interiors")
        interiors.update(path[1:-1])
    require(len(set(branch.values())) == 10 and not set(branch.values()) & interiors,
            "branch/interior collision")
    return tuple(branch[i] for i in range(10))


def flower(k):
    nat(k, 5); require(k % 2 == 1, "odd flower length")
    edges = [(4*i, 4*i+j) for i in range(k) for j in (1,2,3)]
    edges += [(4*i+1, 4*((i+1) % k)+1) for i in range(k)]
    cycle = tuple(4*i+j for j in (2,3) for i in range(k))
    edges += list(zip(cycle, cycle[1:]+cycle[:1]))
    return graph(4*k, edges)


def examples():
    b1 = ((0,5),(1,17),(2,14),(3,8),(4,17),(6,11),(7,17),(9,13),(10,15),(12,16),
          *((i,i+1) for i in range(16)),(0,16))
    b2 = ((0,1),(0,2),(0,14),(1,5),(1,11),(2,3),(2,6),(3,4),(3,9),(4,5),(4,7),
          (5,6),(6,8),(7,8),(7,17),(8,9),(9,15),(10,11),(10,14),(10,16),(11,12),
          (12,13),(12,17),(13,14),(13,15),(15,16),(16,17))
    des = ((0,1),(0,5),(0,19),(1,2),(1,16),(2,3),(2,11),(3,4),(3,14),(4,5),(4,9),
           (5,6),(6,7),(6,15),(7,8),(7,18),(8,9),(8,13),(9,10),(10,11),(10,19),
           (11,12),(12,13),(12,17),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19))
    return (("Petersen",graph(10,PETERSEN_EDGES)),("BlanusaI",graph(18,b1)),
            ("BlanusaII",graph(18,b2)),("J5",flower(5)),("J7",flower(7)),
            ("K4",graph(4,combinations(range(4),2))),("Desargues",graph(20,des)))


def subdivisions():
    paths = {
      "Petersen": PETERSEN_EDGES,
      "J5": ((2,19),(2,6,10,8,9,13),(2,0,3),(19,16),(19,15),(16,17),(16,18),
             (17,13),(17,1,5,4,7),(13,12),(3,18),(3,7),(15,11,7),(15,12),(18,14,12)),
      "BlanusaI": ((7,8),(7,6,11),(7,17,1,2,14),(8,3,4,5,0,16),(8,9),(16,12),
                   (16,15),(12,11),(12,13),(11,10),(14,15),(14,13),(9,13),(9,10),(15,10)),
      "BlanusaII": ((8,7,17),(8,6,2,0,14),(8,9),(17,16),(17,12),(16,10),(16,15),
                    (10,14),(10,11),(14,13),(9,15),(9,3,4,5,1,11),(12,11),(12,13),(15,13)),
      "Desargues": ((7,8),(7,6,15),(7,18),(8,9,10),(8,13),(10,11,2,1,16),(10,19),
                    (16,15),(16,17),(15,14),(18,19),(18,17),(13,12,17),(13,14),(19,0,5,4,3,14))}
    replace = {(17,1):(17,21,25,1),(18,3):(18,22,26,3),(19,2):(19,23,27,2)}
    replace.update({e[::-1]: p[::-1] for e,p in tuple(replace.items())})
    paths["J7"] = tuple(tuple([p[0]] + [x for e in zip(p,p[1:])
                      for x in replace.get(e,e)[1:]]) for p in paths["J5"])
    return paths


def perfect_pairings(items):
    if not items:
        yield ()
    else:
        a = items[0]
        for b in items[1:]:
            for rest in perfect_pairings(tuple(x for x in items[1:] if x != b)):
                yield ((a,b),)+rest


CHECKS = 0
def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def rejects(fn):
    try:
        fn()
    except ValueError:
        check(True, "reject")
    else:
        check(False, "malformed input accepted")


def main():
    records = []
    paths = subdivisions()
    for name, g in examples():
        cubic(g)
        adj = adjacency(g)
        for removed in g.edges:
            seen = {0}; todo = [0]
            while todo:
                v = todo.pop()
                for w in adj[v]:
                    if tuple(sorted((v,w))) != removed and w not in seen:
                        seen.add(w); todo.append(w)
            check(len(seen) == g.n, "bridgeless connected source")
        tait = tuple(colourings(line_graph(g),3))
        check(bool(tait) == (name in ("K4","Desargues")), "exact Tait status")
        for colours in tait:
            audit_colouring(line_graph(g),3,colours)
        p = make_pole(g,0)
        counts, witnesses = profile(p)
        for word in product(range(3),repeat=4):
            check(counts[word] == matching_count(p,word), "independent perfect matching profile")
        row = {"name":name,"vertices":g.n,"edges":len(g.edges),"removed_edge":g.edges[0],
               "ordered_ports":tuple(g.edges[i] for i in p.ports),
               "tait_colourings":len(tait),"four_pole_mass":sum(counts.values()),
               "boundary_profile":[[w,c] for w,c in sorted(counts.items())]}
        if not tait:
            m = counts[(0,0,0,0)]
            check(m > 0 and m % 2 == 0, "positive even inherited profile")
            check(counts == Counter({(a,a,b,b):m for a in range(3) for b in range(3)}),
                  "complete snark boundary relation")
            for word, labels in witnesses.items():
                for c in range(3):
                    if c in word:
                        continue
                    shades = [0] * len(g.edges)
                    shades[p.ports[1]] = shades[p.ports[3]] = 1
                    for flip in (0,1):
                        lifted = lift_two_defects(p,labels,c,tuple(s^flip for s in shades))
                        check(sum(d != 0 for d in defects(g,lifted)) == 2, "exact two defects")
                        check(__import__('functools').reduce(int.__xor__, defects(g,lifted),0)==0,
                              "charge conservation")
            total = sum(counts[w] * (3-len(set(w))) for w in counts) * 2**(len(g.edges)-2)
            check(total == 3*m*2**len(g.edges), "counted defect transfer")
            row.update(m_e=m,two_defect_six_colourings=total,per_charge=m*2**len(g.edges),
                       minimum_defect_vertices=2,palette_six_colourings=1080*m)
            full_palette_profile = Counter()
            for palette in permutations(range(6),3):
                remaining = sorted(set(range(6))-set(palette))
                for word, multiplicity in counts.items():
                    full_palette_profile[tuple(remaining[c] for c in word)] += multiplicity
            check(full_palette_profile == Counter({(a,a,b,b):(60 if a==b else 24)*m
                       for a,b in product(range(6),repeat=2)}), "full weighted palette relation")
            check(sum(full_palette_profile.values()) == 1080*m, "complete palette count")
            row["palette_four_port_multiplicities"] = {"aaaa":60*m,"aabb_a_ne_b":24*m}
            for labels in witnesses.values():
                col = palette_lift(p,labels,(0,2,5))
                i,j = (3+p.active.index(e) for e in p.ports[:2])
                check(col[i] == col[j], "actual abstract six-equality port")
            pg = palette_graph(p)
            check(all(e in pg.edges for e in combinations(range(4),2)), "palette geometry K4 hostile")
        if name in paths:
            row["petersen_branch_vertices"] = audit_subdivision(g,paths[name])
            row["petersen_subdivision"] = paths[name]
            check(True,"authenticated minor")
        if name == "Desargues":
            row["three_edge_colour_witness"] = tait[0]
            row["source_edges"] = g.edges
        # Test the full two-port relation, not just one chosen pair.
        lg = line_graph(g)
        for a,b in combinations(range(lg.n),2):
            linked = (a,b) in lg.edges
            for ca,cb in product(range(6),repeat=2):
                if ca == cb and linked:
                    continue
                col = greedy_extend(lg,6,((a,ca),(b,cb)))
                check((col[a],col[b]) == (ca,cb), "full marked six-boundary support")
        records.append(row)
    # Independent binary enumeration for one Petersen pair assignment.
    p = make_pole(examples()[0][1],0)
    count = 0
    for bits in product((0,1),repeat=len(p.source.edges)):
        count += bits[p.ports[0]] != bits[p.ports[1]] and bits[p.ports[2]] != bits[p.ports[3]]
    check(count == 2**13,"independent shade multiplicity")
    for r,s in product(range(13),repeat=2):
        check(compose_relation(wire_relation(r),wire_relation(s)) == wire_relation(r+s),
              "full shade-sensitive wire composition")
    for r in range(4):
        counts = Counter()
        for cols in product(range(6),repeat=r+1):
            if all(a != b and a//2 == b//2 for a,b in zip(cols,cols[1:])):
                counts[(cols[0],cols[-1])] += 1
        check(tuple(tuple(counts[(a,b)] for b in range(6)) for a in range(6)) == wire_relation(r),
              "literal zero-charge wire")
    check(wire_relation(1)[0][1] == 1 and wire_relation(1)[0][0] == 0,"shade-blind hostile")
    path = graph(3,((0,1),(1,2)))
    path_counts = Counter((c[0],c[2]) for c in colourings(path,6))
    check(path_counts[(0,0)] == 5 and path_counts[(0,1)] == 4,"ordinary vertex path hostile")
    # Six colours as K4 edges, with opposite-edge pairing.
    k4edges = tuple(combinations(range(4),2))
    opposite = ((0,5),(1,4),(2,3))
    labels = {v:i+1 for i,pair in enumerate(opposite) for v in pair}
    shapes = Counter()
    for triple in combinations(range(6),3):
        deg = sorted(Counter(v for i in triple for v in k4edges[i]).values())
        shape = "star" if deg == [1,1,1,3] else "triangle" if deg == [2,2,2] else "P4"
        charge = labels[triple[0]] ^ labels[triple[1]] ^ labels[triple[2]]
        check((charge == 0) == (shape != "P4"),"tetrahedral triple classification")
        shapes[shape] += 1
    pairings = tuple(perfect_pairings(tuple(range(6))))
    check(len(pairings) == 15 and shapes == Counter(star=4,triangle=4,P4=12),"small-object counts")
    for triple in combinations(range(6),3):
        zero = 0
        for pairs in pairings:
            pi = {v:i+1 for i,pair in enumerate(pairs) for v in pair}
            zero += pi[triple[0]] ^ pi[triple[1]] ^ pi[triple[2]] == 0
        check(zero == 6,"pairing-average erases local graph information")
    # Exact pair/shade partition function, independently counted on K4.
    k4 = graph(4,combinations(range(4),2))
    by_pairs = Counter()
    for pairs in product(range(3),repeat=6):
        fibre = shade_fibre(k4,pairs)
        if fibre["count"]:
            col = shade_lift(k4,pairs,(0,)*len(fibre["components"]))
            count_defects = sum(x != 0 for x in defects(k4,col))
            by_pairs[count_defects] += fibre["count"]
            check(tuple(c//2 for c in col) == pairs,"source-bound pair/shade lift")
    literal = Counter(sum(x != 0 for x in defects(k4,col)) for col in colourings(line_graph(k4),6))
    check(by_pairs == literal,"complete K4 defect partition function")
    prism = graph(6,((0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)))
    prism_pairs = tuple(0 if max(e)<3 else 2 if min(e)>=3 else 1 for e in prism.edges)
    blocked = shade_fibre(prism,prism_pairs)
    charges = tuple(__import__('functools').reduce(int.__xor__,
                (prism_pairs[i]+1 for i,e in enumerate(prism.edges) if v in e),0) for v in range(6))
    check(charges == (2,)*6 and __import__('functools').reduce(int.__xor__,charges,0)==0,
          "global pair charge does not see shade holonomy")
    check(blocked["status"] == "odd_shade_cycle" and len(blocked["cycle"])==3,
          "actual cubic shade obstruction")
    # Adversarial exact-type and source-reconstruction controls.
    from dataclasses import replace
    rejects(lambda:audit_graph(Graph(True,())))
    rejects(lambda:audit_graph(Graph(3,((0,1.0),))))
    rejects(lambda:audit_graph(Graph(3,((0,1),(0,1)))))
    rejects(lambda:audit_pole(replace(p,removed_edge=False)))
    rejects(lambda:audit_pole(replace(p,ports=tuple(float(x) for x in p.ports))))
    rejects(lambda:audit_pole(replace(p,ports=p.ports[::-1])))
    rejects(lambda:greedy_extend(path,6,((0,True),)))
    rejects(lambda:greedy_extend(path,6,((0,0),(1,0))))
    rejects(lambda:wire_relation(True))
    rejects(lambda:wire_relation(-1))
    rejects(lambda:audit_subdivision(p.source,((False,1),)+PETERSEN_EDGES[1:]))
    rejects(lambda:audit_subdivision(p.source,((0,2,1),)+PETERSEN_EDGES[1:]))
    witness = profile(p)[1][(0,0,0,0)]
    rejects(lambda:lift_two_defects(p,witness,1,(0,)*15))
    rejects(lambda:palette_lift(p,witness,(0,0,1)))
    rejects(lambda:shade_fibre(k4,(False,0,1,1,2,2)))
    rejects(lambda:shade_lift(prism,prism_pairs,()))
    return {"status":"PROVED elementary transfer; FINITE-EXACT declared graphs; no unit realization",
            "checks":CHECKS,"graphs":records,"tetrahedral_triples":dict(shapes),
            "pairing_counts":{"total":15,"zero_at_fixed_triple":6,"nonzero":9},
            "K4_defect_partition_function":dict(sorted(literal.items())),
            "prism_shade_obstruction":{"edges":prism.edges,"pairs":prism_pairs,
                "charges":charges,"odd_cycle_edge_indices":blocked["cycle"]},
            "wire_parity_relation":"output = input xor (suppressed_vertices mod 2)",
            "ordinary_path_endpoint_counts":{"equal":5,"different":4},
            "scope":"No six-colour equality inference from a Petersen minor or from pair charges alone."}


if __name__ == "__main__":
    print(json.dumps(main(),sort_keys=True,indent=2))
