#!/usr/bin/env python3
"""Exact flower six-port reductions and rooted Petersen boundary certificates.

Universe: every boundary word of a flower claw; strip powers through 16;
twisted/untwisted closures through 41; independent whole-graph matching
enumeration through flower length 15; explicit minor chains through 101;
four marked edge types through J15; both inherited Blanusa graphs; K4 and
a bridged cubic hostile. Standard library only. All checks survive -O.
"""

from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import runpy


CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def edge(a, b):
    return tuple(sorted((a, b)))


def graph(edges):
    edges = frozenset(edge(*e) for e in edges)
    check(all(a != b for a, b in edges), "loop")
    return edges


def adjacency(edges):
    result = {}
    for a, b in edges:
        result.setdefault(a, set()).add(b)
        result.setdefault(b, set()).add(a)
    return result


def cubic(edges):
    adj = adjacency(edges)
    check(all(len(v) == 3 for v in adj.values()), "not cubic")
    return adj


def flower(k, twist=True):
    check(k >= 3, "flower length")
    e = [(4*i, 4*i+j) for i in range(k) for j in (1, 2, 3)]
    e += [(4*i+j, 4*(i+1)+j) for i in range(k-1) for j in (1, 2, 3)]
    terminal = (1, 3, 2) if twist else (1, 2, 3)
    e += [(4*(k-1)+j, terminal[j-1]) for j in (1, 2, 3)]
    result = graph(e)
    cubic(result)
    return result


def matmul(a, b):
    n = len(a)
    out = [[0]*n for _ in range(n)]
    for i in range(n):
        for j, x in enumerate(a[i]):
            if x:
                for k, y in enumerate(b[j]):
                    if y:
                        out[i][k] += x*y
    return out


def support(a):
    return frozenset((i, j) for i, row in enumerate(a)
                     for j, value in enumerate(row) if value)


def components(a):
    todo, result = set(range(len(a))), []
    while todo:
        queue, comp = [min(todo)], set()
        while queue:
            u = queue.pop()
            if u in comp:
                continue
            comp.add(u)
            queue.extend(v for v, x in enumerate(a[u]) if x and v not in comp)
        todo -= comp
        result.append(comp)
    return result


def transfer_audit():
    words = tuple(product(range(3), repeat=3))
    index = {w: i for i, w in enumerate(words)}
    t = [[0]*27 for _ in words]
    # Enumerate the six possible centre-spoke permutations independently.
    for a in words:
        for spokes in permutations(range(3)):
            if all(a[j] != spokes[j] for j in range(3)):
                b = tuple(3-a[j]-spokes[j] for j in range(3))
                t[index[a]][index[b]] += 1
    comps = components(t)
    check(Counter(len(c) for c in comps) == {1: 3, 3: 2, 6: 3},
          "wrong transition components")
    check(all(sum(t[i]) == (0 if len(set(w)) == 1 else 2)
              for i, w in enumerate(words)), "transition degrees")
    check(all(t[i][j] == t[j][i] for i in range(27) for j in range(27)),
          "transition asymmetry")
    t2 = matmul(t, t)
    blocks = [c for c in components(t2) if any(t2[min(c)])]
    check(len(blocks) == 8 and all(len(c) == 3 for c in blocks),
          "T squared block partition")
    for c in blocks:
        for i in c:
            for j in range(27):
                check(t2[i][j] == (int(i == j)+1 if j in c else 0),
                      "T squared differs from I+J")
    powers = [None, t]
    for k in range(2, 44):
        powers.append(matmul(powers[-1], t))
    for k in range(2, 15):
        check(support(powers[k]) == support(powers[k+2]),
              "strip support period two")
    check(len(support(t)) == 48 and len(support(t2)) == 72,
          "support census")
    check(support(t) != support(powers[3]), "T and T cubed hostile lost")
    identity = [[int(i == j) for j in range(27)] for i in range(27)]
    check(support(t2) != support(identity), "bare wires hostile lost")
    check(t2[index[(0, 0, 0)]][index[(0, 0, 0)]] == 0,
          "constant boundary hostile")
    closures = {}
    for k in range(3, 42):
        twisted = sum(powers[k][i][index[(w[0], w[2], w[1])]]
                      for i, w in enumerate(words))
        untwisted = sum(powers[k][i][i] for i in range(27))
        cycled = sum(powers[k][i][index[(w[1], w[2], w[0])]]
                     for i, w in enumerate(words))
        check(twisted == (0 if k % 2 else 3*2**(k+1)),
              "twisted closure count")
        check(untwisted == (2**(k+1)-4 if k % 2 else 2**(k+3)+16),
              "untwisted closure count")
        check(cycled == (2**(k+1)+2 if k % 2 else 2**(k+3)-8),
              "three-cycle closure count")
        closures[k] = twisted, untwisted
    print("transfer: 2 K3 + 3 C6 + 3 isolates; T^2 = 8*(I3+J3) + 3*[0]")
    print("supports: |T|=48, |T^2|=72, supp(T^(r+2))=supp(T^r), r=2..14")
    print("closure formulas checked k=3..41: twisted odd=0, even=3*2^(k+1);"
          " untwisted odd=2^(k+1)-4, even=2^(k+3)+16")
    print("three-cycle strand closure: odd=2^(k+1)+2, even=2^(k+3)-8")
    return closures


def matching_colour_count(edges):
    """Independent enumeration: colour-0 matching plus even two-factor cycles."""
    adj = cubic(edges)
    vertices = tuple(sorted(adj))
    position = {v: i for i, v in enumerate(vertices)}
    neighbors = {v: sum(1 << position[w] for w in adj[v]) for v in vertices}
    count = matching_count = 0

    def visit(mask, selected):
        nonlocal count, matching_count
        if mask == 0:
            matching_count += 1
            other = {v: set(adj[v]) for v in vertices}
            for a, b in selected:
                other[a].remove(b)
                other[b].remove(a)
            unseen, cycles = set(vertices), 0
            while unseen:
                start = min(unseen)
                queue, component = [start], set()
                while queue:
                    v = queue.pop()
                    if v in component:
                        continue
                    component.add(v)
                    queue.extend(other[v] - component)
                if len(component) % 2:
                    return
                unseen -= component
                cycles += 1
            count += 2**cycles
            return
        live = [v for v in vertices if mask >> position[v] & 1]
        a = min(live, key=lambda v: (mask & neighbors[v]).bit_count())
        choices = mask & neighbors[a]
        while choices:
            bit = choices & -choices
            choices ^= bit
            b = vertices[bit.bit_length()-1]
            visit(mask ^ (1 << position[a]) ^ bit, selected+((a, b),))

    visit((1 << len(vertices))-1, ())
    return count, matching_count


def fourpole_count(edges, removed, word):
    """Exact MRV recursion on internal edges with four fixed boundary colours."""
    adj = cubic(edges)
    u, v = removed
    terminals = tuple(sorted(adj[u]-{v}) + sorted(adj[v]-{u}))
    live = sorted(set(adj)-{u, v})
    used = {x: 0 for x in live}
    for x, color in zip(terminals, word):
        bit = 1 << color
        if used[x] & bit:
            return 0
        used[x] |= bit
    internal = tuple(e for e in sorted(edges) if u not in e and v not in e)

    def visit(left):
        if not left:
            return 1
        chosen, choices, size = None, 0, 4
        for i in left:
            a, b = internal[i]
            possible = 7 & ~(used[a] | used[b])
            n = possible.bit_count()
            if n < size:
                chosen, choices, size = i, possible, n
                if n <= 1:
                    break
        if not choices:
            return 0
        a, b = internal[chosen]
        remain, total = left-{chosen}, 0
        while choices:
            bit = choices & -choices
            choices ^= bit
            used[a] |= bit
            used[b] |= bit
            total += visit(remain)
            used[a] ^= bit
            used[b] ^= bit
        return total

    return visit(frozenset(range(len(internal))))


REPS = ((0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0))


def smooth(edges, removed, pairing):
    adj = cubic(edges)
    u, v = removed
    boundary = tuple(sorted(adj[u]-{v})+sorted(adj[v]-{u}))
    result = {e for e in edges if u not in e and v not in e}
    result.update(edge(boundary[a], boundary[b]) for a, b in pairing)
    return graph(result)


def reduce_flower(edges, k):
    """Keep modules 0,3,...,k-1; remove centres 1,2, suppress their leaves."""
    work = set(edges)
    for center in (4, 8):
        work = {e for e in work if center not in e}
    for v in (5, 6, 7, 9, 10, 11):
        adj = adjacency(work)
        check(len(adj[v]) == 2, "suppression vertex degree")
        a, b = sorted(adj[v])
        work = {e for e in work if v not in e}
        check(edge(a, b) not in work, "unwanted parallel edge")
        work.add(edge(a, b))
    relabel = {4*i+j: 4*(i if i == 0 else i-2)+j
               for i in (0, *range(3, k)) for j in range(4)}
    result = graph((relabel[a], relabel[b]) for a, b in work)
    check(result == flower(k-2), "minor reduction has wrong flower wiring")
    return result


PETERSEN = graph(((0, 1), (0, 4), (0, 5), (1, 2), (1, 6),
                  (2, 3), (2, 7), (3, 4), (3, 8), (4, 9),
                  (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)))


def minor_chain_audit():
    # Contract the b0,b1,b2 triangle of J3; the other labels are an explicit
    # isomorphism to the inherited Petersen labelling.
    pmap = {1: 0, 5: 0, 9: 0, 0: 1, 4: 4, 8: 5,
            2: 2, 6: 3, 10: 8, 3: 6, 7: 9, 11: 7}
    contracted = graph((pmap[a], pmap[b]) for a, b in flower(3)
                       if pmap[a] != pmap[b])
    check(contracted == PETERSEN, "J3 triangle is not Petersen")
    for k in (5, 7, 9, 11, 15, 21, 51, 101):
        work = flower(k)
        for j in range(k, 3, -2):
            work = reduce_flower(work, j)
        check(work == flower(3), "minor chain endpoint")
        print(f"minor J{k}: vertices={4*k}, six-port replacements={(k-3)//2},"
              " final triangle contraction=Petersen")


def pole_and_smoothing_audit():
    inherited = runpy.run_path(str(Path(__file__).with_name(
        "snark_boundary_fourpole_thm4116_independent.py")))
    examples = {name: graph(edges) for name, _, edges
                in inherited["graph_definitions"]()}
    # Direct inherited matching engine independently checks the two marked
    # Blanusa four-poles. Chosen source edges are length-one paths in the
    # authenticated Petersen subdivisions inherited below.
    for name, marked in (("BlanusaI", (7, 8)), ("BlanusaII", (9, 15))):
        edges = examples[name]
        sig = tuple(fourpole_count(edges, marked, w) for w in REPS)
        alternate = tuple(inherited["matching_extension_count"](
            tuple(sorted(adjacency(edges))), tuple(sorted(edges)), marked, w)
            for w in REPS)
        check(sig == alternate, "independent Blanusa profile")
        check(sig[0] == sig[1] > 0 and sig[2:] == (0, 0), "Blanusa signature")
        pairs = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
        counts = tuple(matching_colour_count(smooth(edges, marked, p))[0]
                       for p in pairs)
        check(counts == (9*sig[0], 3*sig[0], 3*sig[0]), "smoothing count")
        print(f"{name}: marked_edge={marked}, signature={sig},"
              f" three_smoothing_counts={counts}")
    for k in (3, 5, 7, 9, 11, 13, 15):
        edges = flower(k)
        signatures = {}
        for marked in ((0, 1), (0, 2), (1, 5), (2, 6)):
            sig = tuple(fourpole_count(edges, marked, w) for w in REPS)
            check(sig[0] == sig[1] and (sig[0] > 0 or k == 3)
                  and sig[2:] == (0, 0),
                  "flower four-pole signature")
            signatures[marked] = sig[0]
        expected = {
            (0, 1): (2**k+4)//3,
            (0, 2): (2**k-2)//3,
            (1, 5): (2**(k-1)-4)//3,
            (2, 6): (2**(k-1)+2)//3,
        }
        check(signatures == expected, "flower four-orbit formulas")
        print(f"J{k} four-poles: m_e={signatures}; all (m,m,0,0)")
    # Ten-vertex simple cubic bridge example: two subdivided K4s.
    bridged = set(combinations(range(4), 2)) | set(combinations(range(4, 8), 2))
    bridged -= {(0, 1), (4, 5)}
    bridged |= {(0, 8), (1, 8), (4, 9), (5, 9), (8, 9)}
    bridged = graph(bridged)
    sig = tuple(fourpole_count(bridged, (8, 9), w) for w in REPS)
    check(sig == (4, 4, 0, 0), "bridged hostile signature")
    check(len(adjacency(bridged)) == 10 and len(bridged) == 15,
          "bridged hostile size")
    check(matching_colour_count(bridged)[0] == 0, "bridge should obstruct")
    print("bridged hostile: 10 vertices, 15 edges, signature=(4,4,0,0),"
          " no Petersen minor by equal-order/equal-size rigidity")
    # Authenticate inherited subdivisions; extract a rooted four-pole
    # reduction by retaining an actual branch edge.
    inherited_transfer = runpy.run_path(str(Path(__file__).parent /
        "experiments" / "snark_boundary_transfer_20261018a.py"))
    inherited_examples = dict(inherited_transfer["examples"]())
    paths = inherited_transfer["subdivisions"]()
    for name, marked in (("BlanusaI", (7, 8)), ("BlanusaII", (9, 15))):
        branch = inherited_transfer["audit_subdivision"](
            inherited_examples[name], paths[name])
        check(any(edge(*p) == marked for p in paths[name] if len(p) == 2),
              "marked edge is not a branch edge")
        check(set(marked) <= set(branch), "marked endpoints not branches")
        print(f"{name} rooted Petersen minor: authenticated 15 paths,"
              f" retained branch edge {marked}, complementary B-pole replacement")


def main():
    closures = transfer_audit()
    for k in range(3, 16):
        twisted, matching_count = matching_colour_count(flower(k))
        check(twisted == closures[k][0], "matching/transfer disagreement")
        # Alternating closures are independently tested on both parities.
        if k <= 8:
            untwisted = matching_colour_count(flower(k, False))[0]
            check(untwisted == closures[k][1], "untwisted matching disagreement")
        print(f"matching J{k}: vertices={4*k}, perfect_matchings={matching_count},"
              f" labelled_Tait_colourings={twisted}")
    check(matching_colour_count(graph(combinations(range(4), 2)))[0] == 6,
          "K4 positive control")
    check(matching_colour_count(PETERSEN)[0] == 0, "Petersen negative control")
    minor_chain_audit()
    pole_and_smoothing_audit()
    print(f"PASS: {CHECKS} exact checks")


if __name__ == "__main__":
    main()
