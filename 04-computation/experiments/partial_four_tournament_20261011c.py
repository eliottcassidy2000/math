"""Marked four-arc square: exact recognition, completion and loss controls.

The four-state section is inherited, not a new tournament classification.
Only the stated finite universes are enumerated. No unit K4 is asserted.
"""
from collections import Counter
from itertools import combinations, permutations, product
import json

import tournament_recursive_four_state_20261004 as old
import plane_odd_cycle_channels_20261011b as channels

PAIRS = tuple(combinations(range(4), 2))
TYPES = {(0, 0): "T", (1, 0): "+", (0, 1): "-", (1, 1): "S"}
BASE = frozenset(((0, 1), (1, 2), (2, 3), (0, 3)))


def need(ok, text):
    if not ok:
        raise ValueError(text)


def audit(arcs, count):
    need(type(arcs) is frozenset and len(arcs) == count, "exact arc set and size")
    need(all(type(e) is tuple and len(e) == 2
             and all(type(v) is int and 0 <= v < 4 for v in e)
             and e[0] != e[1] for e in arcs), "four distinct marked vertex labels")
    need(all((v, u) not in arcs for u, v in arcs), "at most one orientation per pair")
    return arcs


def edges(arcs):
    return frozenset(tuple(sorted(e)) for e in arcs)


def degrees(arcs):
    return tuple(sum(u == v for u, _ in arcs) for v in range(4))


def user_pattern(arcs):
    audit(arcs, 4)
    score = degrees(arcs)
    holes = tuple(e for e in PAIRS if e not in edges(arcs))
    return (sorted(score) == [0, 1, 1, 2]
            and len({v for e in holes for v in e}) == 4
            and sorted(tuple(sorted((score[u], score[v]))) for u, v in holes)
            == [(0, 1), (1, 2)])


def frame(arcs):
    """Recover the unique a,b,c,d with ab,bc,cd,ad oriented forward."""
    need(user_pattern(arcs), "specified missing-edge degree pattern required")
    score = degrees(arcs)
    a, d = score.index(2), score.index(0)
    need((a, d) in arcs, "source and sink are adjacent")
    b, = (v for u, v in arcs if u == a and v != d)
    c, = (v for v in range(4) if v not in (a, b, d))
    need(arcs == frozenset(((a, b), (b, c), (c, d), (a, d))),
         "intrinsic Hamiltonian path and forward closing arc")
    return a, b, c, d


def bitpair(value):
    need(type(value) is tuple and len(value) == 2
         and all(type(v) is int and v in (0, 1) for v in value),
         "two exact marked orientation bits")
    return value


def complete(arcs, bits):
    a, b, c, d = frame(arcs)
    bitpair(bits)
    return arcs | frozenset((v, u) if flip else (u, v)
                            for (u, v), flip in zip(((a, c), (b, d)), bits))


def read(arcs, full):
    """Authenticate extension of this same marked partial tournament."""
    a, b, c, d = frame(arcs)
    audit(full, 6)
    need(arcs <= full, "completion must retain the actual four marked arcs")
    answer = (int((c, a) in full), int((d, b) in full))
    need(complete(arcs, answer) == full, "exact complete source reconstruction")
    return answer


def classify(full):
    audit(full, 6)
    matrix = tuple(tuple(int((u, v) in full) for v in range(4)) for u in range(4))
    return old.class4(matrix)


def from_type(arcs, kind):
    need(type(kind) is str and kind in TYPES.values(), "one of four exact types")
    bits, = (b for b, value in TYPES.items() if value == kind)
    return complete(arcs, bits)


def flip(arcs, full, port):
    need(type(port) is int and port in (0, 1), "named diagonal port")
    bits = list(read(arcs, full))
    bits[port] ^= 1
    return complete(arcs, tuple(bits))


def relabel(arcs, labels):
    need(type(labels) is tuple and len(labels) == 4
         and all(type(v) is int for v in labels) and sorted(labels) == list(range(4)),
         "exact permutation")
    return frozenset((labels[u], labels[v]) for u, v in arcs)


def all_orientations(pairs):
    return tuple(frozenset((v, u) if bit else (u, v)
                           for (u, v), bit in zip(pairs, bits))
                 for bits in product((0, 1), repeat=len(pairs)))


def all_completions(arcs):
    audit(arcs, 4)
    missing = tuple(e for e in PAIRS if e not in edges(arcs))
    return tuple(arcs | extra for extra in all_orientations(missing))


def cycle_edges(vertices):
    return frozenset(tuple(sorted((vertices[i], vertices[(i + 1) % len(vertices)])))
                     for i in range(len(vertices)))


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        need(ok, message)
        checks += 1
    def reject(fn):
        nonlocal checks
        try:
            fn()
        except ValueError:
            checks += 1
        else:
            raise ValueError("hostile accepted")

    universe = tuple(arcs for pairs in combinations(PAIRS, 4)
                     for arcs in all_orientations(pairs))
    check(len(universe) == 240, "complete four-arc oriented simple universe")
    selected = tuple(arcs for arcs in universe if user_pattern(arcs))
    check(len(selected) == 24, "all marked user-pattern objects")
    labelled_c4 = tuple(a for a in universe
                        if all(sum(v in e for e in edges(a)) == 2 for v in range(4)))
    check(len(labelled_c4) == 48, "three square graphs times sixteen orientations")
    signatures = Counter()
    for arcs in labelled_c4:
        score = tuple(sorted(degrees(arcs)))
        kinds = Counter(classify(t) for t in all_completions(arcs))
        signatures[(score, tuple(sorted(kinds.items())))] += 1
    expected = {
        ((0, 1, 1, 2), (("+", 1), ("-", 1), ("S", 1), ("T", 1))): 24,
        ((0, 1, 1, 2), (("S", 2), ("T", 2))): 12,
        ((1, 1, 1, 1), (("S", 4),)): 6,
        ((0, 0, 2, 2), (("T", 4),)): 6,
    }
    check(dict(signatures) == expected, "complete square-orientation classification")
    for arcs in universe:
        all_four = len({classify(t) for t in all_completions(arcs)}) == 4
        if arcs in labelled_c4:
            check(all_four == user_pattern(arcs), "four-type transversal iff on squares")
    for arcs in selected:
        f = frame(arcs)
        check(sum(relabel(arcs, p) == arcs for p in permutations(range(4))) == 1,
              "the marked partial orientation is rigid")
        for bits, kind in TYPES.items():
            full = complete(arcs, bits)
            check(read(arcs, full) == bits and classify(full) == kind,
                  "source-owned complete/read/type")
            check(from_type(arcs, kind) == full, "type decoder on this fixed section")
            for port in (0, 1):
                changed = tuple(bit ^ int(i == port) for i, bit in enumerate(bits))
                check(read(arcs, flip(arcs, full, port)) == changed, "named flip action")
            check(flip(arcs, flip(arcs, full, 0), 1)
                  == flip(arcs, flip(arcs, full, 1), 0), "two flips commute")
            for p in permutations(range(4)):
                marked = relabel(arcs, p)
                check(frame(marked) == tuple(p[v] for v in f), "frame relabelling covariance")
                check(complete(marked, bits) == relabel(full, p), "completion covariance")

    # Isomorphism classes are also independently compared through all permutations.
    representatives = {}
    for full in all_orientations(PAIRS):
        matrix = tuple(tuple(int((u, v) in full) for v in range(4)) for u in range(4))
        representatives.setdefault(classify(full), set()).add(old.canonical(matrix))
    check(set(representatives) == set(TYPES.values())
          and all(len(v) == 1 for v in representatives.values()),
          "score classification matches the four full permutation orbits")

    opposite = frozenset(((0, 1), (0, 3), (1, 2), (3, 2)))
    check(sorted(degrees(opposite)) == sorted(degrees(BASE)), "degree-only hostile")
    check(Counter(classify(t) for t in all_completions(opposite)) == {"T": 2, "S": 2},
          "opposite source/sink is not a four-type section")
    reject(lambda: frame(opposite))
    full = complete(BASE, (0, 0))
    moved = relabel(full, (1, 0, 2, 3))
    check(classify(moved) == classify(full), "unmarked type survives moved ports")
    reject(lambda: read(BASE, moved))
    check(old.class4(old.fixed_path(4, 3)) == old.class4(old.fixed_path(4, 4)) == "S"
          and old.class4(old.fixed_path(4, 2)) == "-"
          and old.class4(old.fixed_path(4, 5)) == "S",
          "inherited third-flip class collision")
    for bits, kind in TYPES.items():
        full = complete(BASE, bits)
        converse = frozenset((v, u) for u, v in full)
        reanchored = relabel(converse, (3, 2, 1, 0))
        check(read(BASE, reanchored) == tuple(reversed(bits)), "actual converse swaps ports")
    check(classify(complete(BASE, (1, 1))) != classify(complete(BASE, (0, 0))),
          "both-bit complement is not tournament converse")

    # Same underlying graph and same channel assignment, irrespective of arc bits.
    labels = (1, 0, 0, 0)
    cycle = edges(BASE)
    check(channels.lift(4, cycle, labels, 2)[0] is not None, "square channel lift")
    check(channels.lift(4, cycle | {(0, 2)}, labels, 2)[0] is not None,
          "cross-channel diagonal does not create a monochromatic cycle")
    for arcs in (BASE | {(1, 3)}, BASE | {(3, 1)}):
        colour, bad = channels.lift(4, edges(arcs), labels, 2)
        check(colour is None and set(bad) == {1, 2, 3},
              "either direction of same-channel diagonal gives identical odd-cycle defect")
    for bits in TYPES:
        colour, bad = channels.lift(4, edges(complete(BASE, bits)), labels, 2)
        check(colour is None and set(bad) == {1, 2, 3},
              "every tournament completion has the same graph-colouring obstruction")
    triangles = tuple(cycle_edges(c) for c in ((0, 1, 2), (0, 1, 3), (0, 2, 3)))
    check(triangles[0] ^ triangles[2] == cycle, "the square is a binary cycle-space sum")
    check(triangles[0] ^ triangles[1] ^ triangles[2] == cycle_edges((1, 2, 3)),
          "the omitted residual triangle is the sum of the star basis")
    check(all(len({labels[v] for edge in triangle for v in edge}) > 1 for triangle in triangles),
          "NAE predicates do not respect that linear sum")

    reject(lambda: complete(BASE, (False, 0)))
    reject(lambda: flip(BASE, full, True))
    reject(lambda: audit(frozenset(((False, 1), (1, 2), (2, 3), (0, 3))), 4))
    reject(lambda: read(BASE, BASE))
    reject(lambda: from_type(BASE, 0))
    reject(lambda: relabel(BASE, (0, 1, 2, 3.0)))
    reject(lambda: audit(frozenset(((0, 1), (1, 0), (2, 3), (0, 3))), 4))
    print(json.dumps({
        "scope": "marked partial-square tournament completions; no plane chromatic decision",
        "checks": checks,
        "four_arc_universe": len(universe),
        "labelled_square_orientations": len(labelled_c4),
        "user_pattern_objects": len(selected),
        "completion_table": [{"bits": list(b), "type": k,
                             "scores": sorted(degrees(complete(BASE, b)))}
                            for b, k in TYPES.items()],
        "square_orbit_sizes": {"adjacent_source_sink": 24, "opposite_source_sink": 12,
                              "directed_cycle": 6, "alternating": 6},
        "boundary": "arc direction, edge presence, and colour equality are different ports",
    }, indent=2))


if __name__ == "__main__":
    main()
