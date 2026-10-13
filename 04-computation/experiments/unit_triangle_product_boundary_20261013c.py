"""Colour-preserving deletion repairs the ternary five-frame's geometry.

An exact 81-point model; complete incidence and three-colour support checks.
No new four-/six-chromatic graph is claimed.  All checks survive -O.
"""
from dataclasses import dataclass, replace
from itertools import combinations, product
from collections import Counter
import json


def need(ok, why):
    if not ok:
        raise ValueError(why)


# Q(zeta_15), in its marked power basis. Phi_15 is irreducible over Q.
PHI = (1, -1, 0, 1, -1, 1, 0, -1, 1)
ZERO = (0,) * 8
ONE = (1,) + (0,) * 7


def coefficient_tuple(a):
    need(type(a) is tuple and len(a) == 8 and all(type(c) is int for c in a),
         "eight exact integer cyclotomic coefficients required")


def reduce(a):
    a = list(a) + [0] * max(0, 8 - len(a))
    for j in range(len(a) - 1, 7, -1):
        c = a[j]
        for k in range(9):
            a[j - 8 + k] -= c * PHI[k]
    return tuple(a[:8])


POWERS = tuple(reduce((0,) * j + (1,)) for j in range(15))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def times(a, b):
    out = [0] * 15
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return reduce(out)


def conjugate(a):
    out = ZERO
    for j, c in enumerate(a):
        out = add(out, scale(c, POWERS[-j % 15]))
    return out


def norm(a):
    return times(a, conjugate(a))


def difference(a, b):
    return add(a, scale(-1, b))


TRIANGLE = (ZERO, ONE, add(ONE, POWERS[5]))


def labels(r):
    need(type(r) is int and r >= 1, "positive exact factor count")
    return tuple(product(range(3), repeat=r))


def checked_edges(n, edges):
    need(type(n) is int and n >= 0, "exact vertex count")
    need(type(edges) is frozenset, "canonical immutable edge set")
    for e in edges:
        need(type(e) is tuple and len(e) == 2 and all(type(x) is int for x in e)
             and 0 <= e[0] < e[1] < n, "canonical exact edge endpoints")


def hamming_edges(r):
    v = labels(r)
    return frozenset((i, j) for i, j in combinations(range(len(v)), 2)
                     if sum(a != b for a, b in zip(v[i], v[j])) == 1)


def frame_edges(r):
    v = labels(r)
    return frozenset((i, j) for i, j in combinations(range(len(v)), 2)
                     if len({(a - b) % 3 for a, b in zip(v[i], v[j])}) == 1
                     and v[i] != v[j])


def points_four():
    out = []
    for v in labels(4):
        z = ZERO
        for j, a in enumerate(v):
            z = add(z, times(POWERS[3*j], TRIANGLE[a]))
        out.append(z)
    return tuple(out)


def unit_edges(points):
    need(type(points) is tuple, "marked point tuple")
    for p in points:
        coefficient_tuple(p)
    need(len(set(points)) == len(points), "distinct original vertices required")
    return frozenset((i, j) for i, j in combinations(range(len(points)), 2)
                     if norm(difference(points[i], points[j])) == ONE)


@dataclass(frozen=True)
class Support:
    factors: int
    extra_edges: frozenset
    phases: tuple


def compile_support(r, extra_edges):
    v = labels(r)
    checked_edges(len(v), extra_edges)
    surviving = []
    for signs in product((1, -1), repeat=r):
        if all(sum(s * (v[a][i] - v[b][i]) for i, s in enumerate(signs)) % 3
               for a, b in extra_edges):
            surviving.append(signs)
    return Support(r, extra_edges, tuple(surviving))


def audit_support(packet):
    need(type(packet) is Support, "source-owned support packet required")
    expected = compile_support(packet.factors, packet.extra_edges)
    need(type(packet.phases) is tuple, "exact phase list required")
    for signs in packet.phases:
        need(type(signs) is tuple and len(signs) == packet.factors
             and all(type(s) is int and s in (1, -1) for s in signs),
             "exact signed phase required")
    need(packet.phases == expected.phases, "complete phase support required")
    return True


def lift(packet, origin, signs):
    audit_support(packet)
    need(type(origin) is int and 0 <= origin < 3, "exact labelled palette origin")
    need(type(signs) is tuple and all(type(s) is int for s in signs)
         and signs in packet.phases, "retained source phase required")
    v = labels(packet.factors)
    colours = tuple((origin + sum(s*a for s, a in zip(signs, x))) % 3 for x in v)
    edges = hamming_edges(packet.factors) | packet.extra_edges
    need(all(colours[a] != colours[b] for a, b in edges), "literal source colouring")
    return colours


def colour_safe_edges(packet):
    """Edges true in EVERY proper 3-colouring; empty support is explicit vacuity."""
    audit_support(packet)
    v = labels(packet.factors)
    return frozenset((a, b) for a, b in combinations(range(len(v)), 2)
                     if all(sum(s * (v[a][i] - v[b][i]) for i, s in enumerate(signs)) % 3
                            for signs in packet.phases))


def partition(colours):
    lookup = {}
    return tuple(lookup.setdefault(c, len(lookup)) for c in colours)


def boundary_counts(packet, ports):
    audit_support(packet)
    need(type(ports) is tuple and all(type(p) is int and 0 <= p < 3**packet.factors for p in ports)
         and len(set(ports)) == len(ports), "distinct original boundary vertices")
    out = Counter()
    v = labels(packet.factors)
    for signs in packet.phases:
        for c in range(3):
            colours = tuple((c + sum(s*a for s, a in zip(signs, x))) % 3 for x in v)
            out[partition(tuple(colours[p] for p in ports))] += 1
    return out


def enumerate_three_colours(n, edges):
    """Independent backtracking path used only for complete r<=2 controls."""
    checked_edges(n, edges)
    earlier = [tuple(a for a, b in edges if b == j) for j in range(n)]
    def visit(c):
        if len(c) == n:
            yield c
            return
        for k in range(3):
            if all(c[a] != k for a in earlier[len(c)]):
                yield from visit(c + (k,))
    return tuple(visit(()))


def main():
    count = 0
    def check(ok, why):
        nonlocal count
        count += 1
        if not ok:
            raise RuntimeError(why)
    def rejects(fn, *args):
        try:
            fn(*args)
        except ValueError:
            check(True, "rejected")
        else:
            check(False, "hostile accepted")

    # Field operations cross-checked in the independent inherited Q(zeta_60) implementation.
    from unit_cayley_realization_20261013a import Cyclotomic, make_cayley, classify, Collapse
    other = Cyclotomic(60)
    def convert(a):
        z = other.zero
        for i, c in enumerate(a):
            z = other.add(z, other.scale(c, other.powers[4*i]))
        return z
    for a, b in product(range(8), repeat=2):
        check(convert(times(POWERS[a], POWERS[b])) ==
              other.times(convert(POWERS[a]), convert(POWERS[b])), "independent multiplication")
    for a in POWERS:
        check(convert(conjugate(a)) == other.conjugate(convert(a)), "independent conjugation")
        check(norm(a) == ONE, "roots of unity exact unit norm")

    # Exhaust all 7^4 possible factor-difference words, independently of point pairs.
    differences = tuple(sorted({difference(x, y) for x in TRIANGLE for y in TRIANGLE}))
    nonzero = tuple(d for d in differences if d != ZERO)
    check(len(differences) == 7 and len(nonzero) == 6, "complete triangle differences")
    unit_words = []
    for ds in product(differences, repeat=4):
        z = ZERO
        for j, d in enumerate(ds):
            z = add(z, times(POWERS[3*j], d))
        actual = norm(z) == ONE
        predicted = sum(d != ZERO for d in ds) == 1 or (len(set(ds)) == 1 and ds[0] != ZERO)
        check(actual == predicted, "complete 2401-word unit classification")
        if actual:
            unit_words.append(ds)
    check(len(unit_words) == 30, "24 coordinate plus six diagonal unit differences")

    pts, v = points_four(), labels(4)
    edges = unit_edges(pts)
    old = hamming_edges(4)
    fifth = frame_edges(4)
    ports = tuple(v.index((j,)*4) for j in range(3))
    triangle = frozenset(combinations(ports, 2))
    check(len(pts) == len(set(pts)) == 81, "source injectivity")
    check(len(old) == 324 and len(fifth) == 81 and not old & fifth, "declared source edge counts")
    check(edges == old | triangle and len(edges) == 327, "complete induced incidence")
    check(len(fifth-triangle) == 78, "exact geometry-repair deletion count")
    for a, b in combinations(range(81), 2):
        check((norm(difference(pts[a], pts[b])) == ONE) == ((a, b) in old | triangle),
              "all original point pairs")

    base = compile_support(4, frozenset())
    single = compile_support(4, frozenset({(ports[0], ports[1])}))
    geometric = compile_support(4, triangle)
    full_frame = compile_support(4, fifth)
    check(len(base.phases) == 16 and len(geometric.phases) == 10, "six balanced phases lost")
    check(single.phases == geometric.phases == full_frame.phases,
          "identical complete labelled colouring support after one chord")
    safe = colour_safe_edges(single)
    check(fifth <= safe, "all fifth-direction edges colour-redundant")
    check(not fifth <= edges, "colour redundancy does not prove metric incidence")
    removed_edge = min(fifth - triangle)
    four_colours = list(lift(geometric, 0, (1,1,1,1)))
    for u in removed_edge:
        four_colours[u] = 3
    check(all(four_colours[a] != four_colours[b] for a,b in edges)
          and four_colours[removed_edge[0]] == four_colours[removed_edge[1]],
          "fresh fourth colour refutes cross-palette support transfer")
    for packet in (base, single, geometric, full_frame):
        for signs in packet.phases:
            for c in range(3):
                colouring = lift(packet, c, signs)
                check(len(set(colouring)) == 3, "positive proper three-colouring")

    # Boundary marginals retain multiplicity, not just existence.
    check(boundary_counts(base, ports) == Counter({(0,0,0):18, (0,1,2):30}),
          "complete original diagonal-boundary count")
    check(boundary_counts(geometric, ports) == Counter({(0,1,2):30}), "completed boundary count")
    for nports in range(4):
        for boundary in combinations((0, 1, 2, 3, 7, 13, 40, 80), nports):
            check(boundary_counts(geometric, boundary) == boundary_counts(full_frame, boundary),
                  "all selected boundary marginals equal")

    # Independent complete colouring enumeration for low-dimensional faces.
    for r in (1, 2):
        packet = compile_support(r, frozenset())
        actual = set(enumerate_three_colours(3**r, hamming_edges(r)))
        predicted = {lift(packet, c, signs) for signs in packet.phases for c in range(3)}
        check(actual == predicted and len(actual) == 3*2**r, "independent Latin-face control")
    bad_extra = frozenset({(0, 1), (0, 4), (0, 5)})
    impossible = compile_support(2, bad_extra)
    check(not impossible.phases, "empty support is not a positive colouring certificate")
    check(len(colour_safe_edges(impossible)) == 36, "vacuous universal-edge hostile")

    frame_source = make_cayley((3,)*4,
        tuple(tuple(int(i==j) for j in range(4)) for i in range(4)) + ((1,1,1,1),))
    check(isinstance(classify(frame_source), Collapse), "full frame's independent geometric obstruction")

    # Bounded controls for the all-prime trace and exact phase-count formulas.
    prime_phase_counts = []
    for p in (5, 7, 11, 13):
        r = p-1
        for m in range(2, r):
            check(p*m-(p-1) > m*m, "trace/triangle-inequality exclusion at intermediate support")
        allowed = sum(sum(s) % 3 != 0 for s in product((1,-1), repeat=r))
        check(3*allowed == 2**(r+1)-2, "all-prime labelled-colouring formula")
        prime_phase_counts.append((p, r, 3*allowed))

    rejects(compile_support, True, frozenset())
    rejects(compile_support, 2, frozenset({(False, 1)}))
    rejects(compile_support, 2, frozenset({(0, 1.0)}))
    rejects(audit_support, replace(single, phases=base.phases))
    rejects(audit_support, replace(single, phases=((True,1,1,1),)))
    rejects(lift, single, True, single.phases[0])
    rejects(lift, single, 0, (1,1,-1,-1))
    rejects(boundary_counts, single, (0,0))
    rejects(boundary_counts, single, (False,1))
    rejects(unit_edges, (pts[0], pts[0]))
    rejects(unit_edges, ((True,)+(0,)*7, pts[1]))

    print(json.dumps({"checks":count,"field":"Q(zeta15)","vertices":81,
        "base_edges":324,"complete_frame_edges":405,"actual_induced_edges":327,
        "deleted_frame_edges":78,"exhaustive_difference_words":2401,"unit_words":30,
        "base_three_colourings":48,"repaired_three_colourings":30,
        "full_frame_three_colourings":30,"diagonal_ports":ports,
        "surviving_phases":geometric.phases,
        "boundary_before":{"all_equal":18,"all_distinct":30},
        "boundary_after":{"all_distinct":30},
        "four_colour_hostile_removed_edge":removed_edge,
        "prime_phase_controls":prime_phase_counts,
        "scope":"Exact colour-preserving non-Cayley geometric escape; still chromatic three."},indent=2))


if __name__ == '__main__':
    main()
