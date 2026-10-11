"""Adaptive cycle search using exact signed-boundary fragments as its oracle.

Declared domain: one, two, or three binary channels. No geometric claim is
inferred from an abstract graph, and negative certificates remain complete.
"""

from dataclasses import replace
from itertools import combinations, product

import plane_odd_cycle_channels_20261011b as direct
import plane_six_reframes_parity_20261011b as parity


def split_graph(n, edges):
    """Alternate sorted edges; retain shared endpoints and every isolated vertex."""
    direct.graph(n, edges)
    ordered = tuple(sorted(edges))
    left_edges, right_edges = ordered[::2], ordered[1::2]
    left_vertices = {v for edge in left_edges for v in edge}
    right_vertices = {v for edge in right_edges for v in edge}
    isolates = set(range(n)) - left_vertices - right_vertices
    left_vertices |= isolates
    boundary = tuple(sorted(left_vertices & right_vertices))
    left = parity.Graph(tuple(sorted(left_vertices)), left_edges)
    right = parity.Graph(tuple(sorted(right_vertices)), right_edges)
    parity.audit_graph(left)
    parity.audit_graph(right)
    if set(left.vertices) | set(right.vertices) != set(range(n)):
        raise ValueError('split lost vertices')
    if frozenset(left.edges) | frozenset(right.edges) != edges:
        raise ValueError('split lost edges')
    return left, right, boundary


def split_oracle(n, edges, labels, channels):
    """Root lift API, implemented by boundary equations and witnessed gluing."""
    direct.mono_edges(n, edges, labels, channels)
    if channels > 3:
        raise ValueError('this adapter uses at most three ternary label names')
    left, right, boundary = split_graph(n, edges)
    lp = parity.compile_boundary(left, tuple(labels[v] for v in left.vertices), boundary)
    rp = parity.compile_boundary(right, tuple(labels[v] for v in right.vertices), boundary)
    joined = parity.glue(lp, rp)
    if joined.graph.vertices != tuple(range(n)) or frozenset(joined.graph.edges) != edges:
        raise ValueError('gluing changed the supplied graph')
    if joined.labels != labels:
        raise ValueError('gluing changed supplied labels')
    if joined.odd_cycle:
        # The gluer uses a closed path; SearchCertificate uses no repeated endpoint.
        cycle = direct.chordless(edges, joined.odd_cycle[:-1])
        direct.cycle_checked(n, edges, cycle)
        if len({labels[v] for v in cycle}) != 1:
            raise ValueError('returned cycle is not monochromatic')
        return None, cycle
    colours = joined.colours
    if any(c >= 2 * channels for c in colours):
        raise ValueError('unused channel was introduced by the lift')
    direct.verify(n, edges, direct.SearchCertificate(channels, (), labels, colours))
    return colours, None


def search(n, edges, channels):
    """Same adaptive candidate scheduler; only the rejection/lift oracle changes."""
    direct.graph(n, edges)
    if type(channels) is not int or not 1 <= channels <= 3:
        raise ValueError('exact channel count in 1..3 required')
    cycles = []
    while True:
        candidate = next((q for q in direct.palette.patterns(n, channels)
                          if all(len({q[v] for v in cycle}) > 1 for cycle in cycles)), None)
        if candidate is None:
            cert = direct.SearchCertificate(channels, tuple(cycles), None, None)
            return direct.verify(n, edges, cert)
        colours, cycle = split_oracle(n, edges, candidate, channels)
        if colours is not None:
            cert = direct.SearchCertificate(channels, tuple(cycles), candidate, colours)
            return direct.verify(n, edges, cert)
        if frozenset(cycle) in {frozenset(old) for old in cycles}:
            raise ValueError('repeated learned constraint did not exclude the candidate')
        cycles.append(cycle)


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)
    def rejects(call):
        nonlocal checks
        checks += 1
        try:
            call()
        except ValueError:
            return
        raise RuntimeError('malformed or incomplete certificate accepted')

    graph_count = label_count = searches = positive = negative = 0
    for n in range(1, 5):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            edges = frozenset(edge for i, edge in enumerate(possible) if mask >> i & 1)
            left, right, boundary = split_graph(n, edges)
            check(set(left.vertices) & set(right.vertices) == set(boundary), 'exact intersection')
            check(set(left.edges).isdisjoint(right.edges), 'deterministic edge partition')
            graph_count += 1
            for channels in (1, 2, 3):
                for labels in product(range(channels), repeat=n):
                    got, cycle = split_oracle(n, edges, labels, channels)
                    original, original_cycle = direct.lift(n, edges, labels, channels)
                    check((got is None) == (original is None), 'all-label oracle equivalence')
                    if cycle is not None:
                        direct.cycle_checked(n, edges, cycle)
                        check(len({labels[v] for v in cycle}) == 1, 'actual monochromatic failure')
                    else:
                        check(direct.verify(n, edges, direct.SearchCertificate(channels, (), labels, got)).colours == got,
                              'full authenticated colouring')
                    label_count += 1
                integrated = search(n, edges, channels)
                original = direct.search(n, edges, channels)
                check((integrated.labels is None) == (original.labels is None), 'search decision agreement')
                check(direct.verify(n, edges, integrated) == integrated, 'complete exported certificate')
                positive += integrated.labels is not None
                negative += integrated.labels is None
                searches += 1

    pentagon = direct.unit_pentagon()
    spindle = direct.palette.spindle()
    special = (
        ('unit C5 / one channel', 5, pentagon, 1, False),
        ('unit spindle / two channels', 7, spindle.edges, 2, True),
        ('abstract K7 / three channels', 7, frozenset(combinations(range(7), 2)), 3, False),
    )
    details = []
    for name, n, edges, channels, expected in special:
        integrated = search(n, edges, channels)
        original = direct.search(n, edges, channels)
        check((integrated.labels is not None) == expected, 'declared nontrivial outcome')
        check((integrated.labels is None) == (original.labels is None), 'special oracle agreement')
        check(direct.verify(n, edges, integrated) == integrated, 'special full certificate')
        details.append((name, 'positive' if expected else 'negative', len(integrated.cycles), len(original.cycles)))
        if not expected:
            rejects(lambda n=n, edges=edges, integrated=integrated:
                    direct.verify(n, edges, replace(integrated, cycles=integrated.cycles[:-1])))
    check(details[2][2:] == (35, 35), 'all 35 triples retained for abstract K7')

    # Isolates and empty fragments may not disappear from the retained source.
    for n, edges in ((0, frozenset()), (3, frozenset()), (4, frozenset(((0, 1),)))):
        cert = search(n, edges, 3)
        check(len(cert.colours) == n, 'all isolated vertices survive')
    triangle = frozenset(((0, 1), (0, 2), (1, 2)))
    colours, cycle = split_oracle(3, triangle, (0, 0, 0), 1)
    check(colours is None and set(cycle) == {0, 1, 2}, 'seam-only triangle rejection')
    rejects(lambda: search(True, frozenset(), 3))
    rejects(lambda: search(2, frozenset(((0, 1),)), True))
    rejects(lambda: search(2, frozenset(((0, 1),)), 4))
    rejects(lambda: split_oracle(2, frozenset(((0, 1),)), (False, 0), 1))
    rejects(lambda: split_oracle(2, frozenset(((0, 1),)), (0, 1), 1))
    rejects(lambda: split_graph(2, frozenset(((False, 1),))))
    rejects(lambda: direct.verify(3, triangle, direct.SearchCertificate(3, ((0, 1, 2),), None, None)))

    print('Adaptive odd-cycle search with exact signed-boundary oracle')
    print(f'Graphs orders1..4: {graph_count}; all channel/label oracle controls: {label_count}')
    print(f'Adaptive search comparisons for channels1..3: {searches}; positive={positive}, negative={negative}')
    for name, outcome, integrated_count, direct_count in details:
        print(f'{name}: {outcome}; learned cycles integrated/direct={integrated_count}/{direct_count}')
    print('Negative certificates recheck all channel assignments; deleting a necessary cycle is rejected')
    print('Edges partition deterministically; exactly shared endpoints are ports; isolates retained on the left')
    print(f'Checks: {checks}')
    print('Scope: finite exact oracle integration; no speed claim or plane six/seven decision')


if __name__ == '__main__':
    main()
