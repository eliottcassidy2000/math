#!/usr/bin/env python3
"""Exact eight-port support sieve; no geometric realization is presumed.

Universe: all 2^19 simple supergraphs of the fixed C5 wedge C4 port graph,
all proper six-colour partitions, and all 19 original nonedge terminal pairs.
The independent path uses the K7 criterion on eight vertices, not colouring.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

BASE = ((0, 1), (1, 2), (2, 3), (3, 4), (0, 4),
        (0, 5), (5, 6), (6, 7), (0, 7))
BASE_SET = {tuple(sorted(e)) for e in BASE}
NONEDGES = tuple(e for e in itertools.combinations(range(8), 2)
                 if e not in BASE_SET)


def partitions():
    result = []

    def visit(row):
        if len(row) == 8:
            result.append(tuple(row))
            return
        v = len(row)
        for color in range(min(6, max(row, default=-1) + 2)):
            if all(row[u] != color for u in range(v)
                   if (u, v) in BASE_SET):
                visit(row + [color])
    visit([])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    pats = partitions()
    assert len(pats) == 570 and len(NONEDGES) == 19
    all_states = (1 << len(pats)) - 1
    separators = [sum(1 << i for i, p in enumerate(pats) if p[a] != p[b])
                  for a, b in NONEDGES]
    all_pairs = (1 << len(NONEDGES)) - 1
    incident = [sum(1 << j for j, e in enumerate(NONEDGES) if v in e)
                for v in range(8)]
    supports = [0] * (all_pairs + 1)
    supports[0] = all_states
    minima = [20] * len(NONEDGES)
    minimal_examples = []
    empty_support = 0
    positive_equalities = 0
    comparison_checks = 0
    for added in range(all_pairs + 1):
        if added:
            bit = added & -added
            supports[added] = supports[added ^ bit] & separators[bit.bit_length() - 1]
        states = supports[added]
        missing = all_pairs ^ added
        # Independent structural characterization: Q+xy contains a K7.
        structural = 0
        omitted_witness = {}
        for v in range(8):
            inside = missing & ~incident[v]
            if inside == 0:
                structural = all_pairs
                break
            if inside.bit_count() == 1:
                structural |= inside
                omitted_witness[inside.bit_length() - 1] = v
        measured = 0
        for j, separating in enumerate(separators):
            comparison_checks += 1
            if states & separating == 0:
                measured |= 1 << j
                minima[j] = min(minima[j], added.bit_count())
                if states:
                    positive_equalities += 1
                    if added.bit_count() == 13:
                        first = (states & -states).bit_length() - 1
                        p = pats[first]
                        assert p[NONEDGES[j][0]] == p[NONEDGES[j][1]]
                        assert len(set(p)) == 6
                        minimal_examples.append({
                            'terminals': list(NONEDGES[j]),
                            'omitted_vertex': omitted_witness[j],
                            'added_pairs': [list(e) for k, e in enumerate(NONEDGES)
                                            if added & (1 << k)],
                            'surviving_partition_count': states.bit_count(),
                            'equal_colouring': list(p),
                        })
        assert measured == structural, (added, measured, structural)
        if not states:
            empty_support += 1
    assert set(minima) == {13}
    assert len(minimal_examples) == 98
    assert len({tuple(map(tuple, x['added_pairs'])) for x in minimal_examples}) == 98
    # Hostile/positive controls: no additions retain all support; K7-e forces
    # equality, while adding that missing edge kills all six-colourings.
    assert supports[0].bit_count() == 570
    assert supports[all_pairs] == 0
    exemplar = minimal_examples[0]
    added_mask = sum(1 << NONEDGES.index(tuple(e)) for e in exemplar['added_pairs'])
    terminal_bit = 1 << NONEDGES.index(tuple(exemplar['terminals']))
    assert supports[added_mask] and supports[added_mask | terminal_bit] == 0
    output = {
        'status': 'FINITE-EXACT exhaustive support census; theorem proved in companion note',
        'universe': {'vertices': 8, 'base_edges': list(map(list, BASE)),
                     'augmentation_graphs': all_pairs + 1, 'terminal_pairs': 19,
                     'base_six_colour_partitions': len(pats)},
        'independent_support_vs_K7_comparisons': comparison_checks,
        'empty_six_colour_support_graphs': empty_support,
        'graph_terminal_equalities_with_nonempty_support': positive_equalities,
        'minimum_added_pairs_by_terminal': [
            {'terminals': list(e), 'minimum': minima[j]} for j, e in enumerate(NONEDGES)],
        'minimum_examples': minimal_examples,
        'minimum_example_support_histogram': {
            str(n): sum(x['surviving_partition_count'] == n for x in minimal_examples)
            for n in sorted({x['surviving_partition_count'] for x in minimal_examples})},
        'controls': {'unaugmented_full_support': True, 'K7_minus_edge_equality': True,
                     'closing_missing_edge_empty_support': True},
        'limitations': ['Projected pairs are not asserted to be physical unit edges.',
                        'A noncolourable quotient does not imply a noncolourable lift.',
                        'Loops, incompatible collision labels and extra anchors require separate analysis.'],
        'script_sha256_lf': hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in output.items() if k not in ('minimum_examples', 'minimum_added_pairs_by_terminal')}, indent=2))


if __name__ == '__main__':
    main()
