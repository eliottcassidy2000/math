"""Width-one symbolic checks and saturation corner cases."""
import itertools
import json
import math
import pathlib
import random
import sys
from audit import graph_step, independent_automaton, sccs, literal_count

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'source'))
import transfer as src


def advance_graph(vector, states, letter, width, cap=None):
    out = {s: 0 for s in states}
    for s, x in zip(states, vector):
        if x:
            for nxt, c in graph_step(s, width, letter).items():
                out[nxt] += x * c
    return tuple(min(x, cap) if cap is not None else x for x in out.values())


def main():
    report = {}
    rng_edges = random.Random(47003)
    graph_checks = 0
    for width in range(1, 5):
        for n in range(1, 31):
            for trial in range(10):
                word = [rng_edges.randrange(1 << width) for j in range(n)]
                graph = [[False] * n for i in range(n)]
                for i in range(n):
                    for j in range(i + 1, n):
                        backward = (j - i <= width) and bool(word[j] & (1 << (width - (j - i))))
                        graph[j][i] = backward
                        graph[i][j] = not backward
                assert graph == src.graph(word, width)
                graph_checks += 1
    report['independent_distance_based_graph_checks'] = graph_checks
    reversals = {(1, 4), (2, 5), (3, 5), (4, 6), (4, 7)}
    graph = [[False] * 7 for i in range(7)]
    for i in range(7):
        for j in range(i + 1, 7):
            graph[j][i] = (i + 1, j + 1) in reversals
            graph[i][j] = not graph[j][i]
    assert literal_count(graph) == 47 and len(sccs(graph)) == 1
    assert graph == src.graph([0, 0, 0, 1, 3, 2, 1], 3)
    report['literal_47_from_explicit_reversed_pairs'] = True
    h = [1, 1, 1]
    for n in range(3, 201):
        h.append(sum(h[-3:]))
    direct = []
    for n in range(1, 9):
        a = src.graph([1] * n, 1)
        val = literal_count(a)
        assert val == h[n]
        direct.append(val)
    report['literal_Rn_counts_n1_to_n8'] = direct
    compositions = []
    for n in range(1, 18):
        total = 0
        for cuts in range(1 << (n - 1)):
            breaks = [-1] + [i for i in range(n - 1) if cuts & (1 << i)] + [n - 1]
            parts = [b - a for a, b in zip(breaks, breaks[1:])]
            total += all(a != 1 or b != 1 for a, b in zip(parts, parts[1:]))
        assert total == h[n]
        compositions.append(total)
    report['literal_composition_counts_n1_to_n17'] = compositions
    checked = 0
    for n in range(1, 13):
        for choices in itertools.product(range(2), repeat=n - 1):
            word = (0,) + choices
            breaks = [0] + [j for j in range(1, n) if word[j] == 0] + [n]
            lengths = [b - a for a, b in zip(breaks, breaks[1:])]
            product = math.prod(h[k] for k in lengths)
            assert src.count(word, 1) == product
            strong = len(sccs(src.graph(word, 1))) == 1
            assert strong == (n == 1 or n >= 3 and all(word[1:]))
            checked += 1
    report['width1_words_full_count_and_strong_check_through_n12'] = checked
    generators = [x for x in h if x <= 101]
    monoid = {1}
    changed = True
    while changed:
        new = {x * y for x in monoid for y in generators if x * y <= 101}
        changed = not new.issubset(monoid)
        monoid |= new
    finite = json.loads((ROOT / 'independent-closure-1-101.json').read_text())
    assert monoid == set(map(int, finite['counts']))
    report['width1_monoid_outputs_through101'] = sorted(monoid)

    rng = random.Random(1970047)
    comparisons = 0
    for width in range(1, 5):
        states, matrices = independent_automaton(width)
        for cap in (1, 2, 3, 4, 8, 48, 102):
            for trial in range(5):
                v = tuple(int(s == ()) for s in states)
                clipped = v
                for step in range(25):
                    letter = rng.randrange(1 << width)
                    v = advance_graph(v, states, letter, width)
                    clipped = advance_graph(clipped, states, letter, width, cap)
                    assert clipped == tuple(min(x, cap) for x in v)
                    exact_h = sum(x for s, x in zip(states, v) if len(s) == 1)
                    clipped_h = sum(x for s, x in zip(states, clipped) if len(s) == 1)
                    assert (exact_h >= cap) == (clipped_h >= cap)
                    if exact_h < cap:
                        assert clipped_h == exact_h
                    comparisons += 1
    report['saturation_full_vector_comparisons'] = comparisons
    # A component with both endpoints forgotten must remain extendible.
    report['fully_forgotten_component_step_width2_mask3'] = {
        str(k): v for k, v in graph_step(((-1, -1),), 2, 3).items()
    }
    assert graph_step(((-1, -1),), 2, 3)[((-1, 1),)] == 1
    (ROOT / 'EDGE_CHECKS.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
