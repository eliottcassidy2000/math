"""Exact observer/complement controls. Standard library only; no assert statements.

Universe: every labelled graph/tournament at n <= 7 is canonized by the full
symmetric group; every ordering of every graph isotype is then inspected.
The adjacent-swap theorem is tested on every pair in every graph isotype.
Run from anywhere: python -B check.py [--write]
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path
import sys


CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def setup(n):
    pairs = list(combinations(range(n), 2))
    index = {e: i for i, e in enumerate(pairs)}
    perms = list(permutations(range(n)))
    actions = []
    for p in perms:
        destinations = [index[tuple(sorted((p[u], p[v])))] for u, v in pairs]
        inversion = sum(1 << destinations[k] for k, (u, v) in enumerate(pairs)
                        if p[u] > p[v])
        actions.append((destinations, inversion))
    return pairs, index, perms, actions


def relabel(mask, destinations):
    result = 0
    while mask:
        bit = mask & -mask
        result |= 1 << destinations[bit.bit_length() - 1]
        mask -= bit
    return result


def orbit_table(size, actions, tournament):
    canon = [-1] * size
    reps = []
    autos = {}
    for mask in range(size):
        if canon[mask] >= 0:
            continue
        reps.append(mask)
        stabilizer = 0
        for destinations, inversion in actions:
            image = relabel(mask, destinations) ^ (inversion if tournament else 0)
            canon[image] = mask
            stabilizer += image == mask
        autos[mask] = stabilizer
    return canon, reps, autos


def observer_mask(order, index):
    positions = {v: i for i, v in enumerate(order)}
    return sum(1 << k for (u, v), k in index.items() if positions[u] > positions[v])


def cycle_count(mask, n, index):
    total = 0
    for a, b, c in combinations(range(n), 3):
        x = mask >> index[a, b] & 1
        y = mask >> index[b, c] & 1
        z = mask >> index[a, c] & 1
        total += x == y and x != z
    return total


def edge(mask, u, v, index):
    return mask >> index[tuple(sorted((u, v)))] & 1


def run():
    census = []
    all_incidence = {}
    collision_controls = []
    expected_graphs = [1, 2, 4, 11, 34, 156, 1044]
    expected_tournaments = [1, 1, 2, 4, 12, 56, 456]
    for n in range(1, 8):
        pairs, index, perms, actions = setup(n)
        size = 1 << len(pairs)
        full = size - 1
        gc, greps, ga = orbit_table(size, actions, False)
        tc, treps, ta = orbit_table(size, actions, True)
        check(len(greps) == expected_graphs[n - 1], (n, 'graph census'))
        check(len(treps) == expected_tournaments[n - 1], (n, 'tournament census'))
        invs = [observer_mask(p, index) for p in perms]
        check(len(set(invs)) == len(perms), (n, 'orders distinct'))
        cycles = [cycle_count(t, n, index) for t in range(size)]
        pair_reps = sorted({min(g, gc[full ^ g]) for g in greps})
        tpair = {t: min(t, tc[full ^ t]) for t in treps}
        distributions = {}
        triangle_distributions = {}
        incidence = {}
        independent = []
        width_hist = Counter()
        max_delta = 0
        for g in greps:
            hist = Counter(tc[g ^ inv] for inv in invs)
            cyc_hist = Counter(cycles[g ^ inv] for inv in invs)
            triangle_distributions[g] = tuple(sorted(cyc_hist.items()))
            independent += [g] if len(hist) == 1 else []
            width_hist[max(cyc_hist) - min(cyc_hist)] += 1
            check(all(value % ga[g] == 0 for value in hist.values()),
                  (n, g, 'free ordered graph action'))
            incidence[g] = {t: value // ga[g] for t, value in hist.items()}
            deg = [sum(edge(g, u, v, index) for v in range(n) if v != u)
                   for u in range(n)]
            expected_mean = Fraction(sum(d * (n - 1 - d) for d in deg), 6)
            actual_mean = Fraction(sum(k * v for k, v in cyc_hist.items()), len(perms))
            check(actual_mean == expected_mean, (n, g, 'triangle expectation'))
            for u, v in pairs:
                e = edge(g, u, v, index)
                before, after = [], []
                different = 0
                for w in range(n):
                    if w in (u, v):
                        continue
                    a, b = edge(g, u, w, index), edge(g, v, w, index)
                    different += a != b
                    (before if a != b and a != e else after).append(w)
                p = before + [u, v] + after
                q = before + [v, u] + after
                left = g ^ observer_mask(p, index)
                right = g ^ observer_mask(q, index)
                check(cycles[right] - cycles[left] == different,
                      (n, g, u, v, 'constructive triangle delta'))
                if different == 0:
                    check(tc[left] == tc[right], (n, g, u, v, 'twin safe swap'))
                max_delta = max(max_delta, different)
            merged = Counter()
            for t, value in hist.items():
                merged[tpair[t]] += value
            distributions[g] = tuple(sorted(merged.items()))
        check(set(independent) == {0, full}, (n, 'only empty/complete independent'))
        # An independent direction: fix each tournament and inspect all observers.
        for t in treps:
            hist = Counter(gc[t ^ inv] for inv in invs)
            check(all(value % ta[t] == 0 for value in hist.values()),
                  (n, t, 'free ordered tournament action'))
            for g in greps:
                check(incidence[g].get(t, 0) == hist[g] // ta[t],
                      (n, g, t, 'two-sided incidence'))
        for g in greps:
            cg = gc[full ^ g]
            check(distributions[g] == distributions[cg],
                  (n, g, 'unordered complement-pair distribution'))
            for t, value in incidence[g].items():
                check(incidence[cg].get(tc[full ^ t], 0) == value,
                      (n, g, t, 'complement/converse incidence'))
        classes_by_distribution = defaultdict(list)
        for g in pair_reps:
            classes_by_distribution[distributions[g]].append(g)
        collisions = sorted([values for values in classes_by_distribution.values()
                             if len(values) > 1])
        if collisions:
            collision_controls.append({'n': n, 'complement_pair_reps': collisions,
                                       'edge_lists': [[[list(pairs[k]) for k in range(len(pairs))
                                                        if g >> k & 1] for g in group]
                                                      for group in collisions]})
        census.append({'n': n, 'labelled_graphs': size, 'graph_types': len(greps),
                       'graph_complement_pairs': len(pair_reps),
                       'tournament_types': len(treps),
                       'converse_merged_tournament_types': len(set(tpair.values())),
                       'orders_per_graph_type': len(perms),
                       'order_independent_graph_types': len(independent),
                       'triangle_range_width_histogram': {str(k): v for k, v in sorted(width_hist.items())},
                       'max_adjacent_swap_delta': max_delta,
                       'distinct_converse_merged_distributions': len(classes_by_distribution),
                       'distinct_triangle_count_distributions': len({triangle_distributions[g] for g in pair_reps}),
                       'distribution_collision_groups': len(collisions)})
        if n <= 4:
            all_incidence[str(n)] = {'graph_order': greps, 'tournament_order': treps,
                                    'matrix': [[incidence[g].get(t, 0) for t in treps] for g in greps],
                                    'graph_automorphisms': [ga[g] for g in greps],
                                    'tournament_automorphisms': [ta[t] for t in treps]}
    return {'status': 'FINITE-EXACT', 'universe': 'all labelled simple graphs and tournaments n<=7; all orders of all graph isotypes',
            'checks': CHECKS, 'census': census, 'incidence_small': all_incidence,
            'distribution_collision_controls': collision_controls}


if __name__ == '__main__':
    result = run()
    receipt = Path(__file__).with_name('receipt.json')
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    else:
        check(result == json.loads(receipt.read_text(encoding='utf-8')), 'saved receipt')
    print(json.dumps({'checks': CHECKS, 'census': result['census']}, indent=2))
