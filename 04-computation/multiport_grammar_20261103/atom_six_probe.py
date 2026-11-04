#!/usr/bin/env python3
"""Complete strong-six atom bank; optional inherited all-shifts inverse.

No tournament classes or signature tables are read from a data file. Classes
are generated inductively and canonicalized under every permutation. The
independent literal signature path enumerates vertex permutations, not the
subset dynamic program used by the inherited signature routine.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
from math import factorial
from pathlib import Path

PINS = {
    'component_inverse.py': 'a101041d142c2c5821eca7cd91b1dec145fbc71c671f4315e505e538f29a2dc5',
    'shifted_returns.py': '970794e4810090c17c6aa4ad324b7cb1a6ba0bfc656511ed2ff0800ac2867b34',
}
for filename, expected in PINS.items():
    if hashlib.sha256(Path(__file__).with_name(filename).read_bytes()).hexdigest() != expected:
        raise RuntimeError('inherited source pin mismatch: ' + filename)

from component_inverse import Atom, attach, direct, make_bank_five, strong
from shifted_returns import ShiftInverse, odd_part

GATES = 0


def need(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def permutations(n):
    return tuple(itertools.permutations(range(n)))


def relabel(rows, p):
    return tuple(sum(((rows[p[i]] >> p[j]) & 1) << j for j in range(len(rows)))
                 for i in range(len(rows)))


def canonical(rows):
    return min(relabel(rows, p) for p in permutations(len(rows)))


def reach_strong(rows):
    n = len(rows)
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(row == (1 << n) - 1 for row in reach)


def classes_through_six():
    classes = {1: [(0,)]}
    expected = (0, 1, 1, 2, 4, 12, 56)
    expected_strong = (0, 1, 0, 1, 1, 6, 35)
    census = []
    for n in range(1, 7):
        if n > 1:
            candidates = set()
            for rows in classes[n - 1]:
                for bits in range(1 << (n - 1)):
                    # The new vertex has out-neighbours given by bits.
                    augmented = tuple(rows[i] | (((bits >> i & 1) ^ 1) << (n - 1))
                                      for i in range(n - 1)) + (bits,)
                    candidates.add(canonical(augmented))
            classes[n] = sorted(candidates)
        need(len(classes[n]) == expected[n], 'unlabelled census')
        mass = strong_mass = strong_count = 0
        autos = []
        for rows in classes[n]:
            need(all(not (rows[i] >> i & 1) for i in range(n)), 'loop')
            need(all(((rows[i] >> j & 1) + (rows[j] >> i & 1)) == 1
                     for i in range(n) for j in range(i + 1, n)), 'tournament')
            aut = sum(relabel(rows, p) == rows for p in permutations(n))
            need(aut > 0 and factorial(n) % aut == 0, 'automorphism divisor')
            labelled = factorial(n) // aut
            mass += labelled
            is_strong = reach_strong(rows)
            need(is_strong == strong(rows), 'independent SCC disagreement')
            if is_strong:
                strong_count += 1
                strong_mass += labelled
            autos.append({'rows': rows, 'automorphisms': aut, 'strong': is_strong})
        need(mass == 1 << (n * (n - 1) // 2), 'orbit mass misses labelled tournaments')
        need(strong_count == expected_strong[n], 'strong class count')
        census.append({'n': n, 'classes': len(classes[n]), 'strong_classes': strong_count,
                       'labelled_mass': mass, 'strong_labelled_mass': strong_mass,
                       'representatives': autos})
    return classes, census


def literal_tables(rows):
    """Enumerate every actual path, including all proper vertex subsets."""
    n = len(rows)
    size = 1 << n
    paths = [0] * size
    starts = [[0] * n for _ in range(size)]
    ends = [[0] * n for _ in range(size)]
    paths[0] = 1
    for length in range(1, n + 1):
        for path in itertools.permutations(range(n), length):
            if all(rows[a] >> b & 1 for a, b in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                paths[mask] += 1
                starts[mask][path[0]] += 1
                ends[mask][path[-1]] += 1
    return paths, starts, ends


def literal_sigma(tables, ports):
    paths, starts, ends = tables
    full = len(paths) - 1
    n = len(starts[0])
    st = [sum(row[v] for v in range(n) if ports >> v & 1) for row in starts]
    en = [sum(row[v] for v in range(n) if ports >> v & 1) for row in ends]
    return (paths[full], sum(paths[m] * paths[full ^ m] for m in range(full + 1)),
            sum(en[m] * paths[full ^ m] for m in range(full + 1)),
            sum(paths[m] * st[full ^ m] for m in range(full + 1)),
            sum(en[m] * st[full ^ m] for m in range(full + 1)), en[full], st[full])


def build_bank(classes):
    records, bank, seen = [], [], set()
    by_order = {}
    raw_by_order = Counter()
    for n in range(1, 7):
        for index, rows in enumerate(classes[n]):
            if not reach_strong(rows):
                continue
            tables = literal_tables(rows)
            for ports in range(1 << n):
                sigma = literal_sigma(tables, ports)
                need(sigma == direct(rows, ports), 'literal path vs DP signature')
                raw_by_order[n] += 1
                if sigma not in seen:
                    seen.add(sigma)
                    name = f'S{n}_{index}_{ports:x}'
                    bank.append(Atom(name, rows, ports, sigma))
                    records.append({'name': name, 'rows': rows, 'ports': ports, 'sigma': sigma})
        by_order[n] = len(bank)
    inherited = {a.sigma for a in make_bank_five()}
    need({a.sigma for a in bank if len(a.rows) <= 5} == inherited,
         'full <=5 bank differs from inherited all-labelled bank')
    need(raw_by_order[6] == 35 * 64, 'raw strong-six port universe')
    need(by_order[6] == 2120, 'complete <=6 signature quotient')
    need(by_order[6] - by_order[5] == 1950, 'new six-vertex signatures')
    return tuple(bank), records, by_order, dict(raw_by_order)


class BoundedInverse:
    """Same exact left recurrence, grouped only before expanding actual atoms."""
    def __init__(self, bank):
        self.bank = bank
        self.groups = {}
        for k, atom in enumerate(bank):
            self.groups.setdefault(atom.sigma[:3], []).append(k)

    @lru_cache(None)
    def left(self, h, f):
        if h > f or h < 1 or f < 1:
            return ()
        out = []
        for (H, P, F), indices in self.groups.items():
            if (h, f) == (H, F):
                out.extend((k,) for k in indices if self.bank[k].ports)
            if h % H:
                continue
            tail_h = h // H
            numerator = f - F * tail_h
            if numerator <= 0 or numerator % P:
                continue
            tails = self.left(tail_h, numerator // P)
            out.extend((k,) + word for k in indices for word in tails)
        return tuple(out)


def shifted_targets(bank):
    solver = ShiftInverse(bank)
    grouped = BoundedInverse(bank)
    results = []
    for address, intercept in [('24', -35), ('44', -29)]:
        searches, witnesses = [], []
        for (pre_h, pre_odd_p), prefix in sorted(solver.prefixes(intercept, 63).items()):
            target_c, target_R = intercept // pre_h, 63 // pre_odd_p
            limit = solver.f_bound(target_c)
            B = max(solver.M ** 2, abs(target_c))
            f = target_R
            while f <= limit:
                words = tuple(sorted(word for h in solver.inv.h_values(min(f, B))
                                     for word in grouped.left(h, f)))
                if f <= 1008:
                    expected = tuple(w for w in solver.inv.words_right(f)
                                     if solver.inv.sigma(w)[0] <= B)
                    need(words == expected, 'grouped left vs independent right inverse')
                normalized = [w for w in words if bank[w[0]].ports]
                searches.append({'prefix_h': pre_h, 'prefix_odd_p': pre_odd_p,
                                 'normalized_c': target_c, 'f': f,
                                 'proved_f_bound': limit, 'proved_h_bound': B,
                                 'candidates': len(normalized)})
                for word in normalized:
                    sigma = solver.inv.sigma(word)
                    if sigma[0] - sigma[4] != target_c:
                        continue
                    full_word = prefix + word
                    full_sigma = solver.inv.sigma(full_word)
                    need(full_sigma[0] - full_sigma[4] == intercept, 'prefix intercept')
                    need(odd_part(full_sigma[2]) == 63, 'prefix odd marginal')
                    rows, ports = solver.inv.graph(full_word)
                    shift = -(full_sigma[2] & -full_sigma[2]).bit_length() + 1
                    literal = len(rows) <= 16
                    if literal:
                        need(direct(rows, ports) == full_sigma, 'literal suffix signature')
                    for t in range(max(1, 1 - shift), max(1, 1 - shift) + 3):
                        attached = attach(rows, ports, t + shift)
                        if len(attached) <= 16:
                            need(direct(attached, 0)[0] == intercept + 63 * 2 ** t,
                                 'literal shifted attachment')
                    witnesses.append({'word': [bank[k].name for k in full_word],
                                      'rows': rows, 'ports': ports, 'sigma': full_sigma,
                                      'clock_shift': shift, 'literal_checked': literal})
                print(json.dumps({'target': address, 'prefix': [pre_h, pre_odd_p],
                                  'f': f, 'candidates': len(normalized),
                                  'total_hits': len(witnesses)}), flush=True)
                f *= 2
        results.append({'address': address, 'c': intercept, 'R': 63,
                        'searches': searches, 'witnesses': witnesses})
    need(not results[0]['witnesses'], 'target24 repaired boundary changed')
    need(len(results[1]['witnesses']) == 1, 'target44 representative count changed')
    return results


def physical_schedule():
    """Literal old-name-preserving path from the original4 sibling."""
    target_six = (2, 4, 57, 3, 11, 27)
    rows = [sum(1 << j for j in range(i + 1, 7)) for i in range(7)]
    ports = 84
    reversals = [(i, j) for i in range(6) for j in range(i + 1, 6)
                 if not (target_six[i] >> j & 1)]
    need(len(reversals) == 10, 'internal inversion count')
    actions = [('internal', i, j) for i, j in reversals] + [('remove_port', 2)]
    steps, previous = [], {}
    for step in range(12):
        action = None
        if step:
            action = actions[step - 1]
            if action[0] == 'internal':
                _, i, j = action
                need(rows[i] >> j & 1 and not (rows[j] >> i & 1), 'old internal orientation')
                rows[i] ^= 1 << j
                rows[j] ^= 1 << i
            else:
                need(ports >> 2 & 1, 'old port orientation')
                ports ^= 1 << 2
        sigma = direct(tuple(rows), ports)
        need(rows[6] == 0 and all(rows[i] >> 6 & 1 for i in range(6)), 'terminal singleton')
        need(ports >> 6 & 1, 'terminal vertex remains marked')
        clocks = []
        for original_clock in (3, 4, 5):
            clock = original_clock - 2
            graph = attach(tuple(rows), ports, clock)
            root = len(graph) - 1
            if step == 0:
                old_tt5 = tuple(sum(1 << j for j in range(i + 1, 5)) for i in range(5))
                need(graph == attach(old_tt5, 21, original_clock), 'reblocking changes physical graph')
            else:
                before = previous[original_clock]
                differences = [(i, j) for i in range(len(graph)) for j in range(i + 1, len(graph))
                               if ((graph[i] ^ before[i]) >> j & 1)]
                expected = ((clock + action[1], clock + action[2]) if action[0] == 'internal'
                            else (clock + 2, root))
                need(differences == [expected], 'step is not the declared single physical arc flip')
            previous[original_clock] = graph
            need(reach_strong(graph), 'physical intermediate attachment not strong')
            observed = direct(graph, 1 << root)
            need(observed[0] == sigma[0] - sigma[4] + 2 ** clock * sigma[2], 'physical count law')
            need(observed[2] - observed[0] == sigma[4], 'physical named root defect')
            clocks.append({'original_clock': original_clock, 'rows': graph,
                           'H': observed[0], 'F_root': observed[2],
                           'b_root': observed[2] - observed[0]})
        steps.append({'step': step, 'action': action, 'rows': tuple(rows),
                      'ports': ports, 'sigma': sigma, 'clocks': clocks})
    need(steps[0]['sigma'] == (1, 128, 84, 21, 12, 1, 0), 'TT7 entry signature')
    need(steps[-1]['sigma'] == (17, 448, 252, 105, 46, 17, 4), 'physical final signature')
    for before, after in zip(steps[0]['clocks'], steps[-1]['clocks']):
        need(after['H'] == 3 * before['H'] + 4, 'physical +4 return')
    return steps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--shift', action='store_true', help='run complete inherited 24/44 all-shifts inverse')
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    classes, census = classes_through_six()
    bank, records, bank_sizes, raw = build_bank(classes)
    hits = []
    for atom in bank:
        h, _, f, _, j, _, _ = atom.sigma
        if len(atom.rows) == 6 and f and h - j in (-35, -29) and odd_part(f) == 63:
            hits.append({'name': atom.name, 'rows': atom.rows, 'ports': atom.ports,
                         'sigma': atom.sigma, 'clock_shift': -(f & -f).bit_length() + 1})
    need(not hits, 'recorded negative single-six target changed')
    # Positive atoms and hostile controls are reconstructed as actual graphs.
    need(direct((2, 4, 1), 7) == (3, 12, 9, 9, 6, 3, 3), 'C3 all-port control')
    need(direct((0,), 0) == (1, 2, 0, 0, 0, 0, 0), 'unmarked singleton control')
    need(direct((0,), 1) == (1, 2, 1, 1, 0, 1, 1), 'marked singleton control')
    six = [r for r in records if len(r['rows']) == 6]
    need(max(r['sigma'][0] for r in six) == 45, 'strong-six Hamilton maximum')
    need(min(r['sigma'][0] for r in six) == 15, 'strong-six Hamilton minimum')
    need(min(r['sigma'][2] for r in six if r['ports']) == 15, 'strong-six marked marginal minimum')
    receipt = {'scope': 'complete actual strong-six atoms/all port subsets; no word-length cutoff',
               'inherited_source_pins': PINS, 'class_census': census,
               'signature_bank_sizes_through_order': bank_sizes, 'raw_port_assignments': raw,
               'bank': records, 'single_six_target_hits': hits,
               'physical_schedule': physical_schedule(),
               'single_targets': [{'address': '24', 'c': -35, 'odd_f': 63},
                                  {'address': '44', 'c': -29, 'odd_f': 63}],
               'shift_inverse': None}
    print(json.dumps({'bank_sizes': bank_sizes, 'strong6': 35, 'single_target_hits': hits,
                      'gates': GATES}, sort_keys=True), flush=True)
    if args.shift:
        receipt['shift_inverse'] = shifted_targets(bank)
    receipt['gates'] = GATES
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'gates': GATES, 'output_bytes': len(content.encode()),
                      'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
