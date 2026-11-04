#!/usr/bin/env python3
"""Reproduce all finite menus and all-tail inverse certificates in COMPONENT.md."""
from pathlib import Path
import hashlib
import itertools
import json
import subprocess
import tempfile

import funnel_inverse as F

HERE = Path(__file__).resolve().parent
ROOT = F.ROOT
SOURCE = ROOT / '05-knowledge/results/triangle_return_170807_20261103/zero_ready.json'
TARGETS = (170807, 512423, 512425)


def cached_data(a, cache={}):
    key = tuple(a)
    if key not in cache:
        cache[key] = F.all_data(key)
    return cache[key]


def controls():
    """Hostile empty boundary, direct assignment, and literal graph checks."""
    assignments = literal = 0
    for z in range(1, 256):
        word = F.tail_bits(z)
        count = 0
        for assignment in itertools.product((0, 1), repeat=len(word)):
            before = [i for i, side in enumerate(assignment) if side == 0]
            after = [i for i, side in enumerate(assignment) if side == 1]
            count += ((not before or word[before[-1]] == 1)
                      and (not after or word[after[0]] == 0))
        F.need(F.before_after(z) == count == F.X.tail_coefficients(word, 1)[2],
               ('independent before/after control', word))
        assignments += 1
    seeds = [(0,), (2, 4, 1), F.R.unrank(4, 1)]
    for a in seeds:
        for root, values in enumerate(cached_data(a)):
            for z in range(16):
                word = () if z == 0 else F.tail_bits(z)
                b, newroot = F.X.funnel(a, {root}, word)
                h, _, _, f = F.c.profile(b, newroot)
                F.need(F.forward(values, word) == (h, f),
                       ('literal funnel law', a, root, word))
                expected_strong = len(a) > 1 or bool(word)
                F.need(F.c.strong(b) == expected_strong,
                       ('singleton empty-tail boundary', a, word))
                if expected_strong:
                    out = F.solve_weighted(values, 9, 11, 9*h + 11*f,
                                           allow_empty=len(a) > 1)
                    F.need(any(hit['tail'] == ''.join(map(str, word))
                               for hit in out['hits']), 'inverse recovers literal word')
                literal += 1
    F.need(F.forward((1, 1, 2, 1, 1), ()) == (1, 1),
           'literal singleton boundary values')
    return {'direct_tail_assignment_words': assignments,
            'literal_funnel_and_inverse_controls': literal,
            'singleton_empty_tail_is_not_strong': True}


def nested_grammar(a, names, root):
    records = []
    while len(a) > 1:
        values = cached_data(a)[root]
        F.need(F.R.unique_start_path(a, root) is not None, 'nested unique start')
        keep = [i for i in range(len(a)) if i != root]
        residual = F.R.induced(a, keep)
        blocks = F.R.ordered_components(residual)
        core = blocks[0]
        F.need(all(len(block) == 1 for block in blocks[1:]), 'singleton suffix')
        core_indices = [keep[v] for v in core]
        tail_indices = [keep[block[0]] for block in blocks[1:]]
        entries = [j for j, v in enumerate(core_indices) if a[root] >> v & 1]
        F.need(len(entries) == 1, 'single actual core entry')
        entry = entries[0]
        word = tuple((a[v] >> root) & 1 for v in tail_indices)
        core_rows = F.R.induced(a, core_indices)
        rebuilt, newroot = F.X.funnel(core_rows, {entry}, word)
        order = core_indices + tail_indices + [root]
        F.need(F.R.induced(a, order) == rebuilt, 'literal named reconstruction')
        records.append({'names': names, 'rows': list(a), 'root': names[root],
                        'fields_h_e_P2_F_S2': list(values),
                        'core_names': [names[v] for v in core_indices],
                        'core_rows': list(core_rows), 'entry': names[core_indices[entry]],
                        'tail_names': [names[v] for v in tail_indices],
                        'tail_word': ''.join(map(str, word)),
                        'reconstruction_name_order': [names[v] for v in order]})
        names = [names[v] for v in core_indices]
        a, root = core_rows, entry
    return records


def nested_inverse(records):
    results = []
    for record in records:
        a = tuple(record['rows'])
        for target in TARGETS:
            tested = clocks = 0
            hits = []
            for root, values in enumerate(cached_data(a)):
                out = F.solve_outer(values, target)
                tested += out['candidate_integers']
                clocks += out['clocks']
                hits.extend(dict(entry=record['names'][root], **hit)
                            for hit in out['hits'])
            results.append({'core_order': len(a), 'target': target,
                            'candidate_integers': tested, 'marked_clocks': clocks,
                            'hits': hits})
    recovered = [result for result in results if result['hits']]
    F.need(len(recovered) == 1 and recovered[0]['core_order'] == 8
           and recovered[0]['target'] == 170807
           and recovered[0]['hits'] == [{'entry': 22, 's': 2, 'tail': '10011',
                                        'h': 2195, 'F': 42153}],
           ('actual source positive recovery', recovered))
    return results


def core_arc_inverse(record):
    a = tuple(record['rows'])
    counts = {target: 0 for target in TARGETS[1:]}
    strong = considered = 0
    hits = []
    for edge in [None] + list(itertools.combinations(range(len(a)), 2)):
        b = list(a)
        if edge is not None:
            i, j = edge
            b[i] ^= 1 << j
            b[j] ^= 1 << i
        considered += 1
        if not F.c.strong(b):
            continue
        strong += 1
        for root, values in enumerate(cached_data(b)):
            for target in counts:
                out = F.solve_outer(values, target)
                counts[target] += out['candidate_integers']
                hits.extend(dict(target=target, arc=edge,
                                 root=record['names'][root], **hit)
                            for hit in out['hits'])
    F.need((considered, strong, counts, hits) ==
           (29, 25, {512423: 1077, 512425: 1078}, []), 'core8 arc menu')
    return {'considered_cores': considered, 'strong_cores': strong,
            'core_names': record['names'], 'target_integer_candidates': counts,
            'hits': hits}


def unit_dial_inverse():
    target = 512426
    bank = [(0, (0,))]
    for n in range(3, 9):
        bank.extend((rank, F.R.unrank(n, rank))
                    for rank in range(F.R.class_counts(n)[n]))
    tested = 0
    hits = []
    for rank, a in bank:
        values = cached_data(a)[0]
        out = F.solve_weighted(values, 9, 11, target, allow_empty=len(a) > 1)
        tested += out['candidate_integers']
        hits.extend(dict(core_order=len(a), rank=rank, **hit) for hit in out['hits'])
    F.need((len(bank), tested, hits) == (233, 8944, []), 'unit dial small core bank')
    return {'objective': '9*h + 11*F = 512426', 'core_classes': len(bank),
            'core_orders': [1, 3, 4, 5, 6, 7, 8],
            'candidate_integers': tested, 'hits': hits}


def one_arc_menu(a, names):
    with tempfile.TemporaryDirectory(prefix='component-menu-') as tmp:
        binary = Path(tmp) / 'menu'
        subprocess.run(['c++', '-std=c++20', '-O2', str(HERE/'menu.cpp'),
                        '-o', str(binary)], check=True)
        raw = subprocess.check_output([str(binary)], text=True,
                                      input=str(len(a))+'\n'+' '.join(map(str,a))+'\n')
    data = json.loads(raw)
    F.need(data == json.loads((HERE/'one_arc.json').read_text()), 'C++ receipt replay')
    F.need((data['component_graphs'], data['strong_component_graphs'], data['hits'])
           == (92, 86, []), 'whole B arc menu')
    distance, i, j, mark, clock, h, f = data['nearest_below'][0]
    probe = list(a)
    probe[i] ^= 1 << j
    probe[j] ^= 1 << i
    direct = F.c.profile(probe, mark)
    F.need((direct[0], direct[3], h + (1 << clock)*f, distance)
           == (h, f, 512225, 198), 'independent nearest witness DP')
    data['names_for_all_indices'] = names
    data['nearest_below_columns'] = ['distance_to_one_target', 'flip_i', 'flip_j',
                                    'mark_index', 'prefix_size', 'h_B', 'F']
    return data


def main():
    source = json.loads(SOURCE.read_text())[1]
    a = tuple(source['suffix_rows'])
    names = source['suffix_names']
    root = names.index(31)
    checks = controls()
    F.need(F.c.profile(a, root) == (2195, 1, 53, 42153), 'actual B profile')
    F.need(F.R.rank_class(a, root) == 8868, 'actual unique start class')
    records = nested_grammar(a, names, root)
    result = {'status': 'PASS', 'controls': checks,
              'source_unique_start_rank': 8868,
              'source_unique_path_names': [names[v] for v in F.R.unique_start_path(a, root)],
              'source_nested_grammar': records,
              'all_mark_nested_inverse': nested_inverse(records),
              'core8_one_arc_inverse': core_arc_inverse(records[1]),
              'unit_dial_small_core_inverse': unit_dial_inverse(),
              'whole_B_one_arc_menu': one_arc_menu(a, names)}
    paths = [SOURCE, HERE/'funnel_inverse.py', HERE/'check.py', HERE/'menu.cpp',
             ROOT/'05-knowledge/results/triangle_return_170807_20261103/search.cpp',
             ROOT/'05-knowledge/results/polarized_path_calculus_20261103/explore.py',
             ROOT/'05-knowledge/results/polarized_path_calculus_20261103/recognize.py',
             ROOT/'05-knowledge/results/orbit_spectrum_session_20261103/coverage/check.py',
             ROOT/'05-knowledge/results/session134_constructive_spectrum_frontier_20261101/two_apex/word_recovery.py']
    result['sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in paths}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
