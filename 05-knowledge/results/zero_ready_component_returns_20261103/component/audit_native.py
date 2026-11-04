#!/usr/bin/env python3
"""Independent literal replay; imports no compiler or graph helper."""
from pathlib import Path
import hashlib
import itertools
import json

HERE = Path(__file__).resolve().parent
INVERSE = HERE.parent / 'inverse'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def valid(a):
    n = len(a)
    need(all(0 <= row < 1 << n and not (row >> i & 1)
             for i, row in enumerate(a)), 'loop-free bounded rows')
    need(all(((a[i] >> j) & 1) + ((a[j] >> i) & 1) == 1
             for i, j in itertools.combinations(range(n), 2)), 'tournament pairs')


def add_vertex(a, cut):
    n = len(a)
    need(0 <= cut < 1 << n, 'birth cut bounds')
    b = []
    for i in range(n):
        b.append(a[i] + ((1 - ((cut >> i) & 1)) << n))
    b.append(cut)
    valid(b)
    return b


def native(a, names, word):
    need(len(set(names)) == len(names) and sorted(word) == sorted(names),
         'Hamilton word uses every physical name exactly once')
    where = {name: i for i, name in enumerate(names)}
    path = [where[name] for name in word]
    need(all(a[u] >> v & 1 for u, v in zip(path, path[1:])), 'Hamilton guard')
    protected = {tuple(sorted((u, v))) for u, v in zip(path, path[1:])}
    b = [0] * len(a)
    for u in range(len(a)):
        for v in range(u + 1, len(a)):
            orientation = (a[u] >> v) & 1
            if (u, v) not in protected:
                orientation = 1 - orientation
            winner, loser = (u, v) if orientation else (v, u)
            b[winner] += 1 << loser
    valid(b)
    need(all(b[u] >> v & 1 for u, v in zip(path, path[1:])), 'guard survives')
    return b


def differences(a, b, names):
    need(len(a) == len(b) == len(names), 'same named order')
    return [[names[i], names[j]] for i, j in itertools.combinations(range(len(a)), 2)
            if ((a[i] >> j) & 1) != ((b[i] >> j) & 1)]


def sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run():
    control = [6, 4, 0]
    need(native(control, [0, 1, 2], [0, 1, 2]) == [2, 4, 1],
         'TT3 positive off-path control')
    need(native([2, 4, 1], [0, 1, 2], [0, 1, 2]) == control,
         'same-guard involution')
    for hostile in ([2, 1, 0], [0, 1, 1]):
        rejected = False
        try:
            native(control, [0, 1, 2], hostile)
        except RuntimeError:
            rejected = True
        need(rejected, 'illegal directed/duplicate Hamilton guard rejected')
    targets = json.loads((INVERSE/'named_targets.json').read_text())
    native_report = json.loads((INVERSE/'native.json').read_text())
    original_names, original_rows = targets['source_names'], targets['source_rows']
    valid(original_rows)
    need(len(original_names) == 17 and len(set(original_names)) == 17, '17 original names')
    old_B = set(original_names) - {39, 50, 51}
    inherited_path = HERE.parent.parent / 'triangle_return_170807_20261103/zero_ready.json'
    inherited = json.loads(inherited_path.read_text())[1]
    source_at = {name: i for i, name in enumerate(original_names)}
    source_arrow = lambda u, v: (original_rows[source_at[u]] >> source_at[v]) & 1
    need(old_B == set(inherited['suffix_names']), 'actual inherited suffix names')
    inherited_induced = [sum(1 << j for j, v in enumerate(inherited['suffix_names'])
                             if source_arrow(u, v)) for u in inherited['suffix_names']]
    need(inherited_induced == inherited['suffix_rows'], 'actual inherited B adjacency')
    need(source_arrow(50, 51) and all(source_arrow(u, v) for u in (50, 51) for v in old_B)
         and source_arrow(31, 39)
         and all(source_arrow(39, v) for v in original_names if v not in (31, 39)),
         'actual inherited source wrapper')
    endpoint_by_H = {target['H']: target for target in targets['targets']}
    expected = {512423: 26, 512425: 4}
    need(set(endpoint_by_H) == set(expected), 'two requested endpoint payloads')
    need(sorted(program['target_H'] for program in native_report['programs']) == sorted(expected),
         'one program for each requested target')
    records = []
    for program in native_report['programs']:
        target = endpoint_by_H[program['target_H']]
        names, rows = original_names[:], original_rows[:]
        endpoint_names, endpoint_rows = target['names'], target['rows']
        valid(endpoint_rows)
        need(len(endpoint_names) == len(set(endpoint_names)) == 21, 'order21 named target')
        need(endpoint_names[:17] == original_names, 'original physical names retained in order')
        need(target['new_names'] == endpoint_names[17:], 'four new permanent names')
        old_target = [row & ((1 << 17)-1) for row in endpoint_rows[:17]]
        changed = differences(original_rows, old_target, original_names)
        need(changed == target['old_arc_differences'] and len(changed) == expected[target['H']],
             'actual endpoint old arc difference list')
        need(all(set(pair) <= old_B for pair in changed), 'every changed old arc is inside B')
        endpoint_loc = {name: i for i, name in enumerate(endpoint_names)}
        for receipt in program['permanent_births']:
            name = receipt['name']
            need(name not in names, 'fresh permanent name')
            cut = sum(1 << i for i, old in enumerate(names)
                      if endpoint_rows[endpoint_loc[name]] >> endpoint_loc[old] & 1)
            need(cut == receipt['out_mask'], 'permanent birth exact final attachment cut')
            rows = add_vertex(rows, cut)
            names.append(name)
            need(rows == receipt['rows'], 'permanent birth matrix')
        need(names == endpoint_names, 'permanent births yield exact endpoint order')
        before_gadgets = rows[:]
        helpers = set()
        guarded_paths = 0
        for number, gadget in enumerate(program['changes']):
            need(gadget['changed_pair'] == changed[number], 'ordered actual changed pairs')
            helper = gadget['helper']
            need(helper not in names and helper not in helpers, 'fresh temporary helper')
            helpers.add(helper)
            local = names + [helper]
            expanded = add_vertex(rows, gadget['helper_out_mask'])
            need(expanded == gadget['birth_rows'], 'helper attachment cut')
            need(len(gadget['paths']) == len(gadget['native_rows']) == 2, 'two native moves')
            for word, saved in zip(gadget['paths'], gadget['native_rows']):
                expanded = native(expanded, local, word)
                need(expanded == saved, 'independent off-path complement matrix')
                guarded_paths += 1
            deleted = [row & ((1 << len(names))-1) for row in expanded[:-1]]
            valid(deleted)
            need(deleted == gadget['after_delete'], 'delete only last temporary helper')
            need(differences(rows, deleted, names) == [changed[number]],
                 'completed gadget changes exactly the requested old arc')
            rows = deleted
            changed_since_births = differences(before_gadgets, rows, names)
            need(all(set(pair) <= old_B for pair in changed_since_births),
                 'all cross-B arcs and outside-B stars survive every completed gadget')
            need({39, 50, 51, 31, 6} <= set(names), 'all named ports survive')
        need(rows == endpoint_rows == program['final_rows'], 'exact final adjacency matrix')
        need(names == program['final_names'], 'exact final physical names')
        need(program['deleted_old_names'] == [], 'no original names deleted')
        need((len(program['permanent_births']), len(helpers), guarded_paths,
              program['peak_order'], program['native_moves'], program['temporary_births'])
             == (4, len(changed), 2*len(changed), 22, 2*len(changed), len(changed)),
             'literal operation counts and peak order')
        bnames = set(names) - {39, 50, 51}
        need(len(target['component_names']) == 18 and set(target['component_names']) == bnames,
             'exact18 named suffix vertices')
        at = {name: i for i, name in enumerate(names)}
        arrow = lambda u, v: (rows[at[u]] >> at[v]) & 1
        need(arrow(50, 51) and all(arrow(u, v) for u in (50, 51) for v in bnames),
             'outer TT2 attachment')
        need(arrow(31, 39) and all(arrow(39, v) for v in names if v not in (31, 39)),
             'unchanged outer root with sole in-port31')
        induced = [sum(1 << j for j, v in enumerate(target['component_names'])
                       if arrow(u, v)) for u in target['component_names']]
        need(induced == target['component_rows'], 'literal marked component payload')
        records.append({'H_label': target['H'], 'old_names_retained': 17,
                        'permanent_births': 4, 'old_B_arc_changes': len(changed),
                        'guarded_native_moves': guarded_paths,
                        'temporary_births_and_deletions': len(helpers), 'peak_order': 22,
                        'exact_final_graph_sha256': sha([names, rows]),
                        'macro_cross_B_and_named_outside_stars_preserved': True})
    need(len(records) == 2, 'both programs audited once')
    payload = {'source_names': original_names, 'source_rows': original_rows,
               'targets': [{key: target[key] for key in
                            ('H', 'names', 'rows', 'component_names', 'component_rows',
                             'new_names', 'old_arc_differences')}
                           for target in targets['targets']]}
    return {'status': 'PASS', 'positive_and_hostile_native_controls': 4, 'programs': records,
            'stable_named_graph_payload_sha256': sha(payload),
            'native_program_payload_sha256': sha(native_report['programs']),
            'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Literal native reachability audit. H labels are not independently recounted here; the separate interface audit counts the final graphs. Native intermediates need not preserve H or port stars.'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
