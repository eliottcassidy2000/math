#!/usr/bin/env python3
"""Replay all14 finite-tree edges from their actual already-replayed parents.

Only the independently written elementary auditor is imported. No native
compiler, path finder, tree builder, or inherited graph implementation is used.
"""
from pathlib import Path
import hashlib
import itertools
import json

import audit_native as A

HERE = Path(__file__).resolve().parent
CONTINUATION = HERE.parent / 'continuation'


def strong(rows):
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(row == (1 << len(rows))-1 for row in reach)


def wrapper(target):
    names, rows = target['names'], target['rows']
    at = {name: i for i, name in enumerate(names)}
    arrow = lambda u, v: (rows[at[u]] >> at[v]) & 1
    prefix = target['prefix']
    component = set(names) - {39} - set(prefix)
    A.need((target['root'], target['mark'], target['spare']) == (39, 31, 6)
           and prefix[:2] == [50, 51] and len(set(prefix)) == len(prefix),
           'physical port roles and original prefix names are unchanged')
    A.need(set(target['component_names']) == component
           and len(target['component_names']) == len(component), 'complete component name set')
    A.need({31, 6} <= component, 'retained mark and spare stay in actual component')
    A.need(all(arrow(u, v) for i, u in enumerate(prefix) for v in prefix[i+1:])
           and all(arrow(u, v) for u in prefix for v in component),
           'actual named transitive prefix dominates actual component')
    A.need(arrow(31, 39) and all(arrow(39, v) for v in names if v not in (31, 39)),
           'root39 retains the actual sole in-port31')
    induced = [sum(1 << j for j, v in enumerate(target['component_names']) if arrow(u, v))
               for u in target['component_names']]
    A.need(induced == target['component_rows'], 'exact component named submatrix')
    A.need(strong(rows) and strong(induced), 'both actual target and component are strong')


def run():
    tree = json.loads((CONTINUATION/'named_tree.json').read_text())
    programs = json.loads((CONTINUATION/'native_tree.json').read_text())
    first_targets = json.loads((HERE.parent/'inverse/named_targets.json').read_text())
    source = tree['source']
    A.valid(source['rows'])
    A.need(source['H'] == 170807 and source['branch'] == '', 'actual initial source label')
    A.need(source['names'] == first_targets['source_names']
           and source['rows'] == first_targets['source_rows'], 'exact previously audited source graph')
    A.need(len(source['names']) == 17 and len(set(source['names'])) == 17, '17 original names')
    wanted = {''.join(word) for length in range(1, 4)
              for word in itertools.product('24', repeat=length)}
    targets = {target['branch']: target for target in tree['targets']}
    instructions = {program['branch']: program for program in programs['programs']}
    A.need(len(tree['targets']) == len(programs['programs']) == 14
           and set(targets) == set(instructions) == wanted, 'complete binary tree through depth3')
    receipts = {receipt['target_H']: receipt for receipt in programs['receipts']}
    A.need(len(receipts) == len(programs['receipts']) == 14, 'one inherited receipt per endpoint')
    # This dictionary contains actual replay outputs, not graph representatives
    # looked up by a scalar Hamilton count.
    realized = {'': {'H': source['H'], 'names': source['names'][:], 'rows': source['rows'][:]}}
    records = []
    total_native = total_helpers = total_permanent = 0
    for branch in sorted(wanted, key=lambda word: (len(word), word)):
        parent_branch = branch[:-1]
        actual_parent = realized[parent_branch]
        target, program = targets[branch], instructions[branch]
        names, rows = actual_parent['names'][:], actual_parent['rows'][:]
        A.need(target['parent_H'] == program['parent_H'] == actual_parent['H'],
               'recorded parent label matches actual prior replay output')
        A.need(target['H'] == program['target_H'] == 3*actual_parent['H']+int(branch[-1]),
               'literal ternary offset labels')
        A.need(target['offset'] == int(branch[-1]) and target['depth'] == len(branch),
               'branch depth and last digit')
        A.need(target['names'][:len(names)] == names, 'every immediate-parent name retained in order')
        A.need(target['old_name_identity'] == [[name, name] for name in names],
               'complete physical old-name identity map')
        A.need({39, 31, 6, 50, 51} <= set(names), 'all five physical ports remain present')
        A.valid(target['rows'])
        old_n = len(names)
        restricted_target = [row & ((1 << old_n)-1) for row in target['rows'][:old_n]]
        changed = A.differences(rows, restricted_target, names)
        A.need(changed == target['old_arc_differences'], 'exact immediate-parent old arc differences')
        old_component = set(names) - {39, 50, 51}
        A.need(all(set(pair) <= old_component for pair in changed), 'only old internal B arcs change')
        at = {name: i for i, name in enumerate(target['names'])}
        A.need(target['new_names'] == target['names'][old_n:]
               == [birth['name'] for birth in program['permanent_births']],
               'exact permanent new-name sequence')
        for birth in program['permanent_births']:
            name = birth['name']
            A.need(name not in names, 'permanent birth freshness')
            cut = sum(1 << i for i, old in enumerate(names)
                      if target['rows'][at[name]] >> at[old] & 1)
            A.need(cut == birth['out_mask'], 'permanent birth uses exact target attachment cut')
            rows = A.add_vertex(rows, cut)
            names.append(name)
        A.need(names == target['names'], 'exact permanent endpoint name order')
        start_of_repairs = rows[:]
        helpers = set()
        A.need(len(program['changes']) == len(changed), 'one helper gadget per old arc difference')
        for number, gadget in enumerate(program['changes']):
            A.need(gadget['changed_pair'] == changed[number], 'actual ordered repair pair')
            helper = gadget['helper']
            A.need(helper not in names and helper not in helpers, 'temporary helper freshness')
            helpers.add(helper)
            expanded = A.add_vertex(rows, gadget['helper_out_mask'])
            local_names = names + [helper]
            A.need(len(gadget['paths']) == 2, 'exactly two native path guards')
            for path in gadget['paths']:
                expanded = A.native(expanded, local_names, path)
            deleted = [row & ((1 << len(names))-1) for row in expanded[:-1]]
            A.valid(deleted)
            A.need(A.differences(rows, deleted, names) == [changed[number]],
                   'helper deletion changes exactly the requested old arc')
            A.need(A.sha(deleted) == gadget['after_delete_sha256'],
                   'every recorded completed-gadget matrix hash matches')
            rows = deleted
            A.need(all(set(pair) <= old_component
                       for pair in A.differences(start_of_repairs, rows, names)),
                   'every completed gadget preserves all outside-B and cross-B arcs')
        A.need(rows == target['rows'] and names == target['names'], 'exact actual target matrix')
        wrapper(target)
        receipt = {'target_H': target['H'], 'native_moves': 2*len(changed),
                   'permanent_births': len(program['permanent_births']),
                   'temporary_births': len(changed), 'deleted_old_names': [],
                   'peak_order': len(names)+int(bool(changed))}
        A.need(receipt == receipts[target['H']], 'all program operation-count receipts')
        total_native += receipt['native_moves']
        total_helpers += receipt['temporary_births']
        total_permanent += receipt['permanent_births']
        # The next generation receives these computed rows and physical names.
        realized[branch] = {'H': target['H'], 'names': names[:], 'rows': rows[:]}
        records.append({'branch': branch, 'parent_branch': parent_branch,
                        'parent_H_label': actual_parent['H'], 'target_H_label': target['H'],
                        'old_names_retained': old_n, 'endpoint_order': len(names),
                        'actual_prefix_names': target['prefix'],
                        'permanent_births': receipt['permanent_births'],
                        'old_B_arcs_changed': len(changed),
                        'guarded_native_moves': receipt['native_moves'],
                        'peak_order': receipt['peak_order'],
                        'actual_parent_graph_sha256': A.sha([actual_parent['names'], actual_parent['rows']]),
                        'actual_target_graph_sha256': A.sha([names, rows])})
    A.need(total_native == 944 and total_helpers == 472
           and max(record['peak_order'] for record in records) == 25,
           'complete claimed operation totals')
    graph_keys = ('branch', 'parent_H', 'H', 'depth', 'offset', 'names', 'rows',
                  'component_names', 'component_rows', 'mark', 'root', 'spare',
                  'prefix', 'new_names', 'old_arc_differences', 'old_name_identity')
    payload = {'source': {key: source[key] for key in ('branch', 'H', 'names', 'rows')},
               'targets': [{key: targets[branch][key] for key in graph_keys}
                           for branch in sorted(wanted, key=lambda word: (len(word), word))]}
    return {'status': 'PASS', 'finite_tree_depth': 3, 'edges': len(records),
            'guarded_native_moves': total_native, 'temporary_births_and_deletions': total_helpers,
            'permanent_births_summed_over_tree_edges': total_permanent,
            'peak_order': 25, 'all_immediate_parent_names_retained': True,
            'fixed_physical_ports': [39, 31, 6, 50, 51], 'records': records,
            'stable_named_tree_payload_sha256': A.sha(payload),
            'stable_native_program_payload_sha256': A.sha(programs['programs']),
            'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'elementary_auditor_sha256': hashlib.sha256(Path(A.__file__).read_bytes()).hexdigest(),
            'scope': 'Literal native reachability through depth3 from actual replayed parents. Physical39,31,6,50,51 stay fixed; additional old names may move between the actual prefix and suffix. No H recount or arbitrary-depth closure claim. H labels are independently counted by sibling interface audits.'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
