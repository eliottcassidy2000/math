#!/usr/bin/env python3
"""Authenticate all82 literal records, pins, capacities and physical current map."""
from array import array
from hashlib import sha256
from pathlib import Path
import csv
import gzip
import json
import sys
import native_functions as native

HERE = Path(__file__).resolve().parent
UP = HERE / 'upstream'
CHECKS = 0
U, V = 1965, 2234


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def read(path):
    return json.loads((UP / path).read_text())


def main():
    provenance = json.loads((HERE / 'provenance.json').read_text())
    for record in provenance['upstream_inputs']:
        need(sha256((HERE / record['path']).read_bytes()).hexdigest() == record['sha256'],
             'exact portable upstream input hash')
    supply = read('compatible_charged_helpers_20261023/supply/selected_supply.json')
    old = read(supply['old_domain'])['sources']
    repair = read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
    host = read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
    electric = read('metagraph_order8_electrical_20261023/receipt.json')['order8']
    nodes = list(csv.DictReader((UP / 'metagraph_order8_20261023/merged_nodes.tsv').open(), delimiter='\t'))
    merged = {int(row[key]): int(row['merged_id']) for row in nodes for key in ('class_a', 'class_b')}
    owners = array('H')
    owners.frombytes(gzip.decompress((UP / 'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder != 'little':
        owners.byteswap()
    X, Y = electric['higher_layers'][-1]['lifts'][0], electric['radical_lifts_mod4']
    potentials, moduli = [X] + Y[:5], [256] + [4] * 5
    labels = {s: ((X[s] % 128) << 4) | sum((Y[j][s] % 2) << (j - 1) for j in range(1, 5))
              for s in range(len(nodes))}
    data = {}
    for row in old:
        records = {r['target']: r for r in row['exact_all_shared_clocks_target_domain']}
        data[row['source']] = dict(records=[records[U], records[V]], fixed=row['fixed_clocks'], blocks=[])
    cores = sorted(data)
    for row in supply['helpers']:
        selected = row['selected']
        data[row['source']] = dict(records=[selected['before'], selected['after']],
                                   fixed=selected['common_clocks'], blocks=[])
    records = {r['target']: r for r in repair['records']}
    data[2125] = dict(records=[records[U], records[V]], fixed=[], blocks=repair['retained_blocks'])
    free = sorted(set(data) - set(cores) - {368, 776, 3139})
    frames = {int(s): frame for s, frame in host['source_frames'].items()}
    need(len(data) == 82 and len(cores) == 28 and len(free) == 51, 'unchanged original actual inventory')
    need(host['whole_source_host_order'] == 652, 'same fixed native source host')
    words, prescribed, actual = set(), {}, {}
    costs = {}
    for s, datum in sorted(data.items()):
        matrix, frame = native.matrix(int(nodes[s]['canon_a'])), frames[s]
        costs[s] = sum(a != b for a, b in zip(datum['records'][0]['path'], datum['records'][1]['path']))
        for i, j in native.PAIRS:
            a, b = frame[i], frame[j]
            edge = min(a, b), max(a, b)
            value = matrix[i][j] if a < b else matrix[j][i]
            need(edge not in prescribed or prescribed[edge] == value, 'literal652-owner host restrictions')
            prescribed[edge] = value
        actual[s] = []
        for bit, record in enumerate(datum['records']):
            p, m = record['path'], record['mask']
            need(sorted(p) == list(range(8)) and all(matrix[a][b] for a, b in zip(p, p[1:])), 'literal Hamiltonian source path')
            need(native.word(matrix, p) == m and merged[owners[m]] == s
                 and merged[owners[native.FULL ^ m]] == (U, V)[bit], 'original actual source/target classes')
            need(native.word(native.partner(matrix, p), p) == native.FULL ^ m, 'literal partner square')
            need(native.word([list(row) for row in zip(*matrix)], p[::-1]) == native.reflect(m), 'literal converse square')
            need(all(p[t] == v for t, v in datum['fixed']), 'all inherited indexed owner pins')
            for block in datum['blocks']:
                need(any(p[i:i + len(block)] == block for i in range(9 - len(block))), 'retained2125 consecutive owner blocks')
            need((X[s] - X[(U, V)[bit]]) % 8 == record['phase'], 'inherited selected ordinary scalar phase')
            resource = native.orbit(m)
            need(len(resource) == 4 and not words & resource, 'all164 reserved normalized alternatives have capacity one')
            words |= resource
            actual[s].append(dict(path=p, global_path=[frame[v] for v in p], words=resource))
    need(len(words) == 656, '164 alternatives,656 reserved closed words')
    need(actual[368][0]['global_path'][1] == actual[776][1]['global_path'][4]
         == actual[3139][1]['global_path'][7] == 4, 'immutable actual externalowner4 equality')
    expected = json.loads((HERE / 'inputs.json').read_text())
    need(expected == dict(cores=cores, records={str(s): dict(label=labels[s], cost=costs[s]) for s in data}),
         'all arithmetic inputs derived from authenticated full native records')
    expected_text = [f'{len(cores)} {len(free)}'] + [f'{s} {labels[s]} {costs[s]}' for s in cores + free]
    need((HERE / 'sheet_input.txt').read_text() == '\n'.join(expected_text) + '\n', 'C++ input exactly same labelled inventory')
    # Universal physical-current implication, not selected endpoint testing.
    # Source row changes vanish because2(f_U-f_V)=0; fixed total V count
    # removes target diagonal debt; Q annihilates the remaining target sums.
    for f, modulus in zip(potentials, moduli):
        need(2 * (f[U] - f[V]) % modulus == 0, 'every source-row transition annihilates original dyadic generators')
    for s in data:
        need((Y[0][s] - X[s] - Y[4][s]) % 2 == 0, 'omitted Y0 parity restored by X plus Y4 code identity')
    # Retain the complete native kernel source/partner certificate verbatim;
    # no lift is made into another option or additional resource.
    for s, certificate in host['full_four_cell_kernels'].items():
        need(len(certificate['full_canonical_resource_lifts']) == 2, 'both complete controlled option kernels retained')
        for option, lifts in enumerate(certificate['full_canonical_resource_lifts']):
            need(data[int(s)]['records'][option]['path'] in lifts, 'selected controlled native path occurs in authenticated full kernel')
    def charge(selected):
        q = 0
        for s in selected:
            q = native.plus(q, labels[s])
        return q
    def rows(core, helpers):
        selected = set(core) | set(helpers) | {776, 3139}
        value = {s: [0] * 6 for s in list(data) + [U, V]}
        for s in data:
            target = V if s in selected else U
            for j, f in enumerate(potentials):
                delta = 2 * (f[s] - f[target])
                value[s][j] += delta
                value[target][j] -= delta
        return value
    examples = []
    low = json.loads((HERE / 'low_grade.normal.json').read_text())
    for hostile in low['sharp33_hostiles']:
        states = []
        for state in hostile['states']:
            core = set(state['core']) if hostile['sign'] == 1 else set(cores) - set(state['core'])
            states.append((core, set(state['helpers'])))
        _, i, j = hostile['crossing_witness']
        examples.append(('sharp33', states, (i, j), 33))
    bridges = json.loads((HERE / 'two_charge.normal.json').read_text())
    for bridge in bridges['literal_bridges']:
        examples.append(('two_point_grade' + str(bridge['grade']),
                         [(set(bridge['first_core']), set(bridge['first_helpers'])),
                          (set(bridge['second_core']), set(bridge['second_helpers']))], (0, 1), 27))
    summaries = []
    for name, states, edge, expected_work in examples:
        baseline = rows(*states[0])
        total_charge = charge(states[0][0] | states[0][1] | {776, 3139})
        for core, helpers in states:
            need(len(helpers) == 2 and set(helpers) <= set(free), 'exact helper4 population and fixed tree records')
            need(charge(core | helpers | {776, 3139}) == total_charge, 'exact physical total charge including fixed helper offset')
            current = rows(core, helpers)
            need(all((current[s][j] - baseline[s][j]) % moduli[j] == 0 for s in current for j in range(6)),
                 'full physical source AND target rows, factor-two weight and target diagonal included')
            occupied = set()
            for s in data:
                resource = actual[s][int(s in core | helpers | {776, 3139})]['words']
                need(not occupied & resource, 'every arrival has capacity one at actual normalized resources')
                occupied |= resource
            need(len(occupied) == 328, 'every arrival occupies82 normalized resources/328 closed words')
        a, b = [states[i] for i in edge]
        changed = (a[0] ^ b[0]) | (a[1] ^ b[1])
        work = sum(costs[s] for s in changed)
        need(work == expected_work and len(changed) <= 8, 'literal claimed bridge support and clock ledger')
        summaries.append(dict(name=name, core_grade=len(states[0][0]), full_fibre_states_replayed=len(states),
                              physical_total_charge=total_charge, changed_sources=sorted(changed), work=work,
                              changed_slot_paths={str(s): [actual[s][int(s in c | h | {776,3139})]['global_path']
                                                          for c, h in (a, b)] for s in sorted(changed)}))
    print(json.dumps(dict(checks=CHECKS, source_slots=82, normalized_options=164, occupied_words=328,
                          reserved_words=656, external_owner=4, current_generators_checked=6,
                          universal_current_implication_verified=True, literal_examples=summaries,
                          inherited_full_kernel_certificate_sha256=sha256((UP / 'guarded_component_completion_20261024/native/xor_tree_receipt.json').read_bytes()).hexdigest()),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
