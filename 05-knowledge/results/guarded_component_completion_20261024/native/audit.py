"""Independent native component review; no producer implementation imports."""
from array import array
from collections import defaultdict
from itertools import combinations, product, permutations
from pathlib import Path
import argparse, csv, gzip, json, sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
CHECKS = 0


def require(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    claim = json.loads((HERE/'receipt.json').read_text())
    supplied = json.loads((RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json').read_text())
    old = json.loads((RESULTS/supplied['old_domain']).read_text())['sources']
    nodes = list(csv.DictReader((RESULTS/'metagraph_order8_20261023/merged_nodes.tsv').open(), delimiter='\t'))
    merged = {int(z[k]): int(z['merged_id']) for z in nodes for k in ('class_a', 'class_b')}
    owner = array('H')
    owner.frombytes(gzip.decompress((RESULTS/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    paths, fixed = {}, {}
    for z in old:
        by_target = {r['target']: r for r in z['exact_all_shared_clocks_target_domain']}
        paths[z['source']] = [tuple(by_target[t]['path']) for t in (1965, 2234)]
        fixed[z['source']] = z['fixed_clocks']
    for z in supplied['helpers']:
        m = z['selected']
        paths[z['source']] = [tuple(m[k]['path']) for k in ('before', 'after')]
        fixed[z['source']] = m['common_clocks']
    sources = sorted(paths)
    local_arcs = {}
    pairs = list(combinations(range(8), 2))
    chords = [(i, j) for i in range(8) for j in range(i+2, 8)]
    for s in sources:
        bits = f"{int(nodes[s]['canon_a']):028b}"
        arcs = set()
        for digit, (u, v) in zip(bits, pairs):
            arcs.add((u, v) if digit == '1' else (v, u))
        local_arcs[s] = arcs
    costs = {s: len({i for i in range(8) if paths[s][0][i] != paths[s][1][i]}) for s in sources}

    def encode(s, p):
        return sum(1 << k for k, (i, j) in enumerate(chords) if (p[i], p[j]) in local_arcs[s])

    # Full record resources and literal kernels checked independently from
    # canonical matrix strings and owner table, not producer helper functions.
    resources = set()
    for s in sources:
        for bit, p in enumerate(paths[s]):
            require(all(arc in local_arcs[s] for arc in zip(p, p[1:])), 'Hamiltonian path')
            require(all(p[t] == v for t, v in fixed[s]), 'unchanged original clock pins')
            m = encode(s, p)
            require(merged[owner[m]] == s, 'actual source class')
            require(merged[owner[m ^ ((1<<21)-1)]] == (1965, 2234)[bit], 'actual target class')
            reflect = sum(((m >> k)&1) << chords.index((7-j, 7-i)) for k, (i, j) in enumerate(chords))
            full = frozenset((m, reflect, m ^ ((1<<21)-1), reflect ^ ((1<<21)-1)))
            require(len(full) == 4 and all(full.isdisjoint(r) for r in resources), 'disjoint four-word resources')
            resources.add(full)
    require(len(resources) == 162, 'full option capacity')
    for ss, kernel in claim['original_special_kernels'].items():
        s = int(ss)
        domain = []
        hp = 0
        for p in permutations(range(8)):
            if not all(e in local_arcs[s] for e in zip(p, p[1:])):
                continue
            hp += 1
            if all(p[t] == v for t, v in fixed[s]):
                target = merged[owner[encode(s, p) ^ ((1<<21)-1)]]
                if target in (1965, 2234):
                    domain.append(p)
        require(set(domain) == set(paths[s]), 'complete pinned binary domain')
        require(hp == kernel['H'], 'independent Hamiltonian path count')

    er = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    X, Y = er['higher_layers'][-1]['lifts'][0], er['radical_lifts_mod4']
    def charge(ids):
        a = sum(X[i] for i in ids) % 128
        b = tuple(sum(Y[j][i] for i in ids) % 2 for j in range(1, 5))
        return a, b
    reservoir = claim['reservoir_sources']
    all_charges = {charge(t) for t in combinations(reservoir, 3)}
    require(len(all_charges) == 2048, 'protected bank onto exact product group')
    reports = []
    for family in claim['families']:
        embeddings = {int(s): v for s, v in family['embeddings'].items()}
        owner_count = len(set(x for frame in embeddings.values() for x in frame))
        require(owner_count == family['source_host_order'], 'ambient owner union')
        for frame in embeddings.values():
            require(len(set(frame)) == 8, 'injective cell frame')
        for s in sources:
            if s in reservoir:
                require(all(set(embeddings[s]).isdisjoint(embeddings[t]) for t in sources if t != s),
                        'literal reservoir shield from every other cell')
        roots = family['root_sources']
        guards = family['guards']
        actual_patterns = []
        for bs in product((0, 1), repeat=len(roots)):
            assignment = dict(zip(roots, bs))
            for g in guards:
                assignment[g['helper']] = assignment[g['root']] ^ g['bit_xor']
            v_controls = sum(assignment[g['helper']] for g in guards)
            require(v_controls * 2 == len(guards), 'balanced margin all root states')
            # Compute each oriented host edge from every occurrence. Sources,
            # target partners and their literal converses must all agree.
            for mode in ('source', 'partner', 'converse_source', 'converse_partner'):
                oriented = {}
                for s in sources:
                    bit = assignment.get(s, 0)
                    p = paths[s][bit]
                    positions = {v: i for i, v in enumerate(p)}
                    for u, v in local_arcs[s]:
                        if 'partner' in mode and abs(positions[u]-positions[v]) != 1:
                            u, v = v, u
                        if mode.startswith('converse'):
                            u, v = v, u
                        a, b = embeddings[s][u], embeddings[s][v]
                        key, sign = tuple(sorted((a, b))), a < b
                        require(key not in oriented or oriented[key] == sign, 'all full source/partner/converse overlaps')
                        oriented[key] = sign
            for g in guards:
                root, h = g['root'], g['helper']
                for rb, hb in product((0, 1), repeat=2):
                    lhs = embeddings[root][paths[root][rb][g['root_clock']]]
                    rhs = embeddings[h][paths[h][hb][g['helper_clock']]]
                    require((lhs == rhs) == (hb == rb ^ g['bit_xor']), 'exact literal shared-owner predicate')
            actual_patterns.append(v_controls)
        maxima = (0, 0)
        for a, b in combinations(claim['protected_core_sources'], 2):
            controls = [g['helper'] for g in guards if g['root'] in (a, b)]
            bound = (8+len(controls), costs[a]+costs[b]+sum(costs[s] for s in controls)+sum(sorted([costs[s] for s in reservoir], reverse=True)[:6]))
            maxima = tuple(max(x, y) for x, y in zip(maxima, bound))
        require(maxima == (family['uniform_changed_slot_bound'], family['uniform_primary_clock_bound']), 'exact uniform bound calculation')
        states = family['returning_route']['states']
        require(states[0]['helper_V'] == states[-1]['helper_V'], 'exact helper return')
        target_charge = charge(states[0]['core_V'] + states[0]['helper_V'])
        for state in states:
            require(charge(state['core_V'] + state['helper_V']) == target_charge, 'actual arrival full charge')
            require(len(state['helper_V']) == family['helper_occupancy'], 'actual arrival helper margin')
        for a, b, event in zip(states, states[1:], family['returning_route']['events']):
            diff = set(a['core_V']+a['helper_V']) ^ set(b['core_V']+b['helper_V'])
            require(sorted(diff) == event['changed_sources'], 'exact route changed records')
            require(sum(costs[s] for s in diff) == event['owner_clock_occurrences'], 'actual route clock sum')
        reports.append({'family': family['name'], 'root_patterns': len(actual_patterns), 'bound': list(maxima)})
    result = {'status': 'PASS', 'independent_of_producer_imports': True,
              'checks': CHECKS, 'families': reports, 'full_resources': len(resources)}
    output = HERE/'audit_receipt.json'
    if args.write:
        output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        require(json.loads(output.read_text()) == result, 'frozen independent receipt')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
