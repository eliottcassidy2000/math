#!/usr/bin/env python3
"""Actual shared-owner XOR tree; one missing native charge repaired by2125."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, csv, gzip, json, sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
OLD = RESULTS/'compatible_charged_helpers_20261023/supply'
sys.path.insert(0, str(OLD))
import shared_owner as n

COUNT = 0


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    supply_path = OLD/'selected_supply.json'
    supply = json.loads(supply_path.read_text())
    repair_path = HERE.parent/'sections/native_2125_receipt.json'
    repair = json.loads(repair_path.read_text())
    old_path = RESULTS/supply['old_domain']
    old = json.loads(old_path.read_text())['sources']
    nodes_path = RESULTS/'metagraph_order8_20261023/merged_nodes.tsv'
    nodes = list(csv.DictReader(nodes_path.open(), delimiter='\t'))
    merged = {int(z[k]): int(z['merged_id']) for z in nodes for k in ('class_a', 'class_b')}
    owner_path = RESULTS/'metagraph_order8_20261023/owner_map.u16.gz'
    owner = array('H')
    owner.frombytes(gzip.decompress(owner_path.read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    electric_path = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    electric = json.loads(electric_path.read_text())['order8']
    X, Y = electric['higher_layers'][-1]['lifts'][0], electric['radical_lifts_mod4']
    potentials, moduli = [X]+Y[:5], [256]+[4]*5

    def label(s):
        return ((X[s] % 128) << 4) | sum((Y[j][s] % 2) << (j-1) for j in range(1, 5))

    def charge(ss):
        result = 0
        for s in ss:
            result = n.plus(result, label(s))
        return result

    data = {}
    for z in old:
        by = {r['target']: r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']] = dict(records=[by[1965], by[2234]], fixed=z['fixed_clocks'], blocks=[])
    cores = sorted(data)
    helpers = []
    for z in supply['helpers']:
        ss = z['selected']
        helpers.append(z['source'])
        data[z['source']] = dict(records=[ss['before'], ss['after']], fixed=ss['common_clocks'], blocks=[])
    need(repair['source'] == 2125 and 2125 not in data, 'genuinely new helper source')
    by = {r['target']: r for r in repair['records']}
    data[2125] = dict(records=[by[1965], by[2234]], fixed=[], blocks=repair['retained_blocks'])
    A = {s: n.matrix(int(nodes[s]['canon_a'])) for s in data}
    need(int(nodes[2125]['canon_a']) == repair['canon'] and X[2125] == repair['X'], 'native new source and original electrical coordinate')
    embeddings = {368: list(range(8)), 776: [8, 9, 4, 10, 11, 12, 13, 3],
                  3139: [14, 15, 16, 17, 18, 3, 4, 19], 2125: list(range(20, 28))}
    next_owner = 28
    for s in sorted(set(data)-set(embeddings)):
        embeddings[s] = list(range(next_owner, next_owner+8))
        next_owner += 8
    need(len(data) == 82 and next_owner == 652, 'all82 source cells in coherent652-owner host')
    for a, b in combinations(data, 2):
        expected = 2 if {a, b} <= {368, 776, 3139} else 0
        need(len(set(embeddings[a]) & set(embeddings[b])) == expected, 'only declared common pair3/4 overlaps')

    def glue(matrices):
        prescribed = {}
        for s, aa in matrices.items():
            em = embeddings[s]
            for i, j in n.PAIRS:
                u, v = em[i], em[j]
                edge = min(u, v), max(u, v)
                value = aa[i][j] if u < v else aa[j][i]
                need(edge not in prescribed or prescribed[edge] == value, 'every common source/partner arc agrees')
                prescribed[edge] = value
        return prescribed

    source_edges = glue(A)
    all_words, costs, records = set(), {}, {}
    for s, d in sorted(data.items()):
        records[s] = []
        costs[s] = sum(a != b for a, b in zip(d['records'][0]['path'], d['records'][1]['path']))
        for bit, r in enumerate(d['records']):
            p, m = r['path'], r['mask']
            gp = [embeddings[s][v] for v in p]
            need(all(A[s][a][b] for a, b in zip(p, p[1:])), 'actual cell Hamiltonian path')
            need(n.word(A[s], p) == m and merged[owner[m]] == s and merged[owner[n.FULL ^ m]] == (1965, 2234)[bit],
                 'source/partner normalized classes retained')
            need(n.word(n.partner(A[s], p), p) == n.FULL ^ m, 'literal local partner word')
            need(n.word([list(z) for z in zip(*A[s])], p[::-1]) == n.reflect(m), 'literal converse word')
            need(all(p[t] == v for t, v in d['fixed']), 'every indexed common pin')
            for block in d['blocks']:
                need(any(p[i:i+len(block)] == block for i in range(9-len(block))), 'actual moving consecutive owner block retained')
                reverse = p[::-1]
                need(any(reverse[i:i+len(block)] == block[::-1] for i in range(9-len(block))), 'converse moving owner block retained')
            words = n.orbit(m)
            need(len(words) == 4 and not all_words & words, 'all164 native option resources pairwise disjoint')
            all_words |= words
            need((X[s]-X[(1965, 2234)[bit]]) % 8 == r['phase'], 'chosen original per-slot phase')
            records[s].append(dict(path=p, global_path=gp, source_word=m, target=(1965, 2234)[bit],
                normalized_resource=min(words), closed_words=sorted(words)))
    need(len(all_words) == 656 and costs[2125] == 8, 'new eight-clock menu and exact enlarged reservoir')

    # Independent full source/path and complete menu checks for all four cells.
    kernels = {}
    for s in (368, 776, 3139, 2125):
        d, total, raw, selected, lifts = data[s], 0, 0, [], [[], []]
        for p0 in permutations(range(8)):
            p = list(p0)
            if not all(A[s][a][b] for a, b in zip(p, p[1:])):
                continue
            total += 1
            m = n.word(A[s], p)
            for bit, r in enumerate(d['records']):
                if m in n.orbit(r['mask']):
                    lifts[bit].append(p)
            guard = all(p[t] == v for t, v in d['fixed']) and all(
                any(p[i:i+len(block)] == block for i in range(9-len(block))) for block in d['blocks'])
            if guard:
                raw += 1
                t = merged[owner[n.FULL ^ m]]
                if t in (1965, 2234):
                    selected.append((tuple(p), m, t))
        need(set(selected) == {(tuple(r['path']), r['mask'], r['target']) for r in d['records']},
             'complete guarded U/V menu contains exactly the two full supplied records')
        kernels[s] = dict(H=total, raw_guard_count=raw, full_canonical_resource_lifts=lifts)
    need(kernels[2125]['H'] == 401 and kernels[2125]['raw_guard_count'] == 19, 'independent new moving-block source census')

    truth = []
    for center, left, right in product(range(2), repeat=3):
        pp, qq, rr = records[368][center]['global_path'], records[776][left]['global_path'], records[3139][right]['global_path']
        allowed = pp[1] == qq[4] == rr[7]
        need(allowed == (center != left and center != right), 'exact native two-XOR shared-owner predicate')
        need((pp[::-1][6] == qq[::-1][3] == rr[::-1][0]) == allowed, 'actual reflected clock equation')
        glue({s: n.partner(A[s], data[s]['records'][bit]['path']) for s, bit in ((368, center), (776, left), (3139, right))})
        truth.append(dict(bits=[center, left, right], allowed=allowed, actual_owners=[pp[1], qq[4], rr[7]]))

    def permitted(ids):
        return (368 in ids) != (776 in ids) and (368 in ids) != (3139 in ids)

    censuses, sections = [], []
    for pool in (helpers, sorted(helpers+[2125])):
        counts, section = Counter(), {}
        repairs = []
        for ids in combinations(pool, 4):
            if permitted(ids):
                q = charge(ids)
                counts[q] += 1
                section.setdefault(q, ids)
                if q == 1509:
                    repairs.append(ids)
        missing = sorted(set(range(2048))-set(counts))
        need(missing == ([1509] if len(pool) == 53 else []), 'exact native-guard charge obstruction and completion')
        censuses.append(dict(helper_count=len(pool), occupancy=4, allowed_quadruples=sum(counts.values()),
            covered_charges=len(counts), missing_charges=missing, charge1509_witnesses=[list(z) for z in repairs]))
        sections.append(section)
    need(censuses[-1]['charge1509_witnesses'] == [[368, 676, 1770, 2125], [776, 1596, 2125, 3139]],
         'both and only two native repairs of the formerly missing charge')

    # Critical charge admits a real moving-owner excursion with full return.
    helper_route = [{776, 1596, 2125, 3139}, {368, 676, 1770, 2125}, {776, 1596, 2125, 3139}]
    core = {1394, 2185}

    def physical_rows(active_helpers):
        result = {s: [0]*6 for s in list(data)+[1965, 2234]}
        for s in data:
            t = 2234 if s in core | active_helpers else 1965
            for j, f in enumerate(potentials):
                val = 2*(f[s]-f[t])
                result[s][j] += val
                result[t][j] -= val
        return result

    baseline = physical_rows(helper_route[0])
    route = []
    for hs in helper_route:
        need(len(hs) == 4 and charge(hs) == 1509 and permitted(hs), 'actual charge1509 and coupled pin at each arrival')
        current = physical_rows(hs)
        need(all((current[s][j]-baseline[s][j]) % moduli[j] == 0 for s in current for j in range(6)),
             'every original paired generator physical row, including target diagonal')
        used = set()
        for s in data:
            ww = set(records[s][int(s in core | hs)]['closed_words'])
            need(not used & ww, 'arrival full inverse/converse resource injectivity')
            used |= ww
        need(len(used) == 328, '82 occupied normalized units')
        route.append(dict(core_V=sorted(core), helper_V=sorted(hs), shared_owner=records[368][int(368 in hs)]['global_path'][1],
            physical_delta_rows={str(s): [current[s][j]-baseline[s][j] for j in range(6)] for s in current if current[s] != baseline[s]}))
    events = []
    for aa, bb in zip(helper_route, helper_route[1:]):
        changed = sorted(aa ^ bb)
        events.append(dict(changed_helpers=changed, support=len(changed), owner_clock_occurrences=sum(costs[s] for s in changed)))
    need([z['shared_owner'] for z in route] == [4, 3, 4], 'coupled owner really moves and returns')
    free_costs = sorted((costs[s] for s in helpers+[2125] if s not in (368, 776, 3139)), reverse=True)
    branch_bounds = dict(same_center_U=sum(free_costs[:4]), same_center_V=sum(free_costs[:6]),
                         changed_center=sum(costs[s] for s in (368, 776, 3139))+sum(free_costs[:5]))
    need(branch_bounds == dict(same_center_U=29, same_center_V=43, changed_center=53),
         'exact XOR branch occupancy gives sharper helper clock upper bounds')
    upper = sum(sorted((costs[s] for s in cores), reverse=True)[:2])+max(branch_bounds.values())
    need(upper == 61, 'actual uniform core-plus-coupled-helper clock upper bound61')
    out = dict(status='FINITE-EXACT genuine native shared-owner obstruction and new-source repair',
        input_hashes={p.relative_to(RESULTS).as_posix(): n.digest(p)
                      for p in (supply_path, repair_path, old_path, nodes_path, owner_path, electric_path)},
        source_frames={str(s): embeddings[s] for s in sorted(data)},
        actual_four_cell_records={str(s): records[s] for s in (368, 776, 3139, 2125)},
        full_four_cell_kernels={str(s): kernels[s] for s in kernels}, actual_guard_truth_table=truth,
        whole_source_host_order=652, shared_host_order=20, new_source_private_host_order=28,
        census=censuses, expanded_section=[list(sections[-1][q]) for q in range(2048)],
        critical_charge_route=dict(states=route, events=events),
        fixed_helper_occupancy=4, uniform_changed_slot_bound=10, helper_clock_upper_bounds_by_branch=branch_bounds,
        uniform_primary_clock_bound=61,
        reserved_units=164, reserved_closed_words=656, arrival_units=82, arrival_closed_words=328,
        scope='Original order8 full dyadic aggregate bank; local induced source/partner host frames, not one global Hamiltonian observer; unchanged original53 and old28 pins, new2125 relative blocks.', checks=COUNT)
    target = HERE/'xor_tree_receipt.json'
    if args.write:
        target.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(json.loads(target.read_text()) == out, 'frozen native XOR-tree certificate matches')
    print(f'PASS {COUNT:,}; '+('wrote receipt' if args.write else 'receipt matches'))
    print('missing', censuses[0]['missing_charges'], 'repaired', censuses[1]['covered_charges'], 'events', events)


if __name__ == '__main__':
    main()
