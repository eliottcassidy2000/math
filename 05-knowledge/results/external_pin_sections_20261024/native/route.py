#!/usr/bin/env python3
"""Fixed external owner4; debt paid by actual old-core records, not new ports.

Both exact cardinalities remain fixed: old core k and helpers4.  The helper
tree is permanently368U,776V,3139V. High core grades complement CORE bits
only, retaining the opposite charge sign and unchanged helper labels.
"""
from array import array
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse, csv, gzip, json, random, sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
SUPPLY = RESULTS/'compatible_charged_helpers_20261023/supply'
sys.path.insert(0, str(SUPPLY))
import shared_owner as native

COUNT = 0
U, V = 1965, 2234


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    supply_path = SUPPLY/'selected_supply.json'
    supply = json.loads(supply_path.read_text())
    old_path = RESULTS/supply['old_domain']
    old = json.loads(old_path.read_text())['sources']
    repair_path = RESULTS/'guarded_component_completion_20261024/sections/native_2125_receipt.json'
    repair = json.loads(repair_path.read_text())
    host_path = RESULTS/'guarded_component_completion_20261024/native/xor_tree_receipt.json'
    host = json.loads(host_path.read_text())
    nodes_path = RESULTS/'metagraph_order8_20261023/merged_nodes.tsv'
    nodes = list(csv.DictReader(nodes_path.open(), delimiter='\t'))
    merged = {int(z[k]): int(z['merged_id']) for z in nodes for k in ('class_a', 'class_b')}
    owner_path = RESULTS/'metagraph_order8_20261023/owner_map.u16.gz'
    owners = array('H')
    owners.frombytes(gzip.decompress(owner_path.read_bytes()))
    if sys.byteorder != 'little':
        owners.byteswap()
    electrical_path = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    electrical = json.loads(electrical_path.read_text())['order8']
    X, Y = electrical['higher_layers'][-1]['lifts'][0], electrical['radical_lifts_mod4']
    potentials, moduli = [X]+Y[:5], [256]+[4]*5
    labels = {s: ((X[s] % 128) << 4) | sum((Y[j][s] % 2) << (j-1) for j in range(1, 5)) for s in range(len(nodes))}

    def charge(ss):
        result = 0
        for s in ss:
            result = native.plus(result, labels[s])
        return result

    def signed(q, sign):
        return q if sign == 1 else native.minus(0, q)

    data = {}
    for z in old:
        rr = {r['target']: r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']] = dict(records=[rr[U], rr[V]], fixed=z['fixed_clocks'], blocks=[])
    cores = sorted(data)
    for z in supply['helpers']:
        ss = z['selected']
        data[z['source']] = dict(records=[ss['before'], ss['after']], fixed=ss['common_clocks'], blocks=[])
    rr = {r['target']: r for r in repair['records']}
    data[2125] = dict(records=[rr[U], rr[V]], fixed=[], blocks=repair['retained_blocks'])
    helpers = sorted(set(data)-set(cores))
    fixed_v, fixed_u = {776, 3139}, {368}
    free = sorted(set(helpers)-fixed_v-fixed_u)
    costs = {s: sum(a != b for a, b in zip(d['records'][0]['path'], d['records'][1]['path'])) for s, d in data.items()}
    cheap = [s for s in cores if costs[s] == 2]
    outside = sorted(set(cores)-set(cheap))
    need(len(data) == 82 and len(free) == 51 and len(cheap) == 20 and len(outside) == 8,
         'exact original slots, fixed tree controls, and20+8 core partition')
    need(max(costs[s] for s in outside) == 4 and sum(sorted((costs[s] for s in free), reverse=True)[:2]) == 15,
         'actual cost4 outside core and cost15 two-free-helper bounds')
    embeddings = {int(s): em for s, em in host['source_frames'].items()}
    all_words, literal = set(), {}
    prescribed = {}
    for s, d in data.items():
        A = native.matrix(int(nodes[s]['canon_a']))
        em = embeddings[s]
        for i, j in native.PAIRS:
            a, b = em[i], em[j]
            edge = min(a, b), max(a, b)
            value = A[i][j] if a < b else A[j][i]
            need(edge not in prescribed or prescribed[edge] == value, 'literal host source restrictions remain coherent')
            prescribed[edge] = value
        literal[s] = []
        for bit, r in enumerate(d['records']):
            p, m = r['path'], r['mask']
            need(all(A[a][b] for a, b in zip(p, p[1:])), 'actual source path retained')
            need(native.word(A, p) == m and merged[owners[m]] == s and merged[owners[native.FULL ^ m]] == (U, V)[bit],
                 'actual source/target class and normalized observer word')
            need(native.word(native.partner(A, p), p) == native.FULL ^ m, 'literal K partner')
            need(native.word([list(z) for z in zip(*A)], p[::-1]) == native.reflect(m), 'literal converse path and word')
            need(all(p[t] == v for t, v in d['fixed']), 'all old indexed common pins')
            for block in d['blocks']:
                need(any(p[i:i+len(block)] == block for i in range(9-len(block))), 'new2125 moving blocks retained')
            words = native.orbit(m)
            need(len(words) == 4 and not all_words & words, 'whole native alternative reservoir has unit capacity')
            all_words |= words
            need((X[s]-X[(U, V)[bit]]) % 8 == r['phase'], 'chosen original scalar phase')
            literal[s].append(dict(global_path=[em[v] for v in p], closed_words=words))
    need(len(all_words) == 656, 'unchanged164 alternative normalized resources')
    need(literal[368][0]['global_path'][1] == literal[776][1]['global_path'][4] == literal[3139][1]['global_path'][7] == 4,
         'actual immutable external boundary owner FOUR')

    # Literal complete free-pair table, with minimum selected weight per charge.
    pair_list, pair_by_q, pair_best = [], defaultdict(list), {}
    for pair in combinations(free, 2):
        pair_list.append(pair)
        q, weight = charge(pair), sum(costs[s] for s in pair)
        pair_by_q[q].append(pair)
        candidate = weight, pair
        if q not in pair_best or candidate < pair_best[q]:
            pair_best[q] = candidate
    need(len(pair_list) == 1275 and len(pair_by_q) == 972, 'external owner4 free-pair image exactly972')

    # Core grades0..4: exact minimum selected weight and full witness.
    core_best = []
    for grade in range(5):
        table = {}
        for subset in combinations(cores, grade):
            q = charge(subset)
            candidate = sum(costs[s] for s in subset), subset
            if q not in table or candidate < table[q]:
                table[q] = candidate
        core_best.append(table)
    need(len(core_best[4]) == 2048 and max(w for w, ss in core_best[4].values()) == 10,
         'all core grade4 charges have selected clock weight at most10')

    # Complete low-grade JOINT hubs, with helper labels unchanged at both core signs.
    joint, joint_summary = {}, []
    for sign in (1, -1):
        for grade in range(4):
            table = {}
            for cq, (cw, cs) in core_best[grade].items():
                for hq, (hw, hs) in pair_best.items():
                    q = native.plus(signed(cq, sign), hq)
                    candidate = cw+hw, cs, hs
                    if q not in table or candidate < table[q]:
                        table[q] = candidate
            joint[sign, grade] = table
            worst = max(w for w, cs, hs in table.values())
            endpoint_bound = sum(sorted((costs[s] for s in cores), reverse=True)[:grade])+15
            need(len(table) == (972 if grade == 0 else 2048), 'exact joint low-grade coverage for this core sign')
            need(endpoint_bound+worst <= 42, 'every low-grade source and target reaches its true joint hub within42')
            for q, (weight, cs, hs) in table.items():
                need(len(cs) == grade and len(hs) == 2 and
                     native.plus(signed(charge(cs), sign), charge(hs)) == q and
                     sum(costs[s] for s in cs+hs) == weight, 'literal full joint hub witness retains signs and both margins')
            joint_summary.append(dict(core_sign=sign, virtual_grade=grade, covered_charges=len(table),
                maximum_minimum_selected_weight=worst, endpoint_selected_weight_bound=endpoint_bound,
                source_to_hub_clock_bound=endpoint_bound+worst))

    # Graded20-core reservoir: independent subset count DP and full Gray census.
    counts = {(0, 0): 1}
    for s in cheap:
        nxt = dict(counts)
        for (k, q), count in counts.items():
            key = k+1, native.plus(q, labels[s])
            nxt[key] = nxt.get(key, 0)+count
        counts = nxt
    gray_counts, reservoir_sections = Counter(), {}
    last, q, grade = 0, 0, 0
    for i in range(1 << 20):
        mask = i ^ (i >> 1)
        if i:
            bit = (mask ^ last).bit_length()-1
            on = bool(mask >> bit & 1)
            q = native.plus(q, labels[cheap[bit]]) if on else native.minus(q, labels[cheap[bit]])
            grade += 1 if on else -1
        gray_counts[grade, q] += 1
        reservoir_sections.setdefault((grade, q), mask)
        last = mask
    need(dict(gray_counts) == counts, 'independent all2^20 Gray census equals full graded count recurrence')
    grade_coverage = [sum((k, q) in reservoir_sections for q in range(2048)) for k in range(21)]
    need(grade_coverage == [1,20,187,910,1899]+[2048]*11+[1899,910,187,20,1], 'exact complete grade5..15 interval')

    def reservoir(grade, q):
        mask = reservoir_sections[grade, q]
        return {cheap[i] for i in range(20) if mask >> i & 1}

    def compile_route(source_core, source_pair, target_core, target_pair):
        sc, sh, tc, th = map(set, (source_core, source_pair, target_core, target_pair))
        k = len(sc)
        need(len(tc) == k and len(sh) == len(th) == 2, 'input separate integer margins')
        need(charge(sc | sh) == charge(tc | th), 'input original actual core-plus-pair charge')
        sign = 1 if k <= 14 else -1
        grade = min(k, 28-k)
        def virt(core):
            return set(core) if sign == 1 else set(cores)-set(core)
        def physical(core):
            return set(core) if sign == 1 else set(cores)-set(core)
        sv, tv = virt(sc), virt(tc)
        q = native.plus(signed(charge(sv), sign), charge(sh))
        need(q == native.plus(signed(charge(tv), sign), charge(th)), 'CORE-only complement has correct unchanged-helper sign')
        out = [(sv, sh)]
        def emit(vc, hp):
            state = set(vc), set(hp)
            if state != out[-1]:
                out.append(state)
        def required_core(hp):
            return signed(native.minus(q, charge(hp)), sign)
        if grade <= 3:
            _, cs, hs = joint[sign, grade][q]
            emit(cs, hs)
            emit(tv, th)
        elif grade == 4:
            emit(core_best[4][required_core(sh)][1], sh)
            hp = set(sh)
            exchanges = list(zip(sorted(sh-th), sorted(th-sh)))
            for a, b in exchanges:
                hp.remove(a);hp.add(b)
                emit(core_best[4][required_core(hp)][1], hp)
            emit(tv, th)
        else:
            star = max(5, grade-len(outside))
            need(star in (5, 6), 'canonical reservoir grade at most6 after CORE-only complement')
            def refresh(ec, rg, hp):
                debt = native.minus(required_core(hp), charge(ec))
                return set(ec) | reservoir(rg, debt)
            def normalize(start, hp):
                path = [set(start)]
                ec, rg = set(start) & set(outside), len(set(start) & set(cheap))
                if rg < 5:
                    for s in sorted(ec)[:5-rg]:
                        ec.remove(s)
                    rg = 5
                    path.append(refresh(ec, rg, hp))
                while rg > star:
                    s = min(set(outside)-ec)
                    ec.add(s);rg -= 1
                    path.append(refresh(ec, rg, hp))
                return path
            source_norm, target_norm = normalize(sv, sh), normalize(tv, th)
            for cs in source_norm[1:]:
                emit(cs, sh)
            ec, target_ec = source_norm[-1] & set(outside), target_norm[-1] & set(outside)
            for a, b in zip(sorted(ec-target_ec), sorted(target_ec-ec)):
                ec.remove(a);ec.add(b)
                emit(refresh(ec, star, sh), sh)
            hp = set(sh)
            for a, b in zip(sorted(sh-th), sorted(th-sh)):
                hp.remove(a);hp.add(b)
                emit(refresh(ec, star, hp), hp)
            emit(target_norm[-1], th)
            for cs in target_norm[-2::-1]:
                emit(cs, th)
        physical_route = [(physical(vc), hp) for vc, hp in out]
        need(physical_route[0] == (sc, sh) and physical_route[-1] == (tc, th), 'all prescribed full endpoints retained')
        for cc, hh in physical_route:
            need(len(cc) == k and len(hh) == 2 and charge(cc | hh) == charge(sc | sh), 'every actual route keeps BOTH margins and physical charge')
        return physical_route

    def rows(core, pair):
        active = core | pair | fixed_v
        result = {s: [0]*6 for s in list(data)+[U, V]}
        for s in data:
            target = V if s in active else U
            for j, f in enumerate(potentials):
                v = 2*(f[s]-f[target])
                result[s][j] += v
                result[target][j] -= v
        return result

    def replay(route, save=False):
        baseline = rows(*route[0])
        states, events = [], []
        for core, pair in route:
            active = core | pair | fixed_v
            need(not active & fixed_u and len(pair | fixed_v) == 4, 'external controls and exact helper marginFOUR')
            boundary = [literal[368][0]['global_path'][1], literal[776][1]['global_path'][4], literal[3139][1]['global_path'][7]]
            need(boundary == [4, 4, 4], 'external owner4 never released, even transiently at an arrival')
            used = set()
            for s in data:
                words = literal[s][int(s in active)]['closed_words']
                need(not used & words, 'actual full normalized arrival capacities')
                used |= words
            need(len(used) == 328, 'all82 source slots remain occupied')
            current = rows(core, pair)
            need(all((current[s][j]-baseline[s][j]) % moduli[j] == 0 for s in current for j in range(6)),
                 'every physical source/target paired dyadic current, including diagonals')
            if save:
                states.append(dict(core_V=sorted(core), free_helper_V=sorted(pair), helper_V=sorted(pair | fixed_v),
                    external_owner=4, physical_delta_rows={str(s): [current[s][j]-baseline[s][j] for j in range(6)]
                        for s in current if current[s] != baseline[s]}))
        for (ca, ha), (cb, hb) in zip(route, route[1:]):
            changed = sorted((ca ^ cb) | (ha ^ hb))
            cost = sum(costs[s] for s in changed)
            need(len(changed) <= 20 and cost <= 42, 'uniform native20slot42clock bound')
            events.append(dict(changed_sources=changed, support=len(changed), owner_clock_occurrences=cost))
        return states, events

    rng = random.Random(20261024)
    case_rows, count, event_count = [], 0, 0
    for k in range(29):
        local_max_support = local_max_clock = 0
        for repeat in range(16):
            source = set(rng.sample(cores, k))
            pair = set(rng.choice(pair_list))
            q = charge(source | pair)
            for attempt in range(10000):
                target = set(rng.sample(cores, k))
                needed = native.minus(q, charge(target))
                if needed in pair_by_q:
                    target_pair = set(rng.choice(pair_by_q[needed]));break
            else:
                raise RuntimeError('could not generate an authenticated same-charge test endpoint')
            route = compile_route(source, pair, target, target_pair)
            _, events = replay(route)
            local_max_support = max([local_max_support]+[e['support'] for e in events])
            local_max_clock = max([local_max_clock]+[e['owner_clock_occurrences'] for e in events])
            count += 1;event_count += len(events)
        case_rows.append(dict(core_grade=k, cases=16, maximum_observed_support=local_max_support, maximum_observed_clocks=local_max_clock))
    # Every charge gets a genuinely independent accepted endpoint pair at
    # physical grades1 AND27; do not inflate replay with identity routes.
    all_charge_endpoint_counts = []
    for sign in (1, -1):
        first, last, multiplicity = {}, {}, Counter()
        for c in cores:
            for hp in pair_list:
                q = native.plus(signed(labels[c], sign), charge(hp))
                first.setdefault(q, (c, hp));last[q] = (c, hp);multiplicity[q] += 1
        need(len(first) == 2048 and min(multiplicity.values()) >= 2,
             'every charge admits two independently enumerated grade1 joint endpoints')
        for q in range(2048):
            c0, hp0 = first[q];c1, hp1 = last[q]
            source = {c0} if sign == 1 else set(cores)-{c0}
            target = {c1} if sign == 1 else set(cores)-{c1}
            need((source, hp0) != (target, hp1), 'all-charge source and target are distinct full records')
            _, events = replay(compile_route(source, hp0, target, hp1))
            count += 1;event_count += len(events)
        all_charge_endpoint_counts.append(dict(core_sign=sign, full_assignments=sum(multiplicity.values()),
                                               minimum_assignments_per_charge=min(multiplicity.values())))
    # A literal interior route whose helpers alone cannot supply all charges.
    k = 14
    source = set(cores[:k]);source_pair = set(pair_list[0])
    q = charge(source | source_pair)
    target = set(cores[-k:])
    needed = native.minus(q, charge(target))
    if needed not in pair_by_q:
        for candidate in combinations(free, 2):
            q = charge(source | set(candidate))
            needed = native.minus(q, charge(target))
            if needed in pair_by_q:
                source_pair = set(candidate);break
    target_pair = set(pair_by_q[needed][-1])
    sample_route = compile_route(source, source_pair, target, target_pair)
    sample_states, sample_events = replay(sample_route, True)
    need(all(s['external_owner'] == 4 for s in sample_states), 'demonstrated core-funded route never enters former owner3 hub')
    fixed_charge = charge(fixed_v)
    helper_image = {native.plus(fixed_charge, q) for q in pair_by_q}
    first_missing = min(set(range(2048))-helper_image)
    _, witness_core, witness_pair = joint[1, 1][native.minus(first_missing, fixed_charge)]
    need(charge(set(witness_core) | set(witness_pair) | fixed_v) == first_missing and
         first_missing not in helper_image, 'one actual coreV slot supplies a charge absent at frozen coregradezero')
    out = dict(status='FINITE-EXACT fixed externalowner4 routes with core-funded residual sections',
        fixed_U_helper=368, fixed_V_helpers=[776,3139], free_helpers=free, free_pair_count=1275, free_pair_charge_values=972,
        core_sources=cores, cheap_core_reservoir=cheap, outside_core=outside,
        graded20_coverage=grade_coverage, core_grade4_max_minimum_weight=10, low_joint_tables=joint_summary,
        low_joint_witnesses={str(sign)+':'+str(k): {str(q): dict(weight=w, virtual_core=list(cs), free_helpers=list(hs))
            for q, (w, cs, hs) in joint[sign, k].items()} for sign in (1,-1) for k in range(4)},
        route_cases=count, checked_events=event_count, per_core_grade_controls=case_rows,
        all_charge_distinct_endpoint_controls=all_charge_endpoint_counts,
        literal_route=dict(core_grade=14, states=sample_states, events=sample_events),
        changed_slot_bound=20, primary_clock_bound=42, core_margin_preserved=True, helper_margin=4,
        reserved_units=164, occupied_units=82, source_host_order=652,
        all_charge_core_grades=list(range(1,28)), extreme_core_grade_charge_count=972,
        explicit_boundary_repair=dict(total_charge=first_missing, core_grade=1,
            core_V=list(witness_core), free_helper_V=list(witness_pair), helper_V=sorted(set(witness_pair) | fixed_v),
            external_owner=4, absent_at_core_grade_zero=True),
        input_hashes={p.relative_to(RESULTS).as_posix(): native.digest(p) for p in
            (supply_path, old_path, repair_path, host_path, nodes_path, owner_path, electrical_path)}, checks=COUNT,
        scope='Every nonempty fibre of declared82 actual menus with fixedcoregrade, helpergrade4, externalowner4 and original aggregate dyadic charge. No newports, no alteredhelpergrade, no globalhost observer.')
    target_file = HERE/'receipt.json'
    if args.write:
        target_file.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(json.loads(target_file.read_text()) == out, 'frozen complete external-pin native route receipt')
    print(f'PASS {COUNT:,}; '+('wrote receipt' if args.write else 'receipt matches'))
    print('route cases', count, 'events', event_count, 'sample', [(e['support'], e['owner_clock_occurrences']) for e in sample_events])


if __name__ == '__main__':
    main()
