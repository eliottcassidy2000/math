#!/usr/bin/env python3
"""Exact source-fixed sections and fixed-core-route optimization on native records."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
COUNT = 0
U, V = 1965, 2234


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


def plus(a, b):
    return ((((a >> 4)+(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def minus(a, b):
    return ((((a >> 4)-(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def digest(p):
    return sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    native_path = HERE.parent/'native/receipt.json'
    supply_path = RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
    old_path = RESULTS/'native_packet_decomposition_20261023/native/domains.json'
    electric_path = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    native, supply, old, electric = [json.loads(p.read_text()) for p in
                                    (native_path, supply_path, old_path, electric_path)]
    family = next(f for f in native['families'] if f['name'] == 'twelve_two_components')
    bank = sorted(native['reservoir_sources'])
    guards = family['guards']
    controls = sorted(g['helper'] for g in guards)
    fixed = family['fixed_U_helpers']
    old = old['sources']
    old_ids = sorted(z['source'] for z in old)
    records = {z['source']: {r['target']: r for r in z['exact_all_shared_clocks_target_domain']} for z in old}
    for z in supply['helpers']:
        records[z['source']] = {U: z['selected']['before'], V: z['selected']['after']}
    need(len(bank) == 40 and len(controls) == 12 and fixed == [2948], 'exact protected40 plus12 controls plus fixed-U')
    need(set(bank).isdisjoint(controls) and set(bank).isdisjoint(fixed) and
         set(controls).isdisjoint(fixed), 'no helper counted in two roles')
    need(set(bank+controls+fixed) == {z['source'] for z in supply['helpers']}, 'all original53 slots stay occupied')
    costs = {s: sum(a != b for a, b in zip(records[s][U]['path'], records[s][V]['path'])) for s in records}
    ee = electric['order8']
    x, ys = ee['higher_layers'][-1]['lifts'][0], ee['radical_lifts_mod4']
    potentials, moduli = [x]+ys[:5], [256]+[4]*5
    labels = {s: ((x[s] % 128) << 4) | sum((ys[j][s] % 2) << (j-1) for j in range(1, 5)) for s in records}
    for z in supply['helpers']:
        need(labels[z['source']] == (z['label'][0] << 4) | sum(z['label'][j+1] << j for j in range(4)),
             'labels rebuilt from ORIGINAL electrical potentials')
    def charge(ss):
        q = 0
        for s in ss:
            q = plus(q, labels[s])
        return q
    def selected(mask):
        return {s for i, s in enumerate(bank) if mask >> i & 1}
    def mask(ss):
        return sum(1 << bank.index(s) for s in ss)
    def clock_distance(a, b):
        return sum(costs[s] for s in a ^ b)
    def control_state(core):
        return {g['helper'] for g in guards if int(g['root'] in core) ^ g['bit_xor']}
    triples, by_charge = [], defaultdict(list)
    for t in combinations(bank, 3):
        entry = (mask(t), frozenset(t))
        triples.append(entry)
        by_charge[charge(t)].append(entry)
    need(len(triples) == 9880 and set(by_charge) == set(range(2048)), 'complete40choose3 native charge universe')

    def source_section(source):
        """Empty-separator min-plus recurrence, retaining exact occupancy <=3.

        Truncation is safe because every remaining slot has nonnegative V count.
        No charge class is dropped, and all slots pay their actual record cost.
        """
        dp = {(0, 0): (0, 0, 1)}
        for i, s in enumerate(bank):
            nxt = {}
            for (k, q), (cost, word, count) in dp.items():
                for bit in (0, 1):
                    if k+bit > 3:
                        continue
                    key = k+bit, plus(q, labels[s]) if bit else q
                    value = cost + costs[s]*(bit != int(s in source))
                    witness = word | ((1 << i) if bit else 0)
                    if key not in nxt or value < nxt[key][0]:
                        nxt[key] = value, witness, count
                    elif value == nxt[key][0]:
                        oldcost, oldword, oldcount = nxt[key]
                        nxt[key] = oldcost, min(oldword, witness), oldcount+count
            dp = nxt
        answer = {q: dp[3, q] for q in range(2048)}
        for q, (cost, word, count) in answer.items():
            # Independent literal triple comparison, not the recurrence.
            literal = [(clock_distance(source, set(t)), m) for m, t in by_charge[q]]
            best = min(v for v, m in literal)
            need(cost == best and word == min(m for v, m in literal if v == best),
                 'source-fixed DP attains exact literal minimum at every charge')
            need(count == sum(v == best for v, m in literal), 'exact number of optimal native sections')
            chosen = selected(word)
            need(len(chosen) == 3 and charge(chosen) == q and clock_distance(source, chosen) == cost,
                 'full-record minimum witness has actual margin, charge and clock work')
        return answer

    core_route = [set(s['core_V']) for s in family['returning_route']['states']]
    initial = set(family['returning_route']['states'][0]['reservoir_V'])
    prescribed_final = set(family['returning_route']['states'][-1]['reservoir_V'])
    total_q = family['returning_route']['total_charge']
    need(initial == prescribed_final == {650, 676, 760}, 'literal prescribed reservoir return')
    need([sorted(c) for c in core_route] == [[1394, 2185], [176, 2185], [1532, 2185],
                                            [176, 1532], [176, 1841], [1477, 1841]], 'inherited fixed core route')
    layer_info = []
    for core in core_route:
        control = control_state(core)
        need(len(control) == 6, 'both actual control components preserve sixV')
        q = minus(minus(total_q, charge(core)), charge(control))
        layer_info.append((core, control, q))
    need(charge(initial) == layer_info[0][2] == layer_info[-1][2], 'return endpoints have same residual debt')

    # A complete layered finite graph, retaining each actual triple as a state.
    layers = [by_charge[q] for _, _, q in layer_info]
    layers[0] = [(mask(initial), frozenset(initial))]
    layers[-1] = [(mask(prescribed_final), frozenset(prescribed_final))]
    def event_bill(j, source_mask, target_state):
        return (clock_distance(layer_info[j-1][0], layer_info[j][0]) +
                clock_distance(layer_info[j-1][1], layer_info[j][1]) +
                clock_distance(selected(source_mask), set(target_state)))
    # First minimize only the bottleneck. A locally better (peak,total) pair
    # need not remain better after a later edge raises both path peaks.
    previous = {mask(initial): 0}
    for j in range(1, len(layers)):
        previous = {to_mask: min(max(peak, event_bill(j, from_mask, to_state))
                                 for from_mask, peak in previous.items())
                    for to_mask, to_state in layers[j]}
    best_peak = previous[mask(prescribed_final)]
    # Then solve the additive shortest-path problem with that peak cap fixed.
    previous = {mask(initial): (0, [mask(initial)])}
    for j in range(1, len(layers)):
        nxt = {}
        for to_mask, to_state in layers[j]:
            candidates = [(total+event_bill(j, from_mask, to_state), path+[to_mask])
                          for from_mask, (total, path) in previous.items()
                          if event_bill(j, from_mask, to_state) <= best_peak]
            if candidates:
                nxt[to_mask] = min(candidates)
        previous = nxt
    best_total, best_path = previous[mask(prescribed_final)]
    optimum = best_peak, best_total, best_path
    # Independent complete path enumeration is cheap on these exact native layers.
    brute_values = []
    for trajectory in product(*layers):
        bills = []
        for j in range(1, len(trajectory)):
            bills.append(clock_distance(layer_info[j-1][0], layer_info[j][0]) +
                         clock_distance(layer_info[j-1][1], layer_info[j][1]) +
                         clock_distance(set(trajectory[j-1][1]), set(trajectory[j][1])))
        brute_values.append((max(bills, default=0), sum(bills), [m for m, t in trajectory]))
    need(optimum == min(brute_values), 'complete literal route enumeration authenticates lexicographic bottleneck optimum')
    greedy_masks, greedy_bills = [mask(initial)], []
    for j in range(1, len(layers)):
        bill, chosen = min((event_bill(j, greedy_masks[-1], t), m) for m, t in layers[j])
        greedy_masks.append(chosen)
        greedy_bills.append(bill)
    need((max(greedy_bills), sum(greedy_bills), greedy_masks) in brute_values,
         'source-greedy comparator is one complete legal fixed-core trajectory')

    def physical_rows(core, helpers):
        active = core | helpers
        result = {s: [0]*6 for s in list(records)+[U, V]}
        for s in records:
            t = V if s in active else U
            for k, p in enumerate(potentials):
                flux = 2*(p[s]-p[t])
                result[s][k] += flux
                result[t][k] -= flux
        return result
    baseline = physical_rows(core_route[0], control_state(core_route[0]) | initial)
    sections, states, events = [], [], []
    for j, word in enumerate(optimum[2]):
        core, control, q = layer_info[j]
        reservoir = selected(word)
        helpers = control | reservoir
        need(len(core) == 2 and len(helpers) == 9 and plus(charge(core), charge(helpers)) == total_q,
             'exact separate core and helper margins and fullQ at every optimized arrival')
        for g in guards:
            rootbit, hbit = int(g['root'] in core), int(g['helper'] in helpers)
            rp, hp = records[g['root']][(U, V)[rootbit]]['path'], records[g['helper']][(U, V)[hbit]]['path']
            re, he = family['embeddings'][str(g['root'])], family['embeddings'][str(g['helper'])]
            need(re[rp[g['root_clock']]] == he[hp[g['helper_clock']]], 'actual fixed-owner shared pin at optimized arrival')
        rows = physical_rows(core, helpers)
        need(all((rows[s][k]-baseline[s][k]) % moduli[k] == 0 for s in rows for k in range(6)),
             'all actual source AND target diagonal current rows retained')
        need(not set(fixed) & helpers, 'fixed helper remains in occupiedU record')
        state = dict(core_V=sorted(core), control_V=sorted(control), reservoir_V=sorted(reservoir),
                     helper_V=sorted(helpers), residual_charge=q)
        states.append(state)
        if j:
            last = states[j-1]
            oldcore, oldhelpers = set(last['core_V']), set(last['helper_V'])
            changed = sorted((core ^ oldcore) | (helpers ^ oldhelpers))
            bill = sum(costs[s] for s in changed)
            need(bill <= optimum[0] and len(changed) <= 20, 'actual optimized event respects exact selected bound')
            events.append(dict(changed_sources=changed, support=len(changed), owner_clock_occurrences=bill,
                               core_cost=clock_distance(core, oldcore),
                               forced_control_cost=clock_distance(control, set(last['control_V'])),
                               reservoir_cost=clock_distance(reservoir, set(last['reservoir_V']))))
        if not any(z['source_reservoir'] == sorted(reservoir) for z in sections):
            answers = source_section(reservoir)
            sections.append(dict(source_reservoir=sorted(reservoir),
                exact_costs=[answers[q][0] for q in range(2048)],
                witness_masks=[answers[q][1] for q in range(2048)],
                number_of_minimizers=[answers[q][2] for q in range(2048)],
                worst_required_charge_cost=max(z[0] for z in answers.values())))
    need(states[0]['helper_V'] == states[-1]['helper_V'], 'all literal initial helper records return')
    need(sum(e['owner_clock_occurrences'] for e in events) == optimum[1] and
         max(e['owner_clock_occurrences'] for e in events) == optimum[0], 'literal path attains optimized peak and total')
    old_events = family['returning_route']['events']
    out = dict(status='FINITE-EXACT source-specific native reset and fixed-core-route optimum',
        reservoir_sources=bank, reservoir_occupancy=3, helper_occupancy=9, total_charge=total_q,
        literal_triples=len(triples), complete_layer_sizes=[len(z) for z in layers],
        complete_route_count=len(brute_values), native_family=family['name'],
        source_specific_sections=sections, optimized_states=states, optimized_events=events,
        minimum_peak_clock=optimum[0], minimum_total_at_that_peak=optimum[1],
        greedy_peak=max(greedy_bills), greedy_total=sum(greedy_bills), greedy_event_costs=greedy_bills,
        greedy_reservoir_states=[sorted(selected(m)) for m in greedy_masks],
        previous_peak=max(e['owner_clock_occurrences'] for e in old_events),
        previous_total=sum(e['owner_clock_occurrences'] for e in old_events),
        uniform_native_upper_bound=family['uniform_primary_clock_bound'],
        source_hashes={p.relative_to(RESULTS).as_posix(): digest(p) for p in
                      (native_path, supply_path, old_path, electric_path, HERE/'check_reset.py')}, checks=COUNT,
        scope='Exact for fixed six-arrival core route and supplied full menus/endpoints; not an optimum over arbitrary core routes or extra native records.')
    target = HERE/'reset_receipt.json'
    if args.write:
        target.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(out == json.loads(target.read_text()), 'frozen source-specific reset certificate')
    print(json.dumps(dict(status='PASS', checks=COUNT, peak=optimum[0], total=optimum[1],
                          previous_peak=out['previous_peak'], previous_total=out['previous_total'],
                          greedy_peak=out['greedy_peak'], greedy_total=out['greedy_total'],
                          layers=out['complete_layer_sizes'], paths=out['complete_route_count'])))


if __name__ == '__main__':
    main()
