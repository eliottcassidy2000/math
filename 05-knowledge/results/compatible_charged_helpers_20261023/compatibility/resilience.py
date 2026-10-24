"""Exact one-bit compatibility resilience of authenticated fixed-occupancy banks.

Enumerates every witness in the supplied six banks. Hitting numbers are proved
exact only through three; larger values retain a lower bound and a witness
upper bound. No arbitrary helper exclusion is called a native pin operation.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import argparse
import gzip
import io
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUPPLY = HERE.parent/'supply'/'selected_supply.json'
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def bit_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def charge(labels, indices):
    x = 0
    parity = 0
    for i in indices:
        label = labels[i]
        x += label[0]
        parity ^= sum(v << j for j, v in enumerate(label[1:]))
    return ((x % 128) << 4) | parity


def intersect(masks, full):
    value = full
    for mask in masks:
        value &= mask
    return value


def hit_le_three(witnesses, full):
    common = intersect(witnesses, full)
    if common:
        return 1, common & -common
    first = witnesses[0]
    for i in bit_indices(first):
        bit_i = 1 << i
        residual = [w for w in witnesses if not w & bit_i]
        common = intersect(residual, full)
        if common:
            return 2, bit_i | (common & -common)
    for i in bit_indices(first):
        bit_i = 1 << i
        residual = [w for w in witnesses if not w & bit_i]
        for j in bit_indices(residual[0]):
            bit_j = 1 << j
            common = intersect((w for w in residual if not w & bit_j), full)
            if common:
                return 3, bit_i | bit_j | (common & -common)
    return None, None


def greedy_hit(witnesses):
    remaining = list(witnesses)
    hit = 0
    while remaining:
        counts = Counter(i for w in remaining for i in bit_indices(w))
        i = min(counts, key=lambda v: (-counts[v], v))
        hit |= 1 << i
        remaining = [w for w in remaining if not w & (1 << i)]
    return hit


def single_pin_rows(witnesses, sources):
    full = (1 << len(sources))-1
    forced_u_missing = [[] for _ in sources]
    forced_v_missing = [[] for _ in sources]
    for q, ww in enumerate(witnesses):
        common = intersect(ww, full)
        union = 0
        for w in ww:
            union |= w
        for i in bit_indices(common):
            forced_u_missing[i].append(q)
        for i in bit_indices(full ^ union):
            forced_v_missing[i].append(q)
    return [dict(source=s, forced_U_covered=2048-len(forced_u_missing[i]),
                 forced_U_missing=forced_u_missing[i],
                 forced_V_covered=2048-len(forced_v_missing[i]),
                 forced_V_missing=forced_v_missing[i])
            for i, s in enumerate(sources)]


def robust_to_one_u(witnesses, allowed):
    for ww in witnesses:
        remaining = [w for w in ww if not w & ~allowed]
        if not remaining or intersect(remaining, allowed):
            return False
    return True


def compact_single_u_certificate(witnesses, allowed, k):
    result = {}
    for q, ww in enumerate(witnesses):
        live = [w for w in ww if not w & ~allowed]
        chosen = [live[0]]
        common = live[0]
        while common:
            candidate = min(live, key=lambda w: ((w & common).bit_count(), w))
            need((candidate & common) != common, "Single-U certificate intersection cannot shrink")
            chosen.append(candidate)
            common &= candidate
        need(len(chosen) <= k+1, "Too many witnesses for compact single-U section")
        need(intersect(chosen, allowed) == 0, "Compact single-U certificate has common helper")
        result[str(q)] = chosen
    return result


def replay_all_single_u_sections(sections, allowed, k, sources, labels):
    rows = []
    for forbidden in [None] + list(bit_indices(allowed)):
        forbidden_bit = 0 if forbidden is None else 1 << forbidden
        chosen = []
        for q in range(2048):
            mask = next((w for w in sections[str(q)] if not w & forbidden_bit), None)
            need(mask is not None, "Missing compact single-erasure section")
            need(mask.bit_count() == k and not mask & ~allowed and not mask & forbidden_bit,
                 "Erasure witness violates actual occupancy or forced record")
            need(charge(labels, list(bit_indices(mask))) == q,
                 "Erasure witness failed independent literal charge replay")
            chosen.append(mask)
        rows.append(dict(forced_U_source=None if forbidden is None else sources[forbidden],
                         covered_charges=len(chosen),
                         section_sha256=sha256(json.dumps(chosen, separators=(',', ':')).encode()).hexdigest()))
    return rows


def bounded_hitting_controls():
    cases = 0
    for n in range(1, 6):
        full = (1 << n)-1
        for k in range(1, n+1):
            edges = [sum(1 << i for i in c) for c in combinations(range(n), k)]
            for family in range(1, 1 << len(edges)):
                ww = [w for i, w in enumerate(edges) if family & (1 << i)]
                size, mask = hit_le_three(ww, full)
                literal = next((s for s in range(1, 4)
                                if any(all(sum(1 << i for i in c) & w for w in ww)
                                       for c in combinations(range(n), s))), None)
                need(size == literal, "Bounded transversal algorithm failed exhaustive control")
                if mask is not None:
                    need(mask.bit_count() == size and all(mask & w for w in ww),
                         "Bounded transversal control has invalid witness")
                cases += 1
    return cases


def analyse_bank(bank, helper_map):
    sources = bank['source_ids']
    n = len(sources)
    k = bank['V_occupancy']
    full = (1 << n)-1
    labels = [helper_map[s]['label'] for s in sources]
    witnesses = [[] for _ in range(2048)]
    for indices in combinations(range(n), k):
        q = charge(labels, indices)
        witnesses[q].append(sum(1 << i for i in indices))
    need(all(witnesses), "Inherited bank fails complete fixed-occupancy coverage")
    section = bank['charge_section']
    for q, ww in enumerate(witnesses):
        need(section[str(q)] in ww, "Inherited charge section not a literal current witness")
    pins = single_pin_rows(witnesses, sources)
    # Independent literal all-witness scan for every single forced choice.
    for i, row in enumerate(pins):
        bit = 1 << i
        u = [q for q, ww in enumerate(witnesses) if not any(not w & bit for w in ww)]
        v = [q for q, ww in enumerate(witnesses) if not any(w & bit for w in ww)]
        need(row['forced_U_missing'] == u and row['forced_V_missing'] == v,
             "Intersection/union single-pin gate disagrees with literal scan")
    hist = Counter()
    hits = []
    for q, ww in enumerate(witnesses):
        size, mask = hit_le_three(ww, full)
        if size is None:
            mask = greedy_hit(ww)
            hist['at_least_4'] += 1
        else:
            need(mask.bit_count() == size, "Bounded hitting witness has wrong size")
            hist[str(size)] += 1
        need(all(mask & w for w in ww), "Saved exclusion mask does not hit every witness")
        hits.append(dict(charge=q, exact_minimum=size,
                         lower_bound=4 if size is None else size,
                         upper_bound=mask.bit_count(), exclusion_mask=mask,
                         excluded_sources=[sources[i] for i in bit_indices(mask)]))
    observed = min(h['lower_bound'] for h in hits)
    sharp = next(h for h in hits if h['exact_minimum'] == observed) if observed <= 3 else None
    out = dict(helper_count=n, V_occupancy=k, source_ids=sources,
               subset_count=sum(map(len, witnesses)), minimum_witness_count=min(map(len, witnesses)),
               maximum_witness_count=max(map(len, witnesses)),
               single_pin_rows=pins,
               every_single_U_survives=all(row['forced_U_covered'] == 2048 for row in pins),
               every_single_V_survives=all(row['forced_V_covered'] == 2048 for row in pins),
               universal_forced_U_tolerance_at_least=observed-1,
               universal_forced_U_tolerance_exact=observed-1 if sharp else None,
               global_minimum_adversary=sharp,
               bounded_hitting_histogram=dict(sorted(hist.items())),
               per_charge_hitting_certificates=hits)
    # A deterministic deletion pass searches only within the authenticated bank.
    # It produces an inclusion-minimal robust supply, not a cardinality minimum.
    if out['every_single_U_survives']:
        allowed = full
        deleted = []
        for i in reversed(range(n)):
            candidate = allowed ^ (1 << i)
            if robust_to_one_u(witnesses, candidate):
                allowed = candidate
                deleted.append(sources[i])
        need(robust_to_one_u(witnesses, allowed), "Reduced robust supply fails")
        # Monotonicity means a rejected deletion cannot become valid after
        # further removals; nevertheless verify inclusion-minimality literally.
        for i in bit_indices(allowed):
            need(not robust_to_one_u(witnesses, allowed ^ (1 << i)),
                 "Greedy robust supply is not inclusion-minimal")
        sections = compact_single_u_certificate(witnesses, allowed, k)
        gates = replay_all_single_u_sections(sections, allowed, k, sources, labels)
        out['single_U_robust_subbank'] = dict(
            helper_count=allowed.bit_count(), source_ids=[sources[i] for i in bit_indices(allowed)],
            original_index_mask=allowed, removed_sources=deleted,
            status='inclusion-minimal under deletion, not a global cardinality minimum',
            compact_charge_sections=sections, all_single_U_gate_replays=gates,
            all_gates_covered=all(g['covered_charges'] == 2048 for g in gates),
            inherited_two_erasure_adversary=sharp if sharp and not sharp['exclusion_mask'] & ~allowed else None,
            max_changed_helper_clock_cost=sum(sorted(
                (helper_map[s]['selected']['cost'] for i, s in enumerate(sources) if allowed & (1 << i)),
                reverse=True)[:2*k]),
            controller_changed_slot_bound=2+2*k)
    return out


def four_occupancy_extension(sources, helper_map, core_pair_cost, triple_adversary):
    n = len(sources)
    labels = [helper_map[s]['label'] for s in sources]
    full = (1 << n)-1
    witnesses = [[] for _ in range(2048)]
    for indices in combinations(range(n), 4):
        witnesses[charge(labels, indices)].append(sum(1 << i for i in indices))
    need(all(witnesses), 'Four-occupancy extension misses a charge')
    pins = single_pin_rows(witnesses, sources)
    for i, row in enumerate(pins):
        bit = 1 << i
        literal_u = sum(any(not w & bit for w in ww) for ww in witnesses)
        literal_v = sum(any(w & bit for w in ww) for ww in witnesses)
        need(literal_u == row['forced_U_covered'] and literal_v == row['forced_V_covered'],
             'Four-occupancy single-bit gate failed literal replay')
    hist = Counter()
    for q, ww in enumerate(witnesses):
        size, mask = hit_le_three(ww, full)
        hist[str(size) if size else 'at_least_4'] += 1
        if mask is not None:
            need(mask.bit_count() == size and all(mask & w for w in ww),
                 'Four-occupancy bounded exclusion certificate failed')
    # A two-erasure triple obstruction gives a mixed U/V quadruple obstruction:
    # require j, forbid i; the other three must avoid BOTH original sources.
    i_source, j_source = triple_adversary['excluded_sources']
    i, j = sources.index(i_source), sources.index(j_source)
    q3 = triple_adversary['charge']
    j_charge = charge(labels, [j])
    q4 = ((((q3 >> 4)+(j_charge >> 4)) % 128) << 4) | ((q3 ^ j_charge) & 15)
    need(not any(w & (1 << j) and not w & (1 << i) for w in witnesses[q4]),
         'Translated mixed forced-bit obstruction unexpectedly repaired')
    forced_v_pair = (1 << 0) | (1 << 1)
    vv_missing = [q for q, ww in enumerate(witnesses)
                  if not any(w & forced_v_pair == forced_v_pair for w in ww)]
    need(bool(vv_missing), 'Two forced-V counting obstruction disappeared')
    return dict(helper_count=n, V_occupancy=4, source_ids=sources,
                subset_count=sum(map(len, witnesses)),
                minimum_witness_count=min(map(len, witnesses)),
                single_pin_rows=pins,
                every_single_U_survives=all(r['forced_U_covered'] == 2048 for r in pins),
                every_single_V_survives=all(r['forced_V_covered'] == 2048 for r in pins),
                bounded_forced_U_hitting_histogram=dict(hist),
                mixed_UV_adversary=dict(forced_U_source=i_source, forced_V_source=j_source,
                                        missing_charge=q4, inherited_triple_charge=q3),
                two_V_adversary=dict(forced_V_sources=sources[:2], first_missing_charge=vv_missing[0],
                                     covered_charges=2048-len(vv_missing)),
                controller_changed_slot_bound=10,
                controller_owner_clock_bound=core_pair_cost+sum(sorted(
                    (helper_map[s]['selected']['cost'] for s in sources), reverse=True)[:8]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    cert = json.loads(SUPPLY.read_text())
    hitting_controls = bounded_hitting_controls()
    helper_map = {h['source']: h for h in cert['helpers']}
    resources = []
    for h in helper_map.values():
        selected = h['selected']
        need(selected['before']['target'] == 1965 and selected['after']['target'] == 2234,
             'Actual helper target pair changed')
        need(selected['before']['phase'] == selected['after']['phase'],
             'Actual chosen per-slot phase changed')
        need(selected['cost'] == sum(x != y for x, y in zip(
            selected['before']['path'], selected['after']['path'])),
             'Stored helper clock cost differs from literal paths')
        need(selected['target_pinned_count'] == 2 and
             sorted(selected['target_pinned_masks']) == sorted(
                 [selected['before']['mask'], selected['after']['mask']]),
             'Complete binary target-pinned helper sidecar changed')
        for key in ('before', 'after'):
            rec = selected[key]
            need(rec['resource'] == min(rec['orbit']), 'Normalized resource sidecar changed')
            resources.append(rec['resource'])
    need(len(set(resources)) == len(resources), 'Authenticated option resource union collides')
    banks = [analyse_bank(bank, helper_map) for bank in cert['banks']]
    old_path = HERE.parent.parent/'native_packet_decomposition_20261023/native/domains.json'
    old = json.loads(old_path.read_text())['sources']
    old_costs = []
    for slot in old:
        choices = {r['target']: r['path'] for r in slot['exact_all_shared_clocks_target_domain']}
        old_costs.append(sum(x != y for x, y in zip(choices[1965], choices[2234])))
    need(len(old_costs) == 28, 'Old native core universe changed')
    core_pair_cost = sum(sorted(old_costs, reverse=True)[:2])
    for bank in banks:
        if 'single_U_robust_subbank' in bank:
            robust = bank['single_U_robust_subbank']
            robust['old_28_pair_clock_cost_bound'] = core_pair_cost
            robust['controller_owner_clock_bound'] = core_pair_cost + robust['max_changed_helper_clock_cost']
    b53 = next(bank for bank in banks if bank['helper_count'] == 53)
    four_extensions = [four_occupancy_extension(sources, helper_map, core_pair_cost,
                                               b53['global_minimum_adversary'])
                       for sources in (b53['source_ids'], b53['single_U_robust_subbank']['source_ids'])]
    out = dict(status='FINITE-EXACT compatibility resilience conditional on actual forced-choice realization',
               input_sha256_lf_normalized=sha256(SUPPLY.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
               complete_native_sidecar='../supply/selected_supply.json',
               old_core_sha256_lf_normalized=sha256(old_path.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
               bounded_hitting_exhaustive_controls=hitting_controls,
               charge_encoding='q=16*(sum X mod128) + four-bit xor; source order is bank source_ids',
               bounds_scope='All single forcedU/forcedV tests exact; per-charge minimum hitting size exact only up to3.',
               banks=banks, four_occupancy_extensions=four_extensions, checks_before_replay=CHECKS)
    full_path = HERE/'resilience.full.json.gz'
    if args.write:
        stream = io.BytesIO()
        with gzip.GzipFile(fileobj=stream, mode='wb', filename='', mtime=0) as archive:
            archive.write((json.dumps(out, sort_keys=True, separators=(',', ':'))+'\n').encode())
        full_path.write_bytes(stream.getvalue())
    else:
        need(out == json.loads(gzip.decompress(full_path.read_bytes())),
             'Frozen full resilience certificate changed')
    summary = json.loads(json.dumps(out))
    for bank in summary['banks']:
        del bank['per_charge_hitting_certificates']
        if 'single_U_robust_subbank' in bank:
            del bank['single_U_robust_subbank']['compact_charge_sections']
    for bank in summary['banks']+summary['four_occupancy_extensions']:
        for row in bank['single_pin_rows']:
            for letter in ('U', 'V'):
                missing = row.pop('forced_'+letter+'_missing')
                row['forced_'+letter+'_first_missing'] = missing[0] if missing else None
    summary['full_certificate_file'] = full_path.name
    summary['full_certificate_sha256'] = sha256(full_path.read_bytes()).hexdigest()
    path = HERE/'resilience.json'
    if args.write:
        path.write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    else:
        need(summary == json.loads(path.read_text()), 'Frozen resilience summary changed')
    for b in banks:
        print(json.dumps({k: b[k] for k in ['helper_count','V_occupancy','minimum_witness_count',
            'every_single_U_survives','every_single_V_survives','universal_forced_U_tolerance_exact',
            'bounded_hitting_histogram']}))
        if 'single_U_robust_subbank' in b:
            print('robust_subbank', b['single_U_robust_subbank']['helper_count'])
    for b in four_extensions:
        print('four_extension', json.dumps({k: b[k] for k in [
            'helper_count','every_single_U_survives','every_single_V_survives',
            'bounded_forced_U_hitting_histogram','controller_owner_clock_bound']}))
    print(f'PASS {CHECKS} checks')


if __name__ == '__main__':
    main()
