#!/usr/bin/env python3
"""Exact finite controls for conditional guarded-packet decomposition.

No checkout writes, imports of discovery scripts, or optimization-sensitive
assertions. Original native authentication is re-run read-only in a subprocess.
"""
import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BASE = Path('05-knowledge/results/session45_injective_multisource_moves_20261023')
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def add(a, b, moduli):
    return tuple((x + y) % m for x, y, m in zip(a, b, moduli))


def neg(a, moduli):
    return tuple((-x) % m for x, m in zip(a, moduli))


def total(labels, indices, moduli):
    out = (0,) * len(moduli)
    for i in indices:
        out = add(out, labels[i], moduli)
    return out


def balanced(labels, positive, negative, moduli):
    return bool(positive) and len(positive) == len(negative) and total(labels, positive, moduli) == total(labels, negative, moduli)


def extract_atom(labels, positive, negative, moduli):
    """Minimum-cardinality cross-subset search, not fixed-pair minimality."""
    positive, negative = tuple(positive), tuple(negative)
    for k in range(1, min(len(positive), len(negative)) + 1):
        left = {}
        for subset in itertools.combinations(positive, k):
            left.setdefault(total(labels, subset, moduli), subset)
        for subset in itertools.combinations(negative, k):
            key = total(labels, subset, moduli)
            if key in left:
                return left[key], subset
    raise ValueError('No balanced nonempty subtrade')


def decompose(labels, positive, negative, moduli):
    positive, negative = list(positive), list(negative)
    need(not set(positive).intersection(negative), 'Disjoint endpoint change indices')
    need(len(set(positive)) == len(positive) and len(set(negative)) == len(negative), 'No repeated slot identity')
    need(len(positive) == len(negative) and total(labels, positive, moduli) == total(labels, negative, moduli), 'Accepted endpoint difference')
    out = []
    while positive:
        a, b = extract_atom(labels, positive, negative, moduli)
        out.append((a, b))
        positive = [i for i in positive if i not in a]
        negative = [i for i in negative if i not in b]
    need(not negative, 'Decomposition exhausted both signs')
    return out


def verify_primitive(labels, positive, negative, moduli):
    need(balanced(labels, positive, negative, moduli), 'Atom is balanced')
    for k in range(1, len(positive)):
        left = {total(labels, s, moduli) for s in itertools.combinations(positive, k)}
        need(not any(total(labels, s, moduli) in left for s in itertools.combinations(negative, k)), 'No proper cross-subtrade')


def signature(labels, state, moduli):
    return len(state), total(labels, state, moduli)


def verify_route(labels, start, finish, atoms, moduli, resources=None):
    state = set(start)
    invariant = signature(labels, state, moduli)
    if resources is None:
        resources = [(('slot', i, 0), ('slot', i, 1)) for i in range(len(labels))]
    selected = lambda s: {resources[i][int(i in s)] for i in range(len(labels))}
    endpoint_union = selected(start) | selected(finish)
    states = [sorted(state)]
    changed = set()
    for a, b in atoms:
        verify_primitive(labels, a, b, moduli)
        need(not set(a).intersection(state) and set(b) <= state, 'Conformal event direction')
        need(not changed.intersection(a + b), 'Every changed slot switches once')
        changed.update(a + b)
        state.difference_update(b)
        state.update(a)
        need(signature(labels, state, moduli) == invariant, 'Every intermediate keeps both guards')
        now = selected(state)
        need(len(now) == len(labels), 'Every intermediate has resource injection')
        need(now <= endpoint_union, 'Intermediate resources stay inside endpoint union')
        states.append(sorted(state))
    need(state == set(finish), 'Route reaches exact endpoint')
    return states


def exhaustive_small():
    summaries = []
    for moduli in [(), (2,), (3,), (4,), (2, 2)]:
        elements = list(itertools.product(*(range(m) for m in moduli)))
        arrays = pairs = orders = atoms_total = max_support = 0
        for n in range(6):
            for labels in itertools.product(elements, repeat=n):
                arrays += 1
                fibers = defaultdict(list)
                for mask in range(1 << n):
                    state = frozenset(i for i in range(n) if (mask >> i) & 1)
                    fibers[signature(labels, state, moduli)].append(state)
                for fiber in fibers.values():
                    for start, finish in itertools.combinations(fiber, 2):
                        pairs += 1
                        atoms = decompose(labels, sorted(finish-start), sorted(start-finish), moduli)
                        atoms_total += len(atoms)
                        max_support = max(max_support, *(len(a)+len(b) for a,b in atoms))
                        for perm in itertools.permutations(atoms):
                            verify_route(labels, start, finish, perm, moduli)
                            orders += 1
        summaries.append({'moduli': moduli, 'max_slots': 5, 'ordered_label_arrays': arrays,
                          'distinct_accepted_endpoint_pairs': pairs, 'atom_orders_checked': orders,
                          'extracted_atoms': atoms_total, 'largest_extracted_support': max_support})
    return summaries


def davenport_controls():
    moduli = (2, 2, 2, 2, 128)
    basis = [tuple(int(i == j) for i in range(5)) for j in range(5)]
    sequence = [g for g, m in zip(basis, moduli) for _ in range(m-1)]
    zero = (0,) * 5
    reachable = {zero: {0}}
    algebra_support = {zero}
    for g in sequence:
        following = {h: set(lengths) for h, lengths in reachable.items()}
        for h, lengths in reachable.items():
            following.setdefault(add(h, g, moduli), set()).update(t+1 for t in lengths)
        reachable = following
        algebra_support ^= {add(h, g, moduli) for h in algebra_support}
    need(reachable[zero] == {0}, 'Length-131 lower sequence is zero-sum-free')
    need(len(algebra_support) == 2048, 'Top augmentation monomial has all group coefficients one')
    special = tuple(1 for _ in moduli)
    following = {h: set(lengths) for h, lengths in reachable.items()}
    for h, lengths in reachable.items():
        following.setdefault(add(h, special, moduli), set()).update(t+1 for t in lengths)
    algebra_support ^= {add(h, special, moduli) for h in algebra_support}
    need(following[zero] == {0,132}, 'Appended sequence is minimal zero-sum of length 132')
    need(not algebra_support, 'Length-132 group-algebra product vanishes')
    return {'moduli': moduli, 'group_order': 2048, 'nilpotence_degree': 132,
            'abstract_minimal_zero_sum_length': 132, 'abstract_balanced_trade_support': 264,
            'abstract_trade_construction': 'Positive minimal zero-sum sequence versus 132 distinct negative slots of label zero',
            'conditional_native_support_bound': 264,
            'native_sharpness_claimed': False, 'all_order_constant_claimed': False}


def buffer_route(original_labels, positive, negative, moduli):
    verify_primitive(original_labels, positive, negative, moduli)
    k = len(positive)
    zero = (0,) * len(moduli)
    partial = [zero]
    for a, b in zip(positive, negative):
        d = add(original_labels[a], neg(original_labels[b], moduli), moduli)
        partial.append(add(partial[-1], d, moduli))
    need(partial[-1] == zero and len(set(partial[:-1])) == k, 'Primitive partial sums are distinct')
    labels = list(original_labels) + [neg(g, moduli) for g in partial[:-1]]
    n = len(original_labels)
    state = set(negative) | {n}
    initial = set(state)
    guard = signature(labels, state, moduli)
    states = [sorted(state)]
    event_sizes = []
    projected_debts = []
    for j, (a, b) in enumerate(zip(positive, negative), 1):
        old = set(state)
        state.remove(b)
        state.add(a)
        state.remove(n+j-1)
        state.add(n+(j % k))
        need(signature(labels, state, moduli) == guard, 'Charged buffer pays exact current and margin')
        need(len({(i, int(i in state)) for i in range(len(labels))}) == len(labels), 'Separate synthetic resource per slot choice')
        event_sizes.append(len(old ^ state))
        debt = add(total(labels, state & set(range(n)), moduli), neg(total(labels, negative, moduli), moduli), moduli)
        projected_debts.append(debt)
        if j < k:
            need(debt != zero, 'Original projection is illegal before restoration')
        states.append(sorted(state))
    need(state == set(positive) | {n}, 'Exact original endpoint and buffer return')
    need(max(event_sizes) <= 4, 'At most four changed slots per event')
    return {'scope': 'Abstract supply-conditional buffer control; not native buffer existence',
            'original_slots': len(original_labels), 'buffer_slots': k, 'labels': labels[n:],
            'states': states, 'event_sizes': event_sizes, 'projected_debts': projected_debts}


def hostile_controls():
    reports = {}
    labels = [(0,), (1,), (1,), (0,)]
    paired_differences = [(2,), (1,)]
    need(total(paired_differences, [0,1], (3,)) == (0,) and all(g != (0,) for g in paired_differences), 'Fixed pairing is minimal zero sum')
    atoms = decompose(labels, [0,1], [2,3], (3,))
    need(len(atoms) == 2, 'Fixed-pair minimality is not balanced primitivity')
    reports['pairing_minimality_loss'] = {'paired_differences': paired_differences, 'actual_atoms': atoms}
    need(total([(0,)], [0], (2,)) == (0,) and not balanced([(0,)], [0], [], (2,)), 'Zero charge without target cardinality is illegal')
    reports['cardinality_cannot_be_dropped'] = True
    full = [(1,0), (0,0)]
    need(full[0][1] == full[1][1] and not balanced(full, [0], [1], (2,4)), 'Cyclic coordinate misses binary debt')
    reports['scalar_phase_loses_dyadic_debt'] = True
    need(0 % 3 == 3 % 3 and 0 != 3, 'Quotient agreement loses actual owner identity')
    reports['shared_owner_quotient_loss'] = {'actual_values': [0,3], 'quotient_modulus': 3}
    raw = ['lift-a', 'lift-b']; normalized = ['tile-a', 'tile-a']
    need(len(set(raw)) == 2 and len(set(normalized)) == 1, 'Different lifts are not two normalized resources')
    reports['resource_quotient_loss'] = True
    # Symmetric zero-row-sum F and genuine rowwise group debt, with scalar sum zero.
    F = [[1,-1],[-1,1]]; g = [0,1]
    debts = [sum(a*b for a,b in zip(row,g)) % 4 for row in F]
    need(sum(debts) % 4 == 0 and any(debts), 'Scalar total erases opposing physical-row debts')
    reports['row_debt_not_scalar_total'] = {'matrix': F, 'labels_mod4': g, 'row_debts': debts}
    split = [1,1]  # Two coincident role contributions at one physical row, modulo 2.
    need(any(split) and sum(split) % 2 == 0, 'Overlapping roles must add before the guard test')
    reports['overlapping_roles_need_actual_owner_map'] = {'separate_contributions_mod2': split, 'merged_debt': 0}
    cycle_reports = []
    for n in range(3,11):
        accepted = []
        for mask in range(1 << n):
            bits = [(mask >> i) & 1 for i in range(n)]
            if all((1-bits[j])+bits[(j-1) % n] == 1 for j in range(n)):
                accepted.append(mask)
        need(accepted == [0,(1 << n)-1], 'Multi-target incidence cycle has one full-support atom')
        cycle_reports.append({'slots': n, 'accepted': accepted, 'charge_group_order': 1, 'minimal_support': n})
    reports['multi_target_cycle'] = cycle_reports
    # Neutral auxiliary bits cannot change the two-state projected relation.
    legal = [(x,a) for x in (0,15) for a in range(4)]
    edges = [(v,w) for v in legal for w in legal if (v[0]^w[0]).bit_count() <= 3]
    reached = {legal[0]}
    while True:
        new = reached | {w for v,w in edges if v in reached}
        if new == reached:
            break
        reached = new
    need(all(x == 0 for x,a in reached), 'Neutral workspace cannot cross primitive projected support')
    reports['neutral_projection_barrier'] = {'projected_support': 4, 'allowed_projected_budget': 3, 'reachable_states': len(reached)}
    # Resource constraints must be in the exact encoding, not checked only at endpoints.
    endpoints = [(0,1),(1,0)]
    mixed = (1,1)
    need(all(len(set(s)) == 2 for s in endpoints) and len(set(mixed)) == 1, 'Two legal endpoint injections can have illegal coordinate mixture')
    reports['missing_resource_row'] = {'endpoints': endpoints, 'illegal_mixture': mixed}
    need(3-2 == 1, 'Two nonprimitive moves generate the full integer lattice')
    generators = [(2,-2),(3,-3),(-2,2),(-3,3)]
    feasible = []
    for state in endpoints:
        for move in generators:
            target = tuple(a+b for a,b in zip(state,move))
            if all(0 <= x <= 1 for x in target):
                feasible.append((state,target))
    need(not feasible, 'Integer lattice generation does not connect a bounded fiber')
    reports['lattice_generation_not_markov_connectivity'] = {'generators': generators, 'binary_endpoints': endpoints, 'feasible_edges': feasible}
    return reports


def native_controls(repo):
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for rec in pins['files']:
        data = (repo/rec['path']).read_bytes()
        need(len(data) == rec['bytes'] and hashlib.sha256(data).hexdigest() == rec['sha256'], 'Inherited input pin: '+rec['path'])
    command = [sys.executable, '-B']
    if sys.flags.optimize:
        command += ['-O']
    command += [str(repo/BASE/'balanced-circuits/check.py'), '--repo', str(repo)]
    inherited = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    native = json.loads(inherited.stdout)
    need(native['status'] == 'PASS', 'Inherited complete native, owner, converse and integer-image replay')
    pool = json.loads((repo/BASE/'balanced-circuits/pool.json').read_text())
    sources = sorted(int(s) for s in pool['labels'])
    labels = [tuple((pool['labels'][str(s)][0] >> j) & 1 for j in range(4)) + (pool['labels'][str(s)][1],) for s in sources]
    moduli = (2,2,2,2,128)
    certs = {r['source']: r for r in native['literal_source_options']}
    resources = [tuple(certs[s]['resource_ids']) for s in sources]
    need(len({r for pair in resources for r in pair}) == 40, 'All forty normalized alternative resources are distinct')
    need(len({tuple(g) for g in labels}) == 20, 'Twenty distinct source charge labels')
    # Discover all native support-four relations in the fixed twenty-option supply.
    sums = defaultdict(list)
    rectangles = []
    for left in itertools.combinations(range(len(sources)), 2):
        key = total(labels, left, moduli)
        for right in sums[key]:
            if set(left).isdisjoint(right):
                rectangles.append((left,right))
        sums[key].append(left)
    need(len(rectangles) == 3, 'Exact fixed-supply support-four inventory')
    need(len(set(i for a,b in rectangles for i in a+b)) == 12, 'All three native rectangles are support-disjoint')
    positive = sorted(i for a,b in rectangles for i in a)
    negative = sorted(i for a,b in rectangles for i in b)
    atoms = decompose(labels, positive, negative, moduli)
    need(sorted(len(a)+len(b) for a,b in atoms) == [4,4,4], 'Native twelve-change endpoint pair decomposes into three atoms')
    routes = []
    for order in itertools.permutations(atoms):
        states = verify_route(labels, set(negative), set(positive), order, moduli, resources)
        routes.append([[sources[i] for i in state] for state in states])
    # Recheck exact primitive subtrade census independently of the inherited mask loops.
    candidates = json.loads((repo/BASE/'balanced-circuits/candidates.json').read_text())['found']
    circuit_reports = []
    sixteen = None
    for key, candidate in sorted(candidates.items(), key=lambda p: int(p[0])):
        p = [sources.index(s) for s in candidate['plus']]
        m = [sources.index(s) for s in candidate['minus']]
        verify_primitive(labels, p, m, moduli)
        found = decompose(labels, p, m, moduli)
        need(len(found) == 1 and len(found[0][0])+len(found[0][1]) == len(p)+len(m), 'Primitive native circuit returns one unsplit atom')
        verify_route(labels, set(m), set(p), found, moduli, resources)
        circuit_reports.append({'slots': len(p)+len(m), 'atoms': len(found), 'plus': candidate['plus'], 'minus': candidate['minus']})
        if len(p)+len(m) == 16:
            ordered = p+m
            sixteen = buffer_route([labels[i] for i in ordered], list(range(8)), list(range(8,16)), moduli)
    return {'source_options': 20, 'distinct_normalized_resources': 40,
            'inherited_authenticator_stdout_sha256': hashlib.sha256(inherited.stdout.encode()).hexdigest(),
            'inherited_authenticator_checks': native['checks'],
            'native_four_slot_atoms': [{'plus':[sources[i] for i in a], 'minus':[sources[i] for i in b]} for a,b in rectangles],
            'native_decomposable_twelve_change_routes': routes,
            'native_primitive_circuits': circuit_reports,
            'sixteen_label_abstract_buffer_control': sixteen,
            'source_pins_verified': len(pins['files'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = {'status': 'PASS', 'scope': 'Conditional all-size proofs; exhaustive finite controls; read-only inherited native replay',
              'exhaustive_small': exhaustive_small(), 'davenport': davenport_controls(),
              'hostiles': hostile_controls(), 'native': native_controls(args.repo.resolve())}
    result['checks'] = CHECKS
    rendered = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
