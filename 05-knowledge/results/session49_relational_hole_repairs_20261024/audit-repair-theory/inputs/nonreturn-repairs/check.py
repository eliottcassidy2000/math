#!/usr/bin/env python3
"""Standalone exact non-reacquisition controls. No asserts or parent imports."""
from collections import deque
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def bits(mask, n):
    return tuple((mask >> i) & 1 for i in range(n))


def mask_of(xs):
    return sum(x << i for i, x in enumerate(xs))


def subsets(mask):
    value = mask
    while True:
        yield value
        if value == 0:
            break
        value = (value - 1) & mask


def cost(mask, weights):
    return sum(w for i, w in enumerate(weights) if mask >> i & 1)


def endpoint_dag(F, x, y, weights, budget=None):
    """Minimum weighted bottleneck, or minimum event count at fixed budget."""
    S = x ^ y
    nodes = sorted((I for I in subsets(S) if x ^ I in F),
                   key=lambda I: (I.bit_count(), I))
    best = {0: 0}
    prior = {}
    for I in nodes:
        if I not in best:
            continue
        for J in nodes:
            if I == J or I & ~J:
                continue
            c = cost(J ^ I, weights)
            if budget is not None and c > budget:
                continue
            score = best[I] + 1 if budget is not None else max(best[I], c)
            if J not in best or score < best[J]:
                best[J] = score
                prior[J] = I
    if S not in best:
        return None, []
    route = [y]
    at = S
    while at:
        at = prior[at]
        route.append(x ^ at)
    route.reverse()
    return best[S], route


def literal_record_search(F, x, y, n, weights, budget):
    """Independent BFS retains full per-slot visited-record masks."""
    start = (x, tuple(1 << bit for bit in bits(x, n)))
    queue = deque([(start, 0)])
    seen = {start}
    while queue:
        (z, history), depth = queue.popleft()
        if z == y:
            return depth
        for t in F:
            delta = z ^ t
            if not delta or cost(delta, weights) > budget:
                continue
            next_bits = bits(t, n)
            if any(delta >> i & 1 and history[i] >> next_bits[i] & 1
                   for i in range(n)):
                continue
            nxt = (t, tuple(h | (1 << b) for h, b in zip(history, next_bits)))
            if nxt not in seen:
                seen.add(nxt)
                queue.append((nxt, depth + 1))
    return None


def ordinary_route(F, x, y, weights, budget):
    queue = deque([x])
    previous = {x: None}
    while queue:
        z = queue.popleft()
        if z == y:
            route = []
            while z is not None:
                route.append(z)
                z = previous[z]
            return route[::-1]
        for t in F:
            if t not in previous and cost(t ^ z, weights) <= budget:
                previous[t] = z
                queue.append(t)
    return []


def binary_checks():
    report = []
    for n in range(4):
        families = pairs = cases = 0
        weights = tuple(range(1, n + 1))
        for selected in range(1, 1 << (1 << n)):
            F = {s for s in range(1 << n) if selected >> s & 1}
            families += 1
            for x in sorted(F):
                for y in sorted(F):
                    pairs += 1
                    optimum, route = endpoint_dag(F, x, y, weights)
                    need(route[0] == x and route[-1] == y, 'DAG has literal endpoints')
                    need(all(z in F for z in route), 'DAG arrivals accepted')
                    changes = [a ^ b for a, b in zip(route, route[1:])]
                    need(sum(z.bit_count() for z in changes) == (x ^ y).bit_count(),
                         'strict total slot work equals endpoint difference')
                    need(sum(cost(z, weights) for z in changes) == cost(x ^ y, weights),
                         'strict total literal clock work is endpoint work')
                    first = None
                    for B in range(sum(weights) + 1):
                        expected = literal_record_search(F, x, y, n, weights, B)
                        observed, unused = endpoint_dag(F, x, y, weights, B)
                        need(observed == expected, 'DAG equals independent visited-record BFS')
                        if observed is not None and first is None:
                            first = B
                        cases += 1
                    need(optimum == first, 'weighted bottleneck is exact')
        report.append(dict(slots=n, all_nonempty_families=families,
                           ordered_endpoint_pairs=pairs, weighted_budget_cases=cases))
    return report


def feasible_chain(F, target):
    chain = [target]
    while target:
        candidates = [target ^ (1 << i) for i in range(target.bit_length())
                      if target >> i & 1 and target ^ (1 << i) in F]
        need(bool(candidates), 'accessibility provides predecessor')
        target = min(candidates)
        chain.append(target)
    return chain[::-1]


def antimatroid_checks():
    report = []
    for n in range(5):
        ground = (1 << n) - 1
        count = pairs = 0
        weights = tuple(i + 1 for i in range(n))
        for encoded in range(1 << (1 << n)):
            if not encoded & 1 or not encoded >> ground & 1:
                continue
            F = {a for a in range(1 << n) if encoded >> a & 1}
            if any(not any((a ^ (1 << i)) in F for i in range(n) if a >> i & 1)
                   for a in F if a):
                continue
            if any(a | b not in F for a in F for b in F):
                continue
            count += 1
            chains = {a: feasible_chain(F, a) for a in F}
            for a in F:
                for b in F:
                    pairs += 1
                    raw = [a | z for z in chains[b]] + [b | z for z in chains[a][::-1]]
                    route = []
                    for z in raw:
                        if not route or route[-1] != z:
                            route.append(z)
                    need(route[0] == a and route[-1] == b and all(z in F for z in route),
                         'antimatroid union construction stays accepted')
                    changes = [u ^ v for u, v in zip(route, route[1:])]
                    need(all(d.bit_count() == 1 for d in changes), 'one block at a time')
                    need(len(changes) == (a ^ b).bit_count(), 'antimatroid geodesic is no-return')
                    expected = max((weights[i] for i in range(n) if (a ^ b) >> i & 1), default=0)
                    optimum, unused = endpoint_dag(F, a, b, weights)
                    need(optimum == expected, 'exact nonuniform block bottleneck')
        report.append(dict(blocks=n, antimatroids=count, ordered_endpoint_pairs=pairs))
    return report


def partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in partitions(rest):
        yield [(first,)] + p
        for i in range(len(p)):
            yield p[:i] + [(first,) + p[i]] + p[i + 1:]


def neutral_checks():
    cases = 0
    for n in range(5):
        weights = tuple(i + 1 for i in range(n))
        for ds in product(tuple(product((-1, 1), range(3))), repeat=n):
            if sum(d[0] for d in ds) or sum(d[1] for d in ds) % 3:
                continue
            F = {I for I in range(1 << n)
                 if sum(ds[i][0] for i in range(n) if I >> i & 1) == 0
                 and sum(ds[i][1] for i in range(n) if I >> i & 1) % 3 == 0}
            best = None
            best_atoms = None
            for p in partitions(tuple(range(n))):
                masks = [sum(1 << i for i in block) for block in p]
                if not all(I in F for I in masks):
                    continue
                value = max((cost(I, weights) for I in masks), default=0)
                best = value if best is None else min(best, value)
                if all(not any(J and J != I and J in F for J in subsets(I)) for I in masks):
                    best_atoms = value if best_atoms is None else min(best_atoms, value)
            optimum, unused = endpoint_dag(F, 0, (1 << n) - 1, weights)
            need(optimum == best == best_atoms, 'neutral partition and atom bottlenecks exact')
            cases += 1
    return dict(group='Z x Z/3', all_signed_integer_mod3_controls_through_slots=4,
                balanced_label_cases=cases)


def davenport_checks():
    reports = []
    for orders in ((2,), (4,), (2, 2), (3, 3)):
        group = tuple(product(*(range(a) for a in orders)))
        zero = (0,) * len(orders)
        D = 1 + sum(a - 1 for a in orders)
        cases = 0
        for sequence in product(group, repeat=D):
            reachable = {zero}
            found = False
            for g in sequence:
                if tuple(-v % a for v, a in zip(g, orders)) in reachable:
                    found = True
                    break
                reachable |= {tuple((x + y) % a for x, y, a in zip(h, g, orders))
                              for h in reachable.copy()}
            need(found, 'small p-group D-term sequence has zero-sum subsequence')
            cases += 1
        reports.append(dict(cyclic_orders=list(orders), D=D, exhaustive_D_term_sequences=cases))
    def plus(a, b):
        return (((a & 127) + (b & 127)) & 127) | ((a ^ b) & 1920)
    generators = (1, 128, 256, 512, 1024)
    lower = [1] * 127 + list(generators[1:])
    coefficients = {0: 1}
    augmentation_product = {0}
    for g in lower:
        updated = coefficients.copy()
        for h, count in coefficients.items():
            target = plus(h, g)
            updated[target] = updated.get(target, 0) + count
        coefficients = updated
        augmentation_product ^= {plus(h, g) for h in augmentation_product}
        need(coefficients[0] == 1, 'length131 lower example remains zero-sumfree')
    need(len(lower) == 131 and augmentation_product == set(range(2048)),
         'sharp nonzero augmentation product at degree131')
    degrees = [(v & 127) + (v >> 7).bit_count() for v in range(2048)]
    need(max(degrees) == 131 and degrees.count(131) == 1, 'truncated algebra maximum degree131')
    for g in range(2048):
        need(augmentation_product ^ {plus(h, g) for h in augmentation_product} == set(),
             'all further augmentation factors kill sharp top product')
    final = 1 | 1920
    updated = coefficients.copy()
    for h, count in coefficients.items():
        target = plus(h, final)
        updated[target] = updated.get(target, 0) + count
    need(updated[0] == 2, 'minimal length132 sequence has only empty and full zero-sum subsets')
    # Independent fully enumerated population-pair control in C4 x C2, D=5.
    seq = ((1, 0),) * 3 + ((0, 1), (1, 1))
    accepted = []
    for xs in product((0, 1), repeat=10):
        if any(xs[2 * i] + xs[2 * i + 1] != 1 for i in range(5)):
            continue
        charge = tuple(sum(xs[2 * i] * seq[i][j] for i in range(5)) % modulus
                       for j, modulus in enumerate((4, 2)))
        if charge == (0, 0):
            accepted.append(xs)
    need(set(accepted) == {(0, 1) * 5, (1, 0) * 5}, 'small exact population-only sharp fibre')
    return dict(small_groups=reports, dyadic_group_order=2048, dyadic_D=132,
                augmentation_nilpotence_index=132, lower_zero_sumfree_length=131,
                minimal_length132_sequence_zero_sum_subsets=updated[0],
                conditional_population_pair_slot_bound=264,
                small_sharp_population_control_accepted_states=len(accepted))


P = ((3, 0, 6, 7, 1, 5, 4, 2), (3, 1, 6, 7, 0, 5, 4, 2))
G = (4, 0, 2, 3, 1, 5, 6, 7)


def source(r):
    seed = [[0] * 8 for unused in range(8)]
    for j, (u, v) in enumerate(combinations(range(8), 2)):
        seed[u][v] = (29558785 >> (27 - j)) & 1
        seed[v][u] = 1 - seed[u][v]
    q = 8 * r
    def arc(u, v):
        if u == v:
            return 0
        if u < q and v < q:
            return seed[u % 8][v % 8] if u // 8 == v // 8 else int(u < v)
        if (u < q) != (v < q):
            return int(u < q)
        return int((v - u) % 3 == 1)
    return tuple(tuple(arc(u, v) for v in range(q + 3)) for u in range(q + 3))


def path(r, key, phase):
    return tuple(8 * j + a for j in range(r) for a in P[(key >> j) & 1]) + \
        tuple(8 * r + (phase + t) % 3 for t in range(3))


def transpose(a):
    return tuple(zip(*a))


def partner(a, p):
    pos = {u: j for j, u in enumerate(p)}
    return tuple(tuple(0 if u == v else a[u][v] if abs(pos[u] - pos[v]) == 1
                       else a[v][u] for v in range(len(a))) for u in range(len(a)))


def word(a, p):
    return tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i + 2, len(p)))


def iso(a, b, perm):
    return all(a[u][v] == b[perm[u]][perm[v]] for u in range(len(a)) for v in range(len(a)))


def reservoir(a, p):
    w = word(a, p)
    rw = word(transpose(a), p[::-1])
    return frozenset((w, rw, tuple(1 - z for z in w), tuple(1 - z for z in rw)))


def native_checks():
    report = []
    for r in range(1, 5):
        a = source(r)
        q = 8 * r
        N = q + 3
        paths = {(b, k): path(r, b, k) for b in range(1 << r) for k in range(3)}
        targets = {key: partner(a, p) for key, p in paths.items()}
        banks = {b: reservoir(a, paths[b, 0]) for b in range(1 << r)}
        need(all(len(bank) == 4 for bank in banks.values()), 'every normalized reservoir has four words')
        need(len(set().union(*banks.values())) == 4 * (1 << r), 'all prefix-key reservoirs disjoint')
        for (b, k), p in paths.items():
            need(sorted(p) == list(range(N)) and all(a[u][v] for u, v in zip(p, p[1:])),
                 'literal native Hamiltonian path')
            need(reservoir(a, p) == banks[b], 'three full kernel lifts share one resource')
            need(p[q] - q == k, 'canonical phase subtracts raw offset')
            need((p[q] == q) == (k == 0), 'phase-zero literal owner pin')
            need(p[1] == b % 2, 'Boolean variable is actual owner at clock one')
            need(p[::-1][2] - q == k, 'converse reflected phase pin')
            need(partner(transpose(a), p[::-1]) == transpose(targets[b, k]), 'literal converse partner square')
            for j in range(r):
                g = G if not (b >> j & 1) else tuple(G.index(u) for u in range(8))
                perm = tuple(8 * j + g[u % 8] if u < q and u // 8 == j else u for u in range(N))
                need(iso(targets[b, k], targets[b ^ (1 << j), k], perm), 'actual partner gluing')
            for step in (1, 2):
                perm = tuple(u if u < q else q + (u - q + step) % 3 for u in range(N))
                need(iso(a, a, perm) and iso(targets[b, k], targets[b, (k + step) % 3], perm),
                     'actual phase lift source and partner transport')
        for b in range(1 << r):
            prefix = paths[b, 0][:-3]
            actual = {prefix + suffix for suffix in permutations(range(q, q + 3))
                      if all(a[u][v] for u, v in zip(prefix + suffix, (prefix + suffix)[1:]))}
            need(actual == {paths[b, k] for k in range(3)}, 'all suffix kernel lifts retained')
        for i in range(1 << (r - 1)):
            p, t = paths[2 * i, 0], paths[2 * i + 1, 0]
            need([j for j in range(N) if p[j] != t[j]] == [1, 4], 'paired keys exact two-clock cost')
            need([j for j in range(N) if p[::-1][j] != t[::-1][j]] == sorted([N - 2, N - 5]),
                 'converse preserves exact two-clock cost')
        report.append(dict(prefix_blocks=r, source_order=N, reconstructed_full_records=len(paths),
                           normalized_resources=len(banks), closed_reservation_words=4 * len(banks)))
    seed = tuple(row[:8] for row in source(1)[:8])
    degrees = tuple(map(sum, seed))
    autos = [p for p in permutations(range(8)) if all(degrees[i] == degrees[p[i]] for i in range(8))
             and iso(seed, seed, p)]
    need(autos == [tuple(range(8))], 'independent seed rigidity census')
    for origin in range(8):
        seen = {origin}
        todo = [origin]
        while todo:
            u = todo.pop()
            nxt = {v for v in range(8) if seed[u][v] and v not in seen}
            seen |= nxt
            todo.extend(nxt)
        need(len(seen) == 8, 'seed is strongly connected')
    # Exhaust all Boolean truth tables on three literal owner predicates.
    r, n = 3, 3
    record_options = {xs: tuple(path(r, 2 * i + x, 0) for i, x in enumerate(xs))
                      for xs in product((0, 1), repeat=n)}
    for encoded in range(256):
        symbolic = {xs for xs in record_options if encoded >> mask_of(xs) & 1}
        literal = {xs for xs, records in record_options.items()
                   if encoded >> mask_of(tuple(p[1] for p in records)) & 1
                   and all(p[24] == 24 for p in records)}
        need(symbolic == literal, 'parsimonious literal Boolean-predicate embedding')
    return dict(carriers=report, seed_automorphisms=len(autos), all_three_bit_boolean_relations=256)


def hostile_and_positive_checks():
    sequences = [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 1, 0)]
    F = {mask_of(xs) for xs in sequences}
    x, y = mask_of(sequences[0]), mask_of(sequences[-1])
    ord_route = ordinary_route(F, x, y, (2, 2, 2), 2)
    need([bits(z, 3) for z in ord_route] == sequences, 'unique native cheap hostile path')
    need(literal_record_search(F, x, y, 3, (2, 2, 2), 2) is None, 'no strict cheap native route')
    optimum, direct = endpoint_dag(F, x, y, (2, 2, 2))
    need(optimum == 4 and direct == [x, y], 'strict native endpoint threshold exactly four clocks')
    modes = [{z for z in F if (z & 1) == mode} for mode in (0, 1)]
    need(all(ordinary_route(mode, a, b, (2, 2, 2), 2) for mode in modes for a in mode for b in mode),
         'both unconditioned fibres internally cheap connected')
    need(min(cost(a ^ b, (2, 2, 2)) for a in modes[0] for b in modes[1]) == 2,
         'unconditioned exact fibre boundary is cheap')
    native_arrivals = []
    a = source(3)
    for xs in sequences:
        keys = tuple(2 * i + bit for i, bit in enumerate(xs))
        records = tuple(path(3, key, 0) for key in keys)
        banks = [reservoir(a, p) for p in records]
        need(len(set().union(*banks)) == 12, 'every hostile arrival has three actual units/twelve words')
        native_arrivals.append(dict(bits=list(xs), keys=list(keys)))
    need([entry['keys'][2] for entry in native_arrivals] == [4, 5, 5, 5, 4],
         'literal normalized key reacquisition, not just a phase return')
    # Positive synchronized order-ideal family: blocks of 1,2,3 physical slots.
    sizes = (1, 2, 3)
    ideals = {0, 1, 3, 7}  # chain 0<1<2
    physical = {I: tuple((I >> j) & 1 for j, size in enumerate(sizes) for unused in range(size))
                for I in ideals}
    positive_pairs = 0
    for I in ideals:
        for J in ideals:
            positive_pairs += 1
            Fp = {mask_of(v) for v in physical.values()}
            optimum, route = endpoint_dag(Fp, mask_of(physical[I]), mask_of(physical[J]), (2,) * 6)
            expected = max((2 * size for j, size in enumerate(sizes) if (I ^ J) >> j & 1), default=0)
            need(optimum == expected, 'synchronized native ideal family exact weighted bottleneck')
            previous_keys = None
            histories = [set() for unused in range(6)]
            for z in route:
                keys = tuple(2 * i + bit for i, bit in enumerate(bits(z, 6)))
                records = tuple(path(4, key, 0) for key in keys)
                need(all(p[32] == 32 for p in records), 'positive actual phase pin')
                for i, key in enumerate(keys):
                    if previous_keys is not None and key != previous_keys[i]:
                        need(key not in histories[i], 'positive route never reacquires actual key')
                    histories[i].add(key)
                previous_keys = keys
    return dict(hostile=dict(accepted_states=[list(xs) for xs in sequences],
                            ordinary_slots=1, ordinary_clocks=2, strict_slots=2, strict_clocks=4,
                            literal_arrivals=native_arrivals),
                positive_ideal_block_sizes=list(sizes), positive_ordered_endpoint_pairs=positive_pairs)


def history_checks():
    # Every distinct-record sequence, including retained holds, on three menus.
    menus = (2, 3, 4)
    maximal = sum(m - 1 for m in menus)
    need(maximal == 6, 'finite horizon sum of menu sizes minus one')
    states = {(tuple(0 for unused in menus), tuple(1 for unused in menus))}
    layer = states.copy()
    largest = 0
    while layer:
        following = set()
        for current, seen in layer:
            for target in product(*(range(m) for m in menus)):
                changed = [i for i in range(len(menus)) if current[i] != target[i]]
                if not changed or any(seen[i] >> target[i] & 1 for i in changed):
                    continue
                new_seen = tuple(h | (1 << v) for h, v in zip(seen, target))
                need(sum(h.bit_count() for h in new_seen) > sum(h.bit_count() for h in seen),
                     'every genuine history transition strictly increases rank')
                nxt = (target, new_seen)
                if nxt not in states:
                    states.add(nxt)
                    following.add(nxt)
        largest = max(largest, max((sum(h.bit_count() - 1 for h in seen) for unused, seen in layer), default=0))
        layer = following
    need(largest == maximal, 'finite history rank bound attained')
    return dict(menu_sizes=list(menus), history_states=len(states), maximum_changed_slot_occurrences=maximal)


def source_checks():
    pins = json.loads((HERE / 'SOURCEPINS.json').read_text())
    for name, info in pins['files'].items():
        data = (HERE / name).read_bytes()
        need(len(data) == info['bytes'] and hashlib.sha256(data).hexdigest() == info['sha256'],
             'exact byte source pin: ' + name)
    return {name: entry['sha256'] for name, entry in pins['files'].items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    receipt = dict(schema='session49-nonreturn-repairs-receipt-v1',
                   sourcepins=source_checks(), binary=binary_checks(),
                   antimatroids=antimatroid_checks(), neutral_partitions=neutral_checks(),
                   davenport=davenport_checks(),
                   native=native_checks(), routing=hostile_and_positive_checks(), history=history_checks())
    receipt['checks'] = CHECKS
    payload = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.write_receipt:
        (HERE / 'receipt.json').write_text(payload)
    else:
        need((HERE / 'receipt.json').read_text() == payload, 'deterministic receipt matches frozen result')
    print(payload, end='')


if __name__ == '__main__':
    main()
