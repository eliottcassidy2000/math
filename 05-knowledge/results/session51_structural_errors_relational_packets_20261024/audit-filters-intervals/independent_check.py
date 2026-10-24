#!/usr/bin/env python3
"""Independent mathematical controls. Does not import either producer checker."""
import functools
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COUNTS = {}


def check(value, label):
    COUNTS[label] = COUNTS.get(label, 0) + 1
    if not value:
        raise RuntimeError(label)


def trial_prime(n):
    if n < 2:
        return False
    for k in range(2, math.isqrt(n) + 1):
        if n % k == 0:
            return False
    return True


def prime_audit():
    data = json.loads((ROOT / 'copies/prime-offsets/certificates.json').read_text())
    result = json.loads((ROOT / 'copies/prime-offsets/results.json').read_text())
    known = {2}
    certs = data['primality_certificates']
    for n in sorted(map(int, certs)):
        if n == 2:
            continue
        factors = certs[str(n)]['factors']
        check(len({q for q, e in factors}) == len(factors), 'independent Lucas distinct factors')
        check(all(type(q) is int and type(e) is int and q in known and e > 0 for q, e in factors), 'independent Lucas certified factor ancestry')
        check(math.prod(q ** e for q, e in factors) == n - 1, 'independent Lucas complete factorization')
        a = certs[str(n)]['base']
        check(type(a) is int and 1 < a < n and pow(a, n - 1, n) == 1, 'independent Lucas Fermat')
        for q, _ in factors:
            check(math.gcd(n, pow(a, (n - 1) // q, n) - 1) == 1, 'independent Lucas full order')
        known.add(n)
    check(len(known) == 200, 'all 200 Lucas certificates proved')

    def classification(n, witness):
        if witness.get('prime'):
            check(n in known, 'prime endpoint belongs to certified set')
            return True
        if witness.get('nonpositive_or_one'):
            check(n < 2, 'nonprime small endpoint')
        elif 'divisor' in witness:
            d = witness['divisor']
            check(type(d) is int and 1 < d < n and n // d * d == n, 'independent proper divisor')
        else:
            a = witness['mr_base']
            check(type(a) is int and 1 < a < n and n % 2 == 1, 'independent MR admissible input')
            t = n - 1
            powers = []
            while t % 2 == 0:
                t //= 2
                powers.append(t)
            # A prime must have a^t=1 or one of a^(t*2^j)=-1.
            residues = [pow(a, e, n) for e in powers]
            check(residues[-1] != 1 and n - 1 not in residues, 'independent MR exact failure')
        return False

    primes = [p for p in range(2, 72) if trial_prime(p)]
    possible = [p for p in primes if p < 71]
    check(len(data['rows']) == len(possible) == 19, 'exact finite cutoff indices')
    total = 1
    computed = {'plus': [], 'minus': []}
    endpoint_count = 0
    for i, (p, row) in enumerate(zip(possible, data['rows']), 1):
        total *= p
        check((row['n'], row['p'], row['P']) == (i, p, total), 'independent primorial indexing')
        for side, sign in [('plus', 1), ('minus', -1)]:
            entry = row[side]
            if side == 'minus' and i == 1:
                check(entry is None, 'L1 is undefined')
                computed[side].append(None)
                continue
            d = entry['offset']
            check(len(entry['witnesses']) == d - 1, 'independent complete minimality range')
            answers = [classification(total + sign * k, w) for k, w in enumerate(entry['witnesses'], 2)]
            endpoint_count += len(answers)
            check(answers == [False] * (d - 2) + [True], 'independent first prime endpoint')
            check(d > p or (i == 2 and side == 'minus' and d == 3), 'cutoff and isolated L2 exception')
            computed[side].append(d)
    check(endpoint_count == 1734, '1734 endpoint classifications')
    check(computed['plus'] == result['F_1_through_19'] and computed['minus'] == result['L_1_through_19'], 'sequence outputs agree')
    check(computed['plus'][5] == computed['minus'][5] == 17, 'common minimum at primorial 30030')

    # Independent census: every short Weierstrass model, not one model per j.
    curve_rows = []
    for row in result['arithmetic_predicates']:
        p = row['p']
        if p < 5:
            continue
        square_multiplicity = [0] * p
        for y in range(p):
            square_multiplicity[y * y % p] += 1
        by_j = {}
        models = 0
        for a in range(p):
            for b in range(p):
                disc = (4 * a ** 3 + 27 * b ** 2) % p
                if not disc:
                    continue
                j = 1728 * 4 * a ** 3 * pow(disc, -1, p) % p
                points = 1 + sum(square_multiplicity[(x ** 3 + a * x + b) % p] for x in range(p))
                by_j.setdefault(j, set()).add(points == p + 1)
                models += 1
        check(len(by_j) == p, 'all rational j represented by brute curve census')
        check(all(len(v) == 1 for v in by_j.values()), 'supersingularity invariant across every twist')
        ss = sorted(j for j, values in by_j.items() if values == {True})
        check(ss == row['rational_supersingular_j'], 'all-model census matches producer single-model j census')

        def forms(D):
            # Intentionally loose independent a-bound, then direct reduction tests.
            found = []
            for a in range(1, abs(D) + 1):
                for b in range(-a, a + 1):
                    if (b * b - D) % (4 * a):
                        continue
                    c = (b * b - D) // (4 * a)
                    if a <= c and math.gcd(a, math.gcd(b, c)) == 1 and not (b < 0 and (abs(b) == a or a == c)):
                        found.append((a, b, c))
            return found

        h4 = len(forms(-4 * p))
        h1 = len(forms(-p)) if p % 4 == 3 else 0
        check(h4 == row['h_minus_4p'] and h1 == row['h_minus_p_if_applicable'], 'independent broad-bound primitive reduced forms')
        check(2 * len(ss) == h4 + h1, 'independent class number and all-model curve agreement')
        # Field-class-number version verifies the conversion from order class numbers.
        classical = h4 // 2 if p % 4 == 1 else h1 if p % 8 == 7 else 2 * h1
        check(classical == len(ss), 'published piecewise rational count agrees')
        curve_rows.append({'p': p, 'models': models, 'ss_j': ss, 'h_minus_4p': h4, 'h_minus_p': h1})

    wheel_hostiles = {}
    for d in range(2, 290):
        ends = [trial_prime(30030 - d), trial_prime(30030 + d)]
        if any(ends):
            check(math.gcd(d, 30030) == 1, 'wheel necessity independent trial division')
            if d < 289:
                check(trial_prime(d), 'strict square window independent trial division')
        if d in (17, 19, 23, 289):
            wheel_hostiles[d] = ends
    check(wheel_hostiles == {17: [True, True], 19: [True, False], 23: [False, False], 289: [True, True]}, 'all sharp wheel hostiles')
    return {'certificate_count': len(known), 'minimality_endpoint_classifications': endpoint_count, 'all_curve_models_tested': sum(r['models'] for r in curve_rows), 'curve_rows': curve_rows, 'wheel_hostiles': wheel_hostiles}


def closure(n, edges):
    reach = [[i == j or (i, j) in edges for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    return tuple(tuple(row) for row in reach)


def parts_of(items):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in parts_of(rest):
        yield (frozenset([first]),) + partition
        for i in range(len(partition)):
            yield partition[:i] + (partition[i] | {first},) + partition[i + 1:]


def powerset(items):
    return [frozenset(c) for n in range(len(items) + 1) for c in itertools.combinations(items, n)]


def ideal_set(n, edges):
    return {s for s in powerset(tuple(range(n))) if all(i not in s or j in s for i, j in edges)}


def pareto(points):
    return {p for p in points if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in points)}


def interval_audit():
    # All arbitrary binary relation tables on the three unordered pairs,
    # followed by all product unary restrictions. No clauses are used.
    n = 3
    pairs = list(itertools.combinations(range(n), 2))
    states = list(itertools.product((0, 1), repeat=n))
    families = set()
    for tables in itertools.product(range(16), repeat=len(pairs)):
        binary = frozenset(x for x in states if all(table >> (2 * x[i] + x[j]) & 1 for table, (i, j) in zip(tables, pairs)))
        for unary in itertools.product(range(4), repeat=n):
            families.add(frozenset(x for x in binary if all(t >> x[i] & 1 for i, t in enumerate(unary))))
    endpoint_pairs = 0
    for family in families:
        for x, y in itertools.product(family, repeat=2):
            active = [i for i in range(n) if x[i] != y[i]]
            actual = set()
            for selected in powerset(tuple(range(len(active)))):
                z = list(x)
                for local in selected:
                    z[active[local]] = y[active[local]]
                if tuple(z) in family:
                    actual.add(selected)
            forced_edges = {(i, j) for i in range(len(active)) for j in range(len(active)) if i != j and not any(i in s and j not in s for s in actual)}
            check(actual == ideal_set(len(active), forced_edges), 'all binary-table endpoint intervals equal semantic implication ideals')
            endpoint_pairs += 1

    graph_rows = []
    posets = {}
    for n in range(1, 5):
        arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
        dag_orders = set()
        for bits in range(1 << len(arcs)):
            edges = {e for k, e in enumerate(arcs) if bits >> k & 1}
            reach = closure(n, edges)
            blocks = []
            left = set(range(n))
            while left:
                i = min(left)
                block = frozenset(j for j in left if reach[i][j] and reach[j][i])
                blocks.append(block)
                left -= block
            quotient = {(i, j) for i in range(len(blocks)) for j in range(len(blocks)) if i != j and reach[next(iter(blocks[i]))][next(iter(blocks[j]))]}
            lifted = {frozenset().union(*(blocks[i] for i in s)) for s in ideal_set(len(blocks), quotient)}
            check(lifted == ideal_set(n, edges), 'every labelled digraph equals lifted SCC-poset ideals')
            if len(blocks) == n:
                dag_orders.add(tuple(sorted((i, j) for i in range(n) for j in range(n) if i != j and reach[i][j])))
        posets[n] = dag_orders
        graph_rows.append({'n': n, 'all_directed_graphs': 1 << len(arcs), 'distinct_labelled_posets': len(dag_orders)})

    packet_instances = 0
    neutral_packets = 0
    max_joint_frontier = 0
    for n, orders in posets.items():
        whole = frozenset(range(n))
        partitions = list(parts_of(tuple(range(n))))
        subsets = powerset(tuple(range(n)))
        # Clock totals are at least physical slot totals, hence literally realizable.
        weights = [(i + 1, 2 * (n + 1 - i)) for i in range(n)]
        costs = {s: tuple(sum(weights[i][k] for i in s) for k in (0, 1)) for s in subsets}
        margins = {tuple([0] * n)}
        if n >= 2:
            margins.add(tuple([1, -1] + [0] * (n - 2)))
        if n == 4:
            margins.add((2, -1, 0, -1))
        for edges in sorted(orders):
            relation = ideal_set(n, edges)
            for prefix_q in itertools.product(range(4), repeat=n - 1):
                last = functools.reduce(int.__xor__, prefix_q, 0)
                q = prefix_q + (last,)
                for margin in margins:
                    packet_instances += 1
                    zero = {s for s in subsets if functools.reduce(int.__xor__, (q[i] for i in s), 0) == 0 and sum(margin[i] for i in s) == 0}
                    accepted = relation & zero

                    def path_frontier(covers):
                        values = {frozenset(): {(0, 0)}}
                        for target in sorted(accepted, key=lambda s: (len(s), sorted(s))):
                            if not target:
                                continue
                            options = set()
                            for source in accepted:
                                if not source < target or (covers and any(source < mid < target for mid in accepted)):
                                    continue
                                c = costs[target - source]
                                options.update((max(v[0], c[0]), max(v[1], c[1])) for v in values[source])
                            values[target] = pareto(options)
                        return values[whole]

                    direct = path_frontier(False)
                    check(path_frontier(True) == direct, 'independent cover-only route Pareto equality')
                    max_joint_frontier = max(max_joint_frontier, len(direct))
                    packet_costs, compiled_costs = set(), set()
                    for partition in partitions:
                        if not all(part in zero for part in partition):
                            continue
                        neutral_packets += 1
                        owner = {v: i for i, part in enumerate(partition) for v in part}
                        deps = {(owner[i], owner[j]) for i, j in edges if owner[i] != owner[j]}
                        # Independent criterion via exhaustive packet permutations.
                        legal_orders = []
                        for sequence in itertools.permutations(range(len(partition))):
                            prefix = frozenset()
                            legal = True
                            for index in sequence:
                                prefix |= partition[index]
                                if prefix not in accepted:
                                    legal = False
                                    break
                            if legal:
                                legal_orders.append(sequence)
                        reach = closure(len(partition), deps)
                        acyclic = not any(i != j and reach[i][j] and reach[j][i] for i in range(len(partition)) for j in range(len(partition)))
                        check(bool(legal_orders) == acyclic, 'neutral partition legal ordering iff dependency contraction acyclic')
                        if legal_orders:
                            packet_costs.add(tuple(max(costs[p][k] for p in partition) for k in (0, 1)))
                        pending = set(range(len(partition)))
                        merged = []
                        packet_scc_sizes = []
                        while pending:
                            i = min(pending)
                            cc = {j for j in pending if reach[i][j] and reach[j][i]}
                            merged.append(frozenset().union(*(partition[j] for j in cc)))
                            packet_scc_sizes.append(len(cc))
                            pending -= cc
                        check(all(p in zero for p in merged), 'cycle-compiled packets retain exact defect neutrality')
                        legal_compiled = False
                        for sequence in itertools.permutations(merged):
                            prefix = frozenset()
                            if all((prefix := prefix | p) in accepted for p in sequence):
                                legal_compiled = True
                                break
                        check(legal_compiled, 'cycle-compiled route verified by literal accepted prefixes')
                        compiled_costs.add(tuple(max(costs[p][k] for p in merged) for k in (0, 1)))
                        if all(v == 0 for v in margin) and all(len(p) <= 3 for p in partition):
                            L = max(packet_scc_sizes)
                            check(all(max(costs[p][axis] for p in merged) <= L * 3 * max(v[axis] for v in weights) for axis in (0, 1)), 'conditional L times Davenport times original cost bound')
                    check(pareto(packet_costs) == direct, 'all acyclic packetizations give exact joint Pareto frontier')
                    check(pareto(compiled_costs) == direct, 'all compiled packetizations give exact joint Pareto frontier')
                    accessible = all(any(s - {i} in accepted for i in s) for s in accepted if s)
                    check(accessible == all(q[i] == 0 and margin[i] == 0 for i in range(n)), 'accessibility iff every original block neutral')
                    for left in accepted:
                        for right in accepted:
                            if not left < right:
                                continue
                            diff = right - left
                            relative = {s for s in powerset(tuple(sorted(diff))) if all(i not in s or j not in diff or j in s for i, j in edges)}
                            mids = {left | s for s in relative if s in zero}
                            check(mids == {s for s in accepted if left <= s <= right}, 'relative neutral ideal interval independently exact')

    # A genuine simultaneous-budget tradeoff, with clock totals >= slot totals.
    # The third neutral matching has cyclic contraction and cannot dominate.
    pareto_edges = {(1, 0), (3, 2)}
    special_weights = [(1, 5), (1, 7), (3, 5), (3, 7)]
    special_accepted = {s for s in ideal_set(4, pareto_edges) if len(s) % 2 == 0}
    costs = {s: tuple(sum(special_weights[i][k] for i in s) for k in (0, 1)) for s in powerset(tuple(range(4)))}
    special_frontier = set()
    def walk(current, bottle):
        if len(current) == 4:
            special_frontier.add(bottle)
            return
        for nxt in special_accepted:
            if current < nxt:
                c = costs[nxt - current]
                walk(nxt, (max(bottle[0], c[0]), max(bottle[1], c[1])))
    walk(frozenset(), (0, 0))
    special_frontier = pareto(special_frontier)
    check(special_frontier == {(6, 12), (4, 14)}, 'explicit incomparable joint bottlenecks')
    check((4, 12) not in special_frontier, 'separate coordinate minima are not jointly attainable')

    hostile_rows = []
    for n in range(2, 17):
        prefixes = [0] + [1 if k % 2 else 2 for k in range(1, n)] + [0]
        increments = [prefixes[k] ^ prefixes[k + 1] for k in range(n)]
        check(all(increments), 'hostile chain all logical increments nonzero')
        source = tuple([0, 1] * n)
        physical_defects = []
        for d in increments:
            physical_defects.extend(v for v in (1, 2, 3) if v != d)
        check(all(physical_defects), 'hostile paired lift all physical increments nonzero')
        accepted_count = 0
        for k in range(n + 1):
            target = tuple(1 - v if i < 2 * k else v for i, v in enumerate(source))
            check(sum(target) == n and all(sum(target[2 * i:2 * i + 2]) == 1 for i in range(n)), 'paired lift every prefix retains all exact populations')
            q = functools.reduce(int.__xor__, physical_defects[:2 * k], 0)
            accepted_count += (q == 0)
            check(q == prefixes[k], 'paired physical charge agrees with logical prefix')
        check(accepted_count == 2, 'paired lift exact two-endpoint accepted fibre')
        hostile_rows.append({'logical_blocks': n, 'physical_slots': 2 * n, 'accepted_states': accepted_count, 'original_SCC_size': 2})
    return {'all_binary_table_families_n3': len(families), 'ordered_endpoint_pairs_n3': endpoint_pairs, 'all_graph_rows': graph_rows, 'packet_and_exact_margin_instances': packet_instances, 'neutral_partitions': neutral_packets, 'largest_joint_Pareto_frontier_in_uniform_sweep': max_joint_frontier, 'explicit_joint_Pareto_frontier': sorted(special_frontier), 'explicit_joint_Pareto_hostile_edges': sorted(pareto_edges), 'explicit_joint_Pareto_physical_weights': special_weights, 'paired_hostile_rows': hostile_rows}


def main():
    result = {'status': 'PASS', 'prime': prime_audit(), 'interval': interval_audit(), 'checks': dict(sorted(COUNTS.items()))}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
