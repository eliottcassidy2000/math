"""Independent literal controls for transitive-slot occupancy responses.

Universe: all labeled tournaments on at most five vertices; all roots with at
least two inports; every maximal all-one singleton run and every all-zero run
in the terminal zero-contact residual suffix. No census or prior fit is used.
"""
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import json
from math import factorial, comb

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def strongly_connected(rows):
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for u in range(n):
        for v in range(n):
            if reach[v] >> u & 1:
                reach[v] |= reach[u]
    return all(mask == (1 << n) - 1 for mask in reach)


def components(rows, vertices):
    reach = {v: (rows[v] & sum(1 << x for x in vertices)) | (1 << v)
             for v in vertices}
    for u in vertices:
        for v in vertices:
            if reach[v] >> u & 1:
                reach[v] |= reach[u]
    left = set(vertices)
    parts = []
    while left:
        v = min(left)
        block = {w for w in left if (reach[v] >> w & 1) and (reach[w] >> v & 1)}
        parts.append(tuple(sorted(block)))
        left -= block
    unsorted_parts = tuple(parts)
    parts.sort(key=lambda block: -sum(rows[block[0]] >> other[0] & 1
                                    for other in unsorted_parts if other != block))
    need(all(rows[parts[i][0]] >> parts[j][0] & 1
             for i in range(len(parts)) for j in range(i + 1, len(parts))),
         'intrinsic residual SCC order')
    return parts


def ordinary_slots(rows, root):
    ports = tuple(v for v in range(len(rows)) if rows[v] >> root & 1)
    residual = [v for v in range(len(rows)) if v != root and v not in ports]
    parts = components(rows, residual)
    last_nonzero = max((i for i, part in enumerate(parts)
                        if any(rows[v] >> a & 1 for v in part for a in ports)), default=-1)
    slots, current, previous_type = [], [], None
    for i, part in enumerate(parts):
        slot_type = None
        if len(part) == 1:
            if all(rows[part[0]] >> a & 1 for a in ports):
                slot_type = 1
            elif i > last_nonzero:
                slot_type = 0
        if slot_type is not None:
            if current and slot_type != previous_type:
                slots.append(tuple(current))
                current = []
            current.append(part[0])
        else:
            if current:
                slots.append(tuple(current))
            current = []
        previous_type = slot_type
    if current:
        slots.append(tuple(current))
    return ports, slots


def inflate(rows, root, slot, length):
    kept = [v for v in range(len(rows)) if v not in slot]
    n, anchor = len(kept), slot[0]
    need(all((rows[u] >> v & 1) == (rows[u] >> anchor & 1)
             for u in kept for v in slot), 'slot is an actual homogeneous module')
    answer = [0] * (n + length)
    for i, u in enumerate(kept):
        for j, v in enumerate(kept):
            if rows[u] >> v & 1:
                answer[i] |= 1 << j
        for j in range(n, n + length):
            if rows[u] >> anchor & 1:
                answer[i] |= 1 << j
            else:
                answer[j] |= 1 << i
    for i in range(n, n + length):
        for j in range(i + 1, n + length):
            answer[i] |= 1 << j
    return tuple(answer), kept.index(root)


@lru_cache(maxsize=30000)
def literal_profile(rows, root):
    n = len(rows)
    full = (1 << n) - 1
    ends = [[0] * n for _ in range(full + 1)]
    H = [0] * (full + 1)
    H[0] = 1
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, count in enumerate(ends[mask]):
            if not count:
                continue
            extension = rows[v] & (full ^ mask)
            while extension:
                bit = extension & -extension
                extension ^= bit
                ends[mask | bit][bit.bit_length() - 1] += count
        H[mask] = sum(ends[mask])
    F = sum(ends[mask][root] * H[full ^ mask]
            for mask in range(1, full + 1) if mask >> root & 1)
    return H[full], F - H[full]


def induce(rows, vertices):
    return tuple(sum((rows[u] >> v & 1) << j for j, v in enumerate(vertices))
                 for u in vertices)


def hamilton(rows):
    return literal_profile(tuple(rows), 0)[0] if rows else 1


def q_two(rows):
    n = len(rows)
    return sum(hamilton(induce(rows, [i for i in range(n) if mask >> i & 1]))
               * hamilton(induce(rows, [i for i in range(n) if not (mask >> i & 1)]))
               for mask in range(1 << n))


def rooted_code(rows, root):
    n = len(rows)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return min(sum((rows[order[i]] >> order[j] & 1) << b
                   for b, (i, j) in enumerate(pairs))
               for rest in permutations([v for v in range(n) if v != root])
               for order in [(root,) + rest])


def ordinal_join(blocks):
    n = sum(len(block) for block in blocks)
    result, offset = [], 0
    for block in blocks:
        later = ((1 << n) - 1) ^ ((1 << (offset + len(block))) - 1)
        result.extend((row << offset) | later for row in block)
        offset += len(block)
    return tuple(result)


def transitive(n):
    return tuple(((1 << n) - 1) ^ ((1 << (v + 1)) - 1) for v in range(n))


def growth_witness(v, a, c, m, gap_after=False):
    blocks = [(2, 4, 1)] * a + [(2, 4, 9, 3)] * c
    blocks.insert(len(blocks) if gap_after else 0, transitive(m))
    P = ordinal_join(blocks)
    W = transitive(v + 2)
    result = list(ordinal_join([P, W])) + [0]
    root, first, last = len(result) - 1, len(P), len(P) + len(W) - 1
    for u in range(root):
        if u in (first, last):
            result[u] |= 1 << root
        else:
            result[root] |= 1 << u
    return tuple(result), root


def surjections(k, j):
    if j == 0:
        return int(k == 0)
    table = [1] + [0] * j
    for _ in range(k):
        table = [0] + [r * (table[r] + table[r - 1]) for r in range(1, j + 1)]
    return table[j]


def recover(values):
    result = []
    for k, value in enumerate(values, 1):
        remainder = value - sum(c * surjections(k, j)
                                for j, c in enumerate(result, 1))
        need(remainder % factorial(k) == 0, 'integral occupancy coefficient')
        result.append(remainder // factorial(k))
    return result


def main():
    slots_checked = literal_queries = roots_checked = strong_graphs = 0
    strong_roots = weak_slots = 0
    source_graphs = 0
    constant_defect_slots = 0
    largest_H_run_in_constant_slot = 0
    smallest_nonconstant_example = None
    histogram = {}
    rooted_classes = {(n, b): set() for n in range(1, 6) for b in range(2, 7)}
    scc_factor_checks = preceding_core_bound_checks = 0
    for n in range(1, 6):
        pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
        for mask in range(1 << len(pairs)):
            rows = [0] * n
            for i, (u, v) in enumerate(pairs):
                if mask >> i & 1:
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
            source_graphs += 1
            strong = strongly_connected(rows)
            strong_graphs += int(strong)
            for root in range(n):
                ports, slots = ordinary_slots(rows, root)
                d = len(ports)
                if d < 2:
                    continue
                roots_checked += 1
                strong_roots += int(strong)
                actual_H, actual_b = literal_profile(tuple(rows), root)
                if 2 <= actual_b <= 6:
                    rooted_classes[n, actual_b].add(rooted_code(rows, root))
                parts = components(rows, range(n))
                index = next(i for i, part in enumerate(parts) if root in part)
                A = induce(rows, [v for part in parts[:index] for v in part])
                core_vertices = list(parts[index])
                C = induce(rows, core_vertices)
                B = induce(rows, [v for part in parts[index + 1:] for v in part])
                HC, bC = literal_profile(C, core_vertices.index(root))
                HA, HB, QA = hamilton(A), hamilton(B), q_two(A)
                need(actual_H == HA * HC * HB, 'root-SCC Hamilton factorization')
                need(actual_H + actual_b == QA * (HC + bC) * HB,
                     'root-SCC F factorization')
                if A:
                    need(actual_H <= actual_b, 'nonempty preceding region bounds H')
                    preceding_core_bound_checks += 1
                scc_factor_checks += 1
                for slot in slots:
                    profiles = [literal_profile(*inflate(rows, root, slot, k))
                                for k in range(1, d + 3)]
                    literal_queries += d + 2
                    Hcoeff = recover([p[0] for p in profiles])
                    bcoeff = recover([p[1] for p in profiles])
                    need(all(c >= 0 for c in Hcoeff), 'positive Hamilton occupancy')
                    need(all(c >= 0 for c in bcoeff), 'positive defect occupancy')
                    need(all(c == 0 for c in Hcoeff[d + 1:]), 'Hamilton at most d+1 runs')
                    need(all(c == 0 for c in bcoeff[d:]), 'defect at most d runs')
                    constant = all(c == 0 for c in bcoeff[1:])
                    if constant:
                        constant_defect_slots += 1
                        largest_H_run_in_constant_slot = max(
                            largest_H_run_in_constant_slot,
                            max(j for j, c in enumerate(Hcoeff, 1) if c))
                    elif smallest_nonconstant_example is None:
                        smallest_nonconstant_example = {
                            'rows': rows, 'root': root, 'slot': slot,
                            'Hamilton_occupancy': Hcoeff, 'defect_occupancy': bcoeff}
                    for k, (_, b) in enumerate(profiles, 1):
                        threshold = (b + 2).bit_length()
                        if k >= threshold:
                            need(constant, 'explicit threshold forces one-run defect')
                    slots_checked += 1
                    weak_slots += int(not strong)
                    key = f'd{d}_strong{int(strong)}_constant{int(constant)}'
                    histogram[key] = histogram.get(key, 0) + 1
    surjection_controls = 0
    for k in range(2, 81):
        for j in range(2, k + 1):
            need(surjections(k, j) >= (1 << k) - 2, 'surjection lower bound')
            surjection_controls += 1
    weak_series = {
        2: ({4: 1, 5: 2, 6: 2}, 2),
        3: ({3: 1, 4: 1}, 1),
        4: ({4: 1, 6: 5, 7: 10}, 2),
        5: ({5: 3}, 1),
        6: ({5: 4, 6: -2, 7: 12, 8: 12, 9: -16}, 3),
    }
    census = []
    for (n, b), classes in sorted(rooted_classes.items()):
        numerator, exponent = weak_series[b]
        predicted = sum(c * comb(n - e + exponent - 1, exponent - 1)
                        for e, c in numerator.items() if n >= e)
        need(len(classes) == predicted, 'independent weak rooted class census')
        census.append({'order': n, 'defect': b, 'ordinary_rooted_classes': len(classes)})
    need(q_two((2, 4, 1)) == 12 and q_two((2, 4, 9, 3)) == 32,
         'growth-family core responses')
    growth_controls = []
    for v in range(1, 7):
        for a in range(2):
            for c in range(2):
                for m in range(3):
                    if a + c + m == 0:
                        continue
                    n = m + 3 * a + 4 * c + v + 3
                    if n > 12:
                        continue
                    expected_H = ((1 << (v + 1)) + 1) * 12**a * 32**c * (1 << m) + (1 - (1 << v)) * 3**a * 5**c
                    expected_b = (1 << v) * 3**a * 5**c
                    for after in (False, True):
                        rows, root = growth_witness(v, a, c, m, after)
                        need(strongly_connected(rows), 'growth witness strongness')
                        need(literal_profile(rows, root) == (expected_H, expected_b),
                             'all-parameter growth witness literal control')
                        growth_controls.append({'v': v, 'a': a, 'c': c, 'm': m,
                                                'gap_after': after, 'H': expected_H, 'b': expected_b})
    report = {
        'status': 'FINITE-EXACT independent literal controls; universal proof in RESPONSE_FAMILIES.md',
        'universe': 'All labeled tournaments of orders 1 through 5, every root with d>=2, every maximal all-one singleton slot and every all-zero slot in the terminal residual zero-contact suffix',
        'source_labeled_graphs': source_graphs,
        'strong_labeled_graphs': strong_graphs,
        'branching_rooted_graphs': roots_checked,
        'strong_branching_rooted_graphs': strong_roots,
        'slots_checked': slots_checked,
        'weak_source_slots_checked': weak_slots,
        'literal_inflation_queries': literal_queries,
        'positive_controls': 'nonnegative integral occupancy coefficients; exact d/d+1 run bounds; explicit threshold implication',
        'constant_defect_slots': constant_defect_slots,
        'largest_H_run_in_constant_slot_observed_only': largest_H_run_in_constant_slot,
        'smallest_nonconstant_hostile_control': smallest_nonconstant_example,
        'histogram': histogram,
        'surjection_lower_bound_controls': surjection_controls,
        'SCC_root_factorization_checks': scc_factor_checks,
        'nonempty_preceding_region_H_bound_checks': preceding_core_bound_checks,
        'weak_rooted_census': census,
        'sharp_growth_witness_literal_queries': len(growth_controls),
        'sharp_growth_witness_controls': growth_controls,
        'no_universal_two_run_H_claim': True,
    }
    (HERE / 'response_families.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({key: report[key] for key in (
        'status', 'branching_rooted_graphs', 'slots_checked',
        'literal_inflation_queries', 'SCC_root_factorization_checks',
        'sharp_growth_witness_literal_queries')}, sort_keys=True))


if __name__ == '__main__':
    main()
