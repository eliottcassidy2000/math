#!/usr/bin/env python3
"""Portable exact checks for the feedback-arc interval-coverage budget gate.
No computation below is a proof of the all-order spectrum conjecture.
"""
from itertools import combinations, permutations
from math import comb
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def ilog(x, base):
    if type(x) is not int or x < 1 or type(base) is not int or base < 2:
        raise ValueError('positive integer x and integer base >= 2 required')
    k = 0
    while x >= base:
        x //= base
        k += 1
    return k


def strong_floor(n):
    if type(n) is not int or n < 3:
        raise ValueError('strong tournament order must be at least 3')
    k, r = divmod(n, 3)
    return (3 * 5 ** (k - 1), 5 ** k, 9 * 5 ** (k - 1))[r]


def capped_core_bound(x, b, witness=False):
    """Exact max core order at H<=x and <=b nonsingleton SCCs, NOT fa<=b."""
    if type(x) is not int or x < 1 or type(b) is not int or b < 0:
        raise ValueError('integer x>=1 and integer b>=0 required')
    best, arg = 0, (0, 0, 0)
    for r in range(1, min(b, ilog(x, 3)) + 1):
        power3 = 1
        for u in range(2 * r + 1):
            if power3 > x:
                break
            v = ilog(x // power3, 5)
            if v >= max(0, r - u):
                n = r + 2 * u + 3 * v
                if (n, r, u, v) > (best, *arg):
                    best, arg = n, (r, u, v)
            power3 *= 3
    return (best, arg) if witness else best


def uncapped_core_bound(x):
    # Inherited S133 exact K(x), independently encoded.
    return max((4 * d + 3 * ilog(x // 5 ** d, 3)
                for d in range(3) if 5 ** d <= x), default=0)


def arc_value_bound(x, b):
    n = capped_core_bound(x, b)
    m = n * (n - 1) // 2
    return sum(comb(m, j) for j in range(min(b, m) + 1))


def rows_from_backarcs(n, backwards):
    rows = [((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n)]
    seen = set()
    for i, j in backwards:
        if not 0 <= i < j < n or (i, j) in seen:
            raise ValueError('distinct backarc pairs 0<=i<j<n required')
        seen.add((i, j))
        rows[i] ^= 1 << j
        rows[j] ^= 1 << i
    return rows


def hamilton(rows):
    n = len(rows)
    if not n:
        return 1
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = size - 1
    for s in range(1, size):
        for v, count in enumerate(dp[s]):
            if count:
                z = rows[v] & (full ^ s)
                while z:
                    bit = z & -z
                    z -= bit
                    dp[s | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def literal_hamilton(rows):
    return sum(all(rows[p[i]] >> p[i + 1] & 1 for i in range(len(p) - 1))
               for p in permutations(range(len(rows))))


def fas(rows):
    n = len(rows)
    dp = [0] + [n * n] * ((1 << n) - 1)
    for s in range(1, 1 << n):
        z = s
        while z:
            bit = z & -z
            z -= bit
            v = bit.bit_length() - 1
            old = s ^ bit
            dp[s] = min(dp[s], dp[old] + (rows[v] & old).bit_count())
    return dp[-1]


def literal_fas(rows):
    return min((sum(rows[p[j]] >> p[i] & 1 for j in range(len(p)) for i in range(j))
                for p in permutations(range(len(rows)))), default=0)


def sccs(rows):
    n = len(rows)
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    todo = set(range(n))
    out = []
    while todo:
        v = min(todo)
        part = {u for u in todo if reach[v] >> u & 1 and reach[u] >> v & 1}
        todo -= part
        out.append(sorted(part))
    return out


def induced(rows, selected):
    return [sum((1 << j) for j, v in enumerate(selected) if rows[u] >> v & 1)
            for u in selected]


def pad_sources(rows, total):
    n = len(rows)
    if total < n:
        raise ValueError('padding total below current order')
    pad = total - n
    return [((1 << total) - 1) ^ ((1 << (i + 1)) - 1) for i in range(pad)] + [r << pad for r in rows]


def floor_witness(n):
    # Actual Moon interval-chain construction, not synthetic H-coordinates.
    k, r = divmod(n, 3)
    first = (2, 3, 4)[r]
    spans = [first] + [3] * (k - 1)
    ends, pos = [], 0
    for span in spans:
        ends.append((pos, pos + span))
        pos += span
    need(pos + 1 == n, 'floor witness order')
    return rows_from_backarcs(n, ends)


def capped_bound_independent(max_x=5000, max_b=6):
    # Direct multiset composition of all exact strong-order floors.
    blocks = []
    n = 3
    while strong_floor(n) <= max_x:
        blocks.append((n, strong_floor(n)))
        n += 1
    states = {1: 0}
    by_cap = [{1: 0}]
    for _ in range(max_b):
        nxt = {}
        for value, order in states.items():
            for n, h in blocks:
                if value * h <= max_x:
                    nxt[value * h] = max(nxt.get(value * h, 0), order + n)
        states = nxt
        merged = dict(by_cap[-1])
        for value, order in nxt.items():
            merged[value] = max(merged.get(value, 0), order)
        by_cap.append(merged)
    cases = 0
    for b, states in enumerate(by_cap):
        best = 0
        for x in range(1, max_x + 1):
            best = max(best, states.get(x, 0))
            need(capped_core_bound(x, b) == best, 'capped core formula versus multiset enumeration')
            need(best <= uncapped_core_bound(x), 'capped bound exceeds uncapped kernel')
            cases += 1
    return cases


def odd_cycles(rows):
    cycles = []
    n = len(rows)
    for k in range(3, n + 1, 2):
        for support in combinations(range(n), k):
            v = support[0]
            for tail in permutations(support[1:]):
                word = (v,) + tail
                if all(rows[word[i]] >> word[(i + 1) % k] & 1 for i in range(k)):
                    cycles.append(word)
    return cycles


def main():
    malformed = [(0, 0), (2, -1), (3.0, 2), (True, 2)]
    for x, b in malformed:
        try:
            capped_core_bound(x, b)
        except ValueError:
            pass
        else:
            raise RuntimeError('malformed budget accepted')
    bound_cases = capped_bound_independent()
    full_controls = 0
    host = None
    for n in range(6):
        edges = list(combinations(range(n), 2))
        for code in range(1 << len(edges)):
            rows = rows_from_backarcs(n, [edge for j, edge in enumerate(edges) if code >> j & 1])
            h = hamilton(rows)
            b = fas(rows)
            need(h == literal_hamilton(rows), 'Hamilton subset DP versus literal permutations')
            need(b == literal_fas(rows), 'FAS subset DP versus all orders')
            parts = [s for s in sccs(rows) if len(s) > 1]
            selected = sorted(u for s in parts for u in s)
            core = induced(rows, selected)
            need(hamilton(core) == h, 'singleton SCC deletion preserves H')
            need(fas(core) <= b and len(parts) <= b, 'SCC and induced arc budget')
            need(len(core) <= capped_core_bound(h, b), 'actual tournament capped bound')
            if n <= 4 or code % 64 == 0:
                padded = pad_sources(core, capped_core_bound(h, b))
                need(hamilton(padded) == h and fas(padded) <= b, 'source padding bound carrier')
            if n == 5 and h == 15 and host is None:
                cyc = odd_cycles(rows)
                if len(cyc) == 7 and len(set(frozenset(c) for c in cyc)) == 5:
                    host = {'rows': rows, 'H': h, 'literal_odd_cycles': len(cyc),
                            'distinct_odd_cycle_supports': 5, 'support_collapsed_I_at_2': 11}
            full_controls += 1
    need(host is not None, 'literal-cycle/support-collapse hostile found')
    floor_controls = []
    for n in range(3, 16):
        rows = floor_witness(n)
        need(len(sccs(rows)) == 1, 'floor witness is strong')
        h = hamilton(rows)
        need(h == strong_floor(n), 'actual strong-floor witness')
        floor_controls.append([n, h])
    six = []
    for n in range(9, 16):
        rows = rows_from_backarcs(n, [(0, 6), (1, 4), (2, 5), (0, n - 1), (4, n - 1), (5, n - 1)])
        h = hamilton(rows)
        need(h == 2 ** (n + 1) - 9, 'inherited six-arc deficit-nine family')
        need(fas(rows) == 6, 'inherited exact six-arc family budget')
        six.append([n, h, 6])
    # The transitive twin block is not deletable merely because it is a module.
    twins = []
    for t in range(1, 9):
        rows = rows_from_backarcs(t + 2, [(0, t + 1)])
        h = hamilton(rows)
        need(h == 2 ** t + 1, 'twin module hostile')
        twins.append([t, h])
    gates = []
    for b in (1, 2, 3, 4, 5, 6, 10):
        k = 2
        while arc_value_bound(2 ** k - 1, b) >= 2 ** (k - 2):
            k += 1
        x = 2 ** k - 1
        n, arg = capped_core_bound(x, b, True)
        bound = arc_value_bound(x, b)
        gates.append({'arc_budget': b, 'first_dyadic_exponent_for_this_bound': k,
                      'core_order_bound': n, 'maximizing_r_u_v': arg,
                      'at_most_low_arc_values_below_upper_endpoint': bound,
                      'odd_targets_in_open_dyadic_interval': 2 ** (k - 2),
                      'certified_targets_without_low_arc_realizations_at_least': 2 ** (k - 2) - bound})
    large = []
    for k in (128, 256, 512, 1024):
        n = capped_core_bound(2 ** k - 1, 6)
        bound = arc_value_bound(2 ** k - 1, 6)
        large.append({'dyadic_exponent': k, 'arc_budget': 6, 'core_order_bound': n,
                      'value_upper_bound': bound, 'odd_targets': 2 ** (k - 2)})
    print(json.dumps({'status': 'passed',
                      'scope': 'Exact finite checks; all-order arc-entropy gate is proved in THEORY.md. No universal spectrum closure.',
                      'complete_labelled_universe': 'every tournament of orders 0 through 5',
                      'actual_tournament_controls': full_controls,
                      'capped_core_multiset_comparisons': bound_cases,
                      'strong_floor_controls': floor_controls,
                      'six_arc_deficit_nine_controls': six,
                      'twin_contraction_hostile_t_H': twins,
                      'literal_cycle_hostile': host,
                      'finite_dyadic_count_gates': gates,
                      'larger_six_arc_gates': large}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
