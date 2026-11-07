#!/usr/bin/env python3
"""Standalone exact controls. Python standard library only; stdout is deterministic."""
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb, prod
import hashlib
import json


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def tournament(n, mask):
    out = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if mask >> bit & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


def scores(out):
    return [row.bit_count() for row in out]


def induced(out, vertices):
    return [sum(((out[u] >> v) & 1) << j for j, v in enumerate(vertices))
            for u in vertices]


def hp(out):
    n = len(out)
    if not n:
        return 1
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            todo = out[v] & (full ^ mask) if count else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def hc(out):
    n = len(out)
    if n < 3:
        return 0
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    dp[1][0] = 1
    for mask in range(1, 1 << n, 2):
        for v, count in enumerate(dp[mask]):
            todo = out[v] & (full ^ mask) if count else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(c for v, c in enumerate(dp[-1]) if out[v] & 1)


def brute_hp(out):
    return sum(all(out[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(out))))


def brute_hc(out):
    if len(out) < 3:
        return 0
    return sum(all(out[u] >> v & 1 for u, v in zip(p, p[1:] + p[:1]))
               for tail in permutations(range(1, len(out))) for p in [(0,) + tail])


def gate(out):
    n = len(out)
    s = scores(out)
    return [(u, v) for u in range(n) for v in range(n)
            if out[u] >> v & 1 and s[u] == 1 and s[v] == n - 2]


def closure(out):
    n = len(out)
    return [row | (1 << n) for row in out] + [1 << (n + 1), (1 << n) - 1]


def switch(out, v):
    n = len(out)
    ans = closure(out)
    for i, j in [(n, n + 1), (n + 1, v), (v, n)]:
        check(ans[i] >> j & 1, 'directed triangle before reversal')
        ans[i] ^= 1 << j
        ans[j] ^= 1 << i
    return ans


def scc(out):
    n = len(out)
    reach = [row | (1 << i) for i, row in enumerate(out)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    groups = []
    unseen = set(range(n))
    while unseen:
        i = min(unseen)
        block = {j for j in unseen if reach[i] >> j & 1 and reach[j] >> i & 1}
        groups.append(sorted(block))
        unseen -= block
    # Every pair of SCCs in a tournament is comparable in one direction.
    block_reps = [b[0] for b in groups]
    groups.sort(key=lambda b: -sum(out[b[0]] >> j & 1 for j in block_reps if j != b[0]))
    return groups


def slack(s, mask):
    k = mask.bit_count()
    return sum(x for i, x in enumerate(s) if mask >> i & 1) - comb(k, 2)


def check_faces(out):
    n = len(out)
    s = scores(out)
    blocks = scc(out)
    actual_tight = []
    for mask in range(1, (1 << n) - 1):
        cross = sum((out[i] & (((1 << n) - 1) ^ mask)).bit_count()
                    for i in range(n) if mask >> i & 1)
        check(slack(s, mask) == cross, 'cut slack')
        if not cross:
            actual_tight.append(mask)
    expected = []
    mask = 0
    for block in reversed(blocks[1:]):
        mask |= sum(1 << i for i in block)
        expected.append(mask)
    check(sorted(actual_tight) == sorted(expected), 'SCC tight subsets')
    check(hp(out) == prod(hp(induced(out, b)) for b in blocks), 'SCC H product')
    return blocks


def check_gate_slacks(y):
    n = len(y)
    x = list(v + 1 for v in y) + [1, n]
    vals = []
    for mask in range(1 << n):
        k = mask.bit_count()
        old = slack(y, mask)
        expected = [old + k, old + 1, old + n, old + n - k]
        for flags in range(4):
            target = mask | ((flags & 1) << n) | ((flags >> 1) << (n + 1))
            check(slack(x, target) == expected[flags], 'four affine slack formulas')
            if target and target != (1 << (n + 2)) - 1:
                vals.append(slack(x, target))
    return min(vals)



def four_switch(out, u, v):
    n = len(out)
    check(out[v] >> u & 1, 'four-switch core arc')
    ans = closure(out)
    for i, j in [(n, n + 1), (n + 1, v), (v, u), (u, n)]:
        check(ans[i] >> j & 1, 'directed four-cycle before reversal')
        ans[i] ^= 1 << j
        ans[j] ^= 1 << i
    return ans


def two_path_formula(out, u, v):
    # Independent direct path-word enumeration, not Hamilton-cycle DP.
    n = len(out)
    core = out[:]
    core[v] ^= 1 << u
    core[u] ^= 1 << v
    total = 0
    terms = []
    for mask in range(1 << n):
        if not (mask >> u & 1 and mask >> v & 1):
            continue
        middle = [w for w in range(n) if mask >> w & 1 and w not in (u, v)]
        first = sum(all(core[i] >> j & 1 for i, j in zip(p, p[1:]))
                    for mid in permutations(middle) for p in [(u,) + mid + (v,)])
        second = brute_hp(induced(core, [w for w in range(n) if not(mask >> w & 1)]))
        total += first * second
        if first * second:
            terms.append({'mask': mask, 'endpoint_paths': first, 'complement_paths': second})
    return total, terms

def main():
    report = {'status': 'PASS', 'scope': 'All labelled cores orders 1..5 and triangle choices; all cores orders 2..5 and directed four-cycle choices; two explicit HC=7,21 score-fiber pairs; exact degree-slice vertices through n=6.'}
    rows = []
    digest = hashlib.sha256()
    target_score = (2, 2, 2, 1, 3)
    fiber = []
    direct_cycles = 0
    for n in range(1, 6):
        by_score = defaultdict(lambda: [set(), set()])
        switches = mixed = 0
        total = 1 << comb(n, 2)
        for mask in range(total):
            t = tournament(n, mask)
            h = hp(t)
            check(h == brute_hp(t), 'independent path scan')
            blocks = check_faces(t)
            by_score[tuple(scores(t))][0].add(bool(gate(t)))
            by_score[tuple(scores(t))][1].add(h)
            g = closure(t)
            check(hc(g) == h, 'inherited closure cycle count')
            check(check_gate_slacks(scores(t)) == 1, 'gate exact interior margin')
            if n == 5 and tuple(scores(t)) == target_score:
                fiber.append({'mask': mask, 'rows': t, 'H': h, 'HC': hc(t), 'gates': gate(t)})
            for v in range(n):
                u = switch(t, v)
                delete_h = hp(induced(t, [w for w in range(n) if w != v]))
                check(scores(g) == scores(u), 'score-preserving triangle')
                check(hc(u) == delete_h, 'switched deletion count')
                if n >= 2:
                    check(len(scc(u)) == 1, 'switched strong')
                if n >= 3:
                    internal = 0 < t[v].bit_count() < n - 1
                    check((not gate(u)) == internal, 'exact gate-loss boundary')
                    mixed += internal
                if n <= 4 or (n == 5 and mask in (0, 341, 682, 1023)):
                    check(hc(u) == brute_hc(u), 'independent cycle scan')
                    direct_cycles += 1
                digest.update(json.dumps([n, mask, v, h, delete_h, scores(u), gate(u)], separators=(',', ':')).encode() + b'\n')
                switches += 1
        if n <= 4:
            check(all(len(z[0]) == 1 for z in by_score.values()), 'no lower-order gate ambiguity')
            check(all(len(z[1]) == 1 for z in by_score.values()), 'no lower-order H ambiguity')
        rows.append({'n': n, 'labelled_cores': total, 'chosen_vertex_switches': switches,
                     'nongated_switches': mixed, 'score_vectors': len(by_score),
                     'score_fibers_mixing_gate_status': sum(len(z[0]) > 1 for z in by_score.values()),
                     'score_fibers_mixing_H': sum(len(z[1]) > 1 for z in by_score.values())})
    four_rows = []
    for n in range(2, 6):
        count = 0
        h_support = set()
        for mask in range(1 << comb(n, 2)):
            t = tournament(n, mask)
            g = closure(t)
            for v in range(n):
                for u in range(n):
                    if not (t[v] >> u & 1):
                        continue
                    z = four_switch(t, u, v)
                    check(scores(z) == scores(g), 'four-switch equal score')
                    check(len(scc(z)) == 1, 'four-switch strong')
                    expected_gate = t[u].bit_count() == 0 and t[v].bit_count() == n - 1
                    check(bool(gate(z)) == expected_gate, 'four-switch gate classification')
                    direct, _ = two_path_formula(t, u, v)
                    check(hc(z) == direct, 'four-switch two-path formula')
                    h_support.add(direct)
                    count += 1
        four_rows.append({'core_order': n, 'four_cycle_switches': count, 'HC_support': sorted(h_support)})
    report['four_switch_exhaustive'] = four_rows
    gap_pairs = []
    for t, u, v, expected in [([24, 1, 3, 6, 14], 0, 1, 7),
                               ([30, 44, 56, 48, 34, 1], 1, 4, 21)]:
        g = closure(t)
        z = four_switch(t, u, v)
        direct, terms = two_path_formula(t, u, v)
        check(scores(g) == scores(z), 'gap pair equal labeled score')
        check(len(scc(z)) == 1 and not gate(z), 'gap pair strong nongated')
        check(hc(z) == brute_hc(z) == direct == expected, 'gap pair three independent cycle counts')
        check(hc(g) == brute_hc(g) == hp(t), 'gap pair gated count')
        check(hp(z) == brute_hp(z) and hp(g) == brute_hp(g), 'gap pair independent path counts')
        gap_pairs.append({'core': t, 'u': u, 'v': v,
            'scores': scores(z), 'two_path_terms': terms,
            'gated': {'rows': g, 'HC': hc(g), 'H': hp(g), 'gates': gate(g)},
            'ungated': {'rows': z, 'HC': hc(z), 'H': hp(z), 'gates': gate(z)}})
    report['seven_and_twenty_one_same_score_pairs'] = gap_pairs
    report['exhaustive'] = rows
    report['switch_digest_sha256'] = digest.hexdigest()
    report['independent_direct_cycle_scans'] = direct_cycles
    t = [2, 4, 1]
    g = closure(t)
    u = switch(t, 0)
    check(hp(g) == brute_hp(g) == 15 and hp(u) == brute_hp(u) == 11, 'minimal H hostile')
    check(hc(g) == 3 and hc(u) == 1, 'minimal HC hostile')
    report['smallest_hostile'] = {'core': t, 'chosen_vertex': 0,
        'G': {'rows': g, 'scores': scores(g), 'H': hp(g), 'HC': hc(g), 'gates': gate(g)},
        'U': {'rows': u, 'scores': scores(u), 'H': hp(u), 'HC': hc(u), 'gates': gate(u)}}
    report['complete_labelled_score_fiber'] = {'scores': target_score, 'realizers': fiber,
        'H_support': sorted({r['H'] for r in fiber}), 'HC_support': sorted({r['HC'] for r in fiber})}
    slice_rows = []
    for n in range(2, 7):
        extremal = [-1] + list(range(1, n - 1)) + [n]
        checked = 0
        for y in permutations(extremal):
            check(sum(y) == comb(n, 2), 'enlarged slice total')
            check(min(slack(y, m) for m in range(1, (1 << n) - 1)) == -1, 'relaxed slice slack')
            check(check_gate_slacks(y) == 0, 'slice boundary control')
            checked += 1
        slice_rows.append({'n': n, 'vertex_vector': extremal, 'vertices_checked': checked,
                           'core_min_slack': -1, 'ambient_min_slack': 0})
    report['full_slice_extreme_controls'] = slice_rows
    # Scalar doubling is impossible for tournament H by oddness; order doubling is not count doubling.
    c3 = [2, 4, 1]
    joined = [r | (7 << 3) for r in c3] + [r << 3 for r in c3]
    report['order_doubling_hostile'] = {'rows': joined, 'order': 6, 'H': hp(joined),
        'SCCs': scc(joined), 'face_dimension': 6 - len(scc(joined)), 'odd_index': (hp(joined) - 1) // 2}
    check(hp(joined) == 9, 'doubling vertices squares H under repeated order join')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
