#!/usr/bin/env python3
"""Literal terminal preparation and three zero-channel steps.

Standalone adjacency counting: no producer import, no response-table input.
The full statement for arbitrary clocks and depths is proved in PREPARED_FOUR.
"""
from array import array
from pathlib import Path
import argparse
import hashlib
import json
from math import factorial

GATES = 0


def need(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def count(rows, marked=(), all_end_two=False):
    """Forward Hamilton-path subset DP, retaining actual final vertices.

    Each stored count is at most n!, and n<=19 makes uint64 exact. All
    products and final convolutions use Python integers without a bit cap.
    """
    n = len(rows)
    need(0 < n <= 19 and factorial(n) < 2 ** 64, 'exact uint64 path-count universe')
    size = 1 << n
    full = size - 1
    dp = array('Q', [0]) * (size * n)
    paths = array('Q', [0]) * size
    paths[0] = 1
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, size):
        base = mask * n
        vertices, total = mask, 0
        while vertices:
            bit = vertices & -vertices
            vertices -= bit
            v = bit.bit_length() - 1
            value = dp[base + v]
            total += value
            if value:
                choices = rows[v] & (full ^ mask)
                while choices:
                    nxt = choices & -choices
                    choices -= nxt
                    w = nxt.bit_length() - 1
                    dp[(mask | nxt) * n + w] += value
        paths[mask] = total
    focus = tuple(range(n)) if all_end_two else tuple(marked)
    responses = {v: sum(int(dp[mask * n + v]) * int(paths[full ^ mask])
                        for mask in range(1, size)) for v in focus}
    return {'H': int(paths[full]), 'F': responses,
            'P2': sum(int(paths[mask]) * int(paths[full ^ mask]) for mask in range(size)),
            'ends': [int(dp[full * n + v]) for v in range(n)]}


def strong(rows):
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(row == (1 << len(rows)) - 1 for row in reach)


def induced(rows, vertices):
    return tuple(sum(((rows[v] >> w) & 1) << j for j, w in enumerate(vertices))
                 for v in vertices)


def converse(rows):
    return tuple(sum(((rows[j] >> i) & 1) << j for j in range(len(rows)))
                 for i in range(len(rows)))


def birth(rows, incoming):
    n = len(rows)
    full = (1 << n) - 1
    need(0 <= incoming <= full, 'birth cut')
    out = tuple(row | (((incoming >> i) & 1) << n) for i, row in enumerate(rows))
    out += (full ^ incoming,)
    need(induced(out, tuple(range(n))) == rows, 'birth changes old named relation')
    return out


def c2(rows, aux):
    n = len(rows)
    first = birth(rows, ((1 << n) - 1) ^ (1 << aux))
    final = birth(first, (1 << n) - 1)
    need(strong(first) and strong(final), 'strong permanent births')
    return first, final, n


def source(clock):
    need(clock >= 3, 'legal original clock')
    k = clock - 2
    D = (2, 4, 57, 3, 11, 27)
    suffix = tuple(row | 64 for row in D) + (0,)
    ports = 80
    root = k + 7
    suffix_mask = 127 << k
    rows = [sum(1 << j for j in range(i + 1, k)) | suffix_mask for i in range(k)]
    rows += [(row << k) | ((1 << root) if ports >> i & 1 else 0)
             for i, row in enumerate(suffix)]
    rows += [((1 << k) - 1) | (suffix_mask ^ (ports << k))]
    out = tuple(rows)
    u = k + 6
    need(out[u] == 1 << root, 'terminal u has exactly one out-neighbour')
    need(strong(out), 'prepared source strong')
    return out, root, u


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    D = (2, 4, 57, 3, 11, 27)
    ddata = count(D, (4,))
    need(ddata['H'] == 17 and ddata['ends'][4] == 2, 'literal strong-six ending supplier')
    records = []
    for t in (3, 4, 5):
        rows, root, u = source(t)
        entry = count(rows, (root,))
        dual = count(converse(rows), (u,))
        keep = tuple(v for v in range(len(rows)) if v != u)
        root_deleted = keep.index(root)
        one_deleted = count(induced(rows, keep), (root_deleted,))
        two_deleted = count(induced(rows, tuple(v for v in keep if v != root)))
        H, F, S = entry['H'], entry['F'][root], dual['F'][u]
        L, K = one_deleted['F'][root_deleted], two_deleted['P2']
        need(H == 63 * 2 ** t - 29 and F == 63 * 2 ** t + 17, 'prepared H/F family')
        need(S - H == 2, 'terminal start defect two')
        need(one_deleted['ends'][root_deleted] == 2, 'failed start joins')
        need(L == 28 * 2 ** (t - 2) + 17 and K == 224 * 2 ** (t - 2), 'independent L/K counts')
        need(L + K == F, 'joint preparation G=F')
        initial_rows = rows
        first, rows, aux = c2(rows, u)
        out = count(rows, (root,))
        need(out['H'] == H + 2 * S and out['F'][root] == F + 2 * (L + K), 'direct C2 observer response')
        need(out['H'] == 189 * 2 ** t - 83 and out['F'][root] == 189 * 2 ** t + 51, 'prepared digit four')
        need(out['F'][root] - out['H'] == 134, 'transported root defect134')
        ready = count(converse(rows), all_end_two=True)
        gaps = [ready['F'][v] - ready['H'] for v in range(len(rows))]
        need(gaps[aux] == 0 and all(g not in (1, 2) for g in gaps), 'literal all-auxiliary stopping control')
        Hprime, Fprime = out['H'], out['F'][root]
        channel = [{'level': 0, 'rows': rows, 'aux': aux, 'H': Hprime, 'F': Fprime,
                    'start_defects': gaps}]
        for level in range(1, 4):
            old = rows
            intermediate, rows, aux = c2(rows, aux)
            need(induced(rows, tuple(range(len(initial_rows)))) == initial_rows, 'channel retains all original names')
            need(induced(rows, tuple(range(len(old)))) == old, 'channel retains every earlier birth')
            observed = count(rows, (root,))
            need((observed['H'], observed['F'][root]) == (3 ** level * Hprime, 3 ** level * Fprime),
                 'independent zero-channel observer count')
            channel.append({'level': level, 'rows': rows, 'aux': aux,
                            'H': observed['H'], 'F': observed['F'][root]})
        records.append({'original_clock': t, 'root': root, 'terminal_aux': u,
                        'rows': initial_rows, 'H': H, 'F': F, 'S_u': S,
                        'L': L, 'K': K, 'G': L + K, 'first_birth': first,
                        'prepared_C2_and_zero_channel': channel})
        print(json.dumps({'clock': t, 'entry': [H, F, S, L, K],
                          'channel': [[r['H'], r['F']] for r in channel]}), flush=True)
    result = {'scope': 'finite literal controls for proved all-clock terminal preparation and all-depth zero channel',
              'clocks': [3, 4, 5], 'zero_steps_after_entry': 3,
              'all_auxiliary_control_scope': 'only each first C2 output; all-depth exclusion is inherited analytic theorem',
              'gates': GATES, 'records': records}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    args.output.write_text(text)
    print(json.dumps({'gates': GATES, 'bytes': len(text.encode()),
                      'sha256': hashlib.sha256(text.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
