#!/usr/bin/env python3
"""Independent H51/H25 induced audit: ending DP and literal path DFS."""
from pathlib import Path
import hashlib
import json
import sys
from collections import Counter

BASE = Path(__file__).resolve().parent
BANK = BASE.parent / 'session194-weak-parent-certificate/inherited_tournaments_through8.txt'

def need(b, s):
    if not b:
        raise RuntimeError(s)

def ending_dp(rows):
    n = len(rows)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, size):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            choices = rows[last] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += val
    return [1] + [sum(dp[m]) for m in range(1, size)]

def literal_paths(rows):
    counts = [0] * (1 << len(rows))
    counts[0] = 1
    def visit(mask, last):
        counts[mask] += 1
        choices = rows[last] & ~mask
        while choices:
            bit = choices & -choices
            choices -= bit
            visit(mask | bit, bit.bit_length() - 1)
    for v in range(len(rows)):
        visit(1 << v, v)
    return counts

def sccs(rows):
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for k in range(n):
        for v in range(n):
            if (reach[v] >> k) & 1:
                reach[v] |= reach[k]
    unused = (1 << n) - 1
    result = []
    while unused:
        bit = unused & -unused
        v = bit.bit_length() - 1
        comp = sum(1 << w for w in range(n)
                   if ((reach[v] >> w) & 1) and ((reach[w] >> v) & 1))
        need(comp & unused == comp, 'overlapping SCCs')
        unused ^= comp
        result.append(comp)
    return result

def main():
    raw = BANK.read_bytes()
    lines = raw.decode().splitlines()
    need(int(lines[0]) == len(lines) - 1 == 7412, 'bank length')
    orders = Counter()
    count51 = Counter()
    strong51 = Counter()
    count25 = Counter()
    inventory = []
    witnesses = {}
    direct_graphs = 0
    direct_subsets = 0
    for ix, line in enumerate(lines[1:], 1):
        items = list(map(int, line.split()))
        n, rows = items[0], items[1:]
        need(1 <= n <= 8 and len(rows) == n, 'bad order')
        for i, row in enumerate(rows):
            need(0 <= row < (1 << n) and not row & (1 << i), 'bad row')
            for j in range(i):
                need(((row >> j) & 1) + ((rows[j] >> i) & 1) == 1, 'bad pair')
        orders[n] += 1
        counts = ending_dp(rows)
        h = counts[-1]
        if h == 25:
            count25[n] += 1
            witnesses.setdefault('25', {'order': n, 'rows': rows})
        if h != 51:
            continue
        count51[n] += 1
        comps = sccs(rows)
        isstrong = len(comps) == 1
        if isstrong:
            strong51[n] += 1
        need(literal_paths(rows) == counts, 'DP versus literal path disagreement')
        direct_graphs += 1
        direct_subsets += len(counts)
        hits = [m for m, v in enumerate(counts) if v == 25]
        need(not hits, 'H25 induced subset found')
        compcounts = [counts[m] for m in comps]
        product = 1
        for value in compcounts:
            product *= value
        need(product == h, 'SCC multiplicativity')
        inventory.append({'bank_row': ix, 'order': n, 'rows': rows,
                          'strong': isstrong, 'scc_masks': comps,
                          'scc_H': compcounts, 'H25_subset_masks': hits,
                          'induced_H_histogram': dict(sorted(Counter(counts).items()))})
    need(dict(orders) == {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}, 'order inventory')
    need(witnesses.get('25') is not None, 'positive H25 control missing')
    need(sum(strong51.values()) > 0, 'positive strong H51 control missing')
    need(5 ** 8 > 51 ** 3, 'Busch numerical cutoff')
    result = {'status': 'PASS', 'bank_sha256': hashlib.sha256(raw).hexdigest(),
              'bank_orders': dict(orders), 'H25_orders': dict(count25),
              'H51_orders': dict(count51), 'strong_H51_orders': dict(strong51),
              'direct_path_DFS_graphs': direct_graphs,
              'direct_path_DFS_subset_counts': direct_subsets,
              'positive_H25_witness': witnesses['25'],
              'strong_order_bound_control': {'5^8': 5 ** 8, '51^3': 51 ** 3},
              'H51_inventory': inventory}
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
