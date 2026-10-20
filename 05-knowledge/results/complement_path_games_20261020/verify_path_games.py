#!/usr/bin/env python3
"""Exact bounded verification; stdlib only. Run: python3 verify_path_games.py."""
from itertools import combinations, permutations
from functools import lru_cache
from math import comb
import json
from pathlib import Path


def graph(n, edges):
    a = [0] * n
    for u, v in edges:
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def complement(a):
    full = (1 << len(a)) - 1
    return tuple(full ^ (1 << v) ^ a[v] for v in range(len(a)))


def components(a):
    remaining = (1 << len(a)) - 1
    result = []
    while remaining:
        todo = remaining & -remaining
        seen = 0
        while todo:
            bit = todo & -todo
            todo ^= bit
            v = bit.bit_length() - 1
            seen |= bit
            todo |= a[v] & remaining & ~seen
        result.append(seen)
        remaining &= ~seen
    return result


def endpoints(a):
    """D[S] contains possible endpoints of Hamilton paths of induced G[S]."""
    D = [0] * (1 << len(a))
    for v in range(len(a)):
        D[1 << v] = 1 << v
    for mask in range(1, len(D)):
        ends = D[mask]
        while ends:
            b = ends & -ends
            ends ^= b
            v = b.bit_length() - 1
            nxt = a[v] & ~mask
            while nxt:
                c = nxt & -nxt
                nxt ^= c
                D[mask | c] |= c
    return D


def path_cover_functions(D):
    @lru_cache(None)
    def pc(mask):
        if mask == 0:
            return 0
        if D[mask]:
            return 1
        first = mask & -mask
        best = mask.bit_count()
        sub = mask
        while sub:
            if sub & first and D[sub]:
                best = min(best, 1 + pc(mask ^ sub))
            sub = (sub - 1) & mask
        return best

    def rooted_pc(mask, root):
        bit = 1 << root
        best = mask.bit_count()
        sub = mask
        while sub:
            if sub & bit and D[sub] & bit:
                best = min(best, 1 + pc(mask ^ sub))
            sub = (sub - 1) & mask
        return best
    return pc, rooted_pc


def walk_ham_cost(n):
    return {6: 6, 8: 10}.get(n, 2 * n - 5)


def exhaustive(max_n=6):
    rows = []
    checked_joins = 0
    checked_roots = 0
    for n in range(1, max_n + 1):
        pairs = list(combinations(range(n), 2))
        full = (1 << n) - 1
        count = [0] * 4
        cost = [None] * 4
        root_count = [0] * 4
        root_cost = [None] * 4
        both_connected = 0
        for code in range(1 << len(pairs)):
            R = graph(n, [e for i, e in enumerate(pairs) if code >> i & 1])
            B = complement(R)
            DR, DB = endpoints(R), endpoints(B)
            cr, cb = components(R), components(B)
            hr, hb = bool(DR[-1]), bool(DB[-1])
            wr, wb = len(cr) == 1, len(cb) == 1
            eb = DB[-1] == full
            outcomes = (not hr and not hb, not hr and not wb,
                        not wr and not hb, not wr and not wb)
            rooted = (not hr and not eb, not hr and not wb,
                      not wr and not eb, not wr and not wb)
            edge_count = code.bit_count()
            for flags, counts, costs in [(outcomes, count, cost),
                                         (rooted, root_count, root_cost)]:
                for j, yes in enumerate(flags):
                    if yes:
                        counts[j] += 1
                        costs[j] = edge_count if costs[j] is None else min(costs[j], edge_count)
            both_connected += bool(wr and wb and not hr and not hb)
            if len(cr) > 1:
                pc, rooted_pc = path_cover_functions(DB)
                assert hb == all(pc(C) <= n - C.bit_count() + 1 for C in cr)
                assert eb == all(pc(C) <= n - C.bit_count() for C in cr)
                for s in range(n):
                    Cs = next(C for C in cr if C >> s & 1)
                    predicted = rooted_pc(Cs, s) <= n - Cs.bit_count() + 1
                    predicted &= all(pc(C) <= n - C.bit_count() for C in cr if C != Cs)
                    assert predicted == bool(DB[-1] >> s & 1)
                    checked_roots += 1
                checked_joins += 1
        expected = [n - 1, n - 1, walk_ham_cost(n), None] if n >= 4 else [None] * 4
        expected_root = [n - 2, n - 1 if n >= 4 else None, n - 2, None] if n >= 3 else [None] * 4
        assert cost == expected, (n, cost, expected)
        assert root_cost == expected_root, (n, root_cost, expected_root)
        rows.append(dict(n=n, labeled_graphs=1 << len(pairs),
                         winning_colorings=count, minimum_red_edges=cost,
                         every_start_winning_colorings=root_count,
                         every_start_minimum_red_edges=root_cost,
                         connected_both_nontraceable=both_connected))
    return rows, checked_joins, checked_roots


def witness_checks():
    result = []
    for n in range(4, 17):
        # Blue clique n-2 plus two leaves adjacent to its one distinguished vertex.
        blue_edges = list(combinations(range(n - 2), 2)) + [(0, n - 2), (0, n - 1)]
        B = graph(n, blue_edges)
        R = complement(B)
        assert len(components(R)) > 1 and not endpoints(B)[-1]
        edges_r = sum(x.bit_count() for x in R) // 2
        assert edges_r == 2 * n - 5
        if n in (6, 8):
            q = n // 2
            R = graph(n, combinations(range(q + 1), 2))
            B = complement(R)
            assert len(components(R)) > 1 and not endpoints(B)[-1]
            edges_r = sum(x.bit_count() for x in R) // 2
        assert edges_r == walk_ham_cost(n)
        result.append(dict(n=n, walk_ham_witness_red_edges=edges_r))
    # Parent's sharper connected-both family: blue clique n-2 minus hv,
    # with two leaves at h. Deleting h leaves three blue components;
    # deleting the two leaves in red leaves n-3 components.
    connected_rows = []
    for n in range(7, 17):
        h, v, a, b = 0, 1, n - 2, n - 1
        blue_edges = [e for e in combinations(range(n - 2), 2) if e != (h, v)]
        blue_edges += [(h, a), (h, b)]
        B = graph(n, blue_edges)
        R = complement(B)
        assert len(components(R)) == len(components(B)) == 1
        assert not endpoints(R)[-1] and not endpoints(B)[-1]
        connected_rows.append(dict(n=n, blue_edges=[list(e) for e in blue_edges]))
    return result, connected_rows



def independent_permutation_checks():
    # A structurally different verifier, all graphs through six vertices.
    checked = 0
    for n in range(1, 7):
        pairs = list(combinations(range(n), 2))
        for code in range(1 << len(pairs)):
            a = graph(n, [e for i, e in enumerate(pairs) if code >> i & 1])
            direct = 0
            for p in permutations(range(n)):
                if all(a[u] >> v & 1 for u, v in zip(p, p[1:])):
                    direct |= 1 << p[0]
            assert direct == endpoints(a)[-1]
            checked += 1
    return checked


if __name__ == '__main__':
    rows, joins, roots = exhaustive()
    witnesses, connected_witness = witness_checks()
    independent = independent_permutation_checks()
    report = dict(status='FINITE-EXACT', universe='All labeled simple undirected graphs on n=1..6; all start vertices; n=1..6 permutation cross-check',
                  variant_order=['H/H', 'H/W', 'W/H', 'W/W'],
                  rows=rows, join_criteria_checked=joins,
                  prescribed_root_criteria_checked=roots,
                  independent_permutation_graphs=independent,
                  witness_rows=witnesses, connected_both_witness=connected_witness,
                  caveat='Finite checks validate all-n proofs. Connected-both minimal order seven uses exhaustive exclusion through six plus the proved family for all n>=7.')
    output = Path(__file__).with_name('verification.json')
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
