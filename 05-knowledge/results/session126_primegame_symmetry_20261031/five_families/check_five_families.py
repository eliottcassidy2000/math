#!/usr/bin/env python3
"""Session 126 bounded exact checks. Standard library only; no external inputs.

Universe: square-sum graphs 2<=n<=26, labelled 5-tournaments, first ten primes,
and XOR kernels 0<=t<=63. No unbounded Fermat or square-sum theorem is inferred.
"""
from collections import Counter
from itertools import combinations, permutations
from math import comb, isqrt
import json


def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def constructible_order(n):
    f = factor(n)
    return all(p == 2 or (e == 1 and (p-1) & (p-2) == 0)
               for p, e in f.items())


def square_path(n):
    adj = {i: [j for j in range(1, n+1)
               if i != j and isqrt(i+j)**2 == i+j] for i in range(1, n+1)}
    calls = 0

    def visit(path, remaining):
        nonlocal calls
        calls += 1
        if not remaining:
            return path
        options = sorted(set(adj[path[-1]]) & remaining,
                         key=lambda j: (len(set(adj[j]) & remaining), j))
        for j in options:
            ans = visit(path+[j], remaining-{j})
            if ans is not None:
                return ans
        return None

    # Every path has an endpoint; no pruning, hence exhaustion is a proof of no path.
    for start in sorted(adj, key=lambda x: (len(adj[x]), x)):
        ans = visit([start], set(adj)-{start})
        if ans is not None:
            if not (sorted(ans) == list(range(1, n+1))):
                raise RuntimeError('check_five_families.py: gate originally at line 54')
            if not (all(isqrt(a+b)**2 == a+b for a, b in zip(ans, ans[1:]))):
                raise RuntimeError('check_five_families.py: gate originally at line 55')
            return {'path': ans, 'search_calls': calls}
    return {'path': None, 'search_calls': calls}


def pair_and_tournament_checks():
    n = 5
    pairs = list(combinations(range(n), 2))
    pair_idx = {p: i for i, p in enumerate(pairs)}
    perms = list(permutations(range(n)))
    actions = []
    for perm in perms:
        target = []
        for a, b in pairs:
            x, y = perm[a], perm[b]
            target.append((pair_idx[tuple(sorted((x, y)))], x > y))
        actions.append(target)

    def relabel(word, action):
        out = 0
        for i, (j, reverse) in enumerate(action):
            if ((word >> i) & 1) ^ reverse:
                out |= 1 << j
        return out

    canon = [min(relabel(w, a) for a in actions) for w in range(1 << len(pairs))]
    reps = sorted(set(canon))
    full = (1 << len(pairs))-1
    reverse = {r: canon[r ^ full] for r in reps}
    merged = {r: min(r, reverse[r]) for r in reps}
    m_reps = sorted(set(merged.values()))
    if not (len(reps) == 12 and sum(reverse[r] == r for r in reps) == 8):
        raise RuntimeError('check_five_families.py: gate originally at line 86')
    if not (len(m_reps) == 10):
        raise RuntimeError('check_five_families.py: gate originally at line 87')
    if not (all(canon[relabel(r, a)] == r for r in reps for a in actions)):
        raise RuntimeError('check_five_families.py: gate originally at line 88')
    pair_permutations = {tuple(j for j, _ in a) for a in actions}
    if not (len(pair_permutations) == 120):
        raise RuntimeError('check_five_families.py: gate originally at line 90')
    j_edges = [(i, j) for i, j in combinations(range(10), 2)
               if set(pairs[i]) & set(pairs[j])]
    p_edges = [(i, j) for i, j in combinations(range(10), 2)
               if not set(pairs[i]) & set(pairs[j])]
    if not (len(j_edges) == 30 and len(p_edges) == 15):
        raise RuntimeError('check_five_families.py: gate originally at line 95')

    def graph_data(nodes, edges):
        deg = Counter(v for e in edges for v in e)
        return {'vertices': len(nodes), 'edges': len(edges),
                'degree_multiset': sorted(deg[v] for v in nodes)}

    flip_edges = set()
    merged_flip_edges = set()
    for r in reps:
        for k in range(10):
            s = canon[r ^ (1 << k)]
            if r != s:
                flip_edges.add(tuple(sorted((r, s))))
            if merged[r] != merged[s]:
                merged_flip_edges.add(tuple(sorted((merged[r], merged[s]))))
    if not (graph_data(range(10), j_edges)['degree_multiset'] == [6]*10):
        raise RuntimeError('check_five_families.py: gate originally at line 111')
    if not (graph_data(range(10), p_edges)['degree_multiset'] == [3]*10):
        raise RuntimeError('check_five_families.py: gate originally at line 112')
    return {'arc_positions': 10, 'pair_action_image_order': len(pair_permutations),
            'isoclasses': 12, 'self_converse': 8, 'converse_merged': 10,
            'induced_relabel_action_on_isoclasses': 'trivial',
            'johnson_J52': graph_data(range(10), j_edges),
            'petersen_KG52': graph_data(range(10), p_edges),
            'single_arc_flip_iso_graph': graph_data(reps, flip_edges),
            'single_arc_flip_converse_merged_graph': graph_data(m_reps, merged_flip_edges),
            'single_arc_flip_converse_merged_edges': sorted(merged_flip_edges)}


def multiplier_checks():
    out = []
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29]:
        m = (p-1)//2
        connection_sets = [frozenset(i if w >> (i-1) & 1 else p-i
                                     for i in range(1, m+1)) for w in range(1 << m)]
        fixed_counts = {}
        for a in range(1, p):
            count = sum(frozenset((a*x) % p for x in S) == S for S in connection_sets)
            order, x = 1, a
            while x != 1:
                order += 1
                x = x*a % p
            expected = 0 if order % 2 == 0 else 2**((p-1)//(2*order))
            if not (count == expected):
                raise RuntimeError('check_five_families.py: gate originally at line 137')
            fixed_counts[a] = count
        rigid = not any(count for a, count in fixed_counts.items() if a != 1)
        if not (rigid == ((p-1) & (p-2) == 0)):
            raise RuntimeError('check_five_families.py: gate originally at line 140')
        out.append({'prime': p, 'p_minus_one_factorization': factor(p-1),
                    'connection_set_count': len(connection_sets),
                    'all_multiplier_stabilizers_trivial': rigid,
                    'nonidentity_multipliers_fixing_some_set':
                        {a: c for a, c in fixed_counts.items() if a != 1 and c}})
    return out


def a5_checks():
    def compose(a,b):
        return tuple(a[b[i]] for i in range(5))
    def inverse(a):
        return tuple(a.index(i) for i in range(5))
    def conjugate(a,b):
        return compose(compose(a,b),inverse(a))
    group = [p for p in permutations(range(5))
             if sum(p[i] > p[j] for i,j in combinations(range(5),2)) % 2 == 0]
    if not (len(group) == 60):
        raise RuntimeError('check_five_families.py: gate originally at line 158')
    remaining = set(group)
    conjugacy_sizes = []
    while remaining:
        a = min(remaining)
        cls = {conjugate(g,a) for g in group}
        conjugacy_sizes.append(len(cls))
        remaining -= cls
    if not (sorted(conjugacy_sizes) == [1,12,12,15,20]):
        raise RuntimeError('check_five_families.py: gate originally at line 166')
    edges = list(combinations(range(5),2))
    subgroups = {}
    for edge in edges:
        a,b,c = sorted(set(range(5))-set(edge))
        cyc = list(range(5))
        cyc[a],cyc[b],cyc[c] = b,c,a
        cyc = tuple(cyc)
        subgroups[edge] = frozenset([tuple(range(5)),cyc,compose(cyc,cyc)])
    if not (len(set(subgroups.values())) == 10):
        raise RuntimeError('check_five_families.py: gate originally at line 175')
    for g in group:
        for edge,subgroup in subgroups.items():
            image_edge = tuple(sorted(g[i] for i in edge))
            if not (frozenset(conjugate(g,h) for h in subgroup) == subgroups[image_edge]):
                raise RuntimeError('check_five_families.py: gate originally at line 179')
    return {'order':60, 'prime_factorization':factor(60),
            'conjugacy_class_sizes':sorted(conjugacy_sizes),
            'pair_to_C3_subgroup_equivariance_checks':len(group)*len(edges),
            'C3_subgroups':len(subgroups)}


def main():
    paths = {n: square_path(n) for n in range(2, 27)}
    successful = [n for n, row in paths.items() if row['path'] is not None]
    if not (successful == [15, 16, 17, 23, 25, 26]):
        raise RuntimeError('check_five_families.py: gate originally at line 189')
    platonic = [(p, q) for p in range(3, 9) for q in range(3, 9)
                if (p-2)*(q-2) < 4]
    if not (platonic == [(3, 3), (3, 4), (3, 5), (4, 3), (5, 3)]):
        raise RuntimeError('check_five_families.py: gate originally at line 192')
    fermat = [2**(2**i)+1 for i in range(6)]
    if not (all(factor(p) == {p: 1} for p in fermat[:5])):
        raise RuntimeError('check_five_families.py: gate originally at line 194')
    if not (fermat[5] == 641*6700417):
        raise RuntimeError('check_five_families.py: gate originally at line 195')
    kernels = []
    row = 1
    raw_rule90_support = {0}
    for t in range(64):
        direct = sum((comb(t, j) % 2) << j for j in range(t+1))
        product = 1
        for i in range(t.bit_length()):
            if t >> i & 1:
                product *= fermat[i]
        if not (row == direct == product):
            raise RuntimeError('check_five_families.py: gate originally at line 204')
        if not all((position+t) % 2 == 0 for position in raw_rule90_support):
            raise RuntimeError('Rule90 occupied-position parity')
        compressed90 = sum(1 << ((position+t)//2) for position in raw_rule90_support)
        if compressed90 != row:
            raise RuntimeError('Rule90 parity-compression equals one-sided Rule60')
        if t == 1:
            rebased90 = sum(1 << (position+t) for position in raw_rule90_support)
            if not (rebased90 == 5 and row == 3):
                raise RuntimeError('Raw Rule90 versus compressed row hostile')
        kernels.append(row)
        row ^= row << 1
        raw_rule90_support = ({j-1 for j in raw_rule90_support} ^
                              {j+1 for j in raw_rule90_support})
    out = {'square_sum_path_orders_2_through_26': successful,
           'square_sum_path_checks': paths,
           'constructible_orders_in_first_five':
                [n for n in successful[:5] if constructible_order(n)],
           'platonic_schlafli_pairs': platonic,
           'first_six_fermat_numbers': fermat,
           'F5_factorization': factor(fermat[5]),
           'xor_kernel_first_64_verified': True,
           'xor_convention': 'Rule60 one-sided XOR; Rule90 only after position 2*j-t -> j compression',
           'rule90_parity_compression_first_64_verified': True,
           'rule90_time1_raw_rebased_integer':5,
           'rule60_time1_integer':3,
           'xor_kernel_first_33': kernels[:33],
           'first_ten_prime_multiplier_checks_except_even_prime_2': multiplier_checks(),
           'A5_geometry_bridge': a5_checks(),
           'pair_and_tournament_checks': pair_and_tournament_checks()}
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
