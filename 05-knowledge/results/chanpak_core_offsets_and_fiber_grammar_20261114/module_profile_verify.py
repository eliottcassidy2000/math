#!/usr/bin/env python3
"""Unpublished exact module-response census. Stdlib only; checks survive -O.

Masks list pairs (i,j), i<j, lexicographically: bit 0 means i->j.
Polynomial tuples are ascending ordinary powers. pc[k] counts UNORDERED
spanning directed path covers with k nonempty paths, singleton paths allowed.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial
import json
import time


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def adj(n, mask):
    out = [0] * n
    for b, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (j, i) if mask >> b & 1 else (i, j)
        out[u] |= 1 << v
    return tuple(out)


def poly_from_pc(pc):
    ans = [0] * len(pc)
    fall = [1]
    for k, c in enumerate(pc):
        for j, x in enumerate(fall):
            ans[j] += c * x
        nxt = [0] * (len(fall) + 1)
        for j, x in enumerate(fall):
            nxt[j] -= k * x
            nxt[j + 1] += x
        fall = nxt
    return tuple(ans)


@lru_cache(None)
def pc_from_poly(p):
    vals = [sum(a * x**i for i, a in enumerate(p)) for x in range(len(p))]
    pc = []
    k = 0
    while vals:
        q, r = divmod(vals[0], factorial(k))
        check(r == 0, 'nonintegral path-cover conversion')
        pc.append(q)
        vals = [b - a for a, b in zip(vals, vals[1:])]
        k += 1
    return tuple(pc)


def poly_mul(p, q):
    ans = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            ans[i + j] += x * y
    return tuple(ans)


def direct_pc(out):
    """Independent exhaustive permutations/cuts, not quotient walk kernel."""
    n = len(out)
    totals = [0] * (n + 1)
    for perm in permutations(range(n)):
        bad = sum(not (out[a] >> b & 1) for a, b in zip(perm, perm[1:]))
        good = n - 1 - bad
        for k in range(bad + 1, n + 1):
            totals[k] += comb(good, k - 1 - bad)
    for k in range(1, n + 1):
        check(totals[k] % factorial(k) == 0, 'nonfree ordered-cover action')
        totals[k] //= factorial(k)
    return tuple(totals)


def h_dp(out):
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for i, ways in enumerate(dp[s]):
            if not ways:
                continue
            rest = out[i] & ~s
            while rest:
                bit = rest & -rest
                j = bit.bit_length() - 1
                dp[s | bit][j] += ways
                rest -= bit
    return sum(dp[-1])


def relabel(mask, n, perm):
    pairs = list(combinations(range(n), 2))
    index = {uv: b for b, uv in enumerate(pairs)}
    result = 0
    for b, (i, j) in enumerate(pairs):
        u, v = perm[i], perm[j]
        reverse = (mask >> b & 1) ^ (u > v)
        if reverse:
            result |= 1 << index[tuple(sorted((u, v)))]
    return result


def class_reps(n):
    seen = set()
    perms = list(permutations(range(n)))
    for mask in range(1 << comb(n, 2)):
        if mask not in seen:
            orbit = {relabel(mask, n, p) for p in perms}
            seen.update(orbit)
            yield mask


def is_strong(out):
    n = len(out)
    reach = [out[i] | 1 << i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << n) - 1 for x in reach)


def is_prime(out):
    n = len(out)
    if n < 3:
        return False
    for s in range(1, (1 << n) - 1):
        if s.bit_count() < 2:
            continue
        if all((out[v] & s) in (0, s) for v in range(n) if not (s >> v & 1)):
            return False
    return True


def cycles(perm):
    seen = set()
    result = []
    for i in range(len(perm)):
        if i not in seen:
            cycle = []
            while i not in seen:
                seen.add(i)
                cycle.append(i)
                i = perm[i]
            result.append(tuple(cycle))
    return tuple(result)


@lru_cache(None)
def walk_powers(qmask, sizes):
    """Coefficients of W_Q^d, truncated componentwise, integer arithmetic."""
    q = len(sizes)
    out = adj(q, qmask)
    zero = (0,) * q
    states = {}
    for c in sorted(product(*(range(s + 1) for s in sizes)), key=sum):
        if not any(c):
            continue
        row = [0] * q
        if sum(c) == 1:
            row[c.index(1)] = 1
        else:
            for last in range(q):
                if c[last]:
                    before = list(c)
                    before[last] -= 1
                    prev = states.get(tuple(before), (0,) * q)
                    row[last] = sum(prev[i] for i in range(q) if out[i] >> last & 1)
        states[c] = row
    w = {c: sum(row) for c, row in states.items() if any(row)}
    powers = [{zero: 1}]
    for _ in range(sum(sizes)):
        nxt = defaultdict(int)
        for a, va in powers[-1].items():
            for b, vb in w.items():
                c = tuple(x + y for x, y in zip(a, b))
                if all(x <= bound for x, bound in zip(c, sizes)):
                    nxt[c] += va * vb
        powers.append(dict(nxt))
    return tuple(powers)


@lru_cache(None)
def substitute(qmask, profiles):
    sizes = tuple(len(p) - 1 for p in profiles)
    fs = [tuple(factorial(k) * c for k, c in enumerate(pc_from_poly(p)))
          for p in profiles]
    weights = {}
    for c in product(*(range(1, s + 1) for s in sizes)):
        weight = 1
        for i, k in enumerate(c):
            weight *= fs[i][k]
        weights[c] = weight
    answer = [0] * (sum(sizes) + 1)
    for d, power in enumerate(walk_powers(qmask, sizes)):
        if d:
            total = sum(weight * power.get(c, 0) for c, weight in weights.items())
            answer[d], rem = divmod(total, factorial(d))
            check(rem == 0, 'nonintegral substitution kernel')
    return poly_from_pc(answer)


def substituted_graph(qmask, blocks):
    sizes = [len(b) for b in blocks]
    qout = adj(len(blocks), qmask)
    ans = [0] * sum(sizes)
    offsets = [0]
    for s in sizes:
        offsets.append(offsets[-1] + s)
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            ans[offsets[i] + u] |= row << offsets[i]
            for j in range(len(blocks)):
                if qout[i] >> j & 1:
                    ans[offsets[i] + u] |= ((1 << sizes[j]) - 1) << offsets[j]
    return tuple(ans)


def gate_poly(p):
    ans = [0] * (len(p) + 2)
    for e, c in enumerate(p):
        ans[e + 2] += c
        for odd in range(1, e + 1, 2):
            ans[e - odd + 1] += 2 * c * comb(e, odd)
    return tuple(ans)


def tower_responses(p, count):
    values = []
    for _ in range(count):
        values.append(sum(p))
        p = gate_poly(p)
    return tuple(values)


def prime_contribution(n, mask, q, a):
    autos = [cycles(p) for p in permutations(range(q)) if relabel(mask, q, p) == mask]
    numerator = Counter()
    for cyc in autos:
        assigned = [None] * q
        def visit(i, remaining, weight):
            if i == len(cyc):
                if remaining == 0:
                    numerator[substitute(mask, tuple(assigned))] += weight
                return
            length = len(cyc[i])
            after_min = sum(map(len, cyc[i + 1:]))
            for size in range(1, (remaining - after_min) // length + 1):
                if size >= len(a):
                    continue
                for p, multiplicity in a[size].items():
                    for vertex in cyc[i]:
                        assigned[vertex] = p
                    visit(i + 1, remaining - length * size, weight * multiplicity)
        visit(0, n, 1)
    for p in numerator:
        check(numerator[p] % len(autos) == 0, 'Burnside profile coefficient not integral')
        numerator[p] //= len(autos)
    return numerator


def modular_census(max_n, primes):
    a = [Counter({(1,): 1})]
    strong = [Counter()]
    for n in range(1, max_n + 1):
        s = Counter({(0, 1): 1}) if n == 1 else Counter()
        for q, mask in primes:
            if q <= n:
                s.update(prime_contribution(n, mask, q, a))
        strong.append(s)
        total = Counter()
        for k in range(1, n + 1):
            for p, c in strong[k].items():
                for r, d in a[n - k].items():
                    total[poly_mul(p, r)] += c * d
        a.append(total)
    return a, strong


def main():
    started = time.perf_counter()
    max_n = 6
    direct = [Counter({(1,): 1})]
    direct_strong = [Counter()]
    primes = []
    reps_by_n = []
    rows = []
    for n in range(1, max_n + 1):
        inv = Counter()
        strong = Counter()
        reps = list(class_reps(n))
        reps_by_n.append(reps)
        response_owner = {}
        prime_count = 0
        for mask in reps:
            out = adj(n, mask)
            pc = direct_pc(out)
            p = poly_from_pc(pc)
            check(sum(p) == h_dp(out), 'Hamiltonian permutation/subset-DP mismatch')
            check(pc_from_poly(p) == pc, 'falling-power round trip failed')
            check(p[-1] == 1 and p[0] == 0, 'monic/constant invariant failed')
            check(all(c == 0 or (i % 2 == n % 2 and c >= 0)
                      for i, c in enumerate(p)), 'parity/nonnegative invariant failed')
            d = (n - 1) // 2
            check(all(p[n - 2 * j] % 2 == 0 for j in range(1, d + 1)), 'even cycle weights failed')
            h = sum(p)
            e = min(d, (h - 1) // 2)
            check(all(p[n - 2 * j] == 0 for j in range(e + 1, d + 1)), 'fixed-h support cutoff failed')
            response = tower_responses(p, d)
            check(response not in response_owner or response_owner[response] == p,
                  'tomography collision')
            response_owner[response] = p
            inv[p] += 1
            if is_strong(out):
                strong[p] += 1
            if is_prime(out):
                check(is_strong(out), 'prime quotient not strong')
                primes.append((n, mask))
                prime_count += 1
        direct.append(inv)
        direct_strong.append(strong)
        h_hist = Counter()
        for p, count in inv.items():
            h_hist[sum(p)] += count
        rows.append({'n': n, 'classes': len(reps), 'strong_classes': sum(strong.values()),
                     'prime_classes': prime_count, 'profiles': len(inv),
                     'H_fibers': dict(sorted(h_hist.items()))})
    check([sum(c.values()) for c in direct] == [1, 1, 1, 2, 4, 12, 56], 'A000568 control failed')
    a, s = modular_census(max_n, primes)
    check(a == direct, 'complete polynomial-class census mismatch')
    check(s == direct_strong, 'strong polynomial-class census mismatch')
    # Every binary ordinal join and C3 profile triple with total order <=6.
    controls = 0
    candidates = [(n, mask, poly_from_pc(direct_pc(adj(n, mask))))
                  for n, reps in enumerate(reps_by_n, 1) for mask in reps]
    for q, mask in ((2, 0), (3, 2)):
        for choices in product(candidates, repeat=q):
            if sum(c[0] for c in choices) > 6:
                continue
            predicted = substitute(mask, tuple(c[2] for c in choices))
            graph = substituted_graph(mask, [adj(c[0], c[1]) for c in choices])
            check(predicted == poly_from_pc(direct_pc(graph)), 'direct composition mismatch')
            controls += 1
    witness = []
    for mask in (69, 76):
        out = adj(5, mask)
        p = poly_from_pc(direct_pc(out))
        g = substituted_graph(2, [out, (0,), (0,)])
        check(sum(gate_poly(p)) == h_dp(g), 'gate operator independent count mismatch')
        witness.append({'mask': mask, 'P': p, 'pc': pc_from_poly(p),
                        'H': sum(p), 'H_gate': h_dp(g)})
    check(witness[0]['H'] == witness[1]['H'] == 15, 'witness H changed')
    check([w['H_gate'] for w in witness] == [123, 135], 'witness response changed')
    restricted, _ = modular_census(9, [(3, 2)])
    report = {'status': 'PASS', 'base_commit': '8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90',
              'exact_rows': rows, 'independent_composition_controls': controls,
              'same_H_counterexample': witness,
              'C3_only_counts_n0_to_n9': [sum(c.values()) for c in restricted],
              'C3_only_profile_counts_n0_to_n9': [len(c) for c in restricted],
              'elapsed_seconds': round(time.perf_counter() - started, 3)}
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
