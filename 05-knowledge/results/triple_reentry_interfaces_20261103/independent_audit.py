#!/usr/bin/env python3
"""Independent audit: NEW permutations, isomorphism orbits, incremental covers.
No import of the producer. A second algorithmic path, not a second human audit.
"""
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


def decode(n, code):
    pairs = list(itertools.combinations(range(n), 2))
    return {(a, b) if code & (1 << i) else (b, a)
            for i, (a, b) in enumerate(pairs)}


def strong(arcs, n):
    for v in range(n):
        seen, pending = {v}, [v]
        while pending:
            x = pending.pop()
            for a, b in arcs:
                if a == x and b not in seen:
                    seen.add(b)
                    pending.append(b)
        if len(seen) < n:
            return False
    return True


def cube(arcs, n):
    """Enumerate final paths, count ALL compatible actual cut masks."""
    counts = [0] * (1 << n)
    for word in itertools.permutations(range(n + 1)):
        required_in, required_out = set(), set()
        for a, b in zip(word, word[1:]):
            if b == n:
                required_in.add(a)
            elif a == n:
                required_out.add(b)
            elif (a, b) not in arcs:
                break
        else:
            fixed = sum(1 << a for a in required_in)
            free = [a for a in range(n) if a not in required_in | required_out]
            for bits in itertools.product((0, 1), repeat=len(free)):
                mask = fixed + sum(bit << a for bit, a in zip(bits, free))
                counts[mask] += 1
    return counts


def prime_support(x):
    return {p for p in range(2, x + 1) if x % p == 0
            and all(p % a for a in range(2, math.isqrt(p) + 1))}


def covers(edges):
    out = {frozenset()}
    for edge in sorted(edges, key=lambda e: (len(e), tuple(sorted(e)))):
        nxt = {s | {p} for s in out for p in edge}
        out = {s for s in nxt if not any(t < s for t in nxt)}
    return out


def best_modulus(values, width, either):
    n = (len(values) - 1).bit_length()
    low = {v for mask, v in enumerate(values)
           if (min(mask.bit_count(), n-mask.bit_count()) if either
               else mask.bit_count()) <= width}
    factors = {v: prime_support(v) for v in set(values)}
    best = None
    for target in set(values) - low:
        edges = [factors[v] - factors[target] for v in low]
        if any(not e for e in edges):
            continue
        m = min(map(math.prod, covers(edges)))
        if best is None or m < best:
            best = m
    return best


def orbit(arcs, n):
    positions = {pair: i for i, pair in enumerate(itertools.combinations(range(n), 2))}
    result = set()
    for p in itertools.permutations(range(n)):
        code = 0
        for a, b in arcs:
            x, y = p[a], p[b]
            if x < y:
                code |= 1 << positions[x, y]
        result.add(code)
    return result


def main():
    receipt = json.loads((HERE / 'receipt.json').read_text())
    cubes = 0
    for e in receipt['examples']:
        rows, n = e['rows'], len(e['rows'])
        arcs = {(a, b) for a, r in enumerate(rows) for b in range(n) if r >> b & 1}
        values = cube(arcs, n)
        check(values == e['values'], 'independent new-permutation count')
        check(strong(arcs, n) == e['strong'], 'strong source')
        m = e['modulus']
        check(best_modulus(values, e['width'], e['two_sided']) == m, 'literal-source minimum')
        factors = {v: prime_support(v) for v in set(values)}
        low = {v for s, v in enumerate(values)
               if (min(s.bit_count(), n-s.bit_count()) if e['two_sided']
                   else s.bit_count()) <= e['width']}
        check({frozenset(s) for s in e['low_minimal_blockers']} == covers([factors[v] for v in low]),
              'complete low blocker antichain')
        check({frozenset(s) for s in e['full_minimal_blockers']} == covers(list(factors.values())),
              'complete full blocker antichain')
        # Prime powers do not alter a unit predicate; zero is not a prime.
        for multiplier in (1, 3, 5, 89, 2**9):
            if prime_support(multiplier) - prime_support(2 * 3 * 5 * 89) == set():
                # Only test a same-radical multiplier (2 is harmless for odd counts).
                if all(p == 2 or m % p == 0 for p in prime_support(multiplier)):
                    check([math.gcd(v, m) == 1 for v in values] ==
                          [math.gcd(v, m*multiplier) == 1 for v in values], 'radical lift')
        if e['strong']:
            for mask in e['unit_masks']:
                attached = arcs | {(a, n) if mask >> a & 1 else (n, a) for a in range(n)}
                check(strong(attached, n+1), 'actual unit output strong')
        cubes += len(values)
    recorded = json.loads((HERE / 'census.json').read_text())
    results = []
    for n in range(1, 7):
        remaining = set(range(1 << (n*(n-1)//2)))
        total, strong_total, classes = 0, 0, 0
        minima = {'fixed': None, 'either': None, 'strong_fixed': None}
        while remaining:
            code = min(remaining)
            arcs = decode(n, code)
            group_orbit = orbit(arcs, n)
            check(group_orbit <= remaining, 'isomorphism orbits must partition universe')
            remaining.difference_update(group_orbit)
            is_strong = strong(arcs, n)
            if n == 6 and not is_strong:
                continue
            classes += 1
            total += len(group_orbit)
            if is_strong:
                strong_total += len(group_orbit)
            values = cube(arcs, n)
            for key, either in (('fixed', False), ('either', True), ('strong_fixed', False)):
                if key == 'strong_fixed' and not is_strong:
                    continue
                candidate = best_modulus(values, 2, either)
                if candidate is not None and (minima[key] is None or candidate < minima[key]):
                    minima[key] = candidate
        r = recorded[n-1]
        check((total, strong_total) == (r['total'], r['strong_total']), 'labelled universe sizes')
        check(minima == {k: (v['modulus'] if v else None) for k, v in r['minima'].items()},
              'independent orbit minimum')
        results.append({'n': n, 'classes': classes, 'labelled': total, 'strong_labelled': strong_total,
                        'minima': minima})
    seeds = json.loads((HERE / 'scalar_seeds.json').read_text())
    check(set(map(int, seeds)) == {1, 23, 25, 63, 65, 67}, 'exact scalar root set')
    seed_audit = []
    for h, seed in seeds.items():
        rows, n = seed['rows'], len(seed['rows'])
        arcs = {(i, j) for i in range(n) for j in range(n) if rows[i] >> j & 1}
        direct = sum(all((a, b) in arcs for a, b in zip(word, word[1:]))
                     for word in itertools.permutations(range(n)))
        check(direct == int(h), 'raw scalar seed witness')
        seed_audit.append({'H': direct, 'order': n, 'strong': strong(arcs, n)})
    # Finite decoder control; the proof of universal coverage is in RETURN_FOREST.md.
    root_counts = dict.fromkeys(map(int, seeds), 0)
    for h in range(1, 100002, 2):
        if h in (7, 21):
            continue
        word, current = [], h
        while current not in root_counts:
            j = (current // 2 - 1) % 3
            parent = (current - 2*j) // 3
            check(parent > 0 and parent % 2 == 1 and parent < current,
                  'strict inverse descent')
            check(parent not in (7, 21), 'inverse forest has an unlisted root')
            word.append(j)
            current = parent
        root_counts[current] += 1
        for j in reversed(word):
            current = 3*current + 2*j
        check(current == h, 'exact ternary address replay')
    lifts = json.loads((HERE / 'seed_lifts.json').read_text())
    for item in lifts:
        rows = item['rows']; n = len(rows)
        K = [0]*n
        for word in itertools.permutations(range(n)):
            for split in range(1, n+1):
                P, Q = word[:split], word[split:]
                if all(rows[a] >> b & 1 for path in (P, Q)
                       for a, b in zip(path, path[1:])):
                    K[P[-1]] += 1
        check([k-item['H'] for k in K] == item['J'], 'permutation-cut slack audit')
        def ear(old, q):
            n0 = len(old); a, b = n0, n0+1
            return [r | ((1 << a) if i == q else 0) for i,r in enumerate(old)] + [
                (((1 << n0)-1) ^ (1 << q)) | (1 << b), (1 << n0)-1]
        def ham_dp(r):
            dp = {(1 << i, i): 1 for i in range(len(r))}
            for mask in range(1, 1 << len(r)):
                for last in range(len(r)):
                    val = dp.get((mask,last), 0)
                    if val:
                        for nxt in range(len(r)):
                            if not mask >> nxt & 1 and r[last] >> nxt & 1:
                                key = (mask | (1 << nxt), nxt)
                                dp[key] = dp.get(key,0) + val
            return sum(dp.get(((1 << len(r))-1, i),0) for i in range(len(r)))
        for q in item['best_ports']:
            first = ear(rows, q)
            second = ear(first, n)
            check(ham_dp(first) == 3*item['H']+2*item['minimum_J'], 'literal first ear')
            check(ham_dp(second) == item['minimum_lift_H'], 'literal two-ear canonical lift')
    r5_copies = {}
    for h, seed in seeds.items():
        rows = seed['rows']
        copies = [ids for ids in itertools.combinations(range(len(rows)),5)
                  if all(sum(rows[a] >> b & 1 for b in ids) == 2 for a in ids)]
        r5_copies[h] = copies
        for ids in copies:
            starts = {a:0 for a in ids}
            for word in itertools.permutations(ids):
                if all(rows[a] >> b & 1 for a,b in zip(word,word[1:])):
                    starts[word[0]] += 1
            check(set(starts.values()) == {3}, 'regular-five starting multiplicities')
    check(r5_copies['65'] == [(0,1,2,3,6)], 'literal induced ready obstruction')
    check(all(not copies for h,copies in r5_copies.items() if h != '65'),
          'other seed R5 controls')
    result = {'status': 'PASS', 'raw_seeds': seed_audit, 'independent_rooted_lift_profiles': len(lifts), 'r5_copies_by_seed': r5_copies, 'scalar_roots_through_100001': root_counts, 'independently_counted_literal_cuts': cubes,
              'orbit_census': results}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
