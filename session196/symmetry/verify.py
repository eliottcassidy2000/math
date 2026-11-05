#!/usr/bin/env python3
"""Independent small controls for Hamilton-path orbits and prime rigidity.
No external packages. Explicit checks stay active under python -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
import json
from pathlib import Path


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def rows_from_bits(n, free, bits):
    rows = [0] * n
    for i in range(n - 1):
        rows[i] |= 1 << (i + 1)
    for k, (i, j) in enumerate(free):
        if bits >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def h_dp(rows):
    n = len(rows)
    # Forward subset DP, not the literal permutation implementation below.
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for last, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = rows[last] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def h_reverse(rows):
    incoming = [sum(((rows[i] >> j) & 1) << i for i in range(len(rows))) for j in range(len(rows))]
    @lru_cache(None)
    def end(mask, v):
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        return sum(end(rest, u) for u in range(len(rows)) if rest >> u & 1 and incoming[v] >> u & 1)
    mask = (1 << len(rows)) - 1
    return sum(end(mask, v) for v in range(len(rows)))


def f2_subsets(rows):
    n = len(rows)
    def induced_h(mask):
        vs = [v for v in range(n) if mask >> v & 1]
        small = tuple(sum(((rows[v] >> w) & 1) << j for j, w in enumerate(vs)) for v in vs)
        return h_dp(small)
    full = (1 << n) - 1
    return sum(induced_h(mask) * induced_h(full ^ mask) for mask in range(1, full))


def f2_literal(rows):
    answer = 0
    n = len(rows)
    for p in permutations(range(n)):
        valid = [bool(rows[p[i]] >> p[i+1] & 1) for i in range(n-1)]
        for cut in range(1, n):
            answer += all(valid[:cut-1]) and all(valid[cut:])
    return answer


def h_literal(rows):
    return sum(all(rows[p[i]] >> p[i + 1] & 1 for i in range(len(rows) - 1))
               for p in permutations(range(len(rows))))


def canonical_with_aut(rows, perms, pairs):
    original = sum(((rows[i] >> j) & 1) << k for k, (i, j) in enumerate(pairs))
    minimum = None
    aut = 0
    for p in perms:
        code = sum(((rows[p[i]] >> p[j]) & 1) << k for k, (i, j) in enumerate(pairs))
        minimum = code if minimum is None else min(minimum, code)
        aut += code == original
    return minimum, aut


def prime(p):
    return p >= 2 and all(p % q for q in range(2, int(p ** .5) + 1))


def census(n):
    # Every tournament has a Hamilton path. Normalizing that path to
    # 0,1,...,n-1 leaves exactly binom(n-1,2) free orientations.
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    free = [(i, j) for i, j in pairs if j != i + 1]
    perms = list(permutations(range(n)))
    multiplicities = Counter()
    classes = {}
    for bits in range(1 << len(free)):
        rows = rows_from_bits(n, free, bits)
        code, aut = canonical_with_aut(rows, perms, pairs)
        multiplicities[code] += 1
        if code in classes:
            check(classes[code]['aut'] == aut, 'automorphism class invariance')
        else:
            h = h_dp(rows)
            check(h == h_literal(rows), 'DP versus literal count')
            check(h % 2 == 1 and aut % 2 == 1, 'odd H and odd automorphisms')
            check(h % aut == 0, 'free path action divisibility')
            if prime(h) and h >= 5:
                check(aut == 1, 'prime rigidity')
            classes[code] = {'h': h, 'aut': aut, 'rows': list(rows)}
    for code, data in classes.items():
        check(multiplicities[code] == data['h'] // data['aut'], 'pointed class count')
    expected = [1, 1, 1, 2, 4, 12, 56][n]
    check(len(classes) == expected, 'A000568')
    check(sum(d['h'] // d['aut'] for d in classes.values()) == 2 ** comb(n-1, 2), 'mass identity')
    return {'n': n, 'unlabeled_count': len(classes), 'path_rooted_count': sum(multiplicities.values()),
            'spectrum': sorted({d['h'] for d in classes.values()}),
            'H_aut_histogram': {f"{h},{a}": c for (h, a), c in sorted(Counter((d['h'], d['aut']) for d in classes.values()).items())}}


def circulants(p):
    hist = Counter()
    representatives = {}
    for bits in range(1 << ((p - 1) // 2)):
        steps = [i if bits >> (i - 1) & 1 else p - i for i in range(1, (p+1)//2)]
        rows = tuple(sum(1 << ((v + s) % p) for s in steps) for v in range(p))
        h = h_dp(rows)
        check(h == h_reverse(rows), 'circulant reverse DP')
        if p <= 7:
            check(h == h_literal(rows), 'circulant independent count')
        check(h % p == 0 and h % 2 == 1, 'cyclic free action')
        triangles = sum(1 for a in range(p) for b in range(a+1, p) for c in range(b+1, p)
                        if (rows[a] >> b & 1) == (rows[b] >> c & 1) == (rows[c] >> a & 1))
        check(triangles == p * (p*p-1)//24, 'regular triangle formula')
        check(h >= 1 + 2 * triangles, 'triangle lower bound control')
        hist[h] += 1
        representatives.setdefault(h, list(rows))
        if p <= 7:
            f2 = f2_subsets(rows)
            check(f2 == f2_literal(rows), 'two-cover literal check')
            if p == 7:
                check(f2 == 1638, 'seven-cycle F2 constant')
    return {'p': p, 'labeled_connection_sets': sum(hist.values()), 'H_histogram': dict(sorted(hist.items())),
            'minimum_H': min(hist), 'minimum_rows': representatives[min(hist)]}


def is_strong(rows):
    full = (1 << len(rows)) - 1
    for start in range(len(rows)):
        seen = 1 << start
        while True:
            reached = seen
            for v in range(len(rows)):
                if seen >> v & 1:
                    reached |= rows[v]
            if reached == seen:
                break
            seen = reached
        if seen != full:
            return False
    return True


def all_rows(n, bits):
    rows = [0] * n
    for k, (i,j) in enumerate(( (i,j) for i in range(n) for j in range(i+1,n))):
        if bits >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return rows


def strong4_coefficients():
    profiles = Counter()
    canon = set()
    pairs = [(i,j) for i in range(4) for j in range(i+1,4)]
    perms = list(permutations(range(4)))
    for bits in range(64):
        rows = all_rows(4,bits)
        if not is_strong(rows):
            continue
        canon.add(canonical_with_aut(rows,perms,pairs)[0])
        for root in range(4):
            coeff = []
            for k in range(1,5):
                words = set(permutations([root]*k+[v for v in range(4) if v != root]))
                coeff.append(sum(all(u != v and rows[u] >> v & 1 for u,v in zip(w,w[1:])) for w in words))
            check(tuple(coeff) in ((5,2,0,0),(5,3,0,0)), 'strong4 profile')
            profiles[tuple(coeff)] += 1
    check(len(canon) == 1, 'unique strong4 class')
    check(sum(profiles.values()) == 96, 'all strong4 roots')
    return {str(k):v for k,v in sorted(profiles.items())}


def witness():
    rows = (240,225,195,135,142,156,184,512,255,383)
    h = h_dp(rows)
    check(h == 4151 and h == h_reverse(rows) and h == h_literal(rows), '4151 witness triple count')
    check(prime(593), '593 prime')
    rotation = tuple(list(range(1,7))+[0,7,8,9])
    check(all((rows[i] >> j & 1) == (rows[rotation[i]] >> rotation[j] & 1)
              for i in range(10) for j in range(10)), 'C7 witness rotation')
    # Degree classes restrict every possible automorphism without losing any.
    degree_classes = {}
    for v,r in enumerate(rows):
        degree_classes.setdefault(r.bit_count(), []).append(v)
    candidates = [tuple(range(10))]
    for group in degree_classes.values():
        new = []
        for p in candidates:
            for imgs in permutations(group):
                out = list(p)
                for v,im in zip(group,imgs):
                    out[v] = im
                new.append(tuple(out))
        candidates = new
    aut = sum(all((rows[i] >> j & 1) == (rows[p[i]] >> p[j] & 1)
                  for i in range(10) for j in range(10)) for p in candidates)
    check(aut == 7, 'full witness automorphism order')
    return {'rows':list(rows),'H':h,'prime_factor':593,'automorphism_order':aut,
            'degree_permutation_candidates':len(candidates)}


def main():
    result = {'status': 'PASS', 'census': [census(n) for n in range(1,7)],
              'circulants': [circulants(p) for p in (3,5,7,11,13)],
              'strong4_profiles': strong4_coefficients(), 'witness': witness()}
    text = json.dumps(result, sort_keys=True, indent=2) + '\n'
    Path(__file__).with_name('receipt.json').write_text(text)
    print(text, end='')

if __name__ == '__main__':
    main()
