"""Exact directed path-profile and one-vertex growth probes; no assertions.

Run: python 05-knowledge/results/zero_heavy_ports_20261020g/profile_growth_check.py
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

E = module('growth_enum', HERE.parent / 'metagraph_growth_20261020/enumeration.py')
G = module('chain_growth', HERE.parent / 'metagraph_growth_20261020/growth_check.py')

def need(test, message):
    if not test:
        raise RuntimeError(message)

def path_data(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(dp[mask]):
            if value:
                for w in range(n):
                    if not mask >> w & 1 and rows[v] >> w & 1:
                        dp[mask | (1 << w)][w] += value
    return tuple(map(sum, dp))

def profile(rows):
    H = path_data(rows)
    @lru_cache(None)
    def pc(mask):
        if not mask:
            return (1,)
        ans = [0] * (mask.bit_count() + 1)
        first = mask & -mask
        sub = mask
        while sub:
            if sub & first:
                for d, c in enumerate(pc(mask ^ sub)):
                    ans[d + 1] += H[sub] * c
            sub = (sub - 1) & mask
        return tuple(ans)
    return pc((1 << len(rows)) - 1)[1:]

def port_tensors(rows):
    """Independent explicit unordered path covers and their endpoint tensors."""
    n = len(rows)
    paths = defaultdict(list)
    for mask in range(1, 1 << n):
        for p in permutations(v for v in range(n) if mask >> v & 1):
            if all(rows[a] >> b & 1 for a, b in zip(p, p[1:])):
                paths[mask].append(p)
    @lru_cache(None)
    def covers(mask):
        if not mask:
            return ((),)
        first, sub, result = mask & -mask, mask, []
        while sub:
            if sub & first:
                result.extend((p,) + rest for p in paths[sub] for rest in covers(mask ^ sub))
            sub = (sub - 1) & mask
        return tuple(result)
    c = [0] * (n + 2)
    ends = [[0] * n for _ in c]
    starts = [[0] * n for _ in c]
    B = [[[0] * n for _ in range(n)] for _ in c]
    for cover in covers((1 << n) - 1):
        d = len(cover)
        c[d] += 1
        for p in cover:
            starts[d][p[0]] += 1
            ends[d][p[-1]] += 1
        for i, p in enumerate(cover):
            for j, q in enumerate(cover):
                if i != j:
                    B[d - 1][p[-1]][q[0]] += 1
    return c, ends, starts, B

def port_response(data, mask):
    c, ends, starts, B = data
    n = len(ends[0])
    return tuple(c[d - 1] + d * c[d]
                 + sum((ends[d][u] - starts[d][u]) * (mask >> u & 1) for u in range(n))
                 + sum(B[d][u][v] for u in range(n) for v in range(n)
                       if mask >> u & 1 and not mask >> v & 1)
                 for d in range(1, n + 2))

def feedback(rows):
    n = len(rows)
    dp = [0] * (1 << n)
    for mask in range(1, 1 << n):
        dp[mask] = min(dp[mask ^ (1 << v)] + (rows[v] & (mask ^ (1 << v))).bit_count()
                       for v in range(n) if mask >> v & 1)
    return dp[-1]

def grow(reps):
    out = {}
    for rows in reps:
        for mask in range(1 << len(rows)):
            new = E.extension(rows, mask)
            key = E.canonical(new)
            out[key] = E.decode(len(new), key)
    return tuple(out[k] for k in sorted(out))

def extension_histogram(rows):
    return Counter((profile(new), tuple(sorted(map(len, G.clone_chains(new)))))
                   for mask in range(1 << len(rows))
                   for new in [E.extension(rows, mask)])

def run_partition(rows, mask):
    """Candidate blocks before the only permitted mergers, through the new x."""
    out = []
    for chain in G.clone_chains(rows):
        run = []
        for v in chain:
            if run and ((mask >> v) ^ (mask >> run[-1])) & 1:
                out.append(tuple(run)); run = []
            run.append(v)
        out.append(tuple(run))
    out.append((len(rows),))
    return tuple(out)

def run_quotient(rows, mask):
    new = E.extension(rows, mask)
    runs = run_partition(rows, mask)
    q = G.induced(new, tuple(c[0] for c in runs))
    return tuple(tuple(v for j in c for v in runs[j]) for c in G.clone_chains(q))

def main():
    counts = (1, 1, 2, 4, 12, 56, 456)
    reps = ((0,),)
    run_tests = 0
    witness = None
    for n, expected in enumerate(counts, 1):
        if n > 1:
            reps = grow(reps)
        need(len(reps) == expected, ('class census', n, len(reps)))
        groups = defaultdict(list)
        for rows in reps:
            if n <= 6:
                tensors = port_tensors(rows)
                need(tuple(tensors[0][1:n + 1]) == profile(rows), 'explicit cover audit')
                # Independent permutation cut identity checks the entire profile.
                pc = [0] * n
                from math import comb, factorial
                for p in permutations(range(n)):
                    bad = sum(not (rows[a] >> b & 1) for a, b in zip(p, p[1:]))
                    for d in range(bad + 1, n + 1):
                        pc[d - 1] += comb(n - 1 - bad, d - 1 - bad)
                need(profile(rows) == tuple(pc[d - 1] // factorial(d) for d in range(1, n + 1)),
                     'independent path-cover cut identity')
                for mask in range(1 << n):
                    actual = G.clone_chains(E.extension(rows, mask))
                    proposed = run_quotient(rows, mask)
                    need(set(map(frozenset, actual)) == set(map(frozenset, proposed)),
                         ('run quotient', rows, mask, actual, proposed))
                    need(port_response(tensors, mask) == profile(E.extension(rows, mask)),
                         ('quadratic cover response', rows, mask))
                    cut = sum(rows[u] >> v & 1 for u in range(n) for v in range(n)
                              if not mask >> u & 1 and mask >> v & 1)
                    need((port_response(tensors, mask)[0] == tensors[0][1]) == (cut == 0),
                         'repaired neutrality theorem')
                    run_tests += 1
            if len(G.clone_chains(rows)) == n:
                key = (profile(rows), tuple(sorted(r.bit_count() for r in rows)),
                       E.canonical(rows, aut=True)[1], feedback(rows))
                groups[key].append(rows)
        candidates = 0
        for key, members in groups.items():
            for i, left in enumerate(members):
                for right in members[i + 1:]:
                    if E.canonical(E.converse(left)) == E.canonical(right):
                        continue
                    candidates += 1
                    lh, rh = extension_histogram(left), extension_histogram(right)
                    if lh != rh and witness is None:
                        witness = (left, right, key, lh, rh)
        print('n', n, 'classes', len(reps), 'pairfree', sum(map(len, groups.values())),
              'nonconverse candidate pairs', candidates, flush=True)
        if witness is not None:
            break
    need(witness is not None, 'no witnessed profile-growth split in declared universe')
    left, right, key, lh, rh = witness
    print('WITNESS rows', left, right)
    print('WITNESS signature', key)
    # Independently verify all 128 extension HP counts by literal permutations.
    maxima = []
    for rows in (left, right):
        literal = []
        for mask in range(1 << len(rows)):
            new = E.extension(rows, mask)
            count = sum(all(new[a] >> b & 1 for a, b in zip(p, p[1:]))
                        for p in permutations(range(len(new))))
            need(count == profile(new)[0], 'literal extension HP audit')
            literal.append(count)
        maxima.append(max(literal))
        print('WITNESS HP extension histogram', sorted(Counter(literal).items()))
        print('WITNESS maximum masks', [i for i, value in enumerate(literal) if value == max(literal)])
        data = port_tensors(rows)
        c, ends, starts, B = data
        L = tuple(ends[1][u] - starts[1][u] + sum(B[1][u]) for u in range(6))
        C = tuple(B[1][u][v] + B[1][v][u] for u, v in E.pairs(6))
        print('WITNESS quadratic H linear', L, 'pair coefficients', C)
        need(sum(literal) == (1 << 6) * (c[1]) + (1 << 5) * c[2], 'mean response')
    need(maxima == [159, 147], 'pinned growth witness')
    need(E.slow_canonical(left) != E.slow_canonical(right), 'independent nonisomorphism')
    need(E.slow_canonical(E.converse(left)) != E.slow_canonical(right), 'independent nonconverse')
    for rows in (left, right):
        need(min(sum(rows[p[j]] >> p[i] & 1 for i in range(6) for j in range(i + 1, 6))
                 for p in permutations(range(6))) == 3, 'independent feedback')
        need(sum(E.ordercode(rows, p) == E.rowcode(rows) for p in permutations(range(6))) == 1,
             'independent automorphism')
    print('WITNESS minimum differing response', min(set(lh) | set(rh), key=lambda k: (lh[k] == rh[k], k)))
    print('RUN-WORD exact clone tests', run_tests)
    print('PASS')

if __name__ == '__main__':
    main()
