#!/usr/bin/env python3
"""AUDIT of the transitive-run table (paley_transitive_runs_20261014.out): trans(P_p) and n(p).

Different algorithm from the author's branch and bound from source 0:
  * Aff(P_p) is transitive on ARCS, so a maximum transitive subtournament may be assumed to start 0 -> 1;
    trans(P_p) = 2 + t(N),  N = out(0) & out(1);
  * t(C) = max_{v in C} 1 + t(C & out(v)) computed EXACTLY by memoised recursion on bitmasks (no pruning heuristics),
    with the trivial cut t(C) = |C| when C is itself transitive... (only memo + exhaustive recursion).
Also checks trans >= n(p) (initial run 0..n(p)-1 transitive) and the discrepancy bound trans <= sqrt(4p+1).
Usage: python audit_transitive.py PMAX
"""
import sys, math, time
from functools import lru_cache

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
CLAIM = {3: 2, 7: 3, 11: 4, 19: 5, 23: 5, 31: 7, 43: 7, 47: 7, 59: 9, 67: 8, 71: 9, 79: 9, 83: 8, 103: 11,
         107: 9, 127: 11, 131: 11, 139: 11, 151: 11, 163: 11, 167: 11, 179: 12, 191: 11, 199: 11, 211: 13,
         223: 13, 227: 13, 239: 12, 251: 12, 263: 13, 271: 12, 283: 14}
NCLAIM = {3: 2, 7: 3, 11: 2, 19: 2, 23: 5, 31: 3, 43: 2, 47: 5, 59: 2, 67: 2, 71: 7, 79: 3, 83: 2, 103: 3, 107: 2,
          127: 3, 131: 2, 139: 2, 151: 3, 163: 2, 167: 5, 179: 2, 191: 7, 199: 3, 211: 2, 223: 3, 227: 2, 239: 7,
          251: 2, 263: 5, 271: 3, 283: 2}


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def trans(p):
    qr = {x * x % p for x in range(1, p)}
    out = [sum(1 << ((v + r) % p) for r in qr) for v in range(p)]

    @lru_cache(maxsize=None)
    def t(C):
        if C == 0:
            return 0
        best = 1
        c = C
        # upper bound for early exit: t(C) <= popcount(C)
        pc = bin(C).count('1')
        while c:
            low = c & -c
            v = low.bit_length() - 1
            c ^= low
            sub = C & out[v]
            if 1 + bin(sub).count('1') <= best:
                continue
            val = 1 + t(sub)
            if val > best:
                best = val
                if best == pc:
                    break
        return best
    if p == 3:
        return 2
    N = out[0] & out[1]
    val = 2 + t(N)
    t.cache_clear()
    return val


def least_nonresidue(p):
    for a in range(2, p):
        if pow(a, (p - 1) // 2, p) == p - 1:
            return a


def initial_run_transitive(p, n):
    qr = {x * x % p for x in range(1, p)}
    return all(((j - i) % p) in qr for i in range(n) for j in range(i + 1, n))


def main():
    bad = []
    for p in primes_upto(PMAX):
        if p % 4 != 3:
            continue
        t0 = time.time()
        tr = trans(p)
        n = least_nonresidue(p)
        ok_run = initial_run_transitive(p, n) and not initial_run_transitive(p, n + 1)
        ok = (CLAIM.get(p) == tr) and (NCLAIM.get(p) == n)
        if not ok:
            bad.append(p)
        print(f'p = {p:4d}: n(p) = {n}, trans = {tr} (claimed {CLAIM.get(p)}, n claimed {NCLAIM.get(p)}): {"OK" if ok else "MISMATCH"}; '
              f'initial run 0..n-1 transitive and 0..n not: {ok_run}; trans <= sqrt(4p+1) = {math.sqrt(4 * p + 1):.1f}: {tr <= math.sqrt(4 * p + 1)}; '
              f'{time.time() - t0:.1f} s', flush=True)
    print('mismatches:', bad if bad else 'none')


if __name__ == '__main__':
    main()
