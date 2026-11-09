#!/usr/bin/env python3
"""Audit P1: analyse the independent census (dom_census census 511 16) and compare with the session's claims.

usage: python3 analyze_census.py census_K511.txt ../../../../04-computation/experiments/paleyeven_20261109/dom_K511_summary.txt

Class line format (dom_census.c): n m d parity naut I2 edges.
"""
import sys
from fractions import Fraction
from math import gcd, factorial
from collections import defaultdict, Counter


def partitions_odd(n, maxpart=None):
    if maxpart is None:
        maxpart = n if n % 2 else n - 1
    if n == 0:
        yield []
        return
    for k in range(min(n, maxpart), 0, -1):
        if k % 2 == 1:
            for rest in partitions_odd(n - k, k):
                yield [k] + rest


def tournaments(n):
    """A000568 by Burnside over permutations with only odd cycles (only these fix tournaments)."""
    tot = Fraction(0)
    for lam in partitions_odd(n):
        t = sum((l - 1) // 2 for l in lam)
        for i in range(len(lam)):
            for j in range(i + 1, len(lam)):
                t += gcd(lam[i], lam[j])
        z = 1
        for k, m in Counter(lam).items():
            z *= k ** m * factorial(m)
        tot += Fraction(2 ** t, z)
    assert tot.denominator == 1
    return int(tot)


def euler(a, N):
    """b = Euler transform of a[1..N] (b[0] = 1)."""
    b = [1] + [0] * N
    c = [0] * (N + 1)
    for k in range(1, N + 1):
        c[k] = sum(d * a[d] for d in range(1, k + 1) if k % d == 0)
    for n in range(1, N + 1):
        b[n] = sum(c[k] * b[n - k] for k in range(1, n + 1)) // n
    return b


def inv_euler(b, N):
    a = [0] * (N + 1)
    for n in range(1, N + 1):
        a[n] = 0
        a[n] = b[n] - euler(a, n)[n]
    return a


def closure(vals, limit):
    """all products <= limit of elements of vals (with 1)."""
    S = {1}
    changed = True
    while changed:
        changed = False
        for x in list(S):
            for v in vals:
                p = x * v
                if p <= limit and p not in S:
                    S.add(p)
                    changed = True
    return S


def main():
    path = sys.argv[1]
    per_n = defaultdict(lambda: [0, 0])
    per_d = defaultdict(lambda: [0, 0])
    realisers = defaultdict(list)
    i2_even_conn = defaultdict(list)
    i2_all_conn_le7 = set()
    for line in open(path):
        p = line.split()
        n, m, d, par, naut, i2 = int(p[0]), int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5])
        edges = p[6] if len(p) > 6 else ''
        per_n[n][0] += 1
        per_d[d][0] += 1
        if par == 'E':
            per_n[n][1] += 1
            per_d[d][1] += 1
        if d in (3, 5, 7, 13, 21, 25, 37):
            realisers[d].append((n, m, par, naut, edges))
        if n <= 7:
            i2_all_conn_le7.add(i2)
            if par == 'E' and i2 <= 31:
                i2_even_conn[i2].append((n, m, edges))
    print("== per order: classes / even")
    for n in sorted(per_n):
        print(f"  n={n}: {per_n[n][0]} / {per_n[n][1]}")
    ce = [0] + [per_n[n][1] for n in range(1, 10)]
    A568 = [tournaments(n) for n in range(0, 11)]
    print("A000568 (own Burnside computation), n = 0..10:", A568)
    print("connected-even counts n = 1..9:", ce[1:])
    print("Euler transform of connected-even counts:", euler(ce, 9)[1:])
    print("matches A000568 for n = 1..9:", euler(ce, 9)[1:] == A568[1:10])
    print("inverse Euler transform of A000568 (A334335):", inv_euler(A568, 10)[1:])
    allconn = [0] + [per_n[n][0] for n in range(1, 10)]
    A88 = euler(allconn, 9)
    print("Euler transform of all connected counts (A000088 check):", A88[1:])
    # value sets
    odd = list(range(1, 512, 2))
    Veven = sorted(d for d in per_d if per_d[d][1] > 0)
    Vodd = sorted(d for d in per_d if per_d[d][0] - per_d[d][1] > 0)
    E = closure(Veven, 511)
    print("== value sets below 512")
    print("odd values NOT attained by any Royle-even graph:", [k for k in odd if k not in E])
    print("odd values NOT attained by a CONNECTED Royle-even graph:", [k for k in odd if k not in Veven])
    # odd graphs: some odd component (value a in Vodd) times any graph (every odd value, Lemma B1)
    O = {a * k for a in Vodd for k in odd if a * k <= 511}
    print("odd values NOT attained by any Royle-odd graph:", [k for k in odd if k not in O])
    print("odd values NOT attained by a CONNECTED Royle-odd graph:", [k for k in odd if k not in Vodd])
    print("all graphs: odd values not attained by connected graphs or products:",
          [k for k in odd if k not in closure(sorted(per_d), 511)])
    print("== connected realisers of selected values (n m parity |Aut| edges)")
    for d in sorted(realisers):
        print(f"  d={d}:", realisers[d])
    print("== I(G,2) on connected even graphs with n <= 7 and I <= 31:", {k: len(v) for k, v in sorted(i2_even_conn.items())})
    for k, v in sorted(i2_even_conn.items()):
        print(f"  I={k}:", v)
    Ivals = closure(sorted(i2_even_conn), 31)
    print("I(G,2) values on even graphs below 32 (products of connected-even values, 1 = null graph):", sorted(Ivals))
    print("odd values below 32 missed:", [k for k in range(1, 32, 2) if k not in Ivals])
    # compare with the session summary table
    if len(sys.argv) > 2:
        bad = 0
        rows = 0
        for line in open(sys.argv[2]):
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            d, ncl, nev = int(p[0]), int(p[1]), int(p[2])
            rows += 1
            if per_d[d][0] != ncl or per_d[d][1] != nev:
                bad += 1
                print("  MISMATCH d=", d, (ncl, nev), tuple(per_d[d]))
        print(f"== session summary table: {rows} rows compared, {bad} mismatches; audit has {len(per_d)} d-values")


if __name__ == '__main__':
    main()
