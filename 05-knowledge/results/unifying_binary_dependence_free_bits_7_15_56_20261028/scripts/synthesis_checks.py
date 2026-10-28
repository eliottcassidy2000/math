#!/usr/bin/env python3
"""synthesis_checks.py -- four finite checks behind the unifying note (klein, 2026-10-28). Standard library only.
(a) Hamming [15,11,3] = subgraphs of K6 whose six degrees have one common parity (duad model of PG(3,2));
(b) the Fano lines are the out-neighbourhoods of the Paley tournament P_7 (cyclic (7,3,1) difference set);
(c) arithmetic Hamming banks of rank 4, 5, 6: for every squarefree d over the row primes, the least prime p = d k^2 - 1
    (so the square class of sigma(p) = p + 1 is exactly d); all primes distinct, hence PG(r-1,2) is a square-class
    matroid of primes for r <= 6 and every simple binary matroid of rank <= 6 is a restriction of one of these banks;
(d) chromatic numbers of the line hypergraphs of PG(2,2), PG(3,2), PG(4,2) (= least number of sum-free sets
    partitioning F2^r minus 0), with explicit colourings."""
import itertools
from collections import Counter
out = []
def say(s): out.append(s); print(s)
# (a)
E = list(itertools.combinations(range(6), 2)); wt = Counter(); kinds = Counter()
for mask in range(1 << 15):
    deg = [0] * 6
    for i, (a, b) in enumerate(E):
        if mask >> i & 1: deg[a] += 1; deg[b] += 1
    par = {d % 2 for d in deg}
    if len(par) == 1:
        w = bin(mask).count('1'); wt[w] += 1
        if w in (3, 4): kinds[(w, 'even' if par == {0} else 'odd')] += 1
cols = list(range(1, 16)); hw = Counter()
for mask in range(1 << 15):
    s = 0
    for i in range(15):
        if mask >> i & 1: s ^= cols[i]
    if s == 0: hw[bin(mask).count('1')] += 1
say('(a) K6 equal-parity-degree subgraphs: %d words; weight distribution %s' % (sum(wt.values()), dict(sorted(wt.items()))))
say('    weight 3 = %d triangles + %d perfect matchings; weight 4 = %d four-cycles + %d claw-plus-edge' % (kinds[(3,'even')], kinds[(3,'odd')], kinds[(4,'even')], kinds[(4,'odd')]))
say('    Hamming [15,11,3] weight distribution %s; identical: %s' % (dict(sorted(hw.items())), dict(hw) == dict(wt)))
# (b)
lines = [frozenset((i + d) % 7 for d in (1, 2, 4)) for i in range(7)]
pairs = Counter(frozenset(p) for L in lines for p in itertools.combinations(sorted(L), 2))
say('(b) P_7 out-neighbourhoods {i+1,i+2,i+4}: 21 pairs each on exactly one line: %s; 3 lines per point: %s' % (set(pairs.values()) == {1} and len(pairs) == 21, set(sum(1 for L in lines if v in L) for v in range(7)) == {3}))
# (c)
def is_prime(n):
    if n < 2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
for rows in ((2,3,5,7), (2,3,5,7,11), (2,3,5,7,11,13)):
    table = []
    for m in range(1, 1 << len(rows)):
        d = 1
        for i, q in enumerate(rows):
            if m >> i & 1: d *= q
        k = 1
        while not is_prime(d * k * k - 1): k += 1
        table.append((d, k, d * k * k - 1))
    ps = [t[2] for t in table]
    say('(c) rank %d bank over row primes %s: %d primes, distinct %s, max k %d, max prime %d' % (len(rows), rows, len(table), len(set(ps)) == len(ps), max(t[1] for t in table), max(ps)))
    say('    (d, k, p = d k^2 - 1): ' + ' '.join('(%d,%d,%d)' % t for t in sorted(table)))
# (d)
def chromatic(r):
    pts = list(range(1, 1 << r)); L = [tuple(sorted(l)) for l in {frozenset((a, b, a ^ b)) for a in pts for b in pts if a != b}]
    inc = {v: [l for l in L if v in l] for v in pts}
    def colourable(k):
        col = {}
        def ok(v, c):
            for l in inc[v]:
                o = [u for u in l if u != v]
                if all(u in col and col[u] == c for u in o): return False
            return True
        def rec(i):
            if i == len(pts): return True
            v = pts[i]
            for c in range(k):
                if c > max(col.values(), default=-1) + 1: break
                if ok(v, c):
                    col[v] = c
                    if rec(i + 1): return True
                    del col[v]
            return False
        return (rec(0), dict(col))
    k = 1
    while True:
        okk, col = colourable(k)
        if okk: return k, len(L), col
        k += 1
for r in (3, 4, 5):
    k, nl, col = chromatic(r)
    classes = [sorted(v for v in col if col[v] == c) for c in range(k)]
    say('(d) PG(%d,2): %d lines, chromatic number %d; a colouring (sum-free classes): %s' % (r - 1, nl, k, classes))
open('synthesis_checks.out', 'w').write('\n'.join(out) + '\n')
