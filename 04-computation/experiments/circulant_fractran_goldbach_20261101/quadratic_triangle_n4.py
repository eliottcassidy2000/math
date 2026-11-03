#!/usr/bin/env python3
"""mac-mini-2026-11-02-triangles: the owner's pasted prompt of 2026-11-02 - the map x -> x^2 - 29/16 on odd/4
against the n = 4 tiling metagraph (8 tilings -> 4 classes -> 3 merged nodes), primes modulo even numbers until
they reach +-1, and the jagged triangles (odd pairs summing to an even number; PRIMEGAME's trial-division region).

Inherited (credited, not re-proved): THM-4139 (codex 2025-08-25): the complete rational preperiodic graph of
x -> x^2 - 29/16 is the eight quarter-integers and the 3-cycle -7/4 -> 5/4 -> -1/4 (classical: Poonen 1998).
Integer form used here: x = a/4 with a odd gives f(x) = g(a)/4, g(a) = (a^2 - 29)/4, again odd.

PRE-REGISTERED (before the first run):
  Q1  For every odd a with |a| <= 100001, the orbit of a/4 is bounded iff |a| <= 7, and the eight bounded points
      reach the cycle in at most two steps (depths: cycle 0; -5/4... computed and printed).
  Q2  The image chain of the bounded set B: |B| = 8, |f(B)| = 4, |f^2(B)| = 3 = |f^3(B)|; fibre shapes (2,2,2,2) then
      (2,1,1).
  M1  n = 4: 8 tilings -> 4 classes with fibre sizes {5,1,1,1} -> 3 merged nodes with fibre sizes {5,2,1}; 4 lines.
  M2  Hence no bijection tilings <-> B carries the class partition onto the f-fibre partition (shapes differ); a
      bijection carrying complement-lines onto {x, -x} pairs exists (any two fixed-point-free involutions on 8
      points are conjugate) - GENERIC.
  M3  No bijection between the 4 classes and f(B) carries the class line-multigraph (undirected, loops kept) onto
      the undirected functional graph of f on f(B).
  P1  For every even m with 4 <= m <= 20000, the least prime p = +-1 (mod m) satisfies p < m^2.
  P3  The clique-carrier formula of chart_clique_norm_20261103 (r = 3, forbidden sizes (3,3,5)) gives 10, 20, 39 at
      n = 3, 4, 5.
usage: python3 quadratic_triangle_n4.py
"""
import importlib.util
import itertools
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def f(x):
    return x * x - Fraction(29, 16)


def g(a):
    assert a % 2 == 1 or a % 2 == -1
    return (a * a - 29) // 4


def bounded_int(a, limit=10 ** 6, steps=60):
    seen = set()
    for _ in range(steps):
        if abs(a) > limit:
            return False, None
        if a in seen:
            return True, a
        seen.add(a)
        a = g(a)
    return abs(a) <= limit, None


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return s


def main():
    print('Q  x -> x^2 - 29/16 on odd/4 (integer form g(a) = (a^2 - 29)/4 on odd a)')
    cyc = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]
    check('cycle -7/4 -> 5/4 -> -1/4 -> -7/4', [f(c) for c in cyc] == [cyc[1], cyc[2], cyc[0]])
    bset = []
    ok_q1 = True
    for a in range(-100001, 100002, 2):
        b, _ = bounded_int(a)
        if b:
            bset.append(a)
        if b != (abs(a) <= 7):
            ok_q1 = False
    depth = {}
    for a in bset:
        x, d = Fraction(a, 4), 0
        while x not in cyc:
            x, d = f(x), d + 1
        depth[a] = d
    check('Q1: bounded iff |a| <= 7 (odd a, |a| <= 100001); all enter the cycle within 2 steps',
          ok_q1 and sorted(bset) == [-7, -5, -3, -1, 1, 3, 5, 7] and max(depth.values()) <= 2,
          f'depths {dict(sorted(depth.items()))}')
    B = {Fraction(a, 4) for a in bset}
    fB = {f(x) for x in B}
    f2B = {f(x) for x in fB}
    f3B = {f(x) for x in f2B}
    sh1 = sorted((sum(1 for x in B if f(x) == y) for y in fB), reverse=True)
    sh2 = sorted((sum(1 for x in fB if f(x) == y) for y in f2B), reverse=True)
    check('Q2: |B|, |f(B)|, |f^2 B|, |f^3 B| = 8, 4, 3, 3; fibre shapes (2,2,2,2), (2,1,1)',
          (len(B), len(fB), len(f2B), len(f3B)) == (8, 4, 3, 3) and sh1 == [2, 2, 2, 2] and sh2 == [2, 1, 1],
          f'f(B) = {sorted(fB)}, f^2(B) = {sorted(f2B)}')
    print('M  the n = 4 tiling metagraph (owner explorer gauge)')
    spec = importlib.util.spec_from_file_location(
        'mm', os.path.join(HERE, '..', 'groups_of_five_primegame_20261031', 'merged_metagraph_n5.py'))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    D, classes, lines, disp, ci2d, medges, cen = mm.analyse(4, verbose=False)
    fib = sorted((c['fibre'] for c in classes), reverse=True)
    mfib = sorted((sum(classes[c]['fibre'] for c in node['members']) for node in disp), reverse=True)
    check('M1: 8 tilings -> 4 classes {5,1,1,1} -> 3 merged {5,2,1}; 4 lines',
          len(D['tilings']) == 8 and fib == [5, 1, 1, 1] and mfib == [5, 2, 1] and len(lines) == 4,
          f'classes {[(c["ci"], c["scores"], c["fibre"], c["type"]) for c in classes]}')
    check('M2: class fibre shape {5,1,1,1} differs from the f-fibre shape {2,2,2,2}', fib != sh1)
    # M3: class line multigraph vs functional graph of f on f(B)
    cl_edges = {}
    for L in lines:
        e = tuple(sorted((L['ca'], L['cb'])))
        cl_edges[e] = cl_edges.get(e, 0) + 1
    pts = sorted(fB)
    f_edges = {}
    for x in pts:
        e = tuple(sorted((pts.index(x), pts.index(f(x))))) if f(x) in fB else None
        if e:
            f_edges[e] = f_edges.get(e, 0) + 1
    iso = False
    for perm in itertools.permutations(range(4)):
        mapped = {}
        for (u, v), mlt in cl_edges.items():
            e = tuple(sorted((perm[u], perm[v])))
            mapped[e] = mapped.get(e, 0) + mlt
        if mapped == f_edges:
            iso = True
    print(f'    class line multigraph {cl_edges}; f on f(B) (undirected) {f_edges} on {pts}')
    check('M3: no bijection carries the class line multigraph onto f\'s functional graph on f(B)', not iso)
    print('P  primes modulo even numbers until they reach +-1')
    N = 20000
    S = sieve(N * N if N <= 200 else 4 * 10 ** 7)
    limit = len(S) - 1
    rows = []
    ok_p1 = True
    near = 0
    for m in range(4, N + 1, 2):
        # candidates = +-1 (mod m) in increasing order: m-1, m+1, 2m-1, 2m+1, ... (1 is not prime)
        k, p = 1, None
        while p is None:
            for c in (k * m - 1, k * m + 1):
                if c > limit:
                    ok_p1 = False
                    p = c
                    break
                if S[c]:
                    p = c
                    break
            k += 1
        if p >= m * m:
            ok_p1 = False
        if p in (m - 1, m + 1):
            near += 1
        rows.append((m, p))
    check('P1: least prime = +-1 (mod m) is below m^2 for every even m in [4, 20000]', ok_p1,
          f'P(m) in {{m-1, m+1}} for {near} of {len(rows)} m ({near / len(rows):.1%}); largest P(m)/m = '
          f'{max(p / m for m, p in rows):.1f} at m = {max(rows, key=lambda r: r[1] / r[0])[0]}')
    print('    first rows (m, least prime = +-1 mod m, residues of the odd primes seen before it):')
    for m, p in rows[:9]:
        seen = [q % m for q in range(3, p) if S[q]]
        pairs = sorted({tuple(sorted((r, m - r))) for r in seen if (m - r) in seen})
        gold = [(q, m - q) for q in range(3, m // 2 + 1) if S[q] and S[m - q]]
        print(f'      m={m:>3}: P={p:>3}; residues {seen}; residue pairs r + (m - r) {pairs}; Goldbach pairs {gold}')
    print(f'    {{10, 20}}: P(10) = {dict(rows)[10]}, P(20) = {dict(rows)[20]}')

    def carrier(n, r=3, t=(3, 3, 5)):
        return 1 + (r - 1) * math.comb(n, 2) + sum(math.comb(n, k) for ti in t for k in range(3, min(n, ti - 1) + 1)) \
            + sum(1 for ti in t if n >= ti)
    check('P3: chart_clique_norm carrier formula gives 10, 20, 39 at n = 3, 4, 5',
          [carrier(n) for n in (3, 4, 5)] == [10, 20, 39], f'{[carrier(n) for n in (3, 4, 5, 6)]}')
    print('T  the two jagged triangles')
    print('    odd-pair triangle, rows m = 2k, entry a (odd, a <= k) marked * if a and m - a are both prime:')
    for m in range(6, 41, 2):
        row = ''.join('*' if (S[a] and S[m - a]) else '.' for a in range(1, m // 2 + 1, 2))
        print(f'      m={m:>2}: {row}')
    print('    PRIMEGAME trial-division region {(n, d): D(n) <= d <= n-1}; row n lists the trial divisors tested:')
    for n in range(2, 22):
        spf = next(p for p in range(2, n + 1) if n % p == 0)
        Dn = n // spf
        print(f'      n={n:>2}: d = {n - 1}..{Dn}  ({n - Dn} trials){"  <- prime (reaches d = 1)" if Dn == 1 else ""}')
    print()
    print('ALL PRE-REGISTERED CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
