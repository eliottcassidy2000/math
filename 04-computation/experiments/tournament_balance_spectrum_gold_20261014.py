#!/usr/bin/env python3
"""Balance spectrum of tournaments, the classes attaining 1/3, and a gold-partition test (tournament analogues of the
1/3-2/3 and Peczarski gold partition conjectures).

delta(S) for a set S of permutations = max over pairs x != y of min(p, 1 - p), p = Pr_{sigma in S}[x before y].

Part 1 (all labelled tournaments, n <= 7; a canonical form via all n! relabellings for the classes):
  - the multiset of delta values over strong tournaments (the analogue of posets that are not nontrivial ordinal sums);
  - the least delta > 1/3 overall and over strong tournaments;
  - the number of isomorphism classes with delta = 1/3, compared with A000930(n) - 1 (Narayana's cows,
    compositions of n into parts 1 and 3, minus the all-singletons composition).
Part 2 (closed family): for a tournament T, F(T) = {HP(T) cap Ext(Q)} over posets Q, generated from HP(T) by
  intersecting with the half-spaces A_xy = {x before y}.  Tested on every member S with |S| >= 2:
  (i)  delta(S) >= 1/3                                (tournament + precedence 1/3-2/3)
  (ii) gold partition: there are a first query (a,b) and, for each answer, a second query (c,d) such that every leaf
       satisfies |S| >= |S_1| + |S_12|                 (Peczarski's t0 >= t1 + t2).
  Run for all tournaments up to isomorphism for n <= 6, and for semicomplete digraphs (2-cycles allowed) n <= 5;
  the complete digraph gives back the poset conjectures.
Session opus-2026-10-14-S35.
"""
import itertools, sys, time
from fractions import Fraction
from collections import Counter

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def ham_paths(n, arc):
    """arc[u][v] True if u -> v.  Returns list of tuples."""
    out = []
    def rec(path, used):
        if len(path) == n:
            out.append(tuple(path)); return
        u = path[-1]
        for v in range(n):
            if not used >> v & 1 and arc[u][v]:
                path.append(v); rec(path, used | 1 << v); path.pop()
    for s in range(n):
        rec([s], 1 << s)
    return out


def positions(perm):
    pos = [0] * len(perm)
    for i, v in enumerate(perm):
        pos[v] = i
    return pos


def delta(S, n):
    if len(S) < 2:
        return Fraction(0)
    P = [positions(s) for s in S]
    best = Fraction(0)
    N = len(S)
    for x in range(n):
        for y in range(x + 1, n):
            c = sum(1 for p in P if p[x] < p[y])
            m = Fraction(min(c, N - c), N)
            if m > best:
                best = m
    return best


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        arc = [[False] * n for _ in range(n)]
        for e, (i, j) in enumerate(pairs):
            if code >> e & 1:
                arc[i][j] = True
            else:
                arc[j][i] = True
        yield code, arc


def canon(n, arc):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(arc[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key > best:
            best = key
    return best


def strong(n, arc):
    # strongly connected iff vertex 0 reaches all and all reach 0
    def reach(adjf):
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in range(n):
                if v not in seen and adjf(u, v):
                    seen.add(v); st.append(v)
        return len(seen) == n
    return reach(lambda u, v: arc[u][v]) and reach(lambda u, v: arc[v][u])


def narayana(n):
    c = [1, 1, 1]
    while len(c) <= n:
        c.append(c[-1] + c[-3])
    return c[n]


def part1(nmax):
    print('PART 1: balance spectrum of tournaments')
    for n in range(3, nmax + 1):
        t0 = time.time()
        seen = {}
        for code, arc in tournaments(n):
            k = canon(n, arc) if n <= 6 else None
            if n <= 6:
                if k in seen:
                    continue
            S = ham_paths(n, arc)
            d = delta(S, n)
            st = strong(n, arc)
            if n <= 6:
                seen[k] = (d, st, len(S))
            else:
                seen[code] = (d, st, len(S))
        vals = list(seen.values())
        if n <= 6:
            at13 = sum(1 for d, st, h in vals if d == Fraction(1, 3))
            strong_vals = sorted(set(d for d, st, h in vals if st))
            above = sorted(set(d for d, st, h in vals if d > Fraction(1, 3)))
            print(f' n = {n}: {len(vals)} classes; classes with delta = 1/3: {at13} (A000930(n) - 1 = {narayana(n) - 1}); '
                  f'least delta > 1/3: {above[0] if above else None}; least delta over strong classes: '
                  f'{strong_vals[0] if strong_vals else None}; strong delta values: {[str(v) for v in strong_vals[:8]]}  ({time.time() - t0:.1f} s)', flush=True)
        else:
            dvals = Counter(d for d, st, h in vals)
            strong_min = min(d for d, st, h in vals if st)
            above = min(d for d in dvals if d > Fraction(1, 3))
            print(f' n = {n} (labelled): least delta > 1/3: {above}; least delta over strong: {strong_min}; '
                  f'labelled tournaments at 1/3: {dvals[Fraction(1, 3)]}  ({time.time() - t0:.1f} s)', flush=True)


def closed_family(HP, n):
    """all HP cap Ext(Q): BFS by intersecting with x-before-y half-spaces; sets as frozensets of path indices"""
    P = [positions(s) for s in HP]
    half = {}
    for x in range(n):
        for y in range(n):
            if x != y:
                half[(x, y)] = frozenset(i for i, p in enumerate(P) if p[x] < p[y])
    start = frozenset(range(len(HP)))
    fam = {start}; queue = [start]
    while queue:
        S = queue.pop()
        for h in half.values():
            T = S & h
            if T and T not in fam:
                fam.add(T); queue.append(T)
    return fam, half


def gold_ok(S, half, pairs):
    t0 = len(S)
    for (a, b) in pairs:
        S1a, S1b = S & half[(a, b)], S & half[(b, a)]
        ok = True
        for S1 in (S1a, S1b):
            if not S1:
                continue
            t1 = len(S1)
            if t1 == 1:
                if t0 < 2:
                    ok = False
                continue
            found = False
            for (c, d) in pairs:
                worst = max(len(S1 & half[(c, d)]), len(S1 & half[(d, c)]))
                if t0 >= t1 + worst:
                    found = True; break
            if not found:
                ok = False; break
        if ok:
            return True
    return False


def part2(nmax_t, nmax_s):
    print('PART 2: closed families (precedence constraints), 1/3-2/3 and gold partition')
    for kind, nmax in (('tournament', nmax_t), ('semicomplete', nmax_s)):
        for n in range(3, nmax + 1):
            t0 = time.time()
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            seen = set(); worst = Fraction(1, 2); worst_ex = None; gold_fail = 0; members = 0; classes = 0
            choices = (0, 1) if kind == 'tournament' else (0, 1, 2)
            for code in itertools.product(choices, repeat=len(pairs)):
                arc = [[False] * n for _ in range(n)]
                for e, (i, j) in enumerate(pairs):
                    if code[e] in (0, 2): arc[i][j] = True
                    if code[e] in (1, 2): arc[j][i] = True
                k = canon(n, arc)
                if k in seen:
                    continue
                seen.add(k); classes += 1
                HP = ham_paths(n, arc)
                if len(HP) < 2:
                    continue
                fam, half = closed_family(HP, n)
                for S in fam:
                    if len(S) < 2:
                        continue
                    members += 1
                    Sp = [HP[i] for i in S]
                    d = delta(Sp, n)
                    if d < worst:
                        worst, worst_ex = d, (code, sorted(Sp))
                    if not gold_ok(S, half, pairs):
                        gold_fail += 1
            print(f' {kind} n = {n}: {classes} digraph classes, {members} family members with |S| >= 2; '
                  f'min delta = {worst}; gold-partition failures = {gold_fail}  ({time.time() - t0:.1f} s)', flush=True)
            if worst < Fraction(1, 3):
                print('   counterexample to delta >= 1/3:', worst_ex)


if __name__ == '__main__':
    part1(NMAX)
    part2(int(sys.argv[2]) if len(sys.argv) > 2 else 5, int(sys.argv[3]) if len(sys.argv) > 3 else 4)
