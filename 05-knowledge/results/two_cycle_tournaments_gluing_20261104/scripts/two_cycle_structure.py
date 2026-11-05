#!/usr/bin/env python3
"""two_cycle_structure.py (klein-2026-11-04): structure of every tournament with exactly two Hamiltonian cycles, from the
complete generator's class lists data/cls_{n}.txt ("hc digraph6", canonical).  Usage: two_cycle_structure.py N1 N2.
For each two-cycle class Y of order n in [N1, N2] it enumerates both Hamiltonian cycles C1, C2 (subset DP with
backtracking) and reports, per order:
  swap:  the two cycles differ in exactly three arcs; the three runs of C1 between the changed arcs have sorted
         sizes (j, k, l) with j in {1, 2} (the junction) -- counts by sorted triple; cut arcs lie on both cycles;
         the non-critical vertices form a cyclic interval of exactly one of the two cycles; #cut + #NC = n.
  glue:  for a one-vertex junction s with C1 = P s S2 R, C2 = P S2 s R: Y[S2 + s] and Y[P + R + s] are uniquely
         Hamiltonian with s a sole entry of the first and a sole exit of the second; the only arc from S2 to P u R is
         S2_last -> R_first; the rooted pieces determine Y (distinct pairs = classes) and every pair occurs
         (F(2k-3) F(2l-3) pairs for |S2| = k, |P u R| = l).
  arc:   for a two-vertex junction s -> s2: Y[S2 + s + s2] and Y[P + R + s + s2] are uniquely Hamiltonian, s, s2 are
         consecutive sole entries of the first and consecutive sole exits of the second, and the only arc from S2 to
         P u R is S2_last -> R_first (the mirror reading, with the roles of the two runs exchanged, is reported as such)."""
import sys, os
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
F = [0, 1]
for _ in range(60): F.append(F[-1] + F[-2])
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
def sub_rows(rows, verts):
    return [sum(((rows[v] >> u) & 1) << j for j, u in enumerate(verts)) for v in verts]
def strong_digraph(adj, n):
    for forward in (True, False):
        seen = {0}; fr = [0]
        while fr:
            nf = []
            for v in fr:
                for u in range(n):
                    if u not in seen and ((adj[v] >> u) & 1 if forward else (adj[u] >> v) & 1): seen.add(u); nf.append(u)
            fr = nf
        if len(seen) < n: return False
    return True
def cycles(rows, n):
    N = 1 << n; dp = [[0] * n for _ in range(N)]; dp[1][0] = 1
    for S in range(1, N, 2):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nx = rows[v] & ~S
            while nx:
                u = (nx & -nx).bit_length() - 1; nx &= nx - 1; dp[S | (1 << u)][u] += c
    res = []
    def back(S, v, path):
        if S == 1: res.append(path[::-1]); return
        Sp = S & ~(1 << v)
        for u in range(n):
            if (Sp >> u) & 1 and dp[Sp][u] and (rows[u] >> v) & 1: back(Sp, u, path + [u])
    for v in range(1, n):
        if (rows[v] & 1) and dp[N - 1][v]: back(N - 1, v, [v])
    return res
def runs_of(rows, n, C):
    A = [set((c[i], c[(i + 1) % n]) for i in range(n)) for c in C]
    c = C[0]; diff = A[0] - A[1]; chg = sorted(i for i in range(n) if (c[i], c[(i + 1) % n]) in diff)
    assert len(chg) == 3
    runs = []
    for j in range(3):
        a = chg[j] + 1; b = chg[(j + 1) % 3]; L = (b - a) % n + 1; runs.append([c[(a + t) % n] for t in range(L)])
    return runs, A
def rooted_form(rows, order):
    k = len(order); return tuple(tuple((rows[order[i]] >> order[j]) & 1 for j in range(k)) for i in range(k))
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    triples = Counter(); cutsub = Counter(); ncint = Counter(); ident = Counter()
    gl_hc = Counter(); gl_roles = Counter(); gl_cross = Counter(); pairs = defaultdict(set); ys = Counter()
    ar_hc = Counter(); ar_roles = Counter(); ar_cross = Counter(); ar_pairs = defaultdict(set); ar_ys = Counter()
    total = 0
    for l in open(os.path.join(DATA, 'cls_%d.txt' % n)):
        h, d6 = l.split()
        if h != '2': continue
        total += 1; rows = rows_from_d6(d6); C = cycles(rows, n); assert len(C) == 2
        runs, A = runs_of(rows, n, C); sizes = [len(r) for r in runs]; triples[tuple(sorted(sizes))] += 1
        shared = A[0] & A[1]
        cut = set()
        for u in range(n):
            for w in range(n):
                if (rows[u] >> w) & 1:
                    r2 = rows[:]; r2[u] &= ~(1 << w)
                    if not strong_digraph(r2, n): cut.add((u, w))
        NC = set(v for v in range(n) if strong_digraph(sub_rows(rows, [x for x in range(n) if x != v]), n - 1))
        cutsub[cut <= shared] += 1; ident[len(cut) + len(NC) - n] += 1
        def is_interval(c, S):
            inS = [c[i] in S for i in range(n)]; return sum(1 for i in range(n) if inS[i] != inS[(i + 1) % n]) == 2
        ncint[(is_interval(C[0], NC), is_interval(C[1], NC))] += 1
        j = sizes.index(min(sizes)); seg = runs[j]; S2 = runs[(j + 1) % 3]; PR = runs[(j + 2) % 3]; k, m = len(S2), len(PR)
        if len(seg) == 1:
            s = seg[0]
            hA = len(cycles(sub_rows(rows, [s] + S2), k + 1)); hB = len(cycles(sub_rows(rows, [s] + PR), m + 1))
            gl_hc[(hA, hB)] += 1
            inA = sum(1 for u in S2 if (rows[u] >> s) & 1); outB = sum(1 for u in PR if (rows[s] >> u) & 1)
            gl_roles[('in-degree of s in Y[S2+s] = %d' % inA, 'out-degree of s in Y[PR+s] = %d' % outB)] += 1
            cross = [(S2.index(x) + 1, PR.index(y) + 1) for x in S2 for y in PR if (rows[x] >> y) & 1]
            gl_cross['only S2_%d -> PR_1' % k if cross == [(k, 1)] else str(cross)] += 1
            pairs[(k, m)].add((rooted_form(rows, [s] + S2), rooted_form(rows, [s] + PR))); ys[(k, m)] += 1
        elif len(seg) == 2:
            # the junction is a run of size 2; when two runs have size 2 (triples (2, 2, n-4)) try both, and both readings
            found = None
            for jj in range(3):
                if len(runs[jj]) != 2: continue
                s, s2 = runs[jj]
                for (S2c, PRc) in ((runs[(jj + 1) % 3], runs[(jj + 2) % 3]), (runs[(jj + 2) % 3], runs[(jj + 1) % 3])):
                    Av = S2c + [s, s2]; Bv = PRc + [s, s2]
                    inA = tuple(sum(1 for u in Av if (rows[u] >> v) & 1) for v in (s, s2)); outB = tuple(sum(1 for u in Bv if (rows[v] >> u) & 1) for v in (s, s2))
                    if inA == (1, 1) and outB == (1, 1) and (rows[s] >> s2) & 1:
                        found = (s, s2, S2c, PRc); break
                if found: break
            if not found:
                ar_roles['NO READING WITH s, s2 SOLE ENTRIES / SOLE EXITS'] += 1
            else:
                s, s2, S2c, PRc = found; k, m = len(S2c), len(PRc)
                hA = len(cycles(sub_rows(rows, S2c + [s, s2]), k + 2)); hB = len(cycles(sub_rows(rows, PRc + [s, s2]), m + 2)); ar_hc[(hA, hB)] += 1
                ar_roles['s, s2 consecutive sole entries of Y[S2+s+s2] and consecutive sole exits of Y[PR+s+s2]'] += 1
                cross = [(S2c.index(x) + 1, PRc.index(y) + 1) for x in S2c for y in PRc if (rows[x] >> y) & 1]
                ar_cross['only S2_%d -> PR_1' % k if cross == [(k, 1)] else str(cross)] += 1
                ar_pairs[(k, m)].add((rooted_form(rows, [s, s2] + S2c), rooted_form(rows, [s, s2] + PRc))); ar_ys[(k, m)] += 1
        else:
            triples[('JUNCTION LARGER THAN 2',)] += 1
    print('order %d: %d two-cycle classes' % (n, total))
    print('  swap: sorted run sizes %s' % dict(sorted(triples.items())))
    print('        cut arcs on both cycles: %s; NC a cyclic interval of (C1, C2): %s; #cut + #NC - n: %s' % (dict(cutsub), dict(ncint), dict(ident)))
    print('  glue (vertex junction): (hc of the two pieces) %s; roles %s; cross arcs %s' % (dict(gl_hc), dict(gl_roles), dict(gl_cross)))
    print('        classes by (|S2|, |PR|): %s; distinct rooted pairs: %s; F(2k-3)F(2l-3): %s' % (dict(sorted(ys.items())), {q: len(v) for q, v in sorted(pairs.items())}, {q: F[2 * q[0] - 3] * F[2 * q[1] - 3] for q in sorted(ys)}))
    print('  arc junction: (hc of the two pieces) %s; roles %s; cross arcs %s' % (dict(ar_hc), dict(ar_roles), dict(ar_cross)))
    print('        classes by (|S2|, |PR|): %s; distinct rooted pairs: %s; F(2k-3)F(2l-3): %s' % (dict(sorted(ar_ys.items())), {q: len(v) for q, v in sorted(ar_pairs.items())}, {q: F[2 * q[0] - 3] * F[2 * q[1] - 3] for q in sorted(ar_ys)}))
