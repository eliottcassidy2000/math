#!/usr/bin/env python3
"""cut_vertex_law.py (klein-2026-11-05): (1) verifies the cut-vertex product law hc(Y) = S * M * E (PROOF_CONVERSE.md, Lemma A)
on every strong tournament with a cut vertex of orders 5..8 (nauty gentourng), for every cut vertex; (2) tabulates, for
every two-cycle class of orders 6..12 (data/cls_{n}.txt) and every cut vertex, (S, E, |Q_1|, |Q_r|, r) and checks that
all middle components are singletons.  Usage: cut_vertex_law.py all N1 N2 | two N1 N2."""
import os, subprocess, sys
from collections import Counter
from functools import reduce
from cutvertex_lib import hpaths, strong, components, rows_from_t6, rows_from_d6
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data') + '/'
def law(rows, n, c):
    verts = [v for v in range(n) if v != c]
    if strong(rows, verts): return None
    comps = components(rows, verts); r = len(comps); hs = []; S = E = None
    for i, Q in enumerate(comps):
        h, SE = hpaths(rows, Q); hs.append(h)
        if i == 0: S = sum(cnt for (s, e), cnt in SE.items() if (rows[c] >> s) & 1)
        if i == r - 1: E = sum(cnt for (s, e), cnt in SE.items() if (rows[e] >> c) & 1)
    M = reduce(lambda a, b: a * b, hs[1:-1], 1)
    return S * M * E, S, E, comps
mode = sys.argv[1]
if mode == 'all':
    for n in range(int(sys.argv[2]), int(sys.argv[3]) + 1):
        out = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
        tested = bad = withcut = 0
        for t in out:
            rows = rows_from_t6(t, n)
            if not strong(rows, range(n)): continue
            H, SE = hpaths(rows, range(n)); HC = sum(cnt for (s, e), cnt in SE.items() if (rows[e] >> s) & 1) // n
            cuts = [c for c in range(n) if not strong(rows, [v for v in range(n) if v != c])]
            if cuts: withcut += 1
            for c in cuts:
                tested += 1
                if law(rows, n, c)[0] != HC: bad += 1
        print('order %d: strong classes with a cut vertex %d; (class, cut vertex) pairs tested %d; failures of hc = S*M*E: %d' % (n, withcut, tested, bad), flush=True)
else:
    for n in range(int(sys.argv[2]), int(sys.argv[3]) + 1):
        types = Counter(); ncuts = Counter(); midbig = 0; total = 0
        for l in open(D + 'cls_%d.txt' % n):
            h, d6 = l.split()
            if h != '2': continue
            rows = rows_from_d6(d6); total += 1
            cuts = [c for c in range(n) if not strong(rows, [v for v in range(n) if v != c])]
            ncuts[len(cuts)] += 1
            for c in cuts:
                hcY, S, E, comps = law(rows, n, c)
                assert hcY == 2, (n, d6, c, hcY)
                if any(len(Q) > 1 for Q in comps[1:-1]): midbig += 1
                types[(S, E, len(comps[0]), len(comps[-1]), len(comps))] += 1
        print('order %d: %d two-cycle classes; cut vertices per class %s; cut vertices with a non-singleton middle component: %d' % (n, total, dict(sorted(ncuts.items())), midbig))
        print('    (S, E, |Q_1|, |Q_r|, r) counts:', dict(sorted(types.items())))
