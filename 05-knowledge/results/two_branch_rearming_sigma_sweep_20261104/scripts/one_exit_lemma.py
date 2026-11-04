#!/usr/bin/env python3
"""one_exit_lemma.py (klein-2026-11-04) -- exhaustive check of the one-exit extension law.  For a tournament T on V, a
vertex a, and the extension T + x with x -> a and v -> x for every v != a (the first birth of the C2 operator of
multiport_defect_grammar_20261103/construction/RESULT.md):
    H(T + x) = h(T) - e_a(T) + [P A_a](T),   [P A_a] = sum over S subset of V - a of h(T[S]) * s_a(T[V - S]),
where e_a counts Hamiltonian paths of T ending at a and s_a those starting at a.  With the inherited sigma = [P A_a] - h
this reads H(T + x) = 2h + sigma - e_a.  Checked on every tournament of orders 3..7 (nauty gentourng) and every a."""
import subprocess, itertools, sys
GEN = sys.argv[1] if len(sys.argv) > 1 else '/opt/homebrew/bin/gentourng'
def decode(code, n):
    adj = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if code[idx] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj
def hp(adj, verts, start=None, end=None):
    n = len(verts)
    if n == 0: return 1 if (start is None and end is None) else 0
    f = [[0]*n for _ in range(1 << n)]
    for i, v in enumerate(verts):
        if start is None or v == start: f[1 << i][i] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for i in range(n):
            c = row[i]
            if not c: continue
            for j in range(n):
                if adj[verts[i]][verts[j]] and not (S >> j & 1): f[S | 1 << j][j] += c
    full = (1 << n) - 1
    return sum(f[full][i] for i, v in enumerate(verts) if end is None or v == end)
bad = checks = 0
for n in range(3, 8):
    for code in subprocess.run([GEN, '-q', str(n)], capture_output=True, text=True).stdout.split():
        adj = decode(code, n); V = list(range(n)); h = hp(adj, V)
        for a in range(n):
            e_a = hp(adj, V, end=a); others = [v for v in V if v != a]
            PA = sum(hp(adj, list(S)) * hp(adj, [v for v in V if v not in S], start=a) for r in range(len(others) + 1) for S in itertools.combinations(others, r))
            N = n + 1; A = [row[:] + [0] for row in adj] + [[0]*N]; A[n][a] = 1
            for v in range(n):
                if v != a: A[v][n] = 1
            checks += 1
            if hp(A, list(range(N))) != h - e_a + PA: bad += 1
print('one-exit law H(T+x) = h - e_a + [P A_a]: %d cases (orders 3..7, every a), mismatches %d' % (checks, bad))
