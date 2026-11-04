#!/usr/bin/env python3
"""lemma2_check.py (klein-2026-11-04) -- exhaustive check of Lemma 2 (PROOFS.md section 2): for every tournament U of
orders 3..6 and every r: C(U + a_exit) = s_r(U) and C(U + a_entry) = e_r(U), where C counts directed Hamiltonian cycles."""
import subprocess, sys
GEN = sys.argv[1] if len(sys.argv) > 1 else '/opt/homebrew/bin/gentourng'
def decode(code, n):
    adj = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if code[idx] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj
def paths(adj, verts, start=None, end=None):
    n = len(verts); f = [[0]*n for _ in range(1 << n)]
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
def cycles(adj):
    n = len(adj); V = list(range(1, n)); tot = 0
    for v in V:
        if adj[0][v]:
            # paths of adj[V] starting at v ending at a vertex that beats 0
            m = len(V); f = [[0]*m for _ in range(1 << m)]; f[1 << V.index(v)][V.index(v)] = 1
            for S in range(1, 1 << m):
                row = f[S]
                for i in range(m):
                    c = row[i]
                    if not c: continue
                    for j in range(m):
                        if adj[V[i]][V[j]] and not (S >> j & 1): f[S | 1 << j][j] += c
            tot += sum(f[(1 << m) - 1][i] for i in range(m) if adj[V[i]][0])
    return tot
bad = checks = 0
for n in range(3, 7):
    for code in subprocess.run([GEN, '-q', str(n)], capture_output=True, text=True).stdout.split():
        U = decode(code, n); V = list(range(n))
        for r in range(n):
            for kind in ('exit', 'entry'):
                N = n + 1; A = [row[:] + [0] for row in U] + [[0]*N]
                if kind == 'exit':
                    A[n][r] = 1
                    for v in range(n):
                        if v != r: A[v][n] = 1
                    expect = paths(U, V, start=r)
                else:
                    A[r][n] = 1
                    for v in range(n):
                        if v != r: A[n][v] = 1
                    expect = paths(U, V, end=r)
                checks += 1
                if cycles(A) != expect: bad += 1
print('Lemma 2: C(U + a_exit) = s_r(U), C(U + a_entry) = e_r(U): %d cases (orders 3..6, all r), mismatches %d' % (checks, bad))
