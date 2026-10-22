"""Unique regular completion of almost-regular tournaments: add a vertex v that beats every vertex of score n/2 and loses
to every vertex of score n/2 - 1 (Session30, SELFCONVERSE_POINTING_IDENTITY.md).  Every regular (n+1)-tournament arises,
once per vertex orbit, from the complete almost-regular class list.  Session opus-2026-10-22-S54.
Usage: python complete_regular.py n IN OUT   (IN lines: code H [aut]; OUT lines: code of the (n+1)-tournament)"""
import sys
import numpy as np

n = int(sys.argv[1]); iu = np.triu_indices(n, 1); m = n + 1; ju = np.triu_indices(m, 1)
with open(sys.argv[3], 'w') as out:
    for line in open(sys.argv[2]):
        if not line.strip():
            continue
        code = line.split()[0]
        A = np.zeros((m, m), dtype=int); b = np.array([1 if ch == '1' else 0 for ch in code])
        A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
        s = A[:n, :n].sum(axis=1)
        assert set(s) <= {n // 2, n // 2 - 1}
        for u in range(n):
            if s[u] == n // 2:
                A[n, u] = 1                      # v beats the high-score vertices
            else:
                A[u, n] = 1
        assert all(A.sum(axis=1) == n // 2)
        out.write(''.join('1' if A[i, j] else '0' for i, j in zip(*ju)) + ' 0\n')     # H not computed here
