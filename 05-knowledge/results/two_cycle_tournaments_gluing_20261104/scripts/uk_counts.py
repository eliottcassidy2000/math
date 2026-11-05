#!/usr/bin/env python3
"""uk_counts.py (klein-2026-11-04): t_k(n) = number of tournaments of order n with exactly k Hamiltonian cycles and
u_k(n) = those without a vertex of out-degree 1, from the complete generator's class lists data/clsK3_{n}.txt (k <= 3) and
data/clsK4_{n}.txt (k = 4); tests HYP-9459 (2) u_3(n) = F(2n-9) and (3) u_4(n) = 2 F(2n-11)."""
import os, sys
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
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    out = []
    for K, fn in ((3, 'clsK3_%d.txt'), (4, 'clsK4_%d.txt')):
        p = os.path.join(DATA, fn % n)
        if not os.path.exists(p): continue
        t = {}; u = {}
        for l in open(p):
            h, d6 = l.split(); h = int(h)
            if K == 4 and h != 4: continue
            rows = rows_from_d6(d6); t[h] = t.get(h, 0) + 1
            if all(bin(r).count('1') != 1 for r in rows): u[h] = u.get(h, 0) + 1
        for h in sorted(t): out.append('t_%d = %d, u_%d = %d' % (h, t[h], h, u.get(h, 0)))
    print('order %d: %s; F(2n-9) = %d, 2F(2n-11) = %d' % (n, '; '.join(out), F[2 * n - 9], 2 * F[2 * n - 11] if 2 * n - 11 >= 0 else 0))
