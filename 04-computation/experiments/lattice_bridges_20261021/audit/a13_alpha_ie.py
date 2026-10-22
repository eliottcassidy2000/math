#!/usr/bin/env python3
"""AUDIT pass 2, item 2 (method Z): alpha_k of circulant tournaments by ranked zeta transforms + inclusion-exclusion,
algorithmically independent of alpha_decomp.c (which uses a min-vertex path DP over all masks and the min-element
recursion F[W] over set partitions).
  (1) g(Z) = # directed Hamiltonian cycles of T[Z] for Z containing 0: Held-Karp from vertex 0 (a12_cycles hk, dump);
      f(X) = g(X - min X) for every X (rotation invariance), kept for |X| odd >= 3.
  (2) ranked zeta transform Fz[r](W) = sum_{X subset W, |X| = r} f(X)  (mod p);
  (3) Pk = Fz^k as rank polynomials (pointwise in W, truncated at rank n);
  (4) number of ORDERED k-tuples of pairwise disjoint cycle sets
        = sum_W mu(Pk[|W|])(W) = sum_Y sum_j (-1)^j C(n-|Y|, j) Pk[|Y| + j](Y),
      and alpha_k = that / k!.  Three primes < 2^31, CRT (third prime as a consistency check).
Checks: H = sum_k 2^k alpha_k against the pass-1 class tables; alpha_1 = sum_X f(X).
Usage: python3 a13_alpha_ie.py n [all | R+Emin | "s1 s2 ..."]     Output: a13_alpha_ie_n.out (+ a13_alpha_ie_n.json)
"""
import os, sys, json, subprocess
import numpy as np
from math import gcd, comb, factorial
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'a12_cycles')
PR = [2147483629, 2147483587, 2147483579]
n = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'R+Emin'
rows = [json.loads(l) for l in open(os.path.join(HERE, f'a3_classes_{n}.txt'))]
U = [a for a in range(1, n) if gcd(a, n) == 1]
canon = lambda S: min(tuple(sorted((a * s) % n for s in S)) for a in U)
m = (n - 1) // 2; R = canon(range(1, m + 1))
Emin = min(r['E'] for r in rows)
if mode == 'all': targets = rows
elif mode == 'R+Emin': targets = [r for r in rows if tuple(r['S']) == R or r['E'] == Emin]
else:
    want = canon([int(x) for x in mode.split()]); targets = [r for r in rows if tuple(r['S']) == want]
N = 1 << n; K = n // 3
pc = np.array([bin(i).count('1') for i in range(N)], dtype=np.int64)
low = np.zeros(N, dtype=np.int64)
for i in range(1, N): low[i] = (i & -i).bit_length() - 1
outname = os.path.join(HERE, f'a13_alpha_ie_{n}.out' if mode in ('all', 'R+Emin') else f'a13_alpha_ie_{n}_extra.out')
out = open(outname, 'w'); res = {}
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
def crt(rs):
    x, M = 0, 1
    for r, p in zip(rs, PR):
        t = ((r - x) * pow(M, -1, p)) % p; x += M * t; M *= p
    return x, M
for r in targets:
    S = r['S']
    pref = os.path.join(HERE, f'tmpdump_{n}')
    subprocess.run([BIN, 'hk', pref], input=' '.join(map(str, [n] + S)) + '\n', capture_output=True, text=True, check=True)
    fn = f'{pref}_n{n}.bin'; g = np.fromfile(fn, dtype=np.uint64); os.remove(fn)
    X = np.arange(N, dtype=np.int64)
    idx = np.where(X > 0, X >> (low + 1), 0)
    keep = (pc >= 3) & (pc % 2 == 1)
    alpha_mod = {p: [] for p in PR}
    tot_f = int(sum(int(v) for v in g[idx[keep]]))
    for p in PR:
        fm = np.where(keep, (g[idx] % np.uint64(p)).astype(np.int64), 0)
        Fz = np.zeros((n + 1, N), dtype=np.int64)
        for rr in range(3, n + 1, 2):
            sel = (pc == rr); Fz[rr, sel] = fm[sel]
        for i in range(n):           # zeta transform over each bit, all ranks at once
            v = Fz.reshape(n + 1, N >> (i + 1), 2, 1 << i)
            v[:, :, 1, :] += v[:, :, 0, :]
            v[:, :, 1, :] %= p
        # coefficient table (-1)^(r-|Y|) C(n-|Y|, r-|Y|) for r >= |Y|, built from Cnk below
        Cnk = np.zeros((n + 1, n + 1), dtype=np.int64)
        for a in range(n + 1):
            for b in range(n + 1): Cnk[a, b] = comb(a, b) % p if b <= a else 0
        Pk = Fz.copy(); vals = []
        for k in range(1, K + 1):
            if k > 1:      # in place, ranks from high to low (new rank rr uses only old ranks s <= rr - 3)
                for rr in range(n, 2, -1):
                    if rr < 3 * k: Pk[rr] = 0; continue
                    acc = np.zeros(N, dtype=np.int64)
                    for s in range(3 * (k - 1), rr - 2):
                        acc = (acc + (Pk[s] * Fz[rr - s]) % p) % p
                    Pk[rr] = acc
                Pk[:3] = 0
            tot = 0
            for rr in range(3 * k, n + 1):
                j = rr - pc; ok = j >= 0
                jj = np.where(ok, j, 0); sign = np.where(jj % 2 == 0, 1, p - 1)
                cc = Cnk[n - pc, jj] * ok
                term = (((cc * sign) % p) * Pk[rr]) % p
                tot = (tot + int(term.sum() % p)) % p
            vals.append((tot * pow(factorial(k), -1, p)) % p)
        alpha_mod[p] = vals
    alpha = [1]
    for k in range(K):
        x, M = crt([alpha_mod[p][k] for p in PR[:2]])
        x3 = alpha_mod[PR[2]][k]
        assert x % PR[2] == x3, 'third prime inconsistent'
        alpha.append(x)
    H = sum(a << k for k, a in enumerate(alpha))
    tag = 'R' if tuple(S) == R else ('Emin' if r['E'] == Emin else '')
    res[json.dumps(S)] = {'alpha': alpha, 'H': H, 'E': r['E'], 'tag': tag}
    say(f"n={n} S={S} {tag}: alpha = {alpha}; sum 2^k alpha_k = {H}; equals class-table H: {H == r['H']}; alpha_1 = sum f: {alpha[1] == tot_f}")
json.dump(res, open(outname.replace('.out', '.json'), 'w'))
out.close()
