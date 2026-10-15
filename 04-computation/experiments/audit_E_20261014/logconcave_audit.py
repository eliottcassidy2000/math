#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4655 statement 3: log-concavity of position and gap distributions of uniform Hamiltonian
paths of strong tournaments.  Independent code: digraph6 decoding (gentourng -z -c), own DFS path enumeration.
A sequence a_1..a_k is log-concave iff a_i^2 >= a_(i-1) a_(i+1) for all interior i and there are no internal zeros.
  position distribution of v: N_k = #paths with v at position k (k = 1..n);
  gap distribution of an ordered pair (x, y): N_k = #paths with pos(y) - pos(x) = k (k = 1..n-1)  [Kahn-Saks form]."""
import subprocess, sys

def d6decode(s):
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits.extend((v >> (5 - t)) & 1 for t in range(6))
    out = [[bits[i * n + j] for j in range(n)] for i in range(n)]
    return n, out

def paths(n, out):
    res = []
    nb = [[j for j in range(n) if out[i][j]] for i in range(n)]
    def rec(p, used):
        if len(p) == n:
            res.append(tuple(p)); return
        for v in nb[p[-1]]:
            if not used >> v & 1:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n):
        rec([s], 1 << s)
    return res

def is_lc(a):
    nz = [i for i, x in enumerate(a) if x]
    if not nz:
        return True
    lo, hi = nz[0], nz[-1]
    if any(a[i] == 0 for i in range(lo, hi + 1)):
        return False
    return all(a[i] * a[i] >= a[i - 1] * a[i + 1] for i in range(lo + 1, hi))

for n in range(3, 8):
    lines = subprocess.run(['gentourng', '-q', '-z', '-c', str(n)], capture_output=True, text=True).stdout.split()
    posfail = gapfail = both = 0
    ex_pos = ex_gap = None
    for s in lines:
        n_, out = d6decode(s)
        P = paths(n, out)
        pos = [[0] * n for _ in range(n)]
        for p in P:
            for k, v in enumerate(p):
                pos[v][k] += 1
        pf = [v for v in range(n) if not is_lc(pos[v])]
        gf = []
        for x in range(n):
            for y in range(n):
                if x == y:
                    continue
                g = [0] * (n - 1)
                for p in P:
                    d = p.index(y) - p.index(x)
                    if d > 0:
                        g[d - 1] += 1
                if not is_lc(g):
                    gf.append((x, y, g))
        if pf:
            posfail += 1
            if ex_pos is None:
                ex_pos = (s, pf[0], pos[pf[0]], len(P))
        if gf:
            gapfail += 1
            if ex_gap is None:
                ex_gap = (s, gf[0], len(P))
        if pf and gf:
            both += 1
    print(f"n={n}: strong {len(lines)}; position-distribution failures {posfail}; gap-distribution failures {gapfail}; both {both}")
    if ex_pos:
        print(f"   position example: {ex_pos}")
    if ex_gap:
        print(f"   gap example: {ex_gap}")
