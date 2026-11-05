#!/usr/bin/env python3
"""extension_spectra.py -- one-vertex extension spectra. For a tournament T with Hamiltonian paths P_1..P_h and a cut
(A, B) of V(T) (z -> A, B -> z), H(T + z) = sum_P ([P_1 in A] + [P_last in B] + #{i : P_i in B, P_(i+1) in A}).
For every tournament of order n <= NMAX: is the set of values over all 2^n cuts an interval of odd numbers [h, M]?"""
import sys, subprocess, itertools
import numpy as np

def tournaments(n):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for line in out:
        adj = [0] * n
        for (i, j), b in zip(pairs, line):
            if b == '1': adj[i] |= 1 << j
            else: adj[j] |= 1 << i
        yield adj

def ham_paths_list(n, adj):
    paths = []
    def rec(v, S, p):
        if S == (1 << n) - 1: paths.append(tuple(p)); return
        nxt = adj[v] & ~S
        while nxt:
            u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
            rec(u, S | (1 << u), p + [u])
    for v in range(n): rec(v, 1 << v, [v])
    return paths

def strong(n, adj):
    reach = 1; fr = [0]
    while fr:
        u = fr.pop(); nx = adj[u] & ~reach
        while nx:
            v = (nx & -nx).bit_length() - 1; nx &= nx - 1; reach |= 1 << v; fr.append(v)
    if reach != (1 << n) - 1: return False
    radj = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    reach = 1; fr = [0]
    while fr:
        u = fr.pop(); nx = radj[u] & ~reach
        while nx:
            v = (nx & -nx).bit_length() - 1; nx &= nx - 1; reach |= 1 << v; fr.append(v)
    return reach == (1 << n) - 1

def ext_values(n, paths):
    P = np.array(paths, dtype=np.int64)                       # h x n
    h = len(paths)
    cuts = np.arange(1 << n, dtype=np.int64)
    # membership of P[:, i] in A for every cut: bit P[:, i] of cut
    vals = np.zeros(1 << n, dtype=np.int64)
    inA = [((cuts[:, None] >> P[:, i][None, :]) & 1) for i in range(n)]   # list of (cuts x h) arrays
    first = inA[0].sum(axis=1); last = (1 - inA[n - 1]).sum(axis=1)
    trans = sum(((1 - inA[i]) * inA[i + 1]).sum(axis=1) for i in range(n - 1))
    return first + last + trans

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(3, NMAX + 1):
        tot = 0; interval = 0; strong_tot = 0; strong_int = 0; worst = None; gaps_examples = []
        for adj in tournaments(n):
            paths = ham_paths_list(n, adj); h = len(paths)
            vals = ext_values(n, paths); S = set(int(v) for v in vals)
            M = max(S); isint = all(v in S for v in range(h, M + 1, 2))
            st = strong(n, adj)
            tot += 1; interval += isint
            if st: strong_tot += 1; strong_int += isint
            if not isint and len(gaps_examples) < 3:
                gaps_examples.append((h, M, sorted(v for v in range(h, M + 1, 2) if v not in S)[:6], st))
            # ratio M/h
            r = M / h
            if worst is None or r < worst[0]: worst = (r, h, M, st)
        print("n=%d: tournaments %d, extension spectrum is an odd interval [h, M] for %d; strong: %d of %d; min M/h = %.3f (h=%d, M=%d, strong=%s); gap examples (h, M, missing, strong): %s" % (n, tot, interval, strong_int, strong_tot, worst[0], worst[1], worst[2], worst[3], gaps_examples), flush=True)
