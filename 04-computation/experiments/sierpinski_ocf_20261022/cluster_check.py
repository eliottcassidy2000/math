# Cluster (linear-forest) expansion of kappa(T) = H(T) 2^(n-1)/n!, checked by brute force at low order.
# Exact identity: 2^(n-1) H(T) = sum over directed linear forests F of K_n of (n - |F|)! prod_(e in F) S_e, so
#     kappa(T) = sum_L N_L(S)/(n)_L,    N_L = sum over linear forests with L arcs of prod S_e.
# Write N_L/n^L as a sum over ordered tuples of runs (l_1..l_p)/p! and partitions pi of their positions with
# Moebius weight mu(pi) and the free homomorphism sum Z_pi of Shat = S/n over the quotient graph.  log of the
# generating function = sum over CONNECTED (tuple, pi).  Claim (excess 0 = cycles): the excess-0 part of
# [z^L] log Phi equals [z^L] (1/2) sum_k log((1 + z^2 x_k^2)/(1 - z^2 x_k^2)), x_k = sigma_k/n, i.e. log rho.
# This script enumerates every connected (tuple, pi) with L <= LMAX, sorts by excess, and compares.
import sys, itertools, math
import numpy as np
from math import factorial

def compositions(L):
    if L == 0: yield (); return
    for first in range(1, L + 1):
        for rest in compositions(L - first): yield (first,) + rest

def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p

def mobius(blocks):
    m = 1
    for b in blocks: m *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
    return m

def tournament(kind, n):
    A = np.zeros((n, n))
    if kind == 'paley':
        sq = {(x * x) % n for x in range(1, n)}
        for u in range(n):
            for v in range(n):
                if u != v and (v - u) % n in sq: A[u, v] = 1
    elif kind == 'rot':
        for u in range(n):
            for j in range(1, (n - 1) // 2 + 1): A[u, (u + j) % n] = 1
    return A

def run(kind, n, LMAX):
    A = tournament(kind, n); S = A - A.T; Sh = S / n
    x = np.linalg.eigvals(S).imag / n
    letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    out = {}
    for L in range(1, LMAX + 1):
        by_exc = {}
        for comp in compositions(L):
            p = len(comp)
            pos = [(j, i) for j in range(p) for i in range(comp[j] + 1)]
            ends = {(j, 0) for j in range(p)} | {(j, comp[j]) for j in range(p)}
            for blocks in set_partitions(pos):
                bid = {}
                for t, b in enumerate(blocks):
                    for q in b: bid[q] = t
                # singleton endpoint => Z = 0 (S 1 = 0); loop edge => Z = 0
                if any(len(blocks[bid[e]]) == 1 for e in ends): continue
                edges = [(bid[(j, i)], bid[(j, i + 1)]) for j in range(p) for i in range(comp[j])]
                if any(a == b for a, b in edges): continue
                nb = len(blocks)
                # connectivity of quotient graph
                par = list(range(nb))
                def f(a):
                    while par[a] != a: par[a] = par[par[a]]; a = par[a]
                    return a
                for a, b in edges: par[f(a)] = f(b)
                if len({f(a) for a in range(nb)}) != 1: continue
                exc = len(edges) - nb
                spec = ','.join(letters[a] + letters[b] for a, b in edges) + '->'
                Z = np.einsum(spec, *([Sh] * len(edges)), optimize='greedy')
                by_exc[exc] = by_exc.get(exc, 0.0) + mobius(blocks) * Z / factorial(p)
        pred0 = 0.0
        if L % 2 == 0:
            m = L // 2
            if m % 2 == 1: pred0 = np.sum(x ** L) / m
        out[L] = (by_exc, pred0)
        print(f"{kind} n={n} L={L}: excess-0 sum = {by_exc.get(0, 0):+.12f}   predicted (log rho coeff) = {pred0:+.12f}   "
              + '  '.join(f"exc{e}: {v:+.3e}" for e, v in sorted(by_exc.items()) if e > 0), flush=True)
    return out

if __name__ == '__main__':
    kind = sys.argv[1]; n = int(sys.argv[2]); LMAX = int(sys.argv[3])
    run(kind, n, LMAX)
