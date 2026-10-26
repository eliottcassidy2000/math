"""guard_rank_restored_pairs_20261026.py -- opus-2026-10-26-S62

Named Hamilton-path guards on the marked-path carrier X_n (tournaments containing 0 -> 1 -> ... -> n-1; 2^M states,
M = (n-1)(n-2)/2).  For a permutation p of the vertices, the guard g_p(x) = 1 iff p is a directed Hamiltonian path of x.
Quantities (the pasted summary claims rank 44 and 64 at n = 5):
  r1(n) = rank over Q of the span of the single guards {g_p}            ("single Hamilton-path guards")
  r2(n) = rank over Q of the span of the products {g_p g_q}             ("pairs of restored tests": test p, restore, test q)
  and for comparison the rank of the Boolean guard algebra generated (all products), the number of distinct guard
  vectors, and the support sizes.  Ranks are computed exactly over Q with fractions for n <= 5 and modulo two large
  primes for n = 6 (the modular rank is a lower bound for the rational rank; equality of the two modular ranks is
  reported as the evidence).
Usage: python3 guard_rank_restored_pairs_20261026.py [maxn=6]
"""
import sys, itertools
from fractions import Fraction

def states(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    M = len(pairs); out = []
    for code in range(1 << M):
        adj = [[0] * n for _ in range(n)]
        for i in range(n - 1): adj[i][i + 1] = 1
        for k, (i, j) in enumerate(pairs):
            if code >> k & 1: adj[i][j] = 1
            else: adj[j][i] = 1
        out.append(adj)
    return out

def guard_vectors(n, X):
    vecs = {}
    for p in itertools.permutations(range(n)):
        v = tuple(1 if all(x[p[i]][p[i + 1]] for i in range(n - 1)) else 0 for x in X)
        vecs[p] = v
    return vecs

def rank_mod(rows, mod):
    rows = [list(r) for r in rows]; rank = 0; ncol = len(rows[0]) if rows else 0; r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c] % mod: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], mod - 2, mod)
        rows[r] = [(a * inv) % mod for a in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] % mod:
                f = rows[i][c]; rows[i] = [(a - f * b) % mod for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r

def rank_Q(rows):
    rows = [[Fraction(a) for a in r] for r in rows]; r = 0; ncol = len(rows[0])
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = 1 / rows[r][c]; rows[r] = [a * inv for a in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r

def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(4, maxn + 1):
        X = states(n); N = len(X); G = guard_vectors(n, X)
        singles = list(set(G.values()))
        # products of pairs (distinct permutations), as a set of vectors
        pairs = set()
        for a, b in itertools.combinations(singles, 2):
            v = tuple(x * y for x, y in zip(a, b))
            if any(v): pairs.add(v)
        allvecs = set(singles) | pairs
        if n <= 5:
            r1 = rank_Q(singles); r2 = rank_Q(list(pairs)); r12 = rank_Q(list(allvecs)); tag = "rational"
        else:
            P1, P2 = 1000000007, 998244353
            r1 = (rank_mod(singles, P1), rank_mod(singles, P2)); r2 = (rank_mod(list(pairs), P1), rank_mod(list(pairs), P2)); r12 = (rank_mod(list(allvecs), P1), rank_mod(list(allvecs), P2)); tag = "mod two primes"
        supports = sorted(set(sum(v) for v in singles))
        print(f"n = {n}: states {N}, distinct single guards {len(singles)} (of {len(G)} permutations), pair products {len(pairs)}; "
              f"rank singles {r1}, rank pair products {r2}, rank singles+pairs {r12} ({tag}); guard support sizes {supports[:12]}{'...' if len(supports) > 12 else ''}")

if __name__ == "__main__":
    main()
