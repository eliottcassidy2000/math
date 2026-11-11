#!/usr/bin/env python3
"""Small facts behind the pasted "atoms" text, recomputed: strong tournament classes of orders 3-7, how many are
self-converse, and the edge-usage profile (for every unordered pair, the number of Hamiltonian paths that use it
consecutively) of the classes of orders 3-5.  The profile is the weighted graph W of Session265.
Usage: atoms_small.py"""
import subprocess, itertools
from collections import Counter
def classes(m):
    out = subprocess.run(['gentourng', '-q', str(m)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = [0] * m; k = 0
        for i in range(m):
            for j in range(i + 1, m):
                if line[k] == '1': A[i] |= 1 << j
                else: A[j] |= 1 << i
                k += 1
        yield A
def strong(out, n):
    def reach(adj):
        seen = 1; stack = [0]
        while stack:
            v = stack.pop(); t = adj[v] & ~seen
            while t:
                w = (t & -t).bit_length() - 1; seen |= 1 << w; stack.append(w); t &= t - 1
        return seen
    inn = [0] * n
    for a in range(n):
        for b in range(n):
            if (out[a] >> b) & 1: inn[b] |= 1 << a
    full = (1 << n) - 1
    return reach(out) == full and reach(inn) == full
def self_converse(out, n):
    arc = [[(out[a] >> b) & 1 for b in range(n)] for a in range(n)]
    sc = [bin(o).count('1') for o in out]
    perm = [-1] * n; used = [False] * n
    def rec(i):
        if i == n: return True
        for c in range(n):
            if used[c] or sc[c] != n - 1 - sc[i]: continue
            if all(arc[i][j] == arc[perm[j]][c] for j in range(i)):
                perm[i] = c; used[c] = True
                if rec(i + 1): return True
                perm[i] = -1; used[c] = False
        return False
    return rec(0)
def profile(out, n):
    use = Counter(); H = 0
    for p in itertools.permutations(range(n)):
        if all((out[p[i]] >> p[i + 1]) & 1 for i in range(n - 1)):
            H += 1
            for i in range(n - 1): use[frozenset((p[i], p[i + 1]))] += 1
    return H, sorted(use.values(), reverse=True)
print("order: strong classes, self-converse among them, converse pairs")
for n in range(3, 8):
    st = [A for A in classes(n) if strong(A, n)]; sc = sum(1 for A in st if self_converse(A, n))
    print("  n = %d: %d strong classes, %d self-converse, %d converse pairs" % (n, len(st), sc, (len(st) - sc) // 2))
print("edge-usage profiles (non-zero entries, sorted), orders 3-5:")
for n in range(3, 6):
    for A in classes(n):
        H, prof = profile(A, n)
        print("  n = %d, H = %2d, %s: %s" % (n, H, "strong" if strong(A, n) else "not strong", prof))
