# Cone reduction of a seven-vertex Paley sum.  Choose an independent set C of degree-3 vertices; every c in C has its
# three neighbours in the base B = V - C, and summing x_c over ALL of F_q gives
#   Lambda_c = (-1)^{#neighbours u > c} chi(x_b - x_a) L((x_d - x_a)/(x_b - x_a)),  N(c) = {a < b < d}.
# The cone-reduced sum Y_H = sum over injective x_B (x at two base vertices fixed to 0, 1) of the base arcs times
# prod_c Lambda_c costs O(q^(|B|-2)).  It differs from S_H by merge terms (a cone value meeting another vertex), which
# are sums on fewer free variables; we test S_H - Y_H against the base span numerically.  mac-mini-2026-10-25-pyth.
import itertools
import numpy as np

def tables(q):
    chi = np.full(q, -1, dtype=np.int64); chi[0] = 0; chi[[(x * x) % q for x in range(1, q)]] = 1
    w = np.arange(q); base = chi[(w * (w - 1)) % q]
    L = np.array([int((base * chi[(w - l) % q]).sum()) for l in range(q)], dtype=np.int64)
    inv = np.zeros(q, dtype=np.int64); inv[1:] = [pow(int(x), q - 2, q) for x in range(1, q)]
    return chi, L, inv

def best_cone_set(edges):
    E = {tuple(sorted(e)) for e in edges}
    adj = {v: {u for u in range(7) if tuple(sorted((u, v))) in E} for v in range(7)}
    d3 = [v for v in range(7) if len(adj[v]) == 3]
    best = ()
    for r in range(len(d3), 0, -1):
        for C in itertools.combinations(d3, r):
            if all(b not in adj[a] for a, b in itertools.combinations(C, 2)):
                return C, adj
    return (), adj

def Y_H(q, edges, C=None, T=None):
    chi, L, inv = tables(q) if T is None else T
    E = {tuple(sorted(e)) for e in edges}
    if C is None: C, adj = best_cone_set(edges)
    else: adj = {v: {u for u in range(7) if tuple(sorted((u, v))) in E} for v in range(7)}
    B = [v for v in range(7) if v not in C]
    a0, b0 = B[0], B[1]                       # fixed to 0 and 1
    free = B[2:]
    m = len(free)
    sign = 1
    for c in C: sign *= (-1) ** sum(1 for u in adj[c] if u > c)
    base_arcs = [(i, j) for i, j in itertools.combinations(B, 2) if (i, j) in E]
    cones = [tuple(sorted(adj[c])) for c in C]
    tot = 0
    # enumerate all but the last free variable explicitly, vectorise the last one
    last = np.arange(q)
    for vals in itertools.product(range(q), repeat=m - 1):
        x = {a0: 0, b0: 1}
        for v, t in zip(free[:-1], vals): x[v] = t
        fixedvals = [x[v] for v in x]
        if len(set(fixedvals)) < len(fixedvals): continue
        X = {**x, free[-1]: last}
        mask = np.ones(q, dtype=np.int64)
        for t in fixedvals: mask[t] = 0
        prod = mask.copy()
        for i, j in base_arcs:
            prod = prod * chi[(np.asarray(X[j]) - np.asarray(X[i])) % q]
        for a, b, d in cones:
            xa, xb, xd = np.asarray(X[a]), np.asarray(X[b]), np.asarray(X[d])
            prod = prod * chi[(xb - xa) % q] * L[((xd - xa) % q * inv[(xb - xa) % q]) % q]
        tot += int(prod.sum())
    return sign * tot, C

if __name__ == '__main__':
    import json, re, sys
    G = json.load(open('genuine37.json'))
    NEW = [2, 4, 5, 6, 9, 10, 16, 17, 18, 20, 21, 24, 28, 30, 31, 36]
    for i in NEW:
        C, adj = best_cone_set(G[i]['edges'])
        print(f"H{i:02d}: cone set {C} (base size {7 - len(C)})")
