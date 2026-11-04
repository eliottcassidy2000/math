#!/usr/bin/env python3
"""Independent verification of the (B) certificates in certs.txt.

For every line 'H n row_0 ... row_{n-1}':
  * the rows define a tournament (exactly one arc per pair, no loops);
  * the number of Hamiltonian paths, counted by plain DFS enumeration (no DP), equals H;
  * networkx's VF2 matcher finds exactly one automorphism (the identity).
Also: H is multiplicative over strong components and Aut is trivial iff each strong component is
asymmetric -- spot-checked on random non-strong tournaments built from the certificates.
"""
import sys, random
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def parse(line):
    p = line.split()
    H, n = int(p[0]), int(p[1])
    rows = [int(x, 16) for x in p[2:2 + n]]
    out = [{w for w in range(n) if rows[v] >> w & 1} for v in range(n)]
    return H, n, out

def is_tournament(n, out):
    for v in range(n):
        if v in out[v]:
            return False
        for w in range(v + 1, n):
            if (w in out[v]) == (v in out[w]):
                return False
    return True

def count_hpaths_dfs(n, out):
    total = 0
    sys.setrecursionlimit(10000)
    def dfs(v, visited, k):
        nonlocal total
        if k == n:
            total += 1
            return
        for w in out[v]:
            if not visited >> w & 1:
                dfs(w, visited | 1 << w, k + 1)
    for s in range(n):
        dfs(s, 1 << s, 1)
    return total

def n_automorphisms(n, out, stop=None):
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    G.add_edges_from((v, w) for v in range(n) for w in out[v])
    c = 0
    for _ in DiGraphMatcher(G, G).isomorphisms_iter():
        c += 1
        if stop and c >= stop:
            break
    return c

certs = [parse(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'certs.txt') if l.strip()]
bad = 0
seen = set()
for H, n, out in certs:
    ok = is_tournament(n, out) and count_hpaths_dfs(n, out) == H and n_automorphisms(n, out, stop=2) == 1
    if not ok:
        bad += 1
        print("FAILED certificate", H, n)
    seen.add(H)
odd = set(range(1, 1126, 2))
print("certificates:", len(certs), " failures:", bad)
print("odd H <= 1125 not certified:", sorted(odd - seen))
print("max order used:", max(n for _, n, _ in certs), " order histogram:",
      {k: sum(1 for _, n, _ in certs if n == k) for k in sorted({n for _, n, _ in certs})})

# --- spot-check multiplicativity of H and of Aut over strong components
def compose(parts):
    """transitive arrangement of tournaments: all arcs from earlier parts to later parts"""
    out, base = [], []
    off = 0
    for (n, o) in parts:
        base.append(off); off += n
    N = off
    out = [set() for _ in range(N)]
    for idx, (n, o) in enumerate(parts):
        b = base[idx]
        for v in range(n):
            out[b + v] |= {b + w for w in o[v]}
            out[b + v] |= set(range(b + n, N))
    return N, out

random.seed(7)
small = [(H, n, out) for (H, n, out) in certs if n <= 6]
c3 = (3, [{1}, {2}, {0}])
mult_ok, aut_ok, trials = 0, 0, 0
for _ in range(60):
    k = random.randint(2, 3)
    picks = [random.choice(small) for _ in range(k)]
    parts = [(n, out) for (H, n, out) in picks]
    use_c3 = random.random() < 0.3
    if use_c3:
        parts.insert(random.randint(0, len(parts)), c3)
    N, out = compose(parts)
    if N > 13:
        continue
    trials += 1
    expectH = 1
    for (H, n, o) in picks:
        expectH *= H
    if use_c3:
        expectH *= 3
    # H multiplicative?  (note: certificates need not be strong; H is multiplicative anyway because the
    # composition is itself a transitive arrangement of the certificates' strong components)
    mult_ok += (count_hpaths_dfs(N, out) == expectH)
    a = n_automorphisms(N, out)
    expect_aut = 3 if use_c3 else 1
    aut_ok += (a == expect_aut)
print("composition spot-checks: %d trials, H multiplicative in %d, |Aut| as predicted in %d" % (trials, mult_ok, aut_ok))
