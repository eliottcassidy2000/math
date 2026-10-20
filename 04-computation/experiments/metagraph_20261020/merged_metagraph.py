#!/usr/bin/env python3
"""The merged tournament metagraph (explorer complement-line convention, 01-canon/definitions.md tiling model).
Base path P_0: n -> n-1 -> ... -> 1. Tiles (a,b), a >= b+2; bit 0: a -> b, bit 1: b -> a. Tiling t in {0,1}^m, m = C(n-1,2),
tournament T_t. Mirror: the relabelling i -> n+1-i (tile (a,b) -> (n+1-b, n+1-a)); t is SYMMETRIC iff invariant.
Line: t -- t^c (all tiles flipped). Blue line: both ends symmetric; black: both asymmetric (the mirror commutes with
complement, so lines are monochrome). Node: isomorphism class of T_t, MERGED with the converse class when distinct.
Point of symmetry: a self-converse node all of whose lines are blue (all tilings symmetric).
Output per n: global counts, per-node table (tilings, symmetric tilings, blue/black incident, self-lines, H, |Aut|),
line multigraph between nodes with colour multiplicities, degree parities, and formula checks."""
import sys, subprocess, itertools
from collections import Counter, defaultdict
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
NMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 7

def tiles_of(n):
    return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
def to_d6(n, arcs):
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def labelg(lines, chunk=200000):
    out = []
    for s0 in range(0, len(lines), chunk):
        part = lines[s0:s0 + chunk]
        r = subprocess.run(['labelg', '-q'], input='\n'.join(part) + '\n', capture_output=True, text=True).stdout.split()
        assert len(r) == len(part); out.extend(r)
    return out
def ham_paths(n, adj):
    # count directed Hamiltonian paths by bitmask DP
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])

for n in range(NMIN, NMAX + 1):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]
    N = 1 << m; allmask = N - 1
    # vertex labels 1..n -> indices 0..n-1 (vertex i -> i-1); base path arcs i+1 -> i
    base = [(i, i - 1) for i in range(1, n)]
    def arcs_of(mask):
        arcs = list(base)
        for k, (a, b) in enumerate(tiles):
            arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    d6 = [to_d6(n, arcs_of(mask)) for mask in range(N)]
    can = labelg(d6)
    classes = {}; cls_of = [0] * N
    for mask, c in enumerate(can):
        if c not in classes: classes[c] = len(classes)
        cls_of[mask] = classes[c]
    C = len(classes); rep = {}
    for mask in range(N):
        rep.setdefault(cls_of[mask], mask)
    # converse class of each class: reverse all arcs of the representative and canonicalise
    conv_d6 = [to_d6(n, [(v, u) for u, v in arcs_of(rep[ci])]) for ci in range(C)]
    conv_can = labelg(conv_d6); conv = [classes.get(c) for c in conv_can]
    assert all(x is not None for x in conv)
    # merged nodes
    node_of = [None] * C; nodes = []
    for ci in range(C):
        if node_of[ci] is None:
            node_of[ci] = len(nodes); nodes.append([ci] if conv[ci] == ci else [ci, conv[ci]])
            if conv[ci] != ci: node_of[conv[ci]] = node_of[ci]
    M = len(nodes)
    def sym(mask):
        return all(((mask >> i) & 1) == ((mask >> j) & 1) for i, j in enumerate(mirror))
    issym = [sym(mask) for mask in range(N)]
    # H and |Aut| per class
    H = {}; Aut = {}
    for ci in range(C):
        arcs = arcs_of(rep[ci]); adj = [[0] * n for _ in range(n)]
        for u, v in arcs: adj[u][v] = 1
        H[ci] = ham_paths(n, adj)
        D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
        Aut[ci] = sum(1 for _ in DiGraphMatcher(D, D).isomorphisms_iter())
    # per-node statistics and the line multigraph
    tn = Counter(); sn = Counter(); blue_inc = Counter(); black_inc = Counter(); selfb = Counter(); selfk = Counter()
    LB = Counter(); LK = Counter()
    for mask in range(N):
        u = node_of[cls_of[mask]]; tn[u] += 1
        if issym[mask]: sn[u] += 1
        if mask < (allmask ^ mask):
            v = node_of[cls_of[allmask ^ mask]]
            key = (min(u, v), max(u, v))
            if issym[mask]:
                LB[key] += 1; blue_inc[u] += 1; blue_inc[v] += 1
                if u == v: selfb[u] += 1
            else:
                LK[key] += 1; black_inc[u] += 1; black_inc[v] += 1
                if u == v: selfk[u] += 1
    nsym = sum(issym); nSC = sum(1 for nd in nodes if len(nd) == 1); nmerged = M - nSC
    pos = [u for u in range(M) if len(nodes[u]) == 1 and black_inc[u] == 0]
    print(f"\n===== n = {n}: tiles m = {m}, tilings {N}, symmetric tilings {nsym} (= 2^{(m + (n-1)//2)//2}? {nsym == 2 ** ((m + (n - 1) // 2) // 2)}), "
          f"classes {C}, self-converse classes {nSC}, merged nodes {M} (= (C + SC)/2 {M == (C + nSC)//2}), lines: blue {sum(LB.values())}, black {sum(LK.values())} =====")
    print(f"points of symmetry (self-converse nodes with only blue lines): {len(pos)}; self-lines: blue {sum(selfb.values())}, black {sum(selfk.values())}")
    # formula checks
    ok_tiles = all(tn[u] == sum(H[ci] // Aut[ci] for ci in nodes[u]) for u in range(M))
    print(f"tilings per node = sum H/|Aut| over its classes: {ok_tiles}; every node count odd (SC) / 2 x odd (merged): "
          f"{all((tn[u] % 2 == 1) if len(nodes[u]) == 1 else (tn[u] % 4 == 2) for u in range(M))}; merged nodes have no symmetric tiling: {all(sn[u] == 0 for u in range(M) if len(nodes[u]) == 2)}; "
          f"SC nodes have an odd number of symmetric tilings: {all(sn[u] % 2 == 1 for u in range(M) if len(nodes[u]) == 1)}")
    # degree parity by colour
    print("blue-degree parity histogram over nodes:", sorted(Counter(blue_inc[u] % 2 for u in range(M)).items()),
          "| black-degree parity histogram:", sorted(Counter(black_inc[u] % 2 for u in range(M)).items()))
    print("black degrees (incl. self-lines twice) all even?", all(black_inc[u] % 2 == 0 for u in range(M)),
          "| blue degree odd exactly on SC nodes?", all((blue_inc[u] % 2 == 1) == (len(nodes[u]) == 1) for u in range(M)))
    # per node table (sorted by tilings)
    if n <= 6 or True:
        print("node | classes (H,|Aut|) | tilings | sym | blue inc | black inc | self blue | self black | type")
        for u in sorted(range(M), key=lambda u: (-tn[u], u)):
            typ = 'POS' if u in pos else ('SC' if len(nodes[u]) == 1 else 'merged')
            print(f"  {u:3d} | {[(H[ci], Aut[ci]) for ci in nodes[u]]} | {tn[u]:4d} | {sn[u]:3d} | {blue_inc[u]:4d} | {black_inc[u]:4d} | {selfb[u]:3d} | {selfk[u]:3d} | {typ}")
    # line multigraph summary
    print("line multigraph: node pairs with blue multiplicity:", sorted(Counter(LB.values()).items()), "; with black multiplicity:", sorted(Counter(LK.values()).items()))
    print("node pairs joined by both colours:", sum(1 for k in set(LB) & set(LK)))
