#!/usr/bin/env python3
"""Electrical / Hodge census of the merged tournament metagraph (THM-4739 checks), n = 4..7.
Per class C: scores s, c3 (cyclic triangles), Hodge energies of the +-1 arc flow on K_n: gradient energy E_g = n sum p_i^2 with
p_i = (2 s_i - (n-1))/n, curl energy E_c = C(n,2) - E_g, identities E_g = (n^2-1)/3 - 8 c3/n, E_c = (n-1)(n-2)/6 + 8 c3/n;
the x-axis of THM-790 (x = sum (s_i - mu)^2 scaled: here X = 4 sum(s_i-mu)^2... we use V = sum (s_i - mu)^2 = n(n^2-1)/12 - 2 c3).
Intrinsic currents of a node u: over its tilings t (Hamiltonian paths P mod Aut), J_c3(u) = sum_t [c3(t^c) - c3(t)] where
c3(t^c) - c3(t) = s(start) - s(end) - 1 (THM-785 (4): d_0 - d_(n-1) - 1 with d_0 the score of the path's first vertex = n... here the
base path is n -> ... -> 1 so start = vertex n, end = vertex 1); checked directly by recomputing c3 of the complement tournament.
Metagraph Laplacian identity: sum over lines at u of (c3(v) - c3(u)) = J_c3(u); same for V. Resistance distances in the merged
multigraph (conductance = multiplicity, loops ignored). Exports nodes and directed closed->open multiplicities to JSON."""
import sys, json, subprocess, itertools
from collections import Counter, defaultdict
from fractions import Fraction as F
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
NMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 7
def tiles_of(n): return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
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
def ham_paths_list(n, adj):
    out = []; full = (1 << n) - 1
    def rec(v, mask, path):
        if mask == full: out.append(tuple(path)); return
        for w in adj[v]:
            if not (mask >> w) & 1: path.append(w); rec(w, mask | (1 << w), path); path.pop()
    for v in range(n): rec(v, 1 << v, [v])
    return out
def ham_cycles(n, adj, adjm):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]; dp[1][0] = 1
    for mask in range(1, 1 << n):
        if not mask & 1: continue
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(n) if adjm[v][0])
def c3_of(n, adjm):
    return sum(1 for a, b, c in itertools.combinations(range(n), 3) if (adjm[a][b] and adjm[b][c] and adjm[c][a]) or (adjm[a][c] and adjm[c][b] and adjm[b][a]))
for n in range(NMIN, NMAX + 1):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]; kclose = idx[(n, 1)]; allmask = (1 << m) - 1
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    def is_sym(mask): return all(((mask >> i) & 1) == ((mask >> j) & 1) for i, j in enumerate(mirror))
    can = labelg([to_d6(n, arcs_of(mask)) for mask in range(1 << m)])
    rep = {}
    for mask, c in enumerate(can): rep.setdefault(c, mask)
    cls = {}
    for c, mk in rep.items():
        arcs = arcs_of(mk); adjm = [[0] * n for _ in range(n)]
        for u, v in arcs: adjm[u][v] = 1
        adj = [[v for v in range(n) if adjm[u][v]] for u in range(n)]
        scores = [sum(adjm[u]) for u in range(n)]
        paths = ham_paths_list(n, adj); H = len(paths)
        D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
        A = sum(1 for _ in DiGraphMatcher(D, D).isomorphisms_iter())
        conv = labelg([to_d6(n, [(v, u) for u, v in arcs])])[0]
        c3 = c3_of(n, adjm); mu = (n - 1) / 2
        V = sum((s - mu) ** 2 for s in scores)
        assert abs(V - (n * (n * n - 1) / 12 - 2 * c3)) < 1e-9, 'THM-785 (1)'
        p = [(2 * s - (n - 1)) / n for s in scores]; Eg = n * sum(x * x for x in p); Ec = n * (n - 1) / 2 - Eg
        assert abs(Eg - ((n * n - 1) / 3 - 8 * c3 / n)) < 1e-9 and abs(Ec - ((n - 1) * (n - 2) / 6 + 8 * c3 / n)) < 1e-9, 'Hodge energies'
        # intrinsic currents from Hamiltonian paths (start = first vertex of the path = the base path's vertex n, end = last)
        Jc3 = F(sum(scores[P[0]] - scores[P[-1]] - 1 for P in paths), A)   # per tiling: c3(t^c) - c3(t) = s(start) - s(end) - 1
        # V current: V(t^c) - V(t): scores of the complement: s' = n-1-s except start: n - s(start), end: n-2-s(end)
        def Vof(sc): return sum((s - mu) ** 2 for s in sc)
        JV = F(0)
        for P in paths:
            sc = [n - 1 - scores[v] for v in range(n)]; sc[P[0]] = n - scores[P[0]]; sc[P[-1]] = n - 2 - scores[P[-1]]
            JV += F(Vof(sc) - V).limit_denominator(10 ** 6)
        JV = JV / A
        cls[c] = dict(scores=sorted(scores), H=H, A=A, conv=conv, c3=c3, V=V, Eg=Eg, Ec=Ec, Jc3=Jc3, JV=JV, HC=ham_cycles(n, adj, adjm),
                      ends=Counter((scores[P[0]], scores[P[-1]]) for P in paths))
    node_of = {c: (c if cls[c]['conv'] == c else min(c, cls[c]['conv'])) for c in cls}
    nodes = sorted(set(node_of.values())); nid = {u: i for i, u in enumerate(nodes)}
    N = len(nodes)
    # directed closed->open multiplicities by colour; also check c3(t^c) - c3(t) = s(start)-s(end)-1 directly
    blue = np.zeros((N, N), dtype=int); black = np.zeros((N, N), dtype=int); flux_ok = True
    tcount = Counter(); scount = Counter()
    for t in range(1 << m):
        u = nid[node_of[can[t]]]; tcount[u] += 1
        if is_sym(t): scount[u] += 1
        if not (t >> kclose) & 1: continue   # count each line once, from its closed end
        tc = allmask ^ t; v = nid[node_of[can[tc]]]
        (blue if is_sym(t) else black)[u, v] += 1
        if n <= 6:   # direct flux check on a sample
            adjm = [[0] * n for _ in range(n)]
            for a, b in arcs_of(t): adjm[a][b] = 1
            adjc = [[0] * n for _ in range(n)]
            for a, b in arcs_of(tc): adjc[a][b] = 1
            sc = [sum(r) for r in adjm]
            if c3_of(n, adjc) - c3_of(n, adjm) != sc[n - 1] - sc[0] - 1: flux_ok = False
    print(f"\n=== n = {n}: nodes {N}, flux law c3(t^c)-c3(t) = s(start)-s(end)-1 direct check: {flux_ok}")
    M = blue + black; Lsym = M + M.T   # undirected multiplicities (lines counted once per pair, loops: diagonal counts each loop once... loop contributes 1 to M[u,u])
    # Laplacian identity: sum over lines at u of (c3(v) - c3(u)) [loops contribute 0] = J_c3(u) (node current = sum over its tilings)
    c3v = np.array([cls[u]['c3'] for u in nodes], dtype=float); Vv = np.array([cls[u]['V'] for u in nodes], dtype=float)
    lhs_c3 = np.array([sum(Lsym[i, j] * (c3v[j] - c3v[i]) for j in range(N) if j != i) for i in range(N)])
    lhs_V = np.array([sum(Lsym[i, j] * (Vv[j] - Vv[i]) for j in range(N) if j != i) for i in range(N)])
    rhs_c3 = np.array([float(cls[u]['Jc3']) * (2 if cls[u]['conv'] != u else 1) for u in nodes])
    rhs_V = np.array([float(cls[u]['JV']) * (2 if cls[u]['conv'] != u else 1) for u in nodes])
    print(f"  Laplacian identity L_M c3 = J_c3: max |diff| = {np.max(np.abs(lhs_c3 - rhs_c3)):.3g};  L_M V = J_V: max |diff| = {np.max(np.abs(lhs_V - rhs_V)):.3g}; total current sum J_c3 = {rhs_c3.sum():.3g}")
    # resistance distances (conductance = multiplicity, loops ignored)
    W = Lsym.astype(float); np.fill_diagonal(W, 0); Lap = np.diag(W.sum(1)) - W; Lplus = np.linalg.pinv(Lap)
    R = np.array([[Lplus[i, i] + Lplus[j, j] - 2 * Lplus[i, j] for j in range(N)] for i in range(N)])
    tt = [i for i, u in enumerate(nodes) if cls[u]['H'] == 1][0]  # transitive node
    rows = sorted(range(N), key=lambda i: (cls[nodes[i]]['c3'], -tcount[i]))
    print("  node | H | Aut | SC | t | s | c3 | V | J_c3 | J_V | R(TT,u) | Kirchhoff index")
    Kf = R.sum() / 2
    for i in rows[:min(N, 14)]:
        u = nodes[i]; d = cls[u]
        print(f"  {i:3d} | {d['H']:4d} | {d['A']:2d} | {'SC' if d['conv']==u else 'mg'} | {tcount[i]:4d} | {scount[i]:3d} | {d['c3']:2d} | {d['V']:6.2f} | {float(d['Jc3'])*(2 if d['conv']!=u else 1):7.2f} | {float(d['JV'])*(2 if d['conv']!=u else 1):8.2f} | {R[tt,i]:.4f}")
    print(f"  Kirchhoff index (sum of resistance distances over pairs) = {Kf:.4f}; max resistance from TT = {R[tt].max():.4f}")
    # resistance from TT vs invariants: rank correlation with c3, H, t
    import math
    def spearman(a, b):
        ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b)); return float(np.corrcoef(ra, rb)[0, 1])
    print(f"  Spearman(R(TT,u), c3) = {spearman(R[tt], c3v):.3f}; (R(TT,u), log t_u) = {spearman(R[tt], np.log([tcount[i] for i in range(N)])):.3f}; (R(TT,u), V) = {spearman(R[tt], Vv):.3f}")
    json.dump({'n': n, 'nodes': [dict(id=i, H=cls[u]['H'], Aut=cls[u]['A'], SC=cls[u]['conv'] == u, t=tcount[i], s=scount[i], c3=cls[u]['c3'], V=cls[u]['V'],
                                     HC=cls[u]['HC'], scores=cls[u]['scores'], Jc3=str(cls[u]['Jc3']), JV=str(cls[u]['JV']), rep_mask=rep[u], conv_mask=(rep[cls[u]['conv']] if cls[u]['conv'] in rep else None),
                                     ends={f"{a},{b}": k for (a, b), k in cls[u]['ends'].items()}) for i, u in enumerate(nodes)],
               'blue_closed_to_open': blue.tolist(), 'black_closed_to_open': black.tolist()}, open(f'electric_n{n}.json', 'w'))
