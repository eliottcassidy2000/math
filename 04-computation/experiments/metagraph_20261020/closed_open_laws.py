#!/usr/bin/env python3
"""Closed/open tilings and the forced arrangement of the merged metagraph (mac-mini-2026-10-20-metagraph).
Tiling model of 01-canon/definitions.md: base path n -> n-1 -> ... -> 1; tiles (a,b), a >= b+2; bit 0: a -> b, bit 1: b -> a.
A tiling is CLOSED when the tile (n,1) carries the arc 1 -> n (bit 1), i.e. the base path closes to a Hamiltonian cycle;
OPEN otherwise. The complement line t -- t^c flips every tile, so every line has exactly one closed end.
Per node checks: closed tilings = n*HC(C)/|Aut C| (x2 for a converse pair), non-strong nodes are adjacent to strong nodes
only and carry no self-line, SC nodes have an odd number of non-loop blue lines, connectivity of the metagraph and of the
blue graph, fully-closed nodes (no open tiling), and the n = 4, 5 node tables as worked examples."""
import sys, subprocess
from collections import defaultdict, Counter
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
def ham_paths(n, adj):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full])
def ham_cycles(n, adj, adjm):
    # Hamiltonian cycles = Hamiltonian paths starting at 0 whose end beats 0
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    dp[1][0] = 1
    for mask in range(1, 1 << n):
        if not mask & 1: continue
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(n) if adjm[v][0])
def aut_order(n, arcs):
    D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
    return sum(1 for _ in DiGraphMatcher(D, D).isomorphisms_iter())
for n in range(NMIN, NMAX + 1):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]
    f = sum(1 for i, j in enumerate(mirror) if i == j)
    allmask = (1 << m) - 1; kclose = idx[(n, 1)]
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    def is_sym(mask):
        return all(((mask >> i) & 1) == ((mask >> j) & 1) for i, j in enumerate(mirror))
    print(f"\n=================== n = {n}: m = {m}, f = {f}, tilings {1 << m} ===================", flush=True)
    can = labelg([to_d6(n, arcs_of(mask)) for mask in range(1 << m)])
    classes = {}
    for mask, c in enumerate(can):
        if c not in classes: classes[c] = mask
    # per class data
    cls = {}
    for c, rep in classes.items():
        arcs = arcs_of(rep); adjm = [[0] * n for _ in range(n)]
        for u, v in arcs: adjm[u][v] = 1
        adj = [[v for v in range(n) if adjm[u][v]] for u in range(n)]
        H = ham_paths(n, adj); HC = ham_cycles(n, adj, adjm); A = aut_order(n, arcs)
        conv = labelg([to_d6(n, [(v, u) for u, v in arcs])])[0]
        cls[c] = dict(H=H, HC=HC, A=A, conv=conv, strong=HC > 0)
    node_of = {}
    for c in classes:
        node_of[c] = c if cls[c]['conv'] == c else min(c, cls[c]['conv'])
    nodes = sorted(set(node_of.values()))
    nid = {u: i for i, u in enumerate(nodes)}
    T = Counter(); S = Counter(); CL = Counter(); SCL = Counter()
    sym = [False] * (1 << m)
    for mask in range(1 << m):
        u = nid[node_of[can[mask]]]; T[u] += 1
        s = is_sym(mask); sym[mask] = s
        closed = (mask >> kclose) & 1
        if s: S[u] += 1
        if closed:
            CL[u] += 1
            if s: SCL[u] += 1
    # lines
    blue_mult = Counter(); black_mult = Counter(); blue_loops = Counter(); black_loops = Counter()
    bad_closed = 0
    for t in range(1 << m):
        tc = allmask ^ t
        if t > tc: continue
        if ((t >> kclose) & 1) == ((tc >> kclose) & 1): bad_closed += 1
        u, v = nid[node_of[can[t]]], nid[node_of[can[tc]]]
        key = (min(u, v), max(u, v))
        if sym[t]:
            blue_mult[key] += 1
            if u == v: blue_loops[u] += 1
        else:
            black_mult[key] += 1
            if u == v: black_loops[u] += 1
    print(f"lines with two closed or two open ends: {bad_closed} (must be 0)")
    # checks
    viol_closed = []; viol_nonstrong = []; viol_blueodd = []
    strong_nodes = set(); sc_nodes = set(); merged_nodes = set()
    for u in nodes:
        i = nid[u]; d = cls[u]; merged = d['conv'] != u
        if merged: merged_nodes.add(i)
        else: sc_nodes.add(i)
        if d['strong']: strong_nodes.add(i)
        pred = n * d['HC'] * (2 if merged else 1) // d['A'] if (n * d['HC'] * (2 if merged else 1)) % d['A'] == 0 else None
        if pred != CL[i]: viol_closed.append((u, pred, CL[i]))
    for (u, v), k in list(blue_mult.items()) + list(black_mult.items()):
        if u not in strong_nodes and v not in strong_nodes: viol_nonstrong.append((u, v, k))
    for i in sc_nodes:
        nonloop_blue = S[i] - 2 * blue_loops[i]
        if nonloop_blue % 2 == 0: viol_blueodd.append((i, S[i], blue_loops[i]))
    print(f"nodes {len(nodes)} (SC {len(sc_nodes)}, merged {len(merged_nodes)}); strong nodes {len(strong_nodes)}, non-strong {len(nodes) - len(strong_nodes)}")
    print(f"closed-tiling law c_u = n*HC/|Aut| (x2 merged): violations {len(viol_closed)} {viol_closed[:3]}")
    print(f"lines between two non-strong nodes or loops at non-strong nodes: {len(viol_nonstrong)} (must be 0)")
    print(f"SC nodes with an even number of non-loop blue lines: {len(viol_blueodd)} (must be 0)")
    tot_closed = sum(CL.values()); tot_scl = sum(SCL.values())
    print(f"closed tilings {tot_closed} = 2^(m-1) = {1 << (m - 1)}: {tot_closed == 1 << (m - 1)}; closed symmetric {tot_scl} = 2^((m+f)/2-1) = {1 << ((m + f) // 2 - 1)}: {tot_scl == 1 << ((m + f) // 2 - 1)}")
    fully_closed = [(u, cls[u]['H'], cls[u]['A'], T[nid[u]], S[nid[u]]) for u in nodes if CL[nid[u]] == T[nid[u]]]
    print(f"fully closed nodes (every tiling closed; all Hamiltonian paths close): {len(fully_closed)} -> (H,|Aut|,t,s): {[(h, a, t, s) for _, h, a, t, s in fully_closed]}")
    print(f"fully open nodes (= non-strong): {len(nodes) - len(strong_nodes)}")
    # line type counts
    ss = sum(k for (u, v), k in list(blue_mult.items()) + list(black_mult.items()) if u in strong_nodes and v in strong_nodes)
    sn = sum(k for (u, v), k in list(blue_mult.items()) + list(black_mult.items()) if (u in strong_nodes) != (v in strong_nodes))
    print(f"lines strong-strong {ss}, strong/non-strong {sn}, non-strong/non-strong {sum(blue_mult.values()) + sum(black_mult.values()) - ss - sn}")
    # connectivity
    G = nx.Graph(); G.add_nodes_from(range(len(nodes)))
    for (u, v), k in list(blue_mult.items()) + list(black_mult.items()):
        if u != v: G.add_edge(u, v)
    B = nx.Graph(); B.add_nodes_from(sc_nodes)
    for (u, v), k in blue_mult.items():
        if u != v: B.add_edge(u, v)
    print(f"metagraph connected: {nx.is_connected(G)} (components {nx.number_connected_components(G)}); blue graph on SC nodes connected: {nx.is_connected(B)} (components {nx.number_connected_components(B)})")
    # degree parity in the folded graph
    print(f"blue loops total {sum(blue_loops.values())}, black loops total {sum(black_loops.values())}; nodes with self-lines: {sum(1 for i in range(len(nodes)) if blue_loops[i] + black_loops[i])}, all strong: {all(i in strong_nodes for i in range(len(nodes)) if blue_loops[i] + black_loops[i])}")
    # pure blue nodes
    pb = [(cls[u]['H'], cls[u]['A'], T[nid[u]]) for u in nodes if S[nid[u]] == T[nid[u]]]
    print(f"pure-blue nodes {len(pb)}: {pb}")
    # distribution of nodes by (strong, SC)
    print("node census by (strong, SC):", Counter((nid[u] in strong_nodes, nid[u] in sc_nodes) for u in nodes))
    if n <= 5:
        print("node | type | strong | H | |Aut| | t | s | closed | closed-sym | blue loops | black loops")
        for u in nodes:
            i = nid[u]; d = cls[u]
            print(f"  {i:2d} | {'SC' if i in sc_nodes else 'merged'} | {int(d['strong'])} | {d['H']} | {d['A']} | {T[i]} | {S[i]} | {CL[i]} | {SCL[i]} | {blue_loops[i]} | {black_loops[i]}")
        print("blue lines (node pairs, multiplicity):", dict(blue_mult)); print("black lines:", dict(black_mult))
