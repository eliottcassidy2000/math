#!/usr/bin/env python3
"""The delta = k zero-heavy layer (charge_palettes_boundary_20261020g/zero_heavy.md): pair-free weighted cores Q[w] with
weighted feedback f <= k, all weights <= k, total weight M_k - k (M_k = 2k^2+3k). Deficit identity forces r = k backward
arcs forming a perfect matching on 2k unit endpoints, and either (A) all 2k+1 gaps present with gap deficits summing to k
(each gap weight in [1,k]) or (B) one gap missing, all others of weight k. Records are built as weighted tournaments
(endpoints weight 1, gap vertices weighted), canonicalized under weight-preserving isomorphism, and tested for
pair-freeness (no homogeneous pair) and the exact weighted feedback value k (subset DP over orders of the core).
Prediction: (B) loses exactly the records whose deleted gap separates a matched pair (homogeneous pair), so
Z_{k,k} = (2k-1)!! (C(3k,k) - 1); k = 2 must give 42 (the packet's compiler value). Strong subfamily counted too."""
import sys, itertools
from collections import defaultdict
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
def matchings(points):
    if not points: yield []; return
    a = points[0]
    for i in range(1, len(points)):
        b = points[i]; rest = points[1:i] + points[i+1:]
        for m in matchings(rest): yield [(a, b)] + m
def compositions(total, parts, maxpart):
    if parts == 0:
        if total == 0: yield ()
        return
    for x in range(0, min(total, maxpart) + 1):
        for rest in compositions(total - x, parts - 1, maxpart): yield (x,) + rest
def build(k, matching, gaps):
    """gaps: list of 2k+1 weights (0 = missing). Vertices in order: gap0, e1, gap1, e2, ..., e_{2k}, gap_{2k}.
    Returns (DiGraph with 'w' attribute, strong flag)."""
    order = []; weights = []
    for i in range(2 * k + 1):
        if gaps[i] > 0: order.append(('g', i)); weights.append(gaps[i])
        if i < 2 * k: order.append(('e', i + 1)); weights.append(1)
    idx = {v: j for j, v in enumerate(order)}
    G = nx.DiGraph()
    for j, v in enumerate(order): G.add_node(j, w=weights[j])
    back = {(min(a, b), max(a, b)) for a, b in matching}
    for p in range(len(order)):
        for q in range(p + 1, len(order)):
            u, v = order[p], order[q]
            if u[0] == 'e' and v[0] == 'e' and (u[1], v[1]) in back: G.add_edge(q, p)
            else: G.add_edge(p, q)
    return G
def homogeneous_pair(G):
    nodes = list(G.nodes)
    for a, b in itertools.combinations(nodes, 2):
        ok = True
        for c in nodes:
            if c in (a, b): continue
            if G.has_edge(c, a) != G.has_edge(c, b): ok = False; break
        if ok: return (a, b)
    return None
def weighted_feedback(G):
    nodes = list(G.nodes); n = len(nodes); w = [G.nodes[v]['w'] for v in nodes]
    # min over orders of sum of w_i w_j over backward arcs i->j (j before i): DP over subsets (placed set), cost of placing v last among placed = sum over placed u with arc v->u of w_v w_u
    INF = float('inf'); dp = [INF] * (1 << n); dp[0] = 0
    cost = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if G.has_edge(nodes[i], nodes[j]): cost[i][j] = w[i] * w[j]
    for mask in range(1 << n):
        if dp[mask] == INF: continue
        for v in range(n):
            if mask >> v & 1: continue
            add = sum(cost[v][u] for u in range(n) if mask >> u & 1)  # placing v after all of mask: arcs v->u with u before v are backward
            nm = mask | (1 << v)
            if dp[mask] + add < dp[nm]: dp[nm] = dp[mask] + add
    return dp[(1 << n) - 1]
def is_strong(G): return nx.is_strongly_connected(G)
def canon_key(G):
    return tuple(sorted((G.nodes[v]['w'], G.out_degree(v), sum(G.nodes[u]['w'] for u in G.successors(v))) for v in G.nodes))
def dfact(m): return 1 if m <= 0 else m * dfact(m - 2)
from math import comb
def main():
  for k in [int(a) for a in sys.argv[1:]] or [2, 3]:
      M = 2 * k * k + 3 * k; target = M - k
      endpoints = list(range(1, 2 * k + 1))
      records = []
      for mt in matchings(endpoints):
          # (A) all gaps present, deficits sum to k, each deficit <= k-1 (weight >= 1)
          for defs in compositions(k, 2 * k + 1, k - 1):
              gaps = [k - d for d in defs]; records.append(('A', mt, gaps))
          # (B) one gap missing, others weight k
          for miss in range(2 * k + 1):
              gaps = [k] * (2 * k + 1); gaps[miss] = 0; records.append(('B', mt, gaps))
      print(f"k = {k}: M_k = {M}, weight {target}; raw records: A {sum(1 for r in records if r[0]=='A')}, B {sum(1 for r in records if r[0]=='B')}, total {len(records)}", flush=True)
      buckets = defaultdict(list); n_pairfree = 0; n_bad_feedback = 0; dropped = []
      classes = []  # representatives
      for kind, mt, gaps in records:
          G = build(k, mt, gaps)
          assert sum(G.nodes[v]['w'] for v in G.nodes) == target
          hp = homogeneous_pair(G)
          if hp is not None: dropped.append((kind, mt, gaps, hp)); continue
          n_pairfree += 1
          f = weighted_feedback(G)
          if f != k: n_bad_feedback += 1; continue
          key = canon_key(G)
          found = False
          for H in buckets[key]:
              if DiGraphMatcher(G, H, node_match=lambda a, b: a['w'] == b['w']).is_isomorphic(): found = True; break
          if not found: buckets[key].append(G); classes.append((kind, mt, gaps, is_strong(G)))
      total = len(classes); strong = sum(1 for c in classes if c[3])
      pred = dfact(2 * k - 1) * (comb(3 * k, k) - 1)
      print(f"  dropped (homogeneous pair): {len(dropped)} -> {[(d[0], d[1], d[2]) for d in dropped][:6]}")
      print(f"  pair-free records {n_pairfree}, wrong feedback {n_bad_feedback}, distinct weighted classes {total} (prediction (2k-1)!!(C(3k,k)-1) = {pred}), strong {strong}", flush=True)
      byk = defaultdict(int)
      for c in classes: byk[c[0]] += 1
      print(f"  distinct classes by family: {dict(byk)}; strong by family: { {kk: sum(1 for c in classes if c[0]==kk and c[3]) for kk in byk} }", flush=True)

if __name__ == '__main__':
    main()
