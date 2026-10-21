#!/usr/bin/env python3
"""Strong delta = k zero-heavy layer: nontrivial strong pair-free cores of weight S_k - k (S_k = 2k^2 + k): matching on 2k
unit endpoints crossing every proper initial cut (I_k of them), no exterior gaps, 2k-1 interior gaps; (A') interior gap
deficits summing to k, each in [0, k-1] (weights >= 1); (B') one interior gap missing (weight 0), all others weight k.
Records canonicalized under weight-preserving isomorphism; homogeneous pairs and feedback != k recorded. Also the
all-cores layer at delta = k (control of zero_heavy_delta_k.py) and the exact delta <= k-1 formulas as controls."""
import sys, itertools
from collections import defaultdict
from math import comb
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
sys.path.insert(0, '.')
from zero_heavy_delta_k import matchings, compositions, homogeneous_pair, weighted_feedback, canon_key, dfact
def build_strong(k, matching, gaps):
    """gaps: 2k-1 interior gap weights (0 = missing); vertices: e1, g1, e2, g2, ..., g_{2k-1}, e_{2k}."""
    order = []; weights = []
    for i in range(1, 2 * k + 1):
        order.append(('e', i)); weights.append(1)
        if i < 2 * k and gaps[i - 1] > 0: order.append(('g', i)); weights.append(gaps[i - 1])
    G = nx.DiGraph()
    for j, v in enumerate(order): G.add_node(j, w=weights[j])
    back = {(min(a, b), max(a, b)) for a, b in matching}
    for p in range(len(order)):
        for q in range(p + 1, len(order)):
            u, v = order[p], order[q]
            if u[0] == 'e' and v[0] == 'e' and (u[1], v[1]) in back: G.add_edge(q, p)
            else: G.add_edge(p, q)
    return G
def crossing(matching, k):
    # every proper initial cut {e1..e_i} (1 <= i < 2k) is crossed by some matched pair
    return all(any(a <= i < b for a, b in ((min(x), max(x)) for x in matching)) for i in range(1, 2 * k))
for k in [int(a) for a in sys.argv[1:]] or [2, 3]:
    S = 2 * k * k + k; target = S - k
    endpoints = list(range(1, 2 * k + 1))
    mats = [m for m in matchings(endpoints) if crossing(m, k)]
    I_k = len(mats)
    records = []
    for mt in mats:
        for defs in compositions(k, 2 * k - 1, k - 1):
            records.append(('A', mt, [k - d for d in defs]))
        for miss in range(2 * k - 1):
            gaps = [k] * (2 * k - 1); gaps[miss] = 0; records.append(('B', mt, gaps))
    print(f"k = {k}: I_k = {I_k} crossing matchings, S_k = {S}, weight {target}; raw records A {sum(1 for r in records if r[0]=='A')}, B {sum(1 for r in records if r[0]=='B')}", flush=True)
    buckets = defaultdict(list); classes = []; dropped = []; badf = 0; nonstrong = 0
    for kind, mt, gaps in records:
        G = build_strong(k, mt, gaps)
        assert sum(G.nodes[v]['w'] for v in G.nodes) == target
        if homogeneous_pair(G) is not None: dropped.append((kind, mt, gaps)); continue
        if not nx.is_strongly_connected(G): nonstrong += 1; continue
        if k <= 3 and weighted_feedback(G) != k: badf += 1; continue
        key = canon_key(G); found = False
        for H in buckets[key]:
            if DiGraphMatcher(G, H, node_match=lambda a, b: a['w'] == b['w']).is_isomorphic(): found = True; break
        if not found: buckets[key].append(G); classes.append((kind, mt, gaps))
    adj_drops = len(dropped)
    byk = defaultdict(int)
    for c in classes: byk[c[0]] += 1
    print(f"  dropped (homogeneous pair) {adj_drops}, nonstrong {nonstrong}, wrong feedback {badf}; distinct strong classes {len(classes)} by family {dict(byk)}")
    print(f"  formula check: I_k*(C(3k-2,k)-(2k-1)) = {I_k*(comb(3*k-2,k)-(2*k-1))} for A; B = I_k*(2k-1) - drops = {I_k*(2*k-1)} - {adj_drops} = {I_k*(2*k-1)-adj_drops}; total {I_k*(comb(3*k-2,k)-(2*k-1)) + I_k*(2*k-1) - adj_drops}", flush=True)
