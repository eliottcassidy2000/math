#!/usr/bin/env python3
"""collatz_pentagon_gadget_20260930_selfdouble.py -- a self-contained expansion proof for the hatted icosahedron I_1:
does some iterate C_5^m(I_1) contain two vertex-disjoint, mutually non-adjacent induced copies of I_1?  If so,
C_5^m(I_1) >= I_1 + I_1 (induced), hence |C_5^(mk)(I_1)| >= 13 * 2^k by monotonicity and additivity over disjoint unions,
independently of the edge-removal step in the proof of Theorem 3.5 of Gervacio-Maehara-Ramos.
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, thirteenth note, part A.)
"""
import time
from collatz_pentagon_gadget_20260930 import icosahedron, pentagon_graph, find_icosahedra, hats, is_induced_icosahedron

T0 = time.time()
I = icosahedron()
tad = None
for x in range(12):
    for y in I[x]:
        for z in I[x] & I[y]:
            for t in I[y] - I[x] - I[z] - {x, z}:
                tad = [x, y, z, t]; break
            if tad: break
        if tad: break
    if tad: break
G = [set(a) for a in I] + [set()]
for v in tad:
    G[v].add(12); G[12].add(v)
I1 = G
sizes = [13]
for m in range(1, 5):
    G, _ = pentagon_graph(G)
    sizes.append(len(G))
    icos = find_icosahedra(G, range(len(G)), max_found=400)
    copies = []
    for Ic in icos:
        for u in hats(G, Ic):
            copies.append((frozenset(Ic) | {u}, frozenset(Ic), u))
    # look for two copies that are vertex-disjoint and have no edges between them
    pair = None
    for i in range(len(copies)):
        for j in range(i + 1, len(copies)):
            A, B = copies[i][0], copies[j][0]
            if A & B:
                continue
            if any(b in G[a] for a in A for b in B):
                continue
            pair = (copies[i], copies[j]); break
        if pair:
            break
    print(" C_5^%d(I_1): %d vertices; induced icosahedra found %d; hatted copies of I_1 %d; a vertex-disjoint non-adjacent pair of copies: %s (%.0fs)" % (
        m, len(G), len(icos), len(copies), "YES" if pair else "no", time.time() - T0), flush=True)
    if pair:
        (A, IA, uA), (B, IB, uB) = pair
        # verify: the induced subgraph on A is I_1 (12 + hat with tadpole neighbourhood) and likewise B; no edges between
        okA = is_induced_icosahedron(G, IA) and len(G[uA] & IA) == 4
        okB = is_induced_icosahedron(G, IB) and len(G[uB] & IB) == 4
        print("  copies verified: %s %s; SELF-DOUBLING: C_5^%d(I_1) contains I_1 + I_1 induced, so |C_5^(%dk)(I_1)| >= 13 * 2^k" % (okA, okB, m, m))
        break
    if len(G) > 500:
        break
print(" sizes:", sizes)
