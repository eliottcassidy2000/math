#!/usr/bin/env python3
"""audit G1: small identifications quoted in THM-4701 (G_2 = K_4, G_8 = H(3,4)) with the independent GF code of g1_gq.py."""
import itertools
import networkx as nx
from g1_gq import gq
for q, ref, name in ((2, nx.complete_graph(4), 'K_4'),
                     (8, nx.cartesian_product(nx.cartesian_product(nx.complete_graph(4), nx.complete_graph(4)), nx.complete_graph(4)), 'H(3,4) = K4 x K4 x K4')):
    F, mu, nb = gq(q)
    G = nx.Graph()
    for a in range(F.Q):
        for b in nb[a]:
            G.add_edge(a, b)
    print('G_%d: %d vertices, %d edges ; isomorphic to %s: %s' % (q, G.number_of_nodes(), G.number_of_edges(), name, nx.is_isomorphic(G, ref)))
