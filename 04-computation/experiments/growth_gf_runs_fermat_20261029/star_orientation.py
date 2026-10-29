"""The (n-1)-run observer is a tournament-induced orientation of the star-transposition Cayley graph of S_n.

Identify an ordered (n-1)-tuple s with the permutation (s_1, ..., s_{n-1}, w), w the missing vertex.  A cover by
n-1 paths has one arc s_j -> w, so M[s,t] = 1 iff t = s with s_j replaced by w and s_j -> w, i.e. t = s.(j n).
Hence M + M^T is the adjacency matrix of Cay(S_n, {(1 n), ..., (n-1 n)}) (each edge {s, s(j n)} carries the pair
{s_j, w} and exactly one of its two directions), the Nauru graph at n = 4.  Checked for every class n = 4, 5, 6;
determinants of the strong classes and the kernels of the non-strong ones are recorded."""
import itertools
import networkx as nx
from runs import tournaments, tour_from_bits, run_matrix, gp_batch, gp_matrix

for n in (4, 5, 6):
    tup = list(itertools.permutations(range(n), n - 1))
    idx = {t: i for i, t in enumerate(tup)}
    G = nx.Graph()
    for t in tup:
        w = (set(range(n)) - set(t)).pop()
        for j in range(n - 1):
            u = list(t); u[j] = w
            G.add_edge(idx[t], idx[tuple(u)])
    assert all(d == n - 1 for _, d in G.degree())
    res = []
    dets = []
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        _, M = run_matrix(A, n - 1, tup)
        N = len(tup)
        for i in range(N):
            for j in range(N):
                a, b = M[i][j], M[j][i]
                assert (a + b == 1) == G.has_edge(i, j) and a + b <= 1 and (a + b == 0 or G.has_edge(i, j))
        dets.append('matdet(%s)' % gp_matrix(M))
    D = gp_batch(dets)
    print('n=%d: star graph on %d vertices, %d edges; every class is an orientation of it. det M_{n-1} by class:' % (n, G.number_of_nodes(), G.number_of_edges()), D)
    if n == 4:
        H = nx.LCF_graph(24, [5, -9, 7, -7, 9, -5], 4)       # Nauru graph LCF notation
        print('   star graph ST_4 isomorphic to the Nauru graph (LCF [5,-9,7,-7,9,-5]^4):', nx.is_isomorphic(G, H))
