"""Audit S50 - a2: THM-4744 (2) for every composition with parts >= 2 and m = N + 1 <= 12 (orders 3..12):
  * out-degree of x_t is N-1 (t = 0), N-t (1 <= t <= N-1), 1 (t = N); sorted scores (1,1,2,...,m-3,m-2,m-2);
  * cyclic triples = m - 2 (direct count over triples, not the score formula);
  * strong; |Aut| = 1 for m >= 4 (canon leaf count; networkx VF2 automorphism count for m <= 10);
  * the spine is recovered from the tournament by the stated rule (scores, then the two tie arcs) and the back arcs
    (arcs pointing backward in the recovered order) spell out the composition;
  * all chains of one order are pairwise non-isomorphic (distinct canonical codes; VF2 cross-check m <= 9);
  * converse(IC(p)) is IC(reversed p) under t -> N - t exactly; self-converse iff palindrome."""
import time, itertools
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def to_nx(adj):
    G = nx.DiGraph(); n = len(adj)
    G.add_nodes_from(range(n))
    G.add_edges_from((u, w) for u in range(n) for w in range(n) if (adj[u] >> w) & 1)
    return G


def c3_direct(adj):
    n = len(adj)
    return sum(1 for a, b, c in itertools.combinations(range(n), 3)
               if ((adj[a] >> b) & 1 and (adj[b] >> c) & 1 and (adj[c] >> a) & 1) or
                  ((adj[a] >> c) & 1 and (adj[c] >> b) & 1 and (adj[b] >> a) & 1))


ok_deg = ok_sc = ok_c3 = ok_strong = ok_aut = ok_aut_nx = ok_spine = ok_distinct = ok_vf2 = ok_conv = ok_pal = True
counts = {}
for m in range(3, 13):
    N = m - 1
    comps = list(L.compositions(N))
    counts[m] = len(comps)
    codes = {}
    for p in comps:
        adj = L.chain(p)
        deg = [bin(a).count('1') for a in adj]
        exp = [N - 1] + [N - t for t in range(1, N)] + [1]
        ok_deg &= deg == exp if N >= 2 else True
        if m >= 4:
            ok_sc &= sorted(deg) == [1, 1] + list(range(2, m - 2)) + [m - 2, m - 2]
        ok_c3 &= c3_direct(adj) == m - 2
        ok_strong &= L.is_strong(adj)
        code, aut = L.canon(adj, want_aut=True)
        codes[p] = code
        if m >= 4:
            ok_aut &= aut == 1
            if m <= 10:
                G = to_nx(adj)
                ok_aut_nx &= sum(1 for _ in DiGraphMatcher(G, G).isomorphisms_iter()) == 1
        else:
            ok_aut &= aut == 3
        if N >= 3:
            # recover the spine: descending score; ties at N-1 and at 1 broken by the arc between the tied pair
            order = sorted(range(m), key=lambda v: -deg[v])
            top = [v for v in range(m) if deg[v] == N - 1]; bot = [v for v in range(m) if deg[v] == 1]
            assert len(top) == 2 and len(bot) == 2
            a, b = top
            if not (adj[a] >> b) & 1:
                a, b = b, a
            c, d = bot
            if not (adj[c] >> d) & 1:
                c, d = d, c
            mid = sorted([v for v in range(m) if 1 < deg[v] < N - 1], key=lambda v: -deg[v])
            spine = [a, b] + mid + [c, d]
            pos = {v: i for i, v in enumerate(spine)}
            backs = sorted((pos[w], pos[u]) for u in range(m) for w in range(m) if (adj[u] >> w) & 1 and pos[w] < pos[u])
            # back arcs (lo, hi): reconstruct blocks
            rec = [hi - lo for lo, hi in backs]
            ok_spine &= spine == list(range(m)) and rec == list(p) and all(backs[i][1] == backs[i + 1][0] for i in range(len(backs) - 1)) \
                and backs[0][0] == 0 and backs[-1][1] == N
        conv = L.converse(adj)
        rev = L.chain(tuple(reversed(p)))
        relab = [0] * m
        for u in range(m):
            for w in range(m):
                if (conv[u] >> w) & 1:
                    relab[N - u] |= 1 << (N - w)
        ok_conv &= relab == rev
        ok_pal &= (L.canon(conv) == code) == (tuple(reversed(p)) == p)
    ok_distinct &= len(set(codes.values())) == len(comps)
    if m <= 9:
        G = {p: to_nx(L.chain(p)) for p in comps}
        for p1, p2 in itertools.combinations(comps, 2):
            if nx.is_isomorphic(G[p1], G[p2]):
                ok_vf2 = False
check(ok_deg and ok_sc, 'out-degrees N-1, N-t (1<=t<=N-1), 1; sorted scores (1,1,2,...,m-3,m-2,m-2), every composition, m = 4..12')
check(ok_c3, 'every chain of order m has exactly m - 2 cyclic triples (direct triple count), m = 3..12')
check(ok_strong, 'every chain (parts >= 2) is strong, m = 3..12')
check(ok_aut and ok_aut_nx, '|Aut(IC(p))| = 1 for m >= 4 (canon, m <= 12; networkx VF2, m <= 10); |Aut(C_3)| = 3')
check(ok_spine, 'the spine is recovered by the stated rule and the back arcs spell out the composition (N >= 3, m <= 12)')
check(ok_distinct and ok_vf2, 'chains of different compositions are non-isomorphic (canonical codes m <= 12; VF2 pairwise m <= 9); '
      f'compositions per order {[counts[m] for m in range(3, 13)]}')
check(ok_conv and ok_pal, 'converse(IC(p)) = IC(reversed p) exactly under t -> N - t; self-converse iff palindromic (m <= 12)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
