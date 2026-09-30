#!/usr/bin/env python3
"""collatz_pentagon_gadget_20260930.py -- an expansion certificate for the Paley graph P_13 under the pentagon operator
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, thirteenth note, part A).

 Monotonicity: an induced subgraph H of G gives C_5(H) as an induced subgraph of C_5(G); C_5 distributes over disjoint unions.
 So if some iterate C_5^k(P_13) contains an induced icosahedron with a tadpole hat (I_1 of Gervacio-Maehara-Ramos, Theorem 3.5),
 P_13 is pentagon-expanding.  This script (1) describes C_5(P_13) (39 vertices, Cayley on Z_13 x| Z_3), (2) searches C_5(P_13)
 for induced P_13, induced icosahedra and induced I_1, (3) if needed searches C_5^2(P_13) (4563 vertices) for induced icosahedra
 and hats, from orbit representatives under Aut(P_13) = Z_13 x| Z_6, (4) verifies the hat-doubling mechanism on I_1 directly.
Usage: python3 collatz_pentagon_gadget_20260930.py
"""
import itertools, time
from collections import Counter

T0 = time.time()


def graph_from_edges(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        if a != b:
            adj[a].add(b); adj[b].add(a)
    return adj


def induced_pentagons(adj, limit=10 ** 6):
    n = len(adj); out = set()
    for a in range(n):
        for b in adj[a]:
            if b < a:
                continue
            for c in adj[b]:
                if c <= a or c == b or c in adj[a]:
                    continue
                for d in adj[c]:
                    if d <= a or d in (b, c) or d in adj[a] or d in adj[b]:
                        continue
                    for e in adj[d]:
                        if e <= a or e in (b, c, d) or e not in adj[a] or e in adj[b] or e in adj[c]:
                            continue
                        out.add(min((a, b, c, d, e), (a, e, d, c, b)))
                        if len(out) > limit:
                            return None
    return sorted(out)


def pentagon_graph(adj):
    pents = induced_pentagons(adj)
    edges_of = [set((min(p[i], p[(i + 1) % 5]), max(p[i], p[(i + 1) % 5])) for i in range(5)) for p in pents]
    m = len(pents)
    # adjacency through shared edges: index pentagons by edge
    by_edge = {}
    for i, es in enumerate(edges_of):
        for e in es:
            by_edge.setdefault(e, []).append(i)
    adj2 = [set() for _ in range(m)]
    for e, lst in by_edge.items():
        for i in lst:
            for j in lst:
                if i != j:
                    adj2[i].add(j)
    return adj2, pents


def circulant(n, S):
    return graph_from_edges(n, [(i, (i + s) % n) for i in range(n) for s in S])


def icosahedron():
    e = [(0, i) for i in range(1, 6)] + [(i, 1 + i % 5) for i in range(1, 6)]
    e += [(11, i) for i in range(6, 11)] + [(6 + i, 6 + (i + 1) % 5) for i in range(5)]
    e += [(1 + i, 6 + i) for i in range(5)] + [(1 + i, 6 + (i + 1) % 5) for i in range(5)]
    return graph_from_edges(12, e)


ICO = icosahedron()


def is_induced_icosahedron(adj, verts):
    verts = list(verts)
    if len(verts) != 12:
        return False
    vs = set(verts)
    degs = sorted(len(adj[v] & vs) for v in verts)
    if degs != [5] * 12:
        return False
    # locally pentagonal check: each link is a 5-cycle
    for v in verts:
        N = adj[v] & vs
        if any(len(adj[u] & N) != 2 for u in N):
            return False
    return True


def find_icosahedra(adj, starts, max_found=5, time_budget=None):
    """induced icosahedra containing a start vertex v as a pole: link pentagon a_0..a_4, lower ring b_i, antipode w."""
    found = []
    for v in starts:
        if time_budget is not None and time.time() - T0 > time_budget:
            break
        N = adj[v]
        Nl = sorted(N); pos = {u: i for i, u in enumerate(Nl)}
        sub = [set(pos[u] for u in adj[x] & N) for x in Nl]
        pents = induced_pentagons(sub, limit=300000)
        if pents is None:
            print("  (pole %d: more than 300000 induced pentagons in its link; skipped)" % v, flush=True)
            continue
        for p in pents:
            A = [Nl[i] for i in p]     # cyclic order a_0..a_4
            cands = []
            ok = True
            for i in range(5):
                ai, aj = A[i], A[(i + 1) % 5]
                others = [A[j] for j in range(5) if j not in (i, (i + 1) % 5)]
                c = [b for b in (adj[ai] & adj[aj]) if b != v and b not in N and not any(b in adj[o] for o in others)]
                if not c:
                    ok = False; break
                cands.append(c)
            if not ok:
                continue
            for B in itertools.product(*cands):
                if len(set(B)) != 5:
                    continue
                if not all(B[(i + 1) % 5] in adj[B[i]] for i in range(5)):
                    continue
                if any(B[(i + 2) % 5] in adj[B[i]] for i in range(5)):
                    continue
                W = set.intersection(*[adj[b] for b in B]) - N - {v} - set(B)
                for w in W:
                    if any(w in adj[a] for a in A):
                        continue
                    verts = [v] + A + list(B) + [w]
                    if is_induced_icosahedron(adj, verts):
                        found.append(verts)
                        if len(found) >= max_found:
                            return found
    return found


def hats(adj, ico):
    vs = set(ico); out = []
    for u in range(len(adj)):
        if u in vs:
            continue
        Nu = adj[u] & vs
        if len(Nu) != 4:
            continue
        inner = sum(1 for x in Nu for y in Nu if x < y and y in adj[x])
        degs = sorted(len(adj[x] & Nu) for x in Nu)
        if inner == 4 and degs == [1, 2, 2, 3]:
            out.append(u)
    return out


def induced_subgraph_iso(pattern, host, time_budget=None):
    """is `pattern` an induced subgraph of `host`? backtracking in BFS order of the pattern."""
    n = len(pattern); order = []; seen = set()
    for s in range(n):
        if s in seen:
            continue
        seen.add(s); q = [s]
        while q:
            x = q.pop(0); order.append(x)
            for y in sorted(pattern[x]):
                if y not in seen:
                    seen.add(y); q.append(y)
    mapping = {}; used = set(); H = len(host)
    def bt(i):
        if time_budget is not None and time.time() - T0 > time_budget:
            return None
        if i == n:
            return True
        x = order[i]
        # candidates: neighbours of an already-mapped pattern-neighbour, else all
        prev = [y for y in pattern[x] if y in mapping]
        cand = set.intersection(*[host[mapping[y]] for y in prev]) if prev else set(range(H))
        for w in cand:
            if w in used or len(host[w]) < len(pattern[x]):
                continue
            if all(((y in pattern[x]) == (mapping[y] in host[w])) for y in mapping):
                mapping[x] = w; used.add(w)
                r = bt(i + 1)
                if r:
                    return True
                if r is None:
                    return None
                del mapping[x]; used.discard(w)
        return False
    return bt(0)


def part1():
    print("== (1) C_5(P_13) ==")
    P13 = circulant(13, [1, 3, 4, 9, 10, 12])
    C1, pents = pentagon_graph(P13)
    print(" P_13: %d induced pentagons; C_5(P_13): %d vertices, degrees %s" % (len(pents), len(C1), sorted(Counter(len(a) for a in C1).items())))
    # the affine group x -> ax + b (a a square) acting on pentagons
    sq = [1, 3, 4, 9, 10, 12]
    idx = {p: i for i, p in enumerate(pents)}
    def canon(t):
        # canonical form of a 5-cycle given as a cyclic sequence of vertices
        best = None
        for r in range(5):
            for d in (1, -1):
                seq = tuple(t[(r + d * i) % 5] for i in range(5))
                if best is None or seq < best:
                    best = seq
        return best
    perms = []
    for a in sq:
        for b in range(13):
            perms.append(tuple(idx[canon(tuple((a * x + b) % 13 for x in p))] for p in pents))
    orbits = set()
    for i in range(len(pents)):
        orbits.add(min(pi[i] for pi in perms))
    free39 = all(all(pi[i] != i for i in range(len(pents))) for pi in perms if pi != tuple(range(len(pents))) and True)
    sub39 = [tuple(idx[canon(tuple((a * x + b) % 13 for x in p))] for p in pents) for a in (1, 3, 9) for b in range(13)]
    free = all(all(pi[i] != i for i in range(39)) for pi in sub39 if pi != tuple(range(39)))
    print(" Aut(P_13) (order 78) acting on the 39 pentagons: orbits %d (vertex-transitive: %s); the index-2 subgroup Z_13 x| Z_3 acts freely: %s -> C_5(P_13) is a Cayley graph of the non-abelian group of order 39" % (len(orbits), len(orbits) == 1, free))
    return P13, C1, pents, perms


def part2(P13, C1):
    print("== (2) what C_5(P_13) contains ==")
    t = time.time()
    r = induced_subgraph_iso(P13, C1, time_budget=120)
    print(" induced P_13 inside C_5(P_13): %s (%.0fs)" % (r, time.time() - t))
    icos = find_icosahedra(C1, range(len(C1)), max_found=50)
    print(" induced icosahedra in C_5(P_13) (found, with a pole search from every vertex): %d" % len(icos))
    if icos:
        hs = [hats(C1, I) for I in icos]
        print("  hats on them: %s" % [len(h) for h in hs][:20])
        ok = [I for I, h in zip(icos, hs) if h]
        if ok:
            print("  EXPANSION CERTIFICATE: C_5(P_13) contains an induced icosahedron %s with hat(s) %s" % (ok[0], hats(C1, ok[0])))
            return True
    r = induced_subgraph_iso(ICO, C1, time_budget=120)
    print(" generic induced-subgraph search for the icosahedron in C_5(P_13): %s" % r)
    return False


def part3(P13, C1, pents, perms):
    print("== (3) C_5^2(P_13): icosahedra and hats from orbit representatives ==")
    t = time.time()
    C2, pents2 = pentagon_graph(C1)
    print(" C_5^2(P_13): %d vertices, %d edges, degree range %d..%d (%.0fs)" % (len(C2), sum(len(a) for a in C2) // 2, min(len(a) for a in C2), max(len(a) for a in C2), time.time() - t))
    # orbit representatives: the automorphisms act on pentagons of C_5(P_13) through their action on the 39 vertices
    idx2 = {p: i for i, p in enumerate(pents2)}
    def canon(t):
        best = None
        for r in range(5):
            for d in (1, -1):
                seq = tuple(t[(r + d * i) % 5] for i in range(5))
                if best is None or seq < best:
                    best = seq
        return best
    reps = {}
    for i, p in enumerate(pents2):
        m = min(idx2[canon(tuple(pi[x] for x in p))] for pi in perms)
        reps.setdefault(m, 0); reps[m] += 1
    print(" orbits of Aut(P_13) on the %d vertices of C_5^2(P_13): %d, sizes %s" % (len(pents2), len(reps), sorted(Counter(reps.values()).items())))
    # explicit icosahedra: for every induced icosahedron I of C_5(P_13), its twelve link pentagons are induced pentagons of
    # C_5(P_13), i.e. vertices of C_5^2(P_13), and by monotonicity they induce C_5(I) = I there
    icos1 = find_icosahedra(C1, range(len(C1)), max_found=400)
    print(" induced icosahedra of C_5(P_13) found by the pole search: %d" % len(icos1))
    lifted = []; seen = set()
    for I in icos1:
        vs = set(I)
        links = []
        for v in I:
            N = C1[v] & vs
            # order the link pentagon cyclically
            start = next(iter(N)); cyc = [start]; prev = None
            while len(cyc) < 5:
                nxt = [u for u in C1[cyc[-1]] & N if u != prev and u not in cyc]
                prev = cyc[-1]; cyc.append(nxt[0])
            links.append(canon(tuple(cyc)))
        key = tuple(sorted(idx2[l] for l in links))
        if key in seen:
            continue
        seen.add(key); lifted.append([idx2[l] for l in links])
    print(" distinct lifted icosahedra in C_5^2(P_13) (the pentagon graphs of the icosahedra of C_5(P_13)): %d; all induced icosahedra: %s" % (
        len(lifted), all(is_induced_icosahedron(C2, L) for L in lifted)))
    nh = [len(hats(C2, L)) for L in lifted]
    print(" hats on the lifted icosahedra: %s" % sorted(Counter(nh).items()))
    for L in lifted:
        h = hats(C2, L)
        if h:
            print("  EXPANSION CERTIFICATE: C_5^2(P_13) contains the induced icosahedron %s with tadpole hat(s) %s -> P_13 is pentagon-expanding" % (L, h[:4]))
            return True
    icos = find_icosahedra(C2, sorted(reps), max_found=20, time_budget=420)
    print(" pole search from orbit representatives in C_5^2(P_13): %d icosahedra (%.0fs)" % (len(icos), time.time() - T0))
    for I in icos[:10]:
        h = hats(C2, I)
        print("  icosahedron %s: hats %d" % (I, len(h)))
        if h:
            print("  EXPANSION CERTIFICATE: C_5^2(P_13) contains an induced icosahedron with a tadpole hat -> P_13 is pentagon-expanding")
            return True
    return False


def part4():
    print("== (4) the hat-doubling mechanism on the icosahedron, verified directly ==")
    I = icosahedron()
    # tadpole in I: triangle 0,1,2 with pendant 3 attached to 1? need 3 ~ 1 only: vertices 0,1,2 triangle (0~1, 1~2, 0~2? 0~1 yes, 0~2 yes, 1~2 yes); 3 ~ 2 and 3 ~ 0: not a pendant. use triangle (0,1,2) and pendant 3? 3 ~ 0 and 3 ~ 2 -> no. choose pendant 6: 6 ~ 1 (yes), 6 ~ 0? no, 6 ~ 2? 2 ~ 6 and 2 ~ 7 -> yes. pendant 10: 10 ~ 5 and 10 ~ 1; 10 ~ 0? no; 10 ~ 2? no -> tadpole (0,1,2)+10 attached to 1
    tad = None
    for x in range(12):
        for y in I[x]:
            for z in I[x] & I[y]:
                for t in I[y] - I[x] - I[z] - {x, z}:
                    tad = [x, y, z, t]; break
                if tad: break
            if tad: break
        if tad: break
    print(" tadpole T_(3,1) in the icosahedron: triangle %s with pendant %d attached to %d" % (tad[:3], tad[3], tad[1]))
    def add_hat(adj, verts):
        n = len(adj); new = [set(a) for a in adj] + [set()]
        for x in verts:
            new[x].add(n); new[n].add(x)
        return new
    I1 = add_hat(I, tad)
    G = I1; sizes = [len(G)]
    for k in range(4):
        G, _ = pentagon_graph(G)
        sizes.append(len(G))
        icos = find_icosahedra(G, range(len(G)), max_found=3)
        nh = [len(hats(G, Ic)) for Ic in icos]
        print(" C_5^%d(I_1): %d vertices; induced icosahedra found %d; hats on the first: %s" % (k + 1, len(G), len(icos), nh[:3]))
        if len(G) > 60:
            break
    print(" sizes along the trajectory of I_1:", sizes)


if __name__ == "__main__":
    P13, C1, pents, perms = part1()
    cert = part2(P13, C1)
    if not cert:
        cert = part3(P13, C1, pents, perms)
    part4()
    print("total time %.0fs" % (time.time() - T0))
