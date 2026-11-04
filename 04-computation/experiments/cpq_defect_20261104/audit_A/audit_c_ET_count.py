#!/usr/bin/env python3
"""audit_A (c): THM-4871 statement 3(c), E_T(p) = 4 t_p t_q.
By the proof of 3(c) (which this audit accepts), the tournament types with m_T = 2 at k = p are exactly those whose
fixed-point 4-coloured tournament kappa is circulant on F with a constant colour, and every such type contains
exactly two N-orbits.  So E_T(p) = (number of N-orbits with circulant constant-colour kappa) / 2, and that number
is t_p (S_P) * t_q (S_Q) * 2 (cross) * 4 (colour) * t_p (circulant tournament on F up to isomorphism) = 8 t_p^2 t_q,
NOT the 8 t_p t_q of the theorem (the factor t_p for S_P is dropped there).  The census cannot see the difference
because t_3 = t_5 = 1.  Here the subfamily is enumerated explicitly for (3,5), (3,7), (5,7), (7,11), (7,13) at
k = p and the types are counted with nauty (directed certificates) and cross-checked with networkx VF2.
Also: t_7, t_11, t_13 by brute force, and the datum transformation of the proof checked on an explicit pair.
"""
import itertools, math, random
from collections import defaultdict
import pynauty, networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def units(m): return [a for a in range(1, m) if math.gcd(a, m) == 1]

def skew_classes(m):
    half = list(range(1, (m + 1) // 2)); seen, reps = set(), []
    for bits in itertools.product([0, 1], repeat=len(half)):
        S = frozenset((d if b == 0 else m - d) for d, b in zip(half, bits))
        if S in seen: continue
        orb = {frozenset((a * s) % m for s in S) for a in units(m)}
        seen |= orb; reps.append(S)
    return reps

def build_T(p, q, k, Sp, Sq, cross, Sf, colour):
    """H-fixed tournament with kappa = circulant tournament (Z_k, Sf) on F, constant colour (x, y)."""
    n = p + q + k; adj = {v: [] for v in range(n)}
    x, y = colour
    for i in range(p):
        for s in Sp: adj[i].append((i + s) % p)
    for i in range(q):
        for s in Sq: adj[p + i].append(p + (i + s) % q)
    for i in range(p):
        for j in range(q):
            if cross == 0: adj[i].append(p + j)
            else: adj[p + j].append(i)
    for t in range(k):
        f = p + q + t
        for i in range(p):
            if x: adj[i].append(f)
            else: adj[f].append(i)
        for j in range(q):
            if y: adj[p + j].append(f)
            else: adj[f].append(p + j)
    for i in range(k):
        for s in Sf: adj[p + q + i].append(p + q + (i + s) % k)
    return adj

def to_nx(adj):
    G = nx.DiGraph(); G.add_nodes_from(adj)
    for u in adj:
        for v in adj[u]: G.add_edge(u, v)
    return G

if __name__ == "__main__":
    print("t_m by brute force:", {m: len(skew_classes(m)) for m in (3, 5, 7, 11, 13)})
    for (p, q) in [(3, 5), (3, 7), (5, 7), (7, 11), (7, 13)]:
        k = p; n = p + q + k
        tp, tq = len(skew_classes(p)), len(skew_classes(q))
        types = defaultdict(list); count = 0
        for Sp in skew_classes(p):
            for Sq in skew_classes(q):
                for cross in (0, 1):
                    for Sf in skew_classes(k):
                        for colour in itertools.product((0, 1), repeat=2):
                            adj = build_T(p, q, k, Sp, Sq, cross, Sf, colour)
                            cert = pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adj))
                            types[cert].append((Sp, Sq, cross, Sf, colour, adj)); count += 1
        mults = sorted(len(v) for v in types.values())
        ET = sum(len(v) - 1 for v in types.values())
        print("(p,q)=(%d,%d) k=%d n=%d: t_p=%d t_q=%d; circulant constant-colour N-orbits=%d (8 t_p^2 t_q = %d, 8 t_p t_q = %d); "
              "types=%d; multiplicities=%s; E_T contribution=%d  [theorem: 4 t_p t_q = %d; corrected 4 t_p^2 t_q = %d]" % (
                  p, q, k, n, tp, tq, count, 8 * tp * tp * tq, 8 * tp * tq, len(types), sorted(set(mults)), ET, 4 * tp * tq, 4 * tp * tp * tq))
        # VF2 cross-check: every coincident pair is really isomorphic; a sample of distinct-certificate pairs is not
        ok = True; checked = 0
        for cert, lst in types.items():
            for a, b in itertools.combinations(lst, 2):
                if not DiGraphMatcher(to_nx(a[5]), to_nx(b[5])).is_isomorphic(): ok = False
                checked += 1
        reps = [v[0][5] for v in types.values()]
        random.seed(3); neg = 0
        for _ in range(min(40, len(reps) * (len(reps) - 1) // 2)):
            a, b = random.sample(reps, 2)
            if DiGraphMatcher(to_nx(a), to_nx(b)).is_isomorphic(): ok = False
            neg += 1
        print("   VF2: %d coincident pairs isomorphic, %d distinct-certificate pairs non-isomorphic: %s" % (checked, neg, "OK" if ok else "FAILED"))
        # datum transformation: (Sp, Sq, cross, Sf, (x,y)) ~ (Sf, Sq, y, Sp, (1-x, cross)) up to dilation of Sf, Sp
        bad = 0
        for cert, lst in types.items():
            (Sp, Sq, cross, Sf, (x, y), _) = lst[0]
            partner = build_T(p, q, k, Sf, Sq, y, Sp, (1 - x, cross))
            if pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=partner)) != cert: bad += 1
        print("   datum transformation (cross,(x,y),S_Q,kappa) -> (y,(1-x,cross),S_Q,T[P]) lands in the same type for all %d types: %s" % (len(types), "yes" if bad == 0 else "NO (%d fail)" % bad))
