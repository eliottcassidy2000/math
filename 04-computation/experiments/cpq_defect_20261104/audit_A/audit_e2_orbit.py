#!/usr/bin/env python3
"""audit_A (e), independent orbit-method implementation of the (p, q, 1^k) census (p, q prime), written from the
statement of THM-4871 (1) and sharing no code with cpq_defect_census.py.
  * circulant classes: skew sets / symmetric sets on Z_m up to dilation, multiplier-even flag from the definition;
  * 4-coloured types on k points: (base type, Aut(base)-orbit of colourings), base types by min over S_k for k <= 5
    (nauty certificate for k = 6), automorphisms by direct check; the count is verified against the Burnside value
    b_k(4) (so the enumeration is complete: an invariant that takes exactly #orbits values is a complete invariant);
    colour-sign-even test for graphs from the definition (all colour-preserving automorphisms reverse an even
    number of edges);
  * full objects on n vertices: nauty certificate (directed for tournaments); sign-evenness from nauty generators;
  * optional VF2 cross-check of nauty (networkx): every coincident pair isomorphic, sampled distinct pairs not.
Usage: python3 audit_e2_orbit.py p q k [--vf2]
"""
import sys, math, itertools, time, random
from collections import defaultdict
import pynauty

def units(m): return [a for a in range(1, m) if math.gcd(a, m) == 1]

def dil_classes(m, skew):
    half = list(range(1, (m + 1) // 2)); seen, reps = set(), []
    for bits in itertools.product([0, 1], repeat=len(half)):
        if skew: S = frozenset((d if b == 0 else m - d) for d, b in zip(half, bits))
        else: S = frozenset(x for d, b in zip(half, bits) if b for x in (d, m - d))
        if S in seen: continue
        orb = {frozenset((a * s) % m for s in S) for a in units(m)}
        seen |= orb; reps.append(S)
    return reps

def mult_even(m, C):
    ed = sorted({(min(i, (i + c) % m), max(i, (i + c) % m)) for i in range(m) for c in C})
    for a in units(m):
        if frozenset((a * c) % m for c in C) == C:
            if sum(1 for (u, v) in ed if (a * u) % m > (a * v) % m) & 1: return False
    return True

def base_types(k, directed):
    """isomorphism types of tournaments/graphs on k points with their automorphism groups (explicit)."""
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    perms = list(itertools.permutations(range(k)))
    def relabel(obj, g):   # obj = frozenset of arcs (u,v) (tournament) or edges (u<v) (graph)
        if directed: return frozenset((g[u], g[v]) for (u, v) in obj)
        return frozenset((min(g[u], g[v]), max(g[u], g[v])) for (u, v) in obj)
    seen = set(); out = []
    for bits in itertools.product([0, 1], repeat=len(pairs)):
        if directed: obj = frozenset((i, j) if b else (j, i) for (i, j), b in zip(pairs, bits))
        else: obj = frozenset((i, j) for (i, j), b in zip(pairs, bits) if b)
        if obj in seen: continue
        images = {relabel(obj, g): g for g in perms}
        seen |= set(images.keys())
        aut = [g for g in perms if relabel(obj, g) == obj]
        out.append((obj, aut))
    return out

def coloured_types(k, directed):
    """list of (base obj, colouring tuple, colour-sign-even flag)."""
    if k == 0: return [(frozenset(), (), True)]
    res = []
    for obj, aut in base_types(k, directed):
        edges = sorted(obj) if not directed else None
        seen = set()
        for col in itertools.product(range(4), repeat=k):
            if col in seen: continue
            orbit = set()
            for g in aut:
                nc = [0] * k
                for v in range(k): nc[g[v]] = col[v]
                orbit.add(tuple(nc))
            seen |= orbit
            even = True
            if not directed:
                stab = [g for g in aut if all(col[g[v]] == col[v] for v in range(k))]
                even = all((sum(1 for (u, v) in edges if g[u] > g[v]) & 1) == 0 for g in stab)
            res.append((obj, col, even))
    return res

def burnside_b(k):
    from fractions import Fraction
    def parts(n, mx=None):
        mx = n if mx is None else mx
        if n == 0: yield (); return
        for p in range(min(n, mx), 0, -1):
            for r in parts(n - p, p): yield (p,) + r
    tot = Fraction(0)
    for lam in parts(k):
        if all(l % 2 for l in lam):
            z = 1
            from collections import Counter
            for part, mult in Counter(lam).items(): z *= part ** mult * math.factorial(mult)
            c2 = sum((l - 1) // 2 for l in lam) + sum(math.gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
            tot += Fraction(4 ** len(lam) * 2 ** c2, z)
    return int(tot)

def full_object(p, q, k, SP, SQ, cross, obj, col, directed):
    n = p + q + k; adj = {v: [] for v in range(n)}
    def put(u, v):
        adj[u].append(v)
        if not directed: adj[v].append(u)
    for i in range(p):
        for s in SP:
            j = (i + s) % p
            if directed or i < j: put(i, j)
    for i in range(q):
        for s in SQ:
            j = (i + s) % q
            if directed or i < j: put(p + i, p + j)
    for i in range(p):
        for j in range(q):
            if directed: put(i, p + j) if cross == 0 else put(p + j, i)
            elif cross: put(i, p + j)
    for t in range(k):
        f = p + q + t; x, y = col[t] & 1, col[t] >> 1
        for i in range(p):
            if directed: put(i, f) if x else put(f, i)
            elif x: put(i, f)
        for j in range(q):
            if directed: put(p + j, f) if y else put(f, p + j)
            elif y: put(p + j, f)
    for (u, v) in obj: put(p + q + u, p + q + v)
    return adj

def main(p, q, k, vf2=False):
    t0 = time.time(); n = p + q + k
    tP, tQ = dil_classes(p, True), dil_classes(q, True)
    eP = [C for C in dil_classes(p, False) if mult_even(p, C)]; eQ = [C for C in dil_classes(q, False) if mult_even(q, C)]
    cT = coloured_types(k, True); cG = coloured_types(k, False); cGe = [x for x in cG if x[2]]
    bk = burnside_b(k)
    print("(p,q,k)=(%d,%d,%d) n=%d: t=(%d,%d) e=(%d,%d); coloured tournament types %d, coloured graph types %d, colour-sign-even %d; Burnside b_k(4)=%d  %s  [%.1fs]" % (
        p, q, k, n, len(tP), len(tQ), len(eP), len(eQ), len(cT), len(cG), len(cGe), bk,
        "OK (b_k(4) = coloured tournament types = colour-sign-even coloured graph types)" if len(cT) == bk == len(cGe) else "MISMATCH", time.time() - t0), flush=True)
    # tournaments
    typesT = defaultdict(list); R = 0
    for SP in tP:
        for SQ in tQ:
            for cross in (0, 1):
                for (obj, col, _) in cT:
                    adj = full_object(p, q, k, SP, SQ, cross, obj, col, True)
                    typesT[pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adj))].append((SP, SQ, cross, obj, col)); R += 1
    tau_T = len(typesT); E_T = sum(len(v) - 1 for v in typesT.values())
    print("  tournaments: R=%d tau_T=%d E_T=%d  [%.1fs]" % (R, tau_T, E_T, time.time() - t0), flush=True)
    typesG = {}; RG = 0
    for CP in eP:
        for CQ in eQ:
            for cross in (0, 1):
                for (obj, col, _) in cGe:
                    adj = full_object(p, q, k, CP, CQ, cross, obj, col, False)
                    g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
                    cert = pynauty.certificate(g); RG += 1
                    if cert in typesG:
                        typesG[cert][0] += 1; continue
                    edges = [(u, v) for u in adj for v in adj[u] if u < v]
                    gens = pynauty.autgrp(g)[0]
                    even = all((sum(1 for (u, v) in edges if gg[u] > gg[v]) & 1) == 0 for gg in gens)
                    typesG[cert] = [1, even, (CP, CQ, cross, obj, col)]
    tau_E = sum(1 for m, ev, _ in typesG.values() if ev); E_even = sum(m - 1 for m, ev, _ in typesG.values() if ev)
    W = sum(1 for m, ev, _ in typesG.values() if not ev); E_wit = sum(m - 1 for m, ev, _ in typesG.values() if not ev)
    defect = tau_T - tau_E
    print("  graphs: RG=%d tau_E=%d E_even=%d W=%d E_wit=%d; defect=%d; identity W+E_even+E_wit-E_T=%d %s; Theorem R %s  [%.1fs]" % (
        RG, tau_E, E_even, W, E_wit, defect, W + E_even + E_wit - E_T, "OK" if defect == W + E_even + E_wit - E_T else "FAILS",
        "OK" if R == RG else "FAILS", time.time() - t0), flush=True)
    if vf2:
        import networkx as nx
        from networkx.algorithms.isomorphism import DiGraphMatcher, GraphMatcher
        def nxT(d):
            G = nx.DiGraph(); G.add_nodes_from(range(n))
            for SP, SQ, cross, obj, col in [d]:
                adj = full_object(p, q, k, SP, SQ, cross, obj, col, True)
            for u in adj:
                for v in adj[u]: G.add_edge(u, v)
            return G
        def nxG(d):
            G = nx.Graph(); G.add_nodes_from(range(n))
            CP, CQ, cross, obj, col = d
            adj = full_object(p, q, k, CP, CQ, cross, obj, col, False)
            for u in adj:
                for v in adj[u]: G.add_edge(u, v)
            return G
        bad = 0; pos = 0
        for cert, lst in typesT.items():
            if len(lst) >= 2:
                for a, b in itertools.combinations(lst, 2):
                    pos += 1
                    if not DiGraphMatcher(nxT(a), nxT(b)).is_isomorphic(): bad += 1
        random.seed(5); reps = [v[0] for v in typesT.values()]; neg = 0
        for _ in range(300):
            a, b = random.sample(reps, 2); neg += 1
            if DiGraphMatcher(nxT(a), nxT(b)).is_isomorphic(): bad += 1
        print("  VF2 tournaments: %d coincident pairs checked isomorphic, %d sampled distinct pairs non-isomorphic, failures=%d" % (pos, neg, bad), flush=True)
        # graphs: need all members per type; recompute members for coincident types
        members = defaultdict(list)
        for CP in eP:
            for CQ in eQ:
                for cross in (0, 1):
                    for (obj, col, _) in cGe:
                        adj = full_object(p, q, k, CP, CQ, cross, obj, col, False)
                        cert = pynauty.certificate(pynauty.Graph(n, directed=False, adjacency_dict=adj))
                        if typesG[cert][0] >= 2: members[cert].append((CP, CQ, cross, obj, col))
        bad = 0; pos = 0
        for cert, lst in members.items():
            for a, b in itertools.combinations(lst, 2):
                pos += 1
                if not GraphMatcher(nxG(a), nxG(b)).is_isomorphic(): bad += 1
        repsG = [v[2] for v in typesG.values()]; neg = 0
        for _ in range(300):
            a, b = random.sample(repsG, 2); neg += 1
            if GraphMatcher(nxG(a), nxG(b)).is_isomorphic(): bad += 1
        print("  VF2 graphs: %d coincident pairs checked isomorphic, %d sampled distinct pairs non-isomorphic, failures=%d  [%.1fs]" % (pos, neg, bad, time.time() - t0), flush=True)
        # witnesses listing for small k
        if k <= 3:
            for cert, (m, ev, d) in typesG.items():
                if not ev: print("  witness rep:", d)
    return dict(R=R, tau_T=tau_T, tau_E=tau_E, defect=defect, W=W, E_T=E_T, E_even=E_even, E_wit=E_wit)

if __name__ == "__main__":
    p, q, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    main(p, q, k, vf2=("--vf2" in sys.argv))
