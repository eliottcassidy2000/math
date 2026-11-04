#!/usr/bin/env python3
"""cpq_defect_census.py -- mac-mini-2026-11-04-oddprimes (THM-4871 lane).

HYP-9340 beyond p-groups: for H = C_pq = <(0 .. p-1)(p .. p+q-1)> <= S_n of orbit type (p, q, 1^k), n = p + q + k,
p, q coprime odd prime powers (the orbit sizes; H = C_pq), compute
  R       = number of N-orbits on H-fixed tournaments  (N = N_{S_n}(H) = AGL(1,p) x AGL(1,q) x S_k),
  tau_T   = number of isomorphism types of tournaments whose automorphism group contains a conjugate of H,
  tau_E   = number of isomorphism types of sign-even graphs whose automorphism group contains a conjugate of H,
  W       = number of sign-odd graph types having an N-even H-fixed representative (witnesses),
  E_T, E_even, E_wit = multiplicity excesses (THM-4728 statement 1),
and check the identity tau_T - tau_E = W + E_even + E_wit - E_T and Theorem R (N-orbits on H-fixed tournaments =
N-orbits on N-even H-fixed graphs).

Sign convention (Royle-Praeger-Glasby-Freedman-Devillers; census script of S45): for an automorphism g of a graph X,
sgn_X(g) = (-1)^#{edges {u<v} of X with g(u) > g(v)}; this is a homomorphism Aut X -> {+-1} and is independent of the
vertex order (adjacent transposition argument); X is sign-even iff sgn_X = 1 on Aut X, N-even iff sgn_X = 1 on N cap Aut X.

N-orbit representatives (proved in the note): an H-fixed tournament is (circulant tournament on Z_p, circulant
tournament on Z_q, cross orientation P->Q or Q->P, a 4-coloured tournament on the k fixed points: colour (x, y) with
x = 1 iff P -> f, y = 1 iff Q -> f), and the N-orbits are (multiplier type on Z_p) x (multiplier type on Z_q) x {2} x
(4-coloured tournament types on k points).  An H-fixed graph is N-even iff its Z_p and Z_q circulant parts are
multiplier-even and its 4-coloured graph on the fixed points is colour-sign-even (product criterion: the reversal
parity of (alpha, beta, pi) in N is c(alpha) + c(beta) + c(pi), each computed inside its own block).

Usage: python3 cpq_defect_census.py p q kmax [--brute kb]   (brute-force cross-check of all 2^c H-fixed objects for k <= kb)
Requires pynauty (pip install pynauty).
"""
import sys, itertools, time
from collections import defaultdict
import pynauty

def reversal_parity(g, edges):
    return sum(1 for (u, v) in edges if g[u] > g[v]) & 1

def all_types(k, directed):
    """isomorphism types of tournaments (directed) or graphs on k points, as adjacency dicts."""
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    seen = {}
    for bits in itertools.product([0, 1], repeat=len(pairs)):
        adj = {v: [] for v in range(k)}
        for (i, j), b in zip(pairs, bits):
            if directed:
                if b: adj[i].append(j)
                else: adj[j].append(i)
            elif b:
                adj[i].append(j); adj[j].append(i)
        g = pynauty.Graph(k, directed=directed, adjacency_dict=adj)
        c = pynauty.certificate(g)
        if c not in seen: seen[c] = adj
    return list(seen.values())

def coloured_types(k, directed):
    """4-coloured tournaments (directed) / graphs on k points up to colour-preserving isomorphism.
    Returns list of (adj, col, even) where even = colour-preserving automorphism group is reversal-even (graphs).
    The key must carry the colour-class sizes: nauty's canonical form is canonical only within a fixed ordered
    partition shape (e.g. the transitive 3-tournament with colourings (1,0,0) and (1,1,0) gets the same packed form)."""
    if k == 0:
        return [({}, (), True)]
    out = {}
    for adj in all_types(k, directed):
        edges = sorted((u, v) for u in adj for v in adj[u] if u < v) if not directed else None
        for col in itertools.product(range(4), repeat=k):
            cells = [set(v for v in range(k) if col[v] == c) for c in range(4)]
            present = tuple(c for c in range(4) if cells[c])
            g = pynauty.Graph(k, directed=directed, adjacency_dict=adj,
                              vertex_coloring=[cells[c] for c in present])
            key = (tuple((c, len(cells[c])) for c in present), pynauty.certificate(g))   # cell shape + canonical form
            if key in out: continue
            if directed:
                out[key] = (adj, col, True)
            else:
                gens = pynauty.autgrp(g)[0]
                even = all(reversal_parity(gen, edges) == 0 for gen in gens)
                out[key] = (adj, col, even)
    return list(out.values())

def units(m):
    import math
    return [a for a in range(1, m) if math.gcd(a, m) == 1]

def circulant_tournament_types(p):
    half = list(range(1, (p + 1) // 2))
    seen, reps = set(), []
    for bits in itertools.product([0, 1], repeat=len(half)):
        S = frozenset((d if b == 0 else p - d) for d, b in zip(half, bits))
        if S in seen: continue
        orb = {frozenset((a * s) % p for s in S) for a in units(p)}
        seen |= orb; reps.append(S)
    return reps

def circulant_graph_types(p):
    """symmetric sets C up to multipliers, with the multiplier-even flag."""
    half = list(range(1, (p + 1) // 2))
    seen, reps = set(), []
    for bits in itertools.product([0, 1], repeat=len(half)):
        C = frozenset(x for d, b in zip(half, bits) if b for x in (d, p - d))
        if C in seen: continue
        orb = {frozenset((a * c) % p for c in C) for a in units(p)}
        seen |= orb
        edges = sorted({(min(i, (i + c) % p), max(i, (i + c) % p)) for i in range(p) for c in C})
        even = True
        for a in units(p):
            if frozenset((a * c) % p for c in C) == C:
                g = [(a * i) % p for i in range(p)]
                if reversal_parity(g, edges): even = False; break
        reps.append((C, even))
    return reps

def build(p, q, k, Sp, Sq, cross, adjF, col, directed):
    """full object on n = p+q+k vertices: P = 0..p-1, Q = p..p+q-1, F = p+q..n-1."""
    n = p + q + k
    adj = {v: [] for v in range(n)}
    def arc(u, v):
        adj[u].append(v)
        if not directed: adj[v].append(u)
    for i in range(p):
        for s in Sp:
            j = (i + s) % p
            if directed or i < j: arc(i, j)
    for i in range(q):
        for s in Sq:
            j = (i + s) % q
            if directed or i < j: arc(p + i, p + j)
    for i in range(p):
        for j in range(q):
            if directed:
                if cross == 0: arc(i, p + j)
                else: arc(p + j, i)
            elif cross: arc(i, p + j)
    for t in range(k):
        f = p + q + t
        x, y = col[t] & 1, col[t] >> 1
        for i in range(p):
            if directed:
                if x: arc(i, f)
                else: arc(f, i)
            elif x: arc(i, f)
        for j in range(q):
            if directed:
                if y: arc(p + j, f)
                else: arc(f, p + j)
            elif y: arc(p + j, f)
    for u in adjF:
        for v in adjF[u]:
            if directed or u < v: arc(p + q + u, p + q + v)
    return adj

def census(p, q, k, verbose=True):
    t0 = time.time()
    n = p + q + k
    tT = circulant_tournament_types(p); tQ = circulant_tournament_types(q)
    gP = [C for C, ev in circulant_graph_types(p) if ev]; gQ = [C for C, ev in circulant_graph_types(q) if ev]
    colT = coloured_types(k, True)
    colG = [(adj, col, ev) for (adj, col, ev) in coloured_types(k, False) if ev]
    # tournaments
    types_T = defaultdict(int)
    R = 0
    for Sp in tT:
        for Sq in tQ:
            for cross in (0, 1):
                for adjF, col, _ in colT:
                    adj = build(p, q, k, Sp, Sq, cross, adjF, col, True)
                    cert = pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adj))
                    types_T[cert] += 1; R += 1
    tau_T = len(types_T); E_T = R - tau_T
    # graphs: N-even N-orbits
    types_G = {}   # cert -> [count, sign_even]
    RG = 0; witness_examples = {}
    for Cp in gP:
        for Cq in gQ:
            for cross in (0, 1):
                for adjF, col, _ in colG:
                    adj = build(p, q, k, Cp, Cq, cross, adjF, col, False)
                    g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
                    cert = pynauty.certificate(g)
                    RG += 1
                    if cert in types_G:
                        types_G[cert][0] += 1; continue
                    gens = pynauty.autgrp(g)[0]
                    edges = [(u, v) for u in adj for v in adj[u] if u < v]
                    even = all(reversal_parity(gen, edges) == 0 for gen in gens)
                    types_G[cert] = [1, even]
                    if not even: witness_examples[cert] = (Cp, Cq, cross, adjF, col)
    tau_E = sum(1 for c, (m, ev) in types_G.items() if ev)
    E_even = sum(m - 1 for c, (m, ev) in types_G.items() if ev)
    W = sum(1 for c, (m, ev) in types_G.items() if not ev)
    E_wit = sum(m - 1 for c, (m, ev) in types_G.items() if not ev)
    defect = tau_T - tau_E
    ident = W + E_even + E_wit - E_T
    row = dict(p=p, q=q, k=k, n=n, R=R, RG=RG, tau_T=tau_T, tau_E=tau_E, defect=defect, W=W, E_T=E_T,
               E_even=E_even, E_wit=E_wit, identity_ok=(defect == ident), theoremR_ok=(R == RG),
               circ_T=(len(tT), len(tQ)), circ_Geven=(len(gP), len(gQ)), colT=len(colT), colG=len(colG),
               secs=round(time.time() - t0, 1))
    if verbose:
        print(" ".join("%s=%s" % (a, b) for a, b in row.items()), flush=True)
    return row, witness_examples

# ---------------------------------------------------------------- brute force cross-check (small k)
def perm_group(gens, n):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = tuple(h[g[i]] for i in range(n))
                if x not in seen: seen.add(x); nxt.append(x)
        frontier = nxt
    return seen

def primitive_root(p):
    for a in range(2, p):
        if len({pow(a, e, p) for e in range(1, p)}) == p - 1: return a
    return 1

def brute(p, q, k):
    """all 2^c H-fixed tournaments and graphs; types directly; N-orbits and N-evenness directly."""
    n = p + q + k
    P = list(range(p)); Q = list(range(p, p + q)); F = list(range(p + q, n))
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    # H-orbits on pairs
    def img(pair, g): return tuple(sorted((g[pair[0]], g[pair[1]])))
    orbits = []; seen = set()
    for pr in pairs:
        if pr in seen: continue
        orb = []; x = pr
        while x not in seen: seen.add(x); orb.append(x); x = img(x, h)
        orbits.append(orb)
    # N generators
    gens = [tuple(h)]
    for (blk, m) in ((P, p), (Q, q)):
        for a in units(m):
            g = list(range(n))
            for i in range(m): g[blk[i]] = blk[(a * i) % m]
            gens.append(tuple(g))
    if k >= 2:
        g = list(range(n)); g[F[0]], g[F[1]] = F[1], F[0]; gens.append(tuple(g))
        g = list(range(n))
        for i in range(k): g[F[i]] = F[(i + 1) % k]
        gens.append(tuple(g))
    N = perm_group(gens, n)
    # tournaments: orientation of each orbit: choose orientation of the first pair, propagate by h
    def orient(orb, bit):
        arcs = []
        (u, v) = orb[0]; src, dst = (u, v) if bit == 0 else (v, u)
        for _ in orb:
            arcs.append((src, dst)); src, dst = h[src], h[dst]
        return arcs
    types_T = {}; fixT = []
    for bits in itertools.product([0, 1], repeat=len(orbits)):
        arcs = [a for orb, b in zip(orbits, bits) for a in orient(orb, b)]
        adj = {v: [] for v in range(n)}
        for (u, v) in arcs: adj[u].append(v)
        arcset = frozenset(arcs)
        fixT.append(arcset)
        types_T[arcset] = pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adj))
    # N-orbits on tournaments
    def act_T(g, arcset): return frozenset((g[u], g[v]) for (u, v) in arcset)
    seenT = set(); RT = 0
    for T in fixT:
        if T in seenT: continue
        RT += 1
        for g in N: seenT.add(act_T(g, T))
    tau_T = len(set(types_T.values()))
    # graphs
    typesG = {}; fixG = []; evenG = {}; NevenG = {}
    for bits in itertools.product([0, 1], repeat=len(orbits)):
        edges = sorted(e for orb, b in zip(orbits, bits) if b for e in orb)
        adj = {v: [] for v in range(n)}
        for (u, v) in edges: adj[u].append(v); adj[v].append(u)
        E = frozenset(edges); fixG.append(E)
        g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
        typesG[E] = pynauty.certificate(g)
        gens_ = pynauty.autgrp(g)[0]
        evenG[E] = all(reversal_parity(gg, edges) == 0 for gg in gens_)
        NevenG[E] = all(reversal_parity(gg, edges) == 0 for gg in N if frozenset(tuple(sorted((gg[u], gg[v]))) for (u, v) in edges) == E)
    def act_G(g, E): return frozenset(tuple(sorted((g[u], g[v]))) for (u, v) in E)
    seenG = set(); RGeven = 0; orbit_reps = []
    for E in fixG:
        if E in seenG: continue
        if NevenG[E]: RGeven += 1; orbit_reps.append(E)
        for g in N: seenG.add(act_G(g, E))
    tau_E = len({typesG[E] for E in fixG if evenG[E]})
    # witnesses: sign-odd types with an N-even H-fixed representative
    W = len({typesG[E] for E in fixG if NevenG[E] and not evenG[E]})
    # sanity: sign-even => N-even
    assert all(NevenG[E] for E in fixG if evenG[E])
    return dict(p=p, q=q, k=k, n=n, c=len(orbits), N=len(N), R=RT, RGeven=RGeven, tau_T=tau_T, tau_E=tau_E,
                defect=tau_T - tau_E, W=W)

if __name__ == "__main__":
    p, q, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    kb = int(sys.argv[sys.argv.index("--brute") + 1]) if "--brute" in sys.argv else -1
    print("C_%d of type (%d,%d,1^k): defect census (pynauty %s)" % (p * q, p, q, pynauty.__version__))
    rows = []
    for k in range(0, kmax + 1):
        row, wit = census(p, q, k)
        rows.append(row)
        if k <= kb:
            b = brute(p, q, k)
            print("  brute:", " ".join("%s=%s" % (a, c) for a, c in b.items()),
                  " AGREE" if (b["R"], b["RGeven"], b["tau_T"], b["tau_E"], b["W"]) == (row["R"], row["RG"], row["tau_T"], row["tau_E"], row["W"]) else " DISAGREE", flush=True)
    print("\nk  n   R  tau_T tau_E defect  W  E_T E_even E_wit")
    for r in rows:
        print("%d %2d %6d %6d %6d %5d %4d %4d %5d %4d" % (r["k"], r["n"], r["R"], r["tau_T"], r["tau_E"], r["defect"], r["W"], r["E_T"], r["E_even"], r["E_wit"]))
