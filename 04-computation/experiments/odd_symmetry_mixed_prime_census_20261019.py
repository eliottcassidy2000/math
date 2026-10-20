"""
odd_symmetry_mixed_prime_census_20261019.py -- S45 (owner: prove HYP-9340 beyond Sylow; think in odd primes and their
powers).  HYP-9340(a) is PROVED for every odd p-group (prime_local_walsh_bijection_20261020) and REFUTED for general odd H
(C15 of orbit type (3,5,1,1) at n = 10: 32 tournament types against 31 sign-even graph types).  This census treats the
MIXED-PRIME odd-order subgroup classes of S_n, n = 8..13, with cyclic or Frobenius transitive constituents (every
odd-order transitive group of degree 3, 5, 7, 11, 13 is C_3, C_5, C_7 or F_21, C_11 or F_55, C_13 or F_39), restricted to
classes whose fixed sets have at most 2^15 objects.  For each class H:
  Fix_T(H), Fix_G(H): the 2^c H-fixed tournaments / graphs (c = number of H-orbits on pairs);
  N = N_(S_n)(H); the N-orbit counts on Fix_T(H) and on the N-even H-fixed graphs (THM-4721 Theorem R predicts equality);
  tau_T(H), tau_E(H): isomorphism types among H-fixed tournaments / sign-even H-fixed graphs (= the types whose
  automorphism group contains a conjugate of H); the defect tau_T - tau_E;
  the witnesses: isomorphism types of H-fixed graphs that are sign-odd but admit an N-even H-fixed representative
  (THM-4721 Proposition Q: these carry the weighted defect); for each witness, whether H is pronormal in the sign kernel
  B(X) of Aut X (the pronormal Frattini lemma of the S45 note says a pronormal H forces N-even => even, so every witness
  must have H non-pronormal in B(X)).
Sign convention (Royle-Praeger-Glasby-Freedman-Devillers): sgn_X(g) = (-1)^#{edges {u<v} of X with g(u) > g(v)}; X is
sign-even iff sgn_X = 1 on Aut X.  Automorphism groups are computed by colour refinement with individualisation and VF2
(networkx), returning Schreier generators of a point-stabiliser chain; isomorphism types by invariant buckets + VF2.
"""
import sys, itertools, math, time
from collections import defaultdict
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, DiGraphMatcher
from sympy.combinatorics import Permutation, PermutationGroup

# ----------------------------------------------------------------------------- groups
def cyc(points):
    return {points[i]: points[(i + 1) % len(points)] for i in range(len(points))}
def perm_from_dict(d, n):
    return tuple(d.get(i, i) for i in range(n))
def compose(p, q):            # (p o q)(i) = p[q[i]]
    return tuple(p[q[i]] for i in range(len(p)))
def inverse(p):
    inv = [0] * len(p)
    for i, pi in enumerate(p): inv[pi] = i
    return tuple(inv)
def frob_mult(points, a):     # x -> a x on Z_p labelled by points
    p = len(points); return {points[i]: points[(a * i) % p] for i in range(p)}

def group_elements(gens, n, cap=10**6):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = compose(h, g)
                if x not in seen:
                    seen.add(x); nxt.append(x)
                    if len(seen) > cap: return None
        frontier = nxt
    return seen

def classes(n):
    """Mixed-prime odd-order subgroup classes of S_n (cyclic / Frobenius constituents), as (name, generators)."""
    out = []
    def C(pts): return perm_from_dict(cyc(pts), n)
    def M(pts, a): return perm_from_dict(frob_mult(pts, a), n)
    def prod(*ds):
        d = {}
        for x in ds: d.update(x)
        return perm_from_dict(d, n)
    r = list(range(n))
    if n >= 8:   # (3,5)
        A, B = r[0:3], r[3:8]; out.append(("C15 (3,5,1^%d)" % (n - 8), [prod(cyc(A), cyc(B))]))
    if n >= 11:  # (3,3,5)
        A, A2, B = r[0:3], r[3:6], r[6:11]
        out.append(("C15 (3,3,5,1^%d) synchronised" % (n - 11), [prod(cyc(A), cyc(A2), cyc(B))]))
        out.append(("C3xC3xC5 (3,3,5,1^%d)" % (n - 11), [C(A), C(A2), C(B)]))
    if n >= 13:  # (3,5,5)
        A, B, B2 = r[0:3], r[3:8], r[8:13]
        out.append(("C15 (3,5,5) synchronised", [prod(cyc(A), cyc(B), cyc(B2))]))
        out.append(("C3xC5xC5 (3,5,5)", [C(A), C(B), C(B2)]))
    if n >= 10:  # (3,7) and F21 variants
        A, S = r[0:3], r[3:10]
        out.append(("C21 (3,7,1^%d)" % (n - 10), [prod(cyc(A), cyc(S))]))
        out.append(("F21 (7,1^%d)" % (n - 7), [C(S), M(S, 2)]))
        out.append(("F21 diag (7,3,1^%d)" % (n - 10), [C(S), prod(frob_mult(S, 2), cyc(A))]))
        out.append(("F21xC3 (7,3,1^%d)" % (n - 10), [C(S), M(S, 2), C(A)]))
    elif n >= 7:
        S = r[0:7]; out.append(("F21 (7,1^%d)" % (n - 7), [C(S), M(S, 2)]))
    if n >= 13:  # (7,3,3)
        S, A, A2 = r[0:7], r[7:10], r[10:13]
        out.append(("C21 (7,3,3) synchronised", [prod(cyc(S), cyc(A), cyc(A2))]))
        out.append(("C7xC3xC3 (7,3,3)", [C(S), C(A), C(A2)]))
        out.append(("F21 diag-diag (7,3,3)", [C(S), prod(frob_mult(S, 2), cyc(A), cyc(A2))]))
        out.append(("F21diag x C3 (7,3,3)", [C(S), prod(frob_mult(S, 2), cyc(A)), C(A2)]))
        # the sixth (7,3,3) class, missed in the first run (found by the subdirect-product enumeration after the
        # Codex census packets arrived): F21 x C3 of order 63 whose C3 complement is the diagonal of the two 3-orbits,
        # so that both orbit-pointwise stabilisers are F21 (in "F21diag x C3" one of them is the cyclic C7 x C3)
        out.append(("F21diag x diagC3 (7,3,3)", [C(S), prod(frob_mult(S, 2), cyc(A)), prod(cyc(A), cyc(A2))]))
        # the seventh (7,3,3) class (Goursat: direct product of F21 with the diagonal C3 of the two 3-orbits; order 63,
        # both orbit-pointwise stabilisers equal to F21 acting on the 7-orbit alone)
        out.append(("F21 x diagC3 (7,3,3)", [C(S), M(S, 2), prod(cyc(A), cyc(A2))]))
        out.append(("F21xC3xC3 (7,3,3)", [C(S), M(S, 2), C(A), C(A2)]))
    if n >= 12:  # (5,7)
        B, S = r[0:5], r[5:12]
        out.append(("C35 (5,7,1^%d)" % (n - 12), [prod(cyc(B), cyc(S))]))
        out.append(("F21xC5 (7,5,1^%d)" % (n - 12), [C(S), M(S, 2), C(B)]))
    if n >= 11:  # F55 on 11
        E = r[0:11]; out.append(("F55 (11,1^%d)" % (n - 11), [C(E), M(E, 3)]))
    if n >= 13:
        E = r[0:13]; out.append(("F39 (13)", [C(E), M(E, 3)]))
    return out

# ----------------------------------------------------------------------------- orbits on pairs, fixed sets
def pair_orbits(gens, n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen = set(); orbs = []
    for pr in pairs:
        if pr in seen: continue
        orb = set(); frontier = [pr]; orb.add(pr)
        while frontier:
            nxt = []
            for (a, b) in frontier:
                for g in gens:
                    x = (g[a], g[b]); x = (min(x), max(x))
                    if x not in orb: orb.add(x); nxt.append(x)
            frontier = nxt
        seen |= orb; orbs.append(sorted(orb))
    return orbs

def orient_orbit(orb, gens, n):
    """A consistent orientation of a pair-orbit under odd H: orient (a,b) as a->b and transport by H.  Returns the set
    of ordered pairs; H odd order never reverses a pair, so the transported orientation is well defined."""
    a, b = orb[0]; arcs = {(a, b)}; frontier = [(a, b)]
    while frontier:
        nxt = []
        for (x, y) in frontier:
            for g in gens:
                z = (g[x], g[y])
                if z not in arcs: arcs.add(z); nxt.append(z)
        frontier = nxt
    assert len(arcs) == len(orb)
    return sorted(arcs)

# ----------------------------------------------------------------------------- automorphisms and signs
def refine(adj, col, n):
    """Colour refinement (1-WL) on an undirected or directed adjacency (adj[v] = (out, in) sets)."""
    while True:
        sig = {}
        for v in range(n):
            o = tuple(sorted(col[u] for u in adj[v][0])); i = tuple(sorted(col[u] for u in adj[v][1]))
            sig[v] = (col[v], o, i)
        keys = sorted(set(sig.values())); idx = {k: j for j, k in enumerate(keys)}
        new = {v: idx[sig[v]] for v in range(n)}
        if len(set(new.values())) == len(set(col.values())): return new
        col = new

def aut_generators(G, directed):
    n = G.number_of_nodes()
    adj = {v: (set(G.successors(v)) if directed else set(G.neighbors(v)), set(G.predecessors(v)) if directed else set(G.neighbors(v))) for v in range(n)}
    Matcher = DiGraphMatcher if directed else GraphMatcher
    gens = []
    def find(colA, colB):
        H1 = G.copy(); H2 = G.copy()
        nx.set_node_attributes(H1, colA, "c"); nx.set_node_attributes(H2, colB, "c")
        gm = Matcher(H1, H2, node_match=lambda a, b: a["c"] == b["c"])
        for iso in gm.isomorphisms_iter():
            return tuple(iso[i] for i in range(n))
        return None
    def search(col):
        cells = defaultdict(list)
        for v in range(n): cells[col[v]].append(v)
        big = [c for c in sorted(cells) if len(cells[c]) >= 2]
        if not big: return
        cell = cells[big[0]]; v = cell[0]
        colv = dict(col); colv[v] = max(col.values()) + 1; colv = refine(adj, colv, n)
        for w in cell[1:]:
            colw = dict(col); colw[w] = max(col.values()) + 1; colw = refine(adj, colw, n)
            iso = find(colv, colw)
            if iso is not None: gens.append(iso)
        search(colv)
    search(refine(adj, {v: 0 for v in range(n)}, n))
    return gens

def sign_on_graph(g, edges):
    return -1 if sum(1 for (u, v) in edges if g[u] > g[v]) % 2 else 1

def invariant(G, directed):
    h = nx.weisfeiler_lehman_graph_hash(G, iterations=3)
    deg = tuple(sorted(dict(G.out_degree()).values())) if directed else tuple(sorted(dict(G.degree()).values()))
    return (G.number_of_edges(), deg, h)

# ----------------------------------------------------------------------------- normaliser
def normaliser(gens, n, orbs):
    """N_(S_n)(H): candidates = product over orbits of the normaliser of the constituent in Sym(orbit) (brute force for
    |orbit| <= 7, AGL(1,p) for the Frobenius/cyclic constituents of degree 11, 13) times permutations of equal-size
    orbits; filtered by normalising H."""
    Hset = group_elements(gens, n)
    orbit_lists = sorted(orbs, key=lambda o: (len(o), o))
    per_orbit = []
    for O in orbit_lists:
        k = len(O)
        if k == 1: per_orbit.append([{O[0]: O[0]}]); continue
        idx = {v: i for i, v in enumerate(O)}
        const = {tuple(idx[h[v]] for v in O) for h in Hset}          # the transitive constituent H|_O on O
        if k <= 7:
            cand = []
            for perm in itertools.permutations(range(k)):             # g: O[i] -> O[perm[i]]
                inv = [0] * k
                for i, pi in enumerate(perm): inv[pi] = i
                ok = True
                for c in const:
                    if tuple(perm[c[inv[i]]] for i in range(k)) not in const: ok = False; break
                if ok: cand.append({O[i]: O[perm[i]] for i in range(k)})
            per_orbit.append(cand)
        else:   # degree 11 or 13: the constituent is C_p or F_(p:d), normaliser AGL(1, p) in the labelling of the orbit
            cand = []
            for a in range(1, k):
                for b in range(k):
                    cand.append({O[i]: O[(a * i + b) % k] for i in range(k)})
            per_orbit.append(cand)
    # permutations of equal-size orbits
    sizes = defaultdict(list)
    for idx, O in enumerate(orbit_lists): sizes[len(O)].append(idx)
    swaps = [dict()]
    for s, idxs in sizes.items():
        if s == 1 or len(idxs) == 1: continue
        new = []
        for sw in swaps:
            for perm in itertools.permutations(idxs):
                d = dict(sw)
                for src, dst in zip(idxs, perm):
                    for i, v in enumerate(orbit_lists[src]): d[v] = orbit_lists[dst][i]
                new.append(d)
        swaps = new
    # fixed points may be permuted arbitrarily: handled as size-1 orbits above? no: include them via a symmetric group
    fixed = [O[0] for O in orbit_lists if len(O) == 1]
    fixed_perms = [dict(zip(fixed, p)) for p in itertools.permutations(fixed)] if fixed else [dict()]
    N = set()
    for choice in itertools.product(*[c for c in per_orbit]):
        base = {}
        for d in choice: base.update(d)
        for sw in swaps:
            for fp in fixed_perms:
                d = dict(base)
                # apply swap after base: v -> sw(base(v)) if base(v) in a swapped orbit
                g = {}
                for v in range(n):
                    u = d.get(v, v); u = sw.get(u, u); u = fp.get(u, u); g[v] = u
                gp = perm_from_dict(g, n)
                if len(set(gp)) != n: continue
                ok = True
                for h in gens:
                    if compose(compose(gp, h), inverse(gp)) not in Hset: ok = False; break
                if ok: N.add(gp)
    return N, Hset

# ----------------------------------------------------------------------------- main census
def census(n, name, gens, cap_c=15):
    orbs = pair_orbits(gens, n); c = len(orbs)
    if c > cap_c:
        print(f"  n={n:2d} {name:32s}: c_2 = {c} > {cap_c}, skipped"); return None
    t0 = time.time()
    N, Hset = normaliser(gens, n, [sorted({v for pr in o for v in pr}) for o in []] or orbits_of_points(gens, n))
    orient = [orient_orbit(o, gens, n) for o in orbs]
    # tournaments
    T_objs = []
    for bits in range(2 ** c):
        arcs = []
        for i, o in enumerate(orient):
            arcs += o if (bits >> i) & 1 == 0 else [(b, a) for (a, b) in o]
        T_objs.append(arcs)
    G_objs = []
    for bits in range(2 ** c):
        edges = []
        for i, o in enumerate(orbs):
            if (bits >> i) & 1: edges += o
        G_objs.append(edges)
    def key_T(arcs): return frozenset(arcs)
    def key_G(edges): return frozenset(edges)
    def act(g, pairs, directed):
        return frozenset((g[a], g[b]) if directed else (min(g[a], g[b]), max(g[a], g[b])) for (a, b) in pairs)
    # isomorphism types
    def iso_types(objs, directed):
        buckets = defaultdict(list); reps = []
        for ob in objs:
            X = nx.DiGraph() if directed else nx.Graph(); X.add_nodes_from(range(n)); X.add_edges_from(ob)
            inv = invariant(X, directed)
            found = None
            for (Y, idx) in buckets[inv]:
                M = (DiGraphMatcher if directed else GraphMatcher)(X, Y)
                if M.is_isomorphic(): found = idx; break
            if found is None:
                found = len(reps); reps.append(ob); buckets[inv].append((X, found))
            yield ob, found
    # graphs: automorphism groups, signs, N-evenness
    info = {}
    for edges in G_objs:
        k = key_G(edges)
        X = nx.Graph(); X.add_nodes_from(range(n)); X.add_edges_from(edges)
        gensX = aut_generators(X, False)
        even = all(sign_on_graph(g, edges) == 1 for g in gensX)
        inN = [g for g in N if act(g, edges, False) == k]
        neven = all(sign_on_graph(g, edges) == 1 for g in inN)
        info[k] = dict(edges=edges, gens=gensX, even=even, neven=neven)
    # N-orbits on graphs (N-evenness is constant on N-orbits) and on tournaments
    orbit_id_G = {}; orbit_neven = []
    for edges in G_objs:
        k = key_G(edges)
        if k in orbit_id_G: continue
        oid = len(orbit_neven); orbit_neven.append(info[k]["neven"])
        for g in N: orbit_id_G[act(g, edges, False)] = oid
    nG_neven = sum(orbit_neven)
    orbit_id_T = {}; nT = 0
    for arcs in T_objs:
        k = key_T(arcs)
        if k in orbit_id_T: continue
        oid = nT; nT += 1
        for g in N: orbit_id_T[act(g, arcs, True)] = oid
    type_T = {}
    for arcs, idx in iso_types(T_objs, True): type_T[key_T(arcs)] = idx
    tauT = len(set(type_T.values()))
    type_of = {}; reps_G = {}
    for edges, idx in iso_types(G_objs, False):
        type_of[key_G(edges)] = idx; reps_G.setdefault(idx, edges)
    even_types = {idx for k, idx in type_of.items() if info[k]["even"]}
    tauE = len(even_types)
    # multiplicities: N-orbits per type (= Aut-classes of conjugates of H), and the N-even ones
    mT = defaultdict(set)
    for k, idx in type_T.items(): mT[idx].add(orbit_id_T[k])
    mG = defaultdict(set); mGev = defaultdict(set)
    for k, idx in type_of.items():
        oid = orbit_id_G[k]; mG[idx].add(oid)
        if orbit_neven[oid]: mGev[idx].add(oid)
    witnesses = sorted(idx for idx in mGev if idx not in even_types)
    W = len(witnesses)
    E_T = sum(len(s) - 1 for s in mT.values())
    E_even = sum(len(mG[idx]) - 1 for idx in even_types)
    E_wit = sum(len(mGev[idx]) - 1 for idx in witnesses)
    ident = (tauT - tauE) == W + E_even + E_wit - E_T
    print(f"  n={n:2d} {name:32s}: |H|={len(Hset):4d} |N|={len(N):5d} c_2={c:2d} |Fix|={2**c:6d} | N-orbits: T {nT:5d}, N-even G {nG_neven:5d} ({'R OK' if nT == nG_neven else 'R FAIL'}) | tau_T={tauT:5d} tau_E={tauE:5d} defect={tauT - tauE:3d} | W={W:3d} E_T={E_T:3d} E_even={E_even:3d} E_wit={E_wit:3d} identity {'OK' if ident else 'FAIL'} | {time.time() - t0:5.1f}s")
    for w in witnesses:
        edges = reps_G[w]; k = key_G(edges); X = nx.Graph(); X.add_nodes_from(range(n)); X.add_edges_from(edges)
        gensX = info[k]["gens"]; Aut = PermutationGroup([Permutation(list(g)) for g in gensX]) if gensX else PermutationGroup([Permutation(list(range(n)))])
        order = Aut.order()
        pron = "skipped (|Aut| large)"
        if order <= 60000:
            elems = [tuple(q.array_form) for q in Aut.generate()]
            Hperm = [Permutation(list(h)) for h in gens]
            pron = "pronormal in Aut X (CONTRADICTS the lemma)"
            # the representative edges may belong to a different H-fixed copy: use an H-fixed representative of this type
            for b in elems:
                Hb = [compose(compose(b, h), inverse(b)) for h in gens]
                if all(x in Hset for x in Hb): continue
                J = PermutationGroup(Hperm + [Permutation(list(x)) for x in Hb])
                if J.order() > 200000: pron = "undecided (large join)"; break
                HbSet = group_elements(Hb, n); conj = False
                for j in J.generate():
                    jt = tuple(j.array_form)
                    if all(compose(compose(jt, h), inverse(jt)) in HbSet for h in gens): conj = True; break
                if not conj: pron = "NOT pronormal in Aut X"; break
        degs = sorted(dict(X.degree()).values())
        print(f"      witness type {w}: {len(edges)} edges, degrees {degs}, |Aut| = {order}, m^ev = {len(mGev[w])}, H {pron}")
    return dict(tauT=tauT, tauE=tauE, W=W)

def orbits_of_points(gens, n):
    seen = set(); orbs = []
    for v in range(n):
        if v in seen: continue
        orb = {v}; frontier = [v]
        while frontier:
            nxt = []
            for x in frontier:
                for g in gens:
                    if g[x] not in orb: orb.add(g[x]); nxt.append(g[x])
            frontier = nxt
        seen |= orb; orbs.append(sorted(orb))
    return orbs

if __name__ == "__main__":
    import sys
    only = sys.argv[1:]   # optional substrings: run only the classes whose name contains one of them
    print("# mixed-prime odd symmetry classes: tau_T(H) against tau_E(H), Theorem R orbit counts, defect witnesses"
          + (" (restricted to %s)" % only if only else ""))
    import os
    cap = int(os.environ.get("CENSUS_CAP", "15"))   # fixed sets of at most 2^cap objects (the main run used 15)
    for n in range(7, 14):
        for name, gens in classes(n):
            if only and not any(o in name for o in only):
                continue
            census(n, name, gens, cap_c=cap)
    print("done")
