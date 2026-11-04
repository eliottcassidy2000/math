#!/usr/bin/env python3
"""audit_B_steps_structure.py -- audit B of THM-4871 statement 7: the objects of Steps 1-5 on actual graphs.

For every N-even N-orbit representative X = (C_P, C_Q, cross, gamma) of the rows (p,q,k):
  Step 1: M = Q if C_Q != 0, else the false-twin class of a Q-vertex (= the largest independent module containing Q);
          check M is an independent module (or Q), and every nauty generator of Aut X preserves M.
  Step 2: D = P if C_P != 0 (check P n M = 0, generators preserve P); else D = M if P <= M, else the false-twin class
          of a P-vertex (check D n M = 0, generators preserve D).
  Step 3: F' = F \ (M u D); for every generator g: sgn_X(g) = +1; the factorisation
          sgn_X(g) = sgn_{X[M]}(g|M) sgn_{X[D]}(g|D) sgn_{X[F']}(g|F') with the first two factors +1; the extension g^
          (g on F', identity on F n (M u D)) is a colour-preserving automorphism of gamma = (X[F], (x_f, y_f)), and
          sgn_gamma(g^) = sgn_{X[F']}(g|F').
  Step 5 (k = q-p): the reps other than K_{q,q} and the empty graph have all generators with g(Q) n Q != 0 (sample);
          K_{q,q} and the empty graph are the only reps with a generator moving Q off itself.
  Step 5 gap: the graphs X_1 = (C_P = 0, C_Q a multiplier-even prime circulant, cross, F independent of colour
          (0, cross)) at k = q-p: a q-cycle on P u F IS an automorphism (so the written argument 'with c'' in place of
          g c_Q g^-1' starts from a real configuration), X_1 is neither K_{q,q} nor empty, X[Q] is NOT empty -- so the
          sentence 'by (i) X[Q] is empty' does not follow from c'' alone; what excludes a conjugate of H with q-orbit
          P u F is the p-cycle c' inside Q: Aut X[Q] <= AGL(1,q) contains no p-cycle (checked by enumerating Aut X[Q]).
  Datum uniqueness: by brute force over all 2^c H-fixed graphs at (3,5) k=2 and (5,7) k=2, exactly one is
          isomorphic to K_{q,q}.
"""
import sys, time, itertools, math
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from auditb_common import *

def neighbourhoods(n, edges):
    nb = [set() for _ in range(n)]
    for (u, v) in edges:
        nb[u].add(v); nb[v].add(u)
    return nb

def is_module(S, n, nb):
    S = set(S)
    for v in range(n):
        if v in S: continue
        x = nb[v] & S
        if x and x != S: return False
    return True

def restricted_sign(g, part, edges):
    """sign of g restricted to the induced subgraph on part (g must preserve part)."""
    part = set(part)
    sub = [(u, v) for (u, v) in edges if u in part and v in part]
    return reversal_parity(g, sub)

def check_row(p, q, k):
    n = p + q + k
    P = set(range(p)); Q = set(range(p, p + q)); F = set(range(p + q, n))
    symP = [C for C in dilation_classes(symmetric_sets(p), p) if multiplier_even(C, p)]
    symQ = [C for C in dilation_classes(symmetric_sets(q), q) if multiplier_even(C, q)]
    colG = [(E, col, ev) for (E, col, ev) in coloured_objects(k, False) if ev]
    stats = Counter(); problems = []
    movers = []; nreps = 0
    for CP in symP:
        for CQ in symQ:
            for cross in (0, 1):
                for (EF, col, _) in colG:
                    edges = build_graph(p, q, k, CP, CQ, cross, EF, col)
                    nb = neighbourhoods(n, edges)
                    cert, gens, order = graph_cert_and_gens(n, edges)
                    gens = [list(g) for g in gens]; nreps += 1
                    # Step 1
                    if CQ:
                        M = set(Q); stats["M=Q"] += 1
                    else:
                        u = p
                        M = {v for v in range(n) if nb[v] == nb[u]}
                        assert Q <= M
                        stats["M independent"] += 1
                    if not is_module(M, n, nb): problems.append((p, q, k, "M not a module", CP, CQ, cross, EF, col))
                    if not CQ and any((u in M and v in M) for (u, v) in edges): problems.append((p, q, k, "M not independent", CP, CQ, cross, EF, col))
                    for g in gens:
                        if {g[v] for v in M} != M: problems.append((p, q, k, "generator moves M", CP, CQ, cross, EF, col)); break
                    # Step 2
                    if CP:
                        D = set(P); stats["D=P"] += 1
                        if P & M: problems.append((p, q, k, "P meets M with C_P prime", CP, CQ, cross, EF, col))
                    elif P <= M:
                        D = M; stats["D=M (P<=M)"] += 1
                    else:
                        D = {v for v in range(n) if nb[v] == nb[0]}
                        assert P <= D
                        stats["D independent, disjoint from M"] += 1
                        if D & M: problems.append((p, q, k, "D meets M", CP, CQ, cross, EF, col))
                    if not is_module(D, n, nb): problems.append((p, q, k, "D not a module", CP, CQ, cross, EF, col))
                    for g in gens:
                        if {g[v] for v in D} != D: problems.append((p, q, k, "generator moves D", CP, CQ, cross, EF, col)); break
                    # Step 3
                    Fp = F - (M | D)
                    colour = {f: ((1 if P <= nb[f] else 0), (1 if Q <= nb[f] else 0)) for f in F}
                    for f in F:   # H-fixed: uniform attachment
                        assert (nb[f] & P) in (set(), P) and (nb[f] & Q) in (set(), Q)
                    edgesF = [(u, v) for (u, v) in edges if u in F and v in F]
                    for g in gens:
                        s = reversal_parity(g, edges)
                        if s: problems.append((p, q, k, "ODD automorphism", CP, CQ, cross, EF, col)); break
                        sM = restricted_sign(g, M, edges); sD = restricted_sign(g, D, edges) if D is not M else 0
                        sF = restricted_sign(g, Fp, edges)
                        if (sM + sD + sF) & 1 != s: problems.append((p, q, k, "factorisation fails", CP, CQ, cross, EF, col))
                        if sM or sD: problems.append((p, q, k, "block sign on M or D odd", CP, CQ, cross, EF, col))
                        ghat = {f: (g[f] if f in Fp else f) for f in F}
                        # colour-preserving automorphism of gamma?
                        if any(colour[ghat[f]] != colour[f] for f in F): problems.append((p, q, k, "g^ not colour-preserving", CP, CQ, cross, EF, col))
                        EF_img = {(min(ghat[u], ghat[v]), max(ghat[u], ghat[v])) for (u, v) in edgesF}
                        if EF_img != set(edgesF): problems.append((p, q, k, "g^ not an automorphism of X[F]", CP, CQ, cross, EF, col))
                        gh = list(range(n))
                        for f in F: gh[f] = ghat[f]
                        if reversal_parity(gh, edgesF) != sF: problems.append((p, q, k, "sgn_gamma(g^) != sgn_{X[F']}(g|F')", CP, CQ, cross, EF, col))
                    # Step 5 sample: generators moving Q off itself
                    for g in gens:
                        if not ({g[v] for v in Q} & Q):
                            movers.append((sorted(CP), sorted(CQ), cross, sorted(EF), col, len(edges)))
                            break
    print("(%d,%d,1^%d): %d N-even reps; M-cases / D-cases %s; problems: %d; reps with a generator g, g(Q) n Q = 0: %s"
          % (p, q, k, nreps, dict(stats), len(problems), movers if len(movers) <= 4 else len(movers)), flush=True)
    for pr in problems[:10]:
        print("   PROBLEM:", pr)
    return problems

def step5_gap(p, q):
    """the configurations X_1 at k = q - p."""
    k = q - p; n = p + q + k
    P = list(range(p)); Q = list(range(p, p + q)); F = list(range(p + q, n))
    symQ = [C for C in dilation_classes(symmetric_sets(q), q) if multiplier_even(C, q) and len(C) > 0]
    for CQ in symQ:
        for cross in (0, 1):
            col = tuple(2 * cross for _ in range(k))   # colour (x, y) = (0, cross)
            edges = build_graph(p, q, k, frozenset(), CQ, cross, frozenset(), col)
            cert, gens, order = graph_cert_and_gens(n, edges)
            # c'': a q-cycle on P u F fixing Q pointwise
            PF = P + F
            cpp = list(range(n))
            for i in range(q): cpp[PF[i]] = PF[(i + 1) % q]
            E = set(edges)
            is_aut = all(((min(cpp[u], cpp[v]), max(cpp[u], cpp[v])) in E) for (u, v) in edges)
            # Aut X[Q]: enumerate AGL(1,q) elements preserving C_Q and look for p-cycles; also brute force all of Aut X[Q] by closure of nauty gens
            edgesQ = circulant_edges(CQ, q)
            certQ, gensQ, orderQ = graph_cert_and_gens(q, edgesQ)
            # closure of Aut X[Q]
            elems = {tuple(range(q))}; frontier = [tuple(range(q))]
            gensQ = [tuple(g) for g in gensQ]
            while frontier:
                nxt = []
                for g in frontier:
                    for s in gensQ:
                        x = tuple(s[g[i]] for i in range(q))
                        if x not in elems: elems.add(x); nxt.append(x)
                frontier = nxt
            def cycle_type(g):
                seen = [False] * len(g); ct = []
                for i in range(len(g)):
                    if not seen[i]:
                        j = i; L = 0
                        while not seen[j]: seen[j] = True; j = g[j]; L += 1
                        ct.append(L)
                return tuple(sorted(ct, reverse=True))
            p_cycles = [g for g in elems if cycle_type(g) == tuple([p] + [1] * (q - p))]
            max_fix = max((sum(1 for i in range(q) if g[i] == i) for g in elems if g != tuple(range(q))), default=0)
            print("   (%d,%d) k=%d: X_1 with C_Q=%s cross=%d: |Aut X_1|=%.6g = |Aut X[Q]| * q! ? %s; q-cycle on P u F is an automorphism: %s; sign-even: %s; "
                  "|Aut X[Q]|=%d, p-cycles in Aut X[Q]: %d, max fixed points of a non-identity element of Aut X[Q]: %d"
                  % (p, q, k, sorted(CQ), cross, order, abs(order - len(elems) * math.factorial(q)) <= 1e-9 * order, is_aut, sign_even(gens, edges),
                     len(elems), len(p_cycles), max_fix))

def datum_uniqueness(p, q):
    k = q - p; n = p + q + k
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    orbits = []; seen = set()
    for pr in pairs:
        if pr in seen: continue
        orb = []; x = pr
        while x not in seen:
            seen.add(x); orb.append(x); x = tuple(sorted((h[x[0]], h[x[1]])))
        orbits.append(orb)
    c = len(orbits)
    Kqq = build_graph(p, q, k, frozenset(), frozenset(), 1, frozenset(), tuple(2 for _ in range(k)))
    certK = graph_cert_and_gens(n, Kqq)[0]
    count = 0; found = []
    for X in range(1 << c):
        edges = sorted(e for o in range(c) if (X >> o) & 1 for e in orbits[o])
        if len(edges) != q * q: continue
        if graph_cert_and_gens(n, edges)[0] == certK:
            count += 1; found.append(X)
    print("   (%d,%d) k=%d: H-fixed graphs isomorphic to K_{q,q}: %d (of 2^%d)" % (p, q, k, count, c))

if __name__ == "__main__":
    t0 = time.time()
    print("audit B: Steps 1-5 of statement 7 on the N-even representatives")
    allp = []
    for (p, q, kmax) in ((3, 5, 2), (5, 7, 2), (7, 11, 4), (11, 13, 2), (3, 7, 2), (5, 11, 4), (13, 17, 2)):
        for k in range(kmax + 1):
            allp += check_row(p, q, k)
    print("total problems: %d  (%.1fs)" % (len(allp), time.time() - t0))
    print("\nStep 5 gap: the configurations X_1 (X[Q] a prime circulant, P u F independent) at k = q - p")
    for (p, q) in ((3, 5), (5, 7), (7, 11), (11, 13)):
        step5_gap(p, q)
    print("\nDatum uniqueness of K_{q,q}:")
    for (p, q) in ((3, 5), (5, 7)):
        datum_uniqueness(p, q)
    print("done (%.1fs)" % (time.time() - t0))
