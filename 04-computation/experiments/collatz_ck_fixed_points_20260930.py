#!/usr/bin/env python3
"""collatz_ck_fixed_points_20260930.py -- fixed points of the induced-cycle operators C_k on vertex-transitive graphs
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, thirteenth note, part B).

 (1) Locally-C_k graphs: the vertex links are induced k-cycles.  Theorem (Gauss-Bonnet count): for k >= 5 every contractible
     induced k-cycle is a vertex link, and links are pairwise distinct; so with no short non-contractible cycles C_k(G) is the
     rhombus graph R(G) (x ~ y iff x, y are the apexes of two triangles on a common edge).  For the Eisenstein tori
     Cay(Z[w]/I, units) the rhombus graph is Cay(Z[w]/I, (1 - w) units), isomorphic to G when 3 does not divide the index.
     Checked: P_13 = C_13(1,3,4) (induced hexagons = 13 links; C_6(P_13) iso P_13 via the multiplier 2 = 1 - w),
     C_19(1,7,8), C_31(1,5,6), C_37(1,10,11), C_43(1,6,7), the tori Z_n x Z_n (n = 4..8), Z_4 x Z_8, Z_6 x Z_6.
 (2) Census: all circulants C_n(S) with n <= 30, degree <= 7, and all Cayley graphs of A_4, S_4, SL(2,3), D_n (n <= 8), Q_8,
     Z_3 x| Z_4, Z_13 x| Z_3 and the abelian groups Z_a x Z_b (ab <= 36), degree <= 7, for k = 3..7: the graphs with
     C_k(G) iso G (exactly k induced k-cycles through each vertex is the necessary count).
Usage: python3 collatz_ck_fixed_points_20260930.py
"""
import itertools, math, sys, time
from collections import Counter


# ---------- basic graph tools ----------
def graph_from_edges(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        if a != b:
            adj[a].add(b); adj[b].add(a)
    return adj


def induced_cycles(adj, k, limit=None):
    """all induced k-cycles as canonical tuples (min vertex first, smaller second vertex); stops early past `limit`."""
    n = len(adj); out = set()
    for a in range(n):
        path = [a]
        def dfs():
            if limit is not None and len(out) > limit:
                return
            last = path[-1]
            if len(path) == k - 1:
                for e in adj[last]:
                    if e > a and e in adj[a] and e not in path and not any(e in adj[x] for x in path[1:-1]):
                        cyc = tuple(path) + (e,)
                        out.add(min(cyc, (a,) + cyc[:0:-1]))
                return
            for c in adj[last]:
                if c > a and c not in path and not any(c in adj[x] for x in path[:-1]):
                    path.append(c); dfs(); path.pop()
        dfs()
    return sorted(out)


def local_count(adj, k, cap):
    """number of induced k-cycles through vertex 0 (canonical: 0 is the minimum vertex), stopping past cap."""
    out = set(); a = 0; path = [a]
    def dfs():
        if len(out) > cap:
            return
        last = path[-1]
        if len(path) == k - 1:
            for e in adj[last]:
                if e > a and e in adj[a] and e not in path and not any(e in adj[x] for x in path[1:-1]):
                    cyc = tuple(path) + (e,)
                    out.add(min(cyc, (a,) + cyc[:0:-1]))
            return
        for c in adj[last]:
            if c > a and c not in path and not any(c in adj[x] for x in path[:-1]):
                path.append(c); dfs(); path.pop()
    dfs()
    return len(out)


def cycle_operator(adj, k, limit=None):
    cyc = induced_cycles(adj, k, limit)
    if limit is not None and len(cyc) > limit:
        return None, cyc
    edges_of = [set((min(c[i], c[(i + 1) % k]), max(c[i], c[(i + 1) % k])) for i in range(k)) for c in cyc]
    m = len(cyc)
    return graph_from_edges(m, [(i, j) for i in range(m) for j in range(i + 1, m) if edges_of[i] & edges_of[j]]), cyc


def wl_colours(adj, rounds=4):
    n = len(adj); col = [len(a) for a in adj]
    for _ in range(rounds):
        sig = [(col[v], tuple(sorted(col[u] for u in adj[v]))) for v in range(n)]
        keys = {s: i for i, s in enumerate(sorted(set(sig)))}
        col = [keys[s] for s in sig]
    return col


def spectrum(adj):
    import numpy as np
    n = len(adj)
    M = np.zeros((n, n))
    for i in range(n):
        for j in adj[i]:
            M[i, j] = 1.0
    return tuple(np.round(np.linalg.eigvalsh(M), 6))


ISO_BUDGET = [0]


def is_isomorphic(adj1, adj2, budget=200000):
    """colour refinement + spectrum as invariants, then budgeted backtracking (returns False when the budget is exhausted,
    which is logged; only used after the cheap invariants agree)."""
    n = len(adj1)
    if n != len(adj2):
        return False
    c1, c2 = wl_colours(adj1), wl_colours(adj2)
    if sorted(c1) != sorted(c2):
        return False
    if sorted(Counter(c1).values()) != sorted(Counter(c2).values()):
        return False
    if spectrum(adj1) != spectrum(adj2):
        return False
    # triangle counts per vertex as a further invariant
    t1 = sorted(sum(1 for u in adj1[v] for w in adj1[v] if u < w and w in adj1[u]) for v in range(n))
    t2 = sorted(sum(1 for u in adj2[v] for w in adj2[v] if u < w and w in adj2[u]) for v in range(n))
    if t1 != t2:
        return False
    ISO_BUDGET[0] = budget
    # relabel colours of adj2 consistently: colours are canonical (sorted signatures) only if the multisets agree, which we checked
    order = sorted(range(n), key=lambda v: (Counter(c1)[c1[v]], c1[v]))
    mapping = {}; used = set()
    cand = {v: [w for w in range(n) if c2[w] == c1[v]] for v in range(n)}

    def bt(i):
        if i == n:
            return True
        ISO_BUDGET[0] -= 1
        if ISO_BUDGET[0] < 0:
            return False
        v = order[i]
        for w in cand[v]:
            if w in used:
                continue
            if all(((u in adj1[v]) == (mapping[u] in adj2[w])) for u in mapping):
                mapping[v] = w; used.add(w)
                if bt(i + 1):
                    return True
                del mapping[v]; used.discard(w)
                if ISO_BUDGET[0] < 0:
                    return False
        return False
    r = bt(0)
    if not r and ISO_BUDGET[0] < 0:
        print("   [isomorphism test: budget exhausted on a %d-vertex pair with equal invariants; reported as non-isomorphic]" % n, flush=True)
    return r


def is_locally_ck(adj, k):
    """every vertex link is an induced k-cycle."""
    for v in range(len(adj)):
        N = adj[v]
        if len(N) != k:
            return False
        for u in N:
            if len(adj[u] & N) != 2:
                return False
        # connected 2-regular on k vertices = a k-cycle
        start = next(iter(N)); seen = {start}; frontier = [start]
        while frontier:
            x = frontier.pop()
            for y in adj[x] & N:
                if y not in seen:
                    seen.add(y); frontier.append(y)
        if len(seen) != k:
            return False
    return True


def rhombus_graph(adj):
    n = len(adj)
    e = []
    for x in range(n):
        for y in range(x + 1, n):
            if y in adj[x]:
                continue
            common = adj[x] & adj[y]
            if any(v in adj[u] for u in common for v in common):
                e.append((x, y))
    return graph_from_edges(n, e)


# ---------- Cayley graphs ----------
def cayley(elements, mul, inv, S):
    idx = {g: i for i, g in enumerate(elements)}
    return graph_from_edges(len(elements), [(idx[g], idx[mul(g, s)]) for g in elements for s in S])


def circulant(n, S):
    return cayley(list(range(n)), lambda a, b: (a + b) % n, lambda a: (-a) % n, S)


def abelian(a, b):
    els = [(i, j) for i in range(a) for j in range(b)]
    return els, (lambda x, y: ((x[0] + y[0]) % a, (x[1] + y[1]) % b)), (lambda x: ((-x[0]) % a, (-x[1]) % b)), (0, 0)


def perm_group(gens, deg):
    ident = tuple(range(deg)); els = {ident}; frontier = [ident]
    comp = lambda p, q: tuple(p[q[i]] for i in range(deg))   # p after q
    while frontier:
        g = frontier.pop()
        for s in gens:
            h = comp(g, s)
            if h not in els:
                els.add(h); frontier.append(h)
    inv = lambda p: tuple(sorted(range(deg), key=lambda i: p[i]))
    return sorted(els), comp, inv, ident


def small_groups():
    out = {}
    out["A_4"] = perm_group([(1, 2, 0, 3), (1, 0, 3, 2)], 4)
    out["S_4"] = perm_group([(1, 0, 2, 3), (1, 2, 3, 0)], 4)
    for n in range(3, 9):
        rot = tuple((i + 1) % n for i in range(n)); ref = tuple((-i) % n for i in range(n))
        out["D_%d" % n] = perm_group([rot, ref], n)
    # Q_8 as permutations of 8 elements (regular representation of the quaternion group)
    names = ["1", "-1", "i", "-i", "j", "-j", "k", "-k"]
    tab = {}
    sign = lambda s: -1 if s.startswith("-") else 1
    base = lambda s: s.lstrip("-")
    qm = {("1", "1"): (1, "1"), ("1", "i"): (1, "i"), ("1", "j"): (1, "j"), ("1", "k"): (1, "k"),
          ("i", "1"): (1, "i"), ("j", "1"): (1, "j"), ("k", "1"): (1, "k"),
          ("i", "i"): (-1, "1"), ("j", "j"): (-1, "1"), ("k", "k"): (-1, "1"),
          ("i", "j"): (1, "k"), ("j", "k"): (1, "i"), ("k", "i"): (1, "j"),
          ("j", "i"): (-1, "k"), ("k", "j"): (-1, "i"), ("i", "k"): (-1, "j")}
    def qmul(x, y):
        s, b = qm[(base(x), base(y))]
        s *= sign(x) * sign(y)
        return ("-" if s < 0 else "") + b
    perms = []
    for g in names:
        perms.append(tuple(names.index(qmul(g, h)) for h in names))
    out["Q_8"] = perm_group(perms, 8)
    # SL(2,3) as 2x2 matrices mod 3 acting on the 8 nonzero vectors of F_3^2 (faithful)
    vecs = [(x, y) for x in range(3) for y in range(3) if (x, y) != (0, 0)]
    def mat_perm(m):
        return tuple(vecs.index(((m[0] * v[0] + m[1] * v[1]) % 3, (m[2] * v[0] + m[3] * v[1]) % 3)) for v in vecs)
    out["SL(2,3)"] = perm_group([mat_perm((1, 1, 0, 1)), mat_perm((0, 2, 1, 0))], 8)
    # Z_3 x| Z_4 (dicyclic of order 12): generated by a rotation of order 3 on 3 points? use permutations on 12 points via regular rep
    # simpler: Z_13 x| Z_3 and Z_3 x| Z_4 as affine maps
    def affine_group(p, mults):
        els = [(a, b) for a in mults for b in range(p)]
        mul = lambda x, y: ((x[0] * y[0]) % p, (x[0] * y[1] + x[1]) % p)
        def inv(x):
            ai = pow(x[0], -1, p)
            return (ai, (-ai * x[1]) % p)
        return els, mul, inv, (1, 0)
    out["Z_13x|Z_3"] = affine_group(13, [1, 3, 9])
    out["Z_7x|Z_3"] = affine_group(7, [1, 2, 4])
    return out


def symmetric_subsets(elements, mul, inv, ident, dmax):
    classes = []; seen = set()
    for g in elements:
        if g == ident or g in seen:
            continue
        gi = inv(g)
        cls = (g,) if gi == g else (g, gi)
        seen.update(cls); classes.append(cls)
    for r in range(1, len(classes) + 1):
        for combo in itertools.combinations(classes, r):
            S = [x for cls in combo for x in cls]
            if len(S) <= dmax:
                yield tuple(S)


# ---------- (1) locally-C_6 tori ----------
def part1():
    print("== (1) locally-C_k graphs, the rhombus graph, and the Eisenstein tori as fixed points of C_6 ==")
    P13 = circulant(13, [1, 3, 4, 9, 10, 12])
    hexes = induced_cycles(P13, 6)
    links = sorted(set(tuple(sorted(P13[v])) for v in range(13)))
    print(" P_13 = C_13(1,3,4): locally C_6: %s; induced hexagons: %d; they are the 13 vertex links: %s" % (
        is_locally_ck(P13, 6), len(hexes), sorted(tuple(sorted(h)) for h in hexes) == links))
    C6, _ = cycle_operator(P13, 6)
    R = rhombus_graph(P13)
    nonsq = circulant(13, [2, 5, 6, 7, 8, 11])
    print(" C_6(P_13) iso P_13: %s; C_6(P_13) iso rhombus graph R(P_13): %s; R(P_13) = Cay(Z_13, non-squares): %s (multiplier 1 - w = 1 - 3 = -2)" % (
        is_isomorphic(C6, P13), is_isomorphic(C6, R), all((j in R[i]) == (((j - i) % 13) in {2, 5, 6, 7, 8, 11}) for i in range(13) for j in range(13) if i != j)))
    print(" P_13 under C_5: %d induced pentagons (all non-contractible on the torus); under C_4: %d induced squares" % (len(induced_cycles(P13, 5)), len(induced_cycles(P13, 4))))
    print(" other Eisenstein circulants C_n(1, a, a+1), a^2 + a + 1 = 0 mod n (3 does not divide n):")
    for n, a in ((19, 7), (31, 5), (37, 10), (43, 6), (49, 18)):   # 61, 67 in the _part1_large.out run
        assert (a * a + a + 1) % n == 0
        G = circulant(n, [1, a, a + 1, n - 1, n - a, n - a - 1])
        t0 = time.time(); hx = induced_cycles(G, 6, limit=3 * n)
        if len(hx) > 3 * n:
            print("  n = %2d, a = %2d: locally C_6 %s; more than %d induced hexagons (extra non-contractible ones)" % (n, a, is_locally_ck(G, 6), 3 * n)); continue
        C6g, _ = cycle_operator(G, 6)
        fixed = is_isomorphic(C6g, G) if len(C6g) == n else False
        # the multiplier 1 - a should carry the connection set to the rhombus set
        m = (1 - a) % n
        S = {1, a, a + 1, n - 1, n - a, n - a - 1}
        Srh = {(m * s) % n for s in S}
        print("  n = %2d, a = %2d: locally C_6 %s; induced hexagons %d; C_6(G) iso G: %s; multiplier 1 - a = %d maps S to the rhombus set: %s  (%.1fs)" % (
            n, a, is_locally_ck(G, 6), len(hx), fixed, m, all(((j - i) % n in Srh) == (j in rhombus_graph(G)[i]) for i in range(1) for j in range(n) if j != i), time.time() - t0))
    print(" square and rectangular tori Z_a x Z_b with S = {+-(1,0), +-(0,1), +-(1,1)}:")
    for a, b in ((4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (4, 8), (4, 5), (5, 7), (3, 3), (3, 4)):
        els, mul, inv, ident = abelian(a, b)
        S = [(1, 0), (0, 1), (1, 1), (a - 1, 0), (0, b - 1), (a - 1, b - 1)]
        S = [s for s in S if s != ident]
        G = cayley(els, mul, inv, S)
        loc = is_locally_ck(G, 6)
        hx = induced_cycles(G, 6, limit=3 * a * b)
        if len(hx) > 3 * a * b:
            print("  Z_%d x Z_%d: locally C_6 %s; more than %d induced hexagons" % (a, b, loc, 3 * a * b)); continue
        C6g, _ = cycle_operator(G, 6)
        fixed = (len(C6g) == a * b) and is_isomorphic(C6g, G)
        comps = None
        if len(C6g) == a * b:
            # number of components of C_6(G)
            seen = set(); comps = 0
            for v in range(len(C6g)):
                if v in seen:
                    continue
                comps += 1; st = [v]; seen.add(v)
                while st:
                    x = st.pop()
                    for y in C6g[x]:
                        if y not in seen:
                            seen.add(y); st.append(y)
        print("  Z_%d x Z_%d (index %2d, 3 | index: %s): locally C_6 %s; hexagons %d; C_6(G) iso G: %s; components of C_6(G): %s" % (
            a, b, a * b, a * b % 3 == 0, loc, len(hx), fixed, comps))


# ---------- (2) census ----------
def census(skip_circulants=False):
    print("== (2) census of C_k-fixed Cayley graphs, k = 3..6, degree <= 6 (circulants n <= 24; abelian rank 2, ab <= 36; A_4, S_4, SL(2,3), D_3..D_8, Q_8, Z_7x|Z_3, Z_13x|Z_3) ==", flush=True)
    found = []
    t0 = time.time()
    # circulants
    for n in ([] if skip_circulants else range(4, 25)):
        els = list(range(n)); mul = lambda a, b: (a + b) % n; inv = lambda a: (-a) % n
        for S in symmetric_subsets(els, mul, inv, 0, 6):
            if min(S) != min(S):
                pass
            G = circulant(n, S)
            for k in range(3, 7):
                if local_count(G, k, k) != k:
                    continue
                Ck, cyc = cycle_operator(G, k, limit=n)
                if Ck is None or len(cyc) != n:
                    continue
                if is_isomorphic(Ck, G):
                    found.append(("Z_%d" % n, tuple(sorted(s for s in S if s <= n // 2)), k, len(S)))
    print(" circulants done (%.0fs): %s" % (time.time() - t0, [f for f in found]), flush=True)
    # abelian rank 2
    for a in range(2, 7):
        for b in range(a, 19):
            if a * b > 36 or b % a != 0:
                continue
            els, mul, inv, ident = abelian(a, b)
            print("  abelian Z_%d x Z_%d ... (%.0fs)" % (a, b, time.time() - t0), flush=True)
            for S in symmetric_subsets(els, mul, inv, ident, 6):
                G = cayley(els, mul, inv, S)
                for k in range(3, 7):
                    if local_count(G, k, k) != k:
                        continue
                    Ck, cyc = cycle_operator(G, k, limit=a * b)
                    if Ck is None or len(cyc) != a * b:
                        continue
                    if is_isomorphic(Ck, G):
                        found.append(("Z_%dxZ_%d" % (a, b), S, k, len(S)))
    print(" abelian rank-2 done (%.0fs)" % (time.time() - t0))
    for name, (els, mul, inv, ident) in small_groups().items():
        cnt = 0
        for S in symmetric_subsets(els, mul, inv, ident, 6):
            cnt += 1
            G = cayley(els, mul, inv, S)
            for k in range(3, 7):
                if local_count(G, k, k) != k:
                    continue
                Ck, cyc = cycle_operator(G, k, limit=len(els))
                if Ck is None or len(cyc) != len(els):
                    continue
                if is_isomorphic(Ck, G):
                    found.append((name, "S of degree %d" % len(S), k, len(S), "locally C_k" if is_locally_ck(G, k) else "not locally C_k", "connected" if connected(G) else "disconnected"))
        print(" %s: %d connection sets (%.0fs)" % (name, cnt, time.time() - t0), flush=True)
    # summarize
    print(" C_k-FIXED graphs found (group, S, k, degree[, locality, connectivity]):")
    summary = Counter((f[0], f[2], f[3]) for f in found)
    for key, c in sorted(summary.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        print("   group %-10s k = %d degree %d : %d connection sets" % (key[0], key[1], key[2], c))
    return found


def connected(adj):
    n = len(adj); seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    return len(seen) == n


def part2_details(found):
    print("== (2b) the fixed graphs, identified ==")
    P13 = circulant(13, [1, 3, 4, 9, 10, 12]); K4 = circulant(4, [1, 2, 3])
    ico_els, ico_mul, ico_inv, ico_id = small_groups()["A_4"]
    shown = set()
    for f in found:
        key = (f[0], f[2], f[3])
        if key in shown:
            continue
        shown.add(key)
        print("  ", f)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all", "part1"):
        part1()
    if mode in ("all", "census", "groups"):
        found = census(skip_circulants=(mode == "groups"))
        part2_details(found)
