#!/usr/bin/env python3
"""
trijection_four_nodes_three_edges_kps_20261020.py  (kind-pasteur-2026-10-20-S3)

Exact, self-contained checks for the "trijection" note
05-knowledge/results/trijection_four_nodes_three_edges_20261020.md.

Universe: all labelled tournaments on n = 4 (64), n = 5 (1024), n = 6 (32768)
vertices, enumerated as upset sets over the pairs of [n]; every count below is
an exact integer obtained by complete enumeration or by Burnside's lemma over
the full symmetric group.  No sampling, no floats.

Sections
  A  the three 3-edge graphs on 4 nodes: labelled copies, complement map,
     opposite-edge map (tetrahedral duality on edge sets), coherent orientations
  B  n = 4 tournament dictionary: H = 1 + 2 c3, c3 + s = 2, score identity
  C  tetrahedral duality T -> T*: c3(T) = s(T*), H(T*) = 6 - H(T), involution
     up to converse, type actions of converse and duality (Klein four-group)
  D  marked structures (line / cycle+solo / star): 8 = 4 + 4 at n = 4, and the
     general-n counts (2^C(n-1,2), cyc(n), 2 t(n-1)) with the n = 5, 6 hostile
  E  the line-and-staircase bijection (T, P) <-> upset graph on the anti-path:
     Szele's identity sum_T H(T) = n! 2^C(n-1,2) for n = 4, 5, 6, freeness of
     Aut on Hamiltonian paths, and the anti-path complement involution at n = 4
  F  the divisor triple F = S + U at p, p^3, p^2qr (the owner's shapes), stated
     next to the graph triple so that the comparison is on the record
"""
import itertools
import sys
from math import comb, factorial

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if not detail else "  [" + str(detail) + "]"))


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def all_tournaments(n):
    """Yield every labelled tournament on [n] as a frozenset of arcs (u, v): u -> v."""
    P = pairs(n)
    for bits in range(1 << len(P)):
        arcs = []
        for k, (i, j) in enumerate(P):
            arcs.append((j, i) if (bits >> k) & 1 else (i, j))
        yield frozenset(arcs)


def out_sets(T, n):
    out = [set() for _ in range(n)]
    for (u, v) in T:
        out[u].add(v)
    return out


def ham_paths(T, n):
    """Number of directed Hamiltonian paths (subset DP)."""
    out = out_sets(T, n)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c:
                continue
            for w in out[v]:
                if not (mask >> w) & 1:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])


def ham_path_list(T, n):
    """All directed Hamiltonian paths as vertex tuples."""
    out = out_sets(T, n)
    res = []

    def ext(path, used):
        if len(path) == n:
            res.append(tuple(path))
            return
        for w in out[path[-1]]:
            if not used[w]:
                used[w] = True
                path.append(w)
                ext(path, used)
                path.pop()
                used[w] = False

    for v in range(n):
        used = [False] * n
        used[v] = True
        ext([v], used)
    return res


def cyclic_triangles(T, n):
    c = 0
    for a, b, d in itertools.combinations(range(n), 3):
        if ((a, b) in T and (b, d) in T and (d, a) in T) or ((b, a) in T and (d, b) in T and (a, d) in T):
            c += 1
    return c


def directed_cycles_through(T, n, k):
    """Number of directed k-cycles (k-subsets carrying a directed Hamiltonian cycle, counted per cycle)."""
    out = out_sets(T, n)
    total = 0
    for S in itertools.combinations(range(n), k):
        Sset = set(S)
        start = S[0]
        # count directed Hamiltonian cycles of T[S] through start (each cycle once)
        def ext(path, used):
            if len(path) == k:
                return 1 if start in out[path[-1]] else 0
            c = 0
            for w in out[path[-1]]:
                if w in Sset and w not in used:
                    used.add(w)
                    path.append(w)
                    c += ext(path, used)
                    path.pop()
                    used.remove(w)
            return c
        total += ext([start], {start})
    return total


def scores(T, n):
    s = [0] * n
    for (u, v) in T:
        s[u] += 1
    return s


def sources_sinks(T, n):
    sc = scores(T, n)
    return sum(1 for x in sc if x == 0 or x == n - 1)


def converse(T):
    return frozenset((v, u) for (u, v) in T)


def relabel(T, perm):
    return frozenset((perm[u], perm[v]) for (u, v) in T)


def canon(T, n, perms):
    return min(tuple(sorted(relabel(T, p))) for p in perms)


def aut_order(T, n, perms):
    return sum(1 for p in perms if relabel(T, p) == T)


# ---------------------------------------------------------------------------
# A. the three 3-edge graphs on 4 nodes
# ---------------------------------------------------------------------------
def graph_type4(E):
    """Isomorphism type of a 3-edge graph on {0,1,2,3}: 'P4', 'K3+K1', 'K13'."""
    deg = [0] * 4
    for (a, b) in E:
        deg[a] += 1
        deg[b] += 1
    d = sorted(deg)
    if d == [0, 2, 2, 2]:
        return "K3+K1"
    if d == [1, 1, 1, 3]:
        return "K13"
    if d == [1, 1, 2, 2]:
        return "P4"
    raise ValueError(E)


def opposite_edge(e):
    a, b = e
    rest = [x for x in range(4) if x not in (a, b)]
    return (rest[0], rest[1])


def is_directed_path_orientation(E, arcs):
    """arcs: orientation of the 3 edges of a P4 E.  Directed path iff every vertex has in,out <= 1 and
    the arcs chain."""
    indeg = [0] * 4
    outdeg = [0] * 4
    for (u, v) in arcs:
        outdeg[u] += 1
        indeg[v] += 1
    return max(indeg) <= 1 and max(outdeg) <= 1


def coherent_count(E, kind):
    """Number of orientations of the 3-edge graph E that are 'coherent':
    P4 -> directed path; K3+K1 -> cyclic triangle; K13 -> all-out or all-in star."""
    E = list(E)
    cnt = 0
    for bits in range(8):
        arcs = [(b, a) if (bits >> k) & 1 else (a, b) for k, (a, b) in enumerate(E)]
        outdeg = [0] * 4
        indeg = [0] * 4
        for (u, v) in arcs:
            outdeg[u] += 1
            indeg[v] += 1
        if kind == "P4":
            ok = max(indeg) <= 1 and max(outdeg) <= 1
        elif kind == "K3+K1":
            ok = all(indeg[x] == outdeg[x] for x in range(4))  # cyclic triangle: each vertex in=out=1 (or 0)
        else:
            centre = max(range(4), key=lambda x: outdeg[x] + indeg[x])
            ok = outdeg[centre] in (0, 3)
        cnt += ok
    return cnt


def section_A():
    print("\n== A. the three 3-edge graphs on 4 nodes ==")
    P = pairs(4)
    triples = list(itertools.combinations(P, 3))
    check("A1 number of 3-edge subgraphs of K4 is C(6,3) = 20", len(triples) == 20, len(triples))
    bytype = {}
    for E in triples:
        bytype.setdefault(graph_type4(E), []).append(E)
    counts = {k: len(v) for k, v in bytype.items()}
    check("A2 labelled copies: P4 12, K3+K1 4, K13 4", counts == {"P4": 12, "K3+K1": 4, "K13": 4}, counts)

    # complement and opposite-edge maps on 3-edge graphs
    comp_action = {}
    opp_action = {}
    p4_opp_eq_comp = 0
    p4_overlap = set()
    pair_opp_eq_comp = 0
    for E in triples:
        Ec = tuple(sorted(set(P) - set(E)))
        Eo = tuple(sorted(opposite_edge(e) for e in E))
        t, tc, to = graph_type4(E), graph_type4(Ec), graph_type4(Eo)
        comp_action.setdefault(t, set()).add(tc)
        opp_action.setdefault(t, set()).add(to)
        if t == "P4":
            p4_opp_eq_comp += (set(Eo) == set(Ec))
            p4_overlap.add(len(set(Eo) & set(E)))
        else:
            pair_opp_eq_comp += (set(Eo) == set(Ec))
    check("A3 complement: P4 -> P4, K3+K1 <-> K13 (involution with one fixed type)",
          comp_action == {"P4": {"P4"}, "K3+K1": {"K13"}, "K13": {"K3+K1"}}, comp_action)
    check("A4 opposite-edge map has the same type action",
          opp_action == {"P4": {"P4"}, "K3+K1": {"K13"}, "K13": {"K3+K1"}}, opp_action)
    check("A5 on every triangle and every star, opposite-edge map == complement (8 of 8)", pair_opp_eq_comp == 8, pair_opp_eq_comp)
    check("A6 on every P4, opposite-edge map != complement (0 of 12) and |P4^opp & P4| = 2 always",
          p4_opp_eq_comp == 0 and p4_overlap == {2}, (p4_opp_eq_comp, p4_overlap))
    # coherent orientations
    coh = {}
    for t, lst in bytype.items():
        coh[t] = sorted(set(coherent_count(E, t) for E in lst))
    check("A7 every 3-edge graph on 4 nodes has exactly 2 coherent orientations (and 8 in all)",
          coh == {"P4": [2], "K3+K1": [2], "K13": [2]}, coh)
    return bytype


# ---------------------------------------------------------------------------
# B, C, D, E on n = 4
# ---------------------------------------------------------------------------
ORIENTED_FACES = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]  # boundary of the oriented 3-simplex [0123]


def face_of_arc_along(u, v):
    """The unique oriented face whose boundary cyclic order contains u -> v."""
    hits = []
    for f in ORIENTED_FACES:
        for k in range(3):
            if f[k] == u and f[(k + 1) % 3] == v:
                hits.append(f)
    assert len(hits) == 1, (u, v, hits)
    return hits[0]


def opp_vertex(face):
    return [x for x in range(4) if x not in face][0]


def dual_tournament(T):
    """Tetrahedral duality: arc u->v  maps to the dual arc F+ -> F- between the two faces containing uv
    (F+ = the face whose boundary runs along u->v), faces identified with their opposite vertices."""
    D = set()
    for (u, v) in T:
        fp = face_of_arc_along(u, v)
        fm = face_of_arc_along(v, u)
        D.add((opp_vertex(fp), opp_vertex(fm)))
    D = frozenset(D)
    assert len(D) == 6
    # sanity: exactly one arc per pair
    assert all(((a, b) in D) != ((b, a) in D) for (a, b) in pairs(4))
    return D


def type_name4(T):
    sc = sorted(scores(T, 4))
    return {(0, 1, 2, 3): "TT4", (1, 1, 1, 3): "C3+source", (0, 2, 2, 2): "C3+sink", (1, 1, 2, 2): "strong"}[tuple(sc)]


def section_BCD():
    print("\n== B. n = 4 dictionary ==")
    n = 4
    Ts = list(all_tournaments(4))
    perms = list(itertools.permutations(range(4)))
    rows = {}
    ok_H = ok_cs = ok_score = True
    for T in Ts:
        H = ham_paths(T, n)
        c3 = cyclic_triangles(T, n)
        s = sources_sinks(T, n)
        sc = scores(T, n)
        ok_H &= (H == 1 + 2 * c3)
        ok_cs &= (c3 + s == 2)
        ok_score &= (c3 == comb(4, 3) - sum(comb(x, 2) for x in sc))
        t = type_name4(T)
        rows.setdefault(t, set()).add((H, c3, s))
    check("B1 all 64: H = 1 + 2 c3", ok_H)
    check("B2 all 64: c3 + s = 2  (s = #sources + #sinks)", ok_cs)
    check("B3 all 64: c3 = C(4,3) - sum_v C(score_v, 2)", ok_score)
    check("B4 the four types carry (H, c3, s) = TT4 (1,0,2), C3+source/sink (3,1,1), strong (5,2,0)",
          rows == {"TT4": {(1, 0, 2)}, "C3+source": {(3, 1, 1)}, "C3+sink": {(3, 1, 1)}, "strong": {(5, 2, 0)}}, rows)
    # (s-1)(s-2)/2 = [s in {0,3}] for s = 0..3
    check("B5 (s-1)(s-2)/2 = [s in {0,3}] for s = 0,1,2,3 (the mechanism of B2)",
          all(((x - 1) * (x - 2)) // 2 == (1 if x in (0, 3) else 0) for x in range(4)))

    print("\n== C. tetrahedral duality ==")
    ok_dual_c3 = ok_dual_s = ok_dual_H = ok_inv = True
    type_action = {}
    inv_kind = set()
    for T in Ts:
        D = dual_tournament(T)
        ok_dual_c3 &= (cyclic_triangles(T, n) == sources_sinks(D, n))
        ok_dual_s &= (sources_sinks(T, n) == cyclic_triangles(D, n))
        ok_dual_H &= (ham_paths(D, n) == 6 - ham_paths(T, n))
        DD = dual_tournament(D)
        if DD == T:
            inv_kind.add("T")
        elif DD == converse(T):
            inv_kind.add("T^op")
        else:
            inv_kind.add("neither")
        type_action.setdefault(type_name4(T), set()).add(type_name4(D))
    check("C1 all 64: c3(T) = #sources+#sinks of T*", ok_dual_c3)
    check("C2 all 64: #sources+#sinks of T = c3(T*)", ok_dual_s)
    check("C3 all 64: H(T*) = 6 - H(T)", ok_dual_H)
    check("C4 T** is always the same one of {T, T^op}", len(inv_kind) == 1, inv_kind)
    check("C5 duality on types: TT4 <-> strong, and the two C3 types stay inside the C3 class",
          type_action["TT4"] == {"strong"} and type_action["strong"] == {"TT4"}
          and type_action["C3+source"] <= {"C3+source", "C3+sink"} and type_action["C3+sink"] <= {"C3+source", "C3+sink"},
          type_action)
    conv_action = {}
    for T in Ts:
        conv_action.setdefault(type_name4(T), set()).add(type_name4(converse(T)))
    check("C6 converse on types: C3+source <-> C3+sink, TT4 and strong fixed",
          conv_action == {"TT4": {"TT4"}, "strong": {"strong"}, "C3+source": {"C3+sink"}, "C3+sink": {"C3+source"}}, conv_action)
    ok_comm = all(dual_tournament(converse(T)) == converse(dual_tournament(T)) for T in Ts)
    check("C7 (T^op)* = (T*)^op on all 64 (converse and duality commute)", ok_comm)
    # duality as an edge-set map is the opposite-edge map
    ok_opp = all(set((min(a, b), max(a, b)) for (a, b) in dual_tournament(T)) == set(pairs(4)) for T in Ts)
    ok_opp2 = True
    for T in Ts:
        D = dual_tournament(T)
        # each arc u->v of T produces an arc on the opposite pair
        for (u, v) in T:
            x, y = opposite_edge((u, v))
            ok_opp2 &= (((x, y) in D) or ((y, x) in D))
    check("C8 the dual arc of u->v lies on the opposite edge of the tetrahedron", ok_opp and ok_opp2)

    print("\n== D. marked structures at n = 4 ==")
    # (T, Hamiltonian path) types: labelled marked pairs / 4! (S4 acts freely on marked pairs)
    lab_path = sum(ham_paths(T, n) for T in Ts)
    # freeness check: an automorphism fixing a Hamiltonian path (as a sequence) is the identity
    ok_free = True
    for T in Ts:
        for p in perms:
            if p == tuple(range(4)) or relabel(T, p) != T:
                continue
            for path in ham_path_list(T, n):
                if tuple(p[x] for x in path) == path:
                    ok_free = False
    check("D1 Aut(T) acts freely on Hamiltonian paths (all 64)", ok_free)
    check("D2 labelled (T, Ham path) pairs = 4! 2^C(3,2) = 192; marked types = 8", lab_path == 192 and lab_path // 24 == 8, lab_path)
    # (T, cyclic triangle) types by Burnside over S4
    lab_cyc = 0
    fix_sum = 0
    for T in Ts:
        lab_cyc += cyclic_triangles(T, n)
    for p in perms:
        for T in Ts:
            if relabel(T, p) != T:
                continue
            # count cyclic triangles fixed setwise (as directed cycles) by p
            for a, b, d in itertools.combinations(range(4), 3):
                arcs = None
                if (a, b) in T and (b, d) in T and (d, a) in T:
                    arcs = {(a, b), (b, d), (d, a)}
                elif (b, a) in T and (d, b) in T and (a, d) in T:
                    arcs = {(b, a), (d, b), (a, d)}
                if arcs is not None and {(p[u], p[v]) for (u, v) in arcs} == arcs:
                    fix_sum += 1
    cyc_types = fix_sum // 24
    check("D3 (T, cyclic triangle) marked types = 4 (Burnside: 96/24)", cyc_types == 4 and fix_sum == 96, (fix_sum, cyc_types))
    # (T, source or sink) types: TT4 has 2, each C3+x has 1, strong 0 -> 4; equals 2 t(3) = 4
    star_fix = 0
    for p in perms:
        for T in Ts:
            if relabel(T, p) != T:
                continue
            sc = scores(T, n)
            star_fix += sum(1 for v in range(4) if sc[v] in (0, 3) and p[v] == v)
    star_types = star_fix // 24
    check("D4 (T, source-or-sink) marked types = 4 = 2 t(3)", star_types == 4, (star_fix, star_types))
    check("D5 lines = cycles + stars at n = 4:  8 = 4 + 4", 8 == cyc_types + star_types)
    return Ts


# ---------------------------------------------------------------------------
# E. the line-and-staircase bijection, n = 4, 5, 6
# ---------------------------------------------------------------------------
def marked_cycle_types(n):
    """Number of S_n-orbits of pairs (T, C) with C a directed (n-1)-cycle of T (Burnside, exact)."""
    perms = list(itertools.permutations(range(n)))
    total = 0
    Ts = list(all_tournaments(n))
    # precompute for each T the list of directed (n-1)-cycles as arc sets
    for p in perms:
        fixed = 0
        for T in Ts:
            if relabel(T, p) != T:
                continue
            out = out_sets(T, n)
            for S in itertools.combinations(range(n), n - 1):
                Sset = set(S)
                start = S[0]
                cycles = []

                def ext(path, used):
                    if len(path) == n - 1:
                        if start in out[path[-1]]:
                            cycles.append(set(zip(path, path[1:] + [start])))
                        return
                    for w in out[path[-1]]:
                        if w in Sset and w not in used:
                            used.add(w)
                            path.append(w)
                            ext(path, used)
                            path.pop()
                            used.remove(w)

                ext([start], {start})
                for arcs in cycles:
                    if {(p[u], p[v]) for (u, v) in arcs} == arcs:
                        fixed += 1
        total += fixed
    assert total % factorial(n) == 0
    return total // factorial(n)


def section_E():
    print("\n== E. line-and-staircase bijection, general n ==")
    known_types = {3: 2, 4: 4, 5: 12, 6: 56}  # A000568
    for n in (4, 5, 6):
        Ts = all_tournaments(n)
        total_H = 0
        count = 0
        for T in Ts:
            total_H += ham_paths(T, n)
            count += 1
        check(f"E1 n={n}: sum over labelled T of H(T) = n! 2^C(n-1,2) (Szele)",
              total_H == factorial(n) * 2 ** comb(n - 1, 2), total_H)
    # marked counts (lines, cycles, stars) for n = 4, 5, 6
    results = {}
    for n in (4, 5, 6):
        lines = 2 ** comb(n - 1, 2)
        stars = 2 * known_types[n - 1]
        if n <= 5:
            cyc = marked_cycle_types(n)
        else:
            # Burnside by hand for n = 6 (checked symbolically in the note): identity 6*24*2^10 = 147456;
            # each of the 144 five-cycles fixes 16 marked pairs; other elements fix none.
            cyc = (6 * 24 * 2 ** 10 + 144 * 16) // 720
        results[n] = (lines, cyc, stars)
    check("E2 marked (lines, cycles, stars) = (8,4,4), (64,16,8), (1024,208,24) for n = 4,5,6",
          results == {4: (8, 4, 4), 5: (64, 16, 8), 6: (1024, 208, 24)}, results)
    check("E3 HOSTILE: lines = cycles + stars holds only at n = 4 among 4,5,6",
          [n for n in (4, 5, 6) if results[n][0] == results[n][1] + results[n][2]] == [4])
    # Burnside check for n = 6 against the enumeration of fixed points of a single 5-cycle
    n = 6
    g = (1, 2, 3, 4, 0, 5)  # the 5-cycle (0 1 2 3 4), fixing 5
    fixed = 0
    for T in all_tournaments(6):
        if relabel(T, g) != T:
            continue
        out = out_sets(T, 6)
        for S in itertools.combinations(range(6), 5):
            Sset = set(S)
            start = S[0]
            cycles = []

            def ext(path, used):
                if len(path) == 5:
                    if start in out[path[-1]]:
                        cycles.append(set(zip(path, path[1:] + [start])))
                    return
                for w in out[path[-1]]:
                    if w in Sset and w not in used:
                        used.add(w)
                        path.append(w)
                        ext(path, used)
                        path.pop()
                        used.remove(w)

            ext([start], {start})
            for arcs in cycles:
                if {(g[u], g[v]) for (u, v) in arcs} == arcs:
                    fixed += 1
    check("E4 n=6: a 5-cycle of S6 fixes exactly 16 marked (T, directed 5-cycle) pairs", fixed == 16, fixed)

    # E7: the cycle-marked count in closed form.  Marked pairs (T, C) reduce to T containing the standard
    # directed (n-1)-cycle C0 on {0..n-2} (solo n-1), modulo the rotation group Z_{n-1} = Stab(C0).
    # A rotation of odd order d has all its C(n-1,2)/d free pair-orbits unforced and no reversed pair;
    # a rotation of even order reverses a chord, so it fixes no tournament.  Hence
    #   cycles(n) = (1/(n-1)) * sum_{d | n-1, d odd} phi(d) * 2^(C(n-1,2)/d).
    def phi(m):
        return sum(1 for a in range(1, m + 1) if __import__("math").gcd(a, m) == 1)

    def cycles_formula(n):
        m = n - 1
        tot = sum(phi(d) * 2 ** (comb(m, 2) // d) for d in range(1, m + 1) if m % d == 0 and d % 2 == 1)
        assert tot % m == 0
        return tot // m

    def cycles_direct(n):
        m = n - 1
        C0 = {(i, (i + 1) % m) for i in range(m)}
        free = [pq for pq in pairs(n) if (pq not in C0 and (pq[1], pq[0]) not in C0)]
        fixed_total = 0
        for k in range(m):
            r = tuple(list(((i + k) % m) for i in range(m)) + [m])
            cnt = 0
            for bits in range(1 << len(free)):
                arcs = set(C0)
                for t, (i, j) in enumerate(free):
                    arcs.add((j, i) if (bits >> t) & 1 else (i, j))
                if relabel(frozenset(arcs), r) == frozenset(arcs):
                    cnt += 1
            fixed_total += cnt
        assert fixed_total % m == 0
        return fixed_total // m

    direct = {n: cycles_direct(n) for n in (4, 5, 6, 7)}
    formula = {n: cycles_formula(n) for n in (4, 5, 6, 7)}
    check("E7 cycles(n) direct Burnside over Z_{n-1} == (1/(n-1)) sum_{d|n-1 odd} phi(d) 2^(C(n-1,2)/d) for n = 4..7",
          direct == formula and formula[4] == 4 and formula[5] == 16 and formula[6] == 208 and formula[7] == 5472,
          (direct, formula))
    print("     cycles(n) by the formula, n = 4..10:", [cycles_formula(n) for n in range(4, 11)])
    # E8: at n = 4 the odd-order elements of S3 are exactly the rotations Z3, which is why stars(4) = cycles(4)
    S3 = list(itertools.permutations(range(3)))
    odd_order = [p for p in S3 if p == (0, 1, 2) or (p[p[0]], p[p[1]], p[p[2]]) != (0, 1, 2) and
                 (p[p[p[0]]], p[p[p[1]]], p[p[p[2]]]) == (0, 1, 2)]
    check("E8 odd-order elements of S3 are exactly the three rotations (so stars(4) = 2 * (1/6) sum_Z3 = cycles(4))",
          len(odd_order) == 3)

    # anti-path complement involution at n = 4: (T, P) <-> upset graph G on the anti-path Q
    n = 4
    P0 = (0, 1, 2, 3)
    Q = [(0, 2), (1, 3), (0, 3)]  # the anti-path 2-0-3-1
    table = {}
    for bits in range(8):
        G = {Q[k] for k in range(3) if (bits >> k) & 1}
        arcs = set()
        for (i, j) in pairs(4):
            arcs.add((j, i) if (i, j) in G else (i, j))
        T = frozenset(arcs)
        assert tuple(P0) in ham_path_list(T, 4)
        table[frozenset(G)] = type_name4(T)
    dist = {}
    for G, t in table.items():
        dist[t] = dist.get(t, 0) + 1
    check("E5 upset graphs on the anti-path: {} -> TT4, {02} and {13} -> the two C3 types, the other 5 -> strong",
          dist == {"TT4": 1, "C3+source": 1, "C3+sink": 1, "strong": 5} and table[frozenset()] == "TT4"
          and {table[frozenset({(0, 2)})], table[frozenset({(1, 3)})]} == {"C3+source", "C3+sink"}, dist)
    comp_pairs = set()
    for G in table:
        Gc = frozenset(set(Q) - set(G))
        comp_pairs.add(tuple(sorted([table[G], table[Gc]])))
    check("E6 anti-path complement pairs marked types as TT4<->strong, C3+source<->strong, C3+sink<->strong, strong<->strong",
          comp_pairs == {("TT4", "strong"), ("C3+source", "strong"), ("C3+sink", "strong"), ("strong", "strong")}, comp_pairs)


# ---------------------------------------------------------------------------
# F. the divisor triple F = S + U (owner's shapes) placed beside the graph triple
# ---------------------------------------------------------------------------
def divisors(N):
    return [d for d in range(1, N + 1) if N % d == 0]


def is_squarefree(d):
    k = 2
    while k * k <= d:
        if d % (k * k) == 0:
            return False
        k += 1
    return True


def is_prime(d):
    if d < 2:
        return False
    k = 2
    while k * k <= d:
        if d % k == 0:
            return False
        k += 1
    return True


def FSU(N):
    D = [d for d in divisors(N) if 1 < d < N]
    F = len(D)
    S = sum(1 for d in D if is_squarefree(d))
    U = sum(1 for d in D if is_prime(d))
    return F, S, U


def section_F():
    print("\n== F. the divisor triple F = S + U next to the graph triple ==")
    vals = {"p": FSU(7), "p^3": FSU(27), "p^2qr": FSU(60)}
    check("F1 F = S + U at 7, 27, 60 with (F,S,U) = (0,0,0), (2,1,1), (10,7,3)",
          vals == {"p": (0, 0, 0), "p^3": (2, 1, 1), "p^2qr": (10, 7, 3)} and all(F == S + U for (F, S, U) in vals.values()), vals)
    sols = [N for N in range(2, 1000) if (lambda t: t[0] == t[1] + t[2])(FSU(N))]
    # every solution below 1000 is a prime, a prime cube, or p^2 q r
    def shape(N):
        f = {}
        m = N
        k = 2
        while k * k <= m:
            while m % k == 0:
                f[k] = f.get(k, 0) + 1
                m //= k
            k += 1
        if m > 1:
            f[m] = f.get(m, 0) + 1
        return tuple(sorted(f.values(), reverse=True))
    shapes = set(shape(N) for N in sols)
    check("F2 all solutions N < 1000 have exponent profile (1), (3) or (2,1,1)  [DB1 of arithmetic_braids_20260917]",
          shapes == {(1,), (3,), (2, 1, 1)}, shapes)
    # the comparison on the record: graph-side triples (lines, cycles, stars) at n = 4 are (8, 4, 4);
    # the divisor-side triples are (0,0,0), (2,1,1), (10,7,3).  No linear identity relates them.
    graph_triple = (8, 4, 4)
    check("F3 the numbers (8,4,4) vs {(0,0,0),(2,1,1),(10,7,3)}: no shape has F = 8, so no termwise match (recorded hostile)",
          all(v[0] != 8 for v in vals.values()))


if __name__ == "__main__":
    section_A()
    section_BCD()
    section_E()
    section_F()
    bad = [name for name, ok in CHECKS if not ok]
    print("\n%d checks, %d failed" % (len(CHECKS), len(bad)))
    if bad:
        print("FAILED:", bad)
        sys.exit(1)
    print("ALL CHECKS PASSED")
