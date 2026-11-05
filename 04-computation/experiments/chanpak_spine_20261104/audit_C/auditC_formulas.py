#!/usr/bin/env python3
"""auditC_formulas.py -- audit C: independent checks of the formulas in THM-4876 (1), (4), (6) and HYP-9469 (1)/(4).

(a) path-pair formula H(X+v) = sum_{S nonempty proper} e_B(S) s_A(V-S) + s_A(V) + e_B(V) against brute force,
    all tournaments of orders 2..6 (gentourng) and all labeled tournaments of order <= 4, all cuts;
    Claim A at every (T, v), orders 3..8; the odd-cycle formula H = I(Omega, 2), orders 3..7.
(d) flip formula H(T') - H(T) = U(L_j, R_i) - U(L_i, R_j), all tournaments of orders 3..6 (+ all labeled <= 5), all arcs;
    U(L_i, R_j) = #paths of T with i immediately before j (brute force); THM-017's signed B(L_i,R_j) = B(L_j,R_i).
(f) one-centre family: H(TT_m + u) = 1 + 2 N(A) brute force for m <= 7, all A; closed form (2^i-1)(2^l-1)+2;
    N(A) = 3 impossible; coverage below 200,000; the reversed-extreme-arc tournament R_n and its end counts;
    the sole-exit increment versus (h - e_a)/2 and (h + sigma - e_a)/2.
"""
import sys, itertools
from math import comb
from auditC_lib import *


def path_pair_H(n, adj, A, E, St):
    """my implementation of the path-pair formula for the cut A (bitmask; v -> A, V \\ A -> v)."""
    full = (1 << n) - 1
    B = full & ~A

    def eB(S):
        return sum(E[S][u] for u in range(n) if (S >> u) & 1 and (B >> u) & 1)

    def sA(S):
        return sum(St[S][w] for w in range(n) if (S >> w) & 1 and (A >> w) & 1)
    tot = 0
    for S in range(1, full):
        tot += eB(S) * sA(full & ~S)
    return tot + sA(full) + eB(full)


def check_a():
    print("== (a) path-pair formula, Claim A, odd-cycle formula ==")
    bad = tot = 0
    for n in range(2, 7):
        for adj in tournaments(n):
            E = dp_end(n, adj)
            St = dp_start(n, adj)
            for A in range(1 << n):
                n2, adj2 = add_vertex(n, adj, A)
                tot += 1
                if H(n2, adj2) != path_pair_H(n, adj, A, E, St):
                    bad += 1
    print("path-pair formula vs brute force: %d (X, cut) pairs over gentourng classes of orders 2..6: mismatches %d" % (tot, bad))
    bad = tot = 0
    for n in range(1, 5):
        for adj in labeled_tournaments(n):
            E = dp_end(n, adj)
            St = dp_start(n, adj)
            for A in range(1 << n):
                n2, adj2 = add_vertex(n, adj, A)
                tot += 1
                if H(n2, adj2) != path_pair_H(n, adj, A, E, St):
                    bad += 1
    print("path-pair formula vs brute force: %d labeled (X, cut) pairs of orders 1..4: mismatches %d" % (tot, bad))
    # Claim A at every (T, v), orders 3..8
    for n in range(3, 9):
        bad = tot = 0
        for adj in tournaments(n):
            E = dp_end(n, adj)
            hs = hsub_table(n, E)
            full = (1 << n) - 1
            Ht = hs[full]
            for v in range(n):
                Hv = hs[full & ~(1 << v)]
                rhs = 2 * sum(hs[full & ~S] for S, L in odd_cycles_through(n, adj, v))
                tot += 1
                if Ht - Hv != rhs:
                    bad += 1
        print("Claim A (H(T) - H(T-v) = 2 sum_C H(T - V(C))) at order %d: %d (T, v) pairs, mismatches %d" % (n, tot, bad))
    # odd-cycle formula, orders 3..7
    for n in range(3, 8):
        bad = tot = 0
        for adj in tournaments(n):
            tot += 1
            if H(n, adj) != I_omega_2(all_odd_cycles(n, adj)):
                bad += 1
        print("odd-cycle formula H = I(Omega, 2) at order %d: %d tournaments, mismatches %d" % (n, tot, bad))
    # lemma: a vertex on an odd cycle of length >= 5 lies on a shorter odd cycle (orders 5..8)
    for n in range(5, 9):
        bad = tot = 0
        for adj in tournaments(n):
            for v in range(n):
                cyc = odd_cycles_through(n, adj, v)
                if not cyc:
                    continue
                tot += 1
                if min(L for S, L in cyc) != 3:
                    bad += 1
        print("lemma (vertex on an odd cycle lies on a triangle) at order %d: %d (T, v) with an odd cycle, failures %d" % (n, tot, bad))


def U_conv(n, E, St, W, x, y):
    """U(L_x, R_y) = sum_{S subset W} L_x(S u {x}) R_y({y} u (W \\ S)); W a bitmask."""
    tot = 0
    sub = W
    # iterate over all subsets S of W
    S = W
    while True:
        tot += E[S | (1 << x)][x] * St[(W & ~S) | (1 << y)][y]
        if S == 0:
            break
        S = (S - 1) & W
    return tot


def B_conv(n, E, St, W, x, y):
    tot = 0
    S = W
    while True:
        tot += (-1) ** bin(S).count('1') * E[S | (1 << x)][x] * St[(W & ~S) | (1 << y)][y]
        if S == 0:
            break
        S = (S - 1) & W
    return tot


def check_d():
    print("== (d) flip formula ==")
    for label, gen in (("gentourng classes orders 3..6", [(n, adj) for n in range(3, 7) for adj in tournaments(n)]),
                       ("labeled tournaments orders 3..5", [(n, adj) for n in range(3, 6) for adj in labeled_tournaments(n)])):
        bad = badU = badB = badT = tot = 0
        for n, adj in gen:
            full = (1 << n) - 1
            E = dp_end(n, adj)
            St = dp_start(n, adj)
            Ht = sum(E[full])
            paths = ham_paths_list(n, adj)
            for i in range(n):
                for j in range(n):
                    if not (adj[i] >> j) & 1:
                        continue
                    W = full & ~(1 << i) & ~(1 << j)
                    adj2 = list(adj)
                    adj2[i] &= ~(1 << j)
                    adj2[j] |= 1 << i
                    E2 = dp_end(n, adj2)
                    St2 = dp_start(n, adj2)
                    Ht2 = sum(E2[full])
                    # U(L_i, R_j) in T, U(L_j, R_i) in T' and in T
                    Uij = U_conv(n, E, St, W, i, j)
                    Uji_T2 = U_conv(n, E2, St2, W, j, i)
                    Uji_T = U_conv(n, E, St, W, j, i)
                    tot += 1
                    if Ht2 - Ht != Uji_T2 - Uij:
                        bad += 1
                    if Uji_T2 != Uji_T:
                        badT += 1
                    # U(L_i,R_j) = number of Hamiltonian paths of T with i immediately before j
                    cnt = sum(1 for p in paths if any(p[t] == i and p[t + 1] == j for t in range(n - 1)))
                    if cnt != Uij:
                        badU += 1
                    # THM-017
                    if B_conv(n, E, St, W, i, j) != B_conv(n, E, St, W, j, i):
                        badB += 1
        print("%s: %d arcs; flip formula mismatches %d; U(L_j,R_i) differs between T and T': %d; U(L_i,R_j) != #paths with i->j consecutive: %d; THM-017 B(L_i,R_j) != B(L_j,R_i): %d" % (label, tot, bad, badT, badU, badB))


def c_of(d):
    return 1 if d == 1 else 2 ** (d - 2)


def N_of(A, m):
    return sum(c_of(b - a) for a in A for b in range(1, m + 1) if b not in A and a < b)


def one_centre(m, A):
    """X = TT_m + u on vertices 0..m-1 (vertex i+1 of the theorem is index i), u = m; u -> A (1-based), rest -> u."""
    n = m + 1
    adj = [0] * n
    for i in range(m):
        for j in range(i + 1, m):
            adj[i] |= 1 << j
        if (i + 1) in A:
            adj[m] |= 1 << i
        else:
            adj[i] |= 1 << m
    return n, adj


def check_f():
    print("== (f) one-centre family ==")
    bad = tot = 0
    vals3 = 0
    for m in range(1, 8):
        for r in range(0, m + 1):
            for A in itertools.combinations(range(1, m + 1), r):
                n, adj = one_centre(m, set(A))
                tot += 1
                if H(n, adj) != 1 + 2 * N_of(set(A), m):
                    bad += 1
    print("H(TT_m + u) = 1 + 2 N(A): %d (m, A) pairs, m <= 7: mismatches %d" % (tot, bad))
    # closed form (2^i - 1)(2^l - 1) + 2 for A = {1..i}, m = i + l
    bad = 0
    for i in range(1, 9):
        for l in range(1, 9):
            if 1 + 2 * N_of(set(range(1, i + 1)), i + l) != (2 ** i - 1) * (2 ** l - 1) + 2:
                bad += 1
    print("closed form (2^i-1)(2^l-1)+2 for A = {1..i}, m = i+l, i,l <= 8: mismatches %d" % bad)
    # symbolic
    import sympy as sp
    i, l = sp.symbols('i l', positive=True, integer=True)
    N_formula = 2 ** (l - 1) + (2 ** (i - 1) - 1) * (2 ** l - 1)
    print("algebra: 1 + 2 N simplifies to (2^i-1)(2^l-1)+2:", sp.simplify(1 + 2 * N_formula - ((2 ** i - 1) * (2 ** l - 1) + 2)) == 0)
    # Fermat case A = {1}: 2^(m-1) + 1 ; and which vertex ends (h-1)/2 paths
    for m in range(2, 9):
        n, adj = one_centre(m, {1})
        E = dp_end(n, adj)
        full = (1 << n) - 1
        h = sum(E[full])
        ends = E[full]
        print("A={1}, m=%d (order %d): H = %d (2^(m-1)+1 = %d); end counts by vertex (1..m, u): %s; max over non-sinks %d = (h-1)/2? %s at vertex %s" % (
            m, n, h, 2 ** (m - 1) + 1, ends, max(ends[v] for v in range(n) if adj[v]), max(ends[v] for v in range(n) if adj[v]) == (h - 1) // 2,
            ["%d" % (v + 1) if v < m else "u" for v in range(n) if adj[v] and ends[v] == max(ends[w] for w in range(n) if adj[w])]))
    # N(A) = 3 impossible: all words of length <= 16
    seen = set()
    for m in range(1, 17):
        for bits in itertools.product((0, 1), repeat=m):
            A = {i + 1 for i in range(m) if bits[i]}
            seen.add(N_of(A, m))
    print("N values attained by all words of length <= 16 below 12:", sorted(x for x in seen if x < 12), " 3 attained:", 3 in seen)
    # coverage below 200,000: canonical words start with 1 and end with 0; every 1 at position a pairs with the last 0 at b:
    # c(b - a) <= N <= 10^5 forces b - a <= 18, so m <= 19 suffices; enumerate m <= 21 to be safe.
    LIMIT = 200000
    Nmax = (LIMIT - 1) // 2
    vals = set([1])
    count_words = 0
    for m in range(2, 22):
        # positions 1 and m fixed (1 ... 0); interior free, prune by N > Nmax impossible in simple product loop -> use recursion
        def rec(pos, A, N):
            nonlocal count_words
            if N > Nmax:
                return
            if pos == m:
                # last letter 0: add pairs (a, m) for a in A
                N2 = N + sum(c_of(m - a) for a in A)
                if N2 <= Nmax:
                    vals.add(1 + 2 * N2)
                    count_words += 1
                return
            # letter 1 at pos
            rec(pos + 1, A + [pos], N)
            # letter 0 at pos
            rec(pos + 1, A, N + sum(c_of(pos - a) for a in A))
        rec(2, [1], 0)
    odd_below = (LIMIT) // 2
    attained = sorted(v for v in vals if v < LIMIT)
    missed = [v for v in range(1, LIMIT, 2) if v not in vals]
    print("one-centre values below %d: %d of %d odd numbers attained (words enumerated: %d); missed (first 30): %s" % (LIMIT, len(attained), odd_below, count_words, missed[:30]))
    print("   7 attained:", 7 in vals, "; 13 attained:", 13 in vals, "; (2^i-1)(2^l-1)+2 samples 11,23,47,51,95,107,227 attained:", all(x in vals for x in (11, 23, 47, 51, 95, 107, 227)))
    # sole-exit increment versus (h - e_a)/2 and (h + sigma - e_a)/2
    eq_simple = eq_sigma = tot = 0
    for n in range(2, 7):
        for adj in tournaments(n):
            E = dp_end(n, adj)
            St = dp_start(n, adj)
            hs = hsub_table(n, E)
            full = (1 << n) - 1
            h = hs[full]
            for a in range(n):
                n2, adj2 = add_vertex(n, adj, 1 << a)
                inc2 = H(n2, adj2) - h
                e_a = E[full][a]
                Ba = sum(hs[S] * St[full & ~S][a] for S in range(1 << n) if not (S >> a) & 1)
                sigma = Ba - h
                tot += 1
                if inc2 == h - e_a:
                    eq_simple += 1
                if inc2 == h + sigma - e_a:
                    eq_sigma += 1
    print("sole-exit cut (v -> a only), orders 2..6, %d (T, a): increment equals (h - e_a)/2 in %d cases; equals (h + sigma - e_a)/2 (THM-4862 (a), sigma = B_a - h >= 0) in %d cases" % (tot, eq_simple, eq_sigma))
    # example: the strong 4-tournament, vertex of out-degree 1
    n, adj = 4, [0b0010, 0b1100, 0b1000, 0b0001 >> 0]
    # build explicitly: 0->1, 1->2, 2->0, 0->3? choose the strong 4-tournament: 0->1->2->0 triangle, 3: 0->3, 3->1, 3->2? check strong
    adj = [0] * 4
    arcs = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 1), (2, 3)]
    for x, y in arcs:
        adj[x] |= 1 << y
    assert is_strong(4, adj) and H(4, adj) == 5
    E = dp_end(4, adj)
    St = dp_start(4, adj)
    hs = hsub_table(4, E)
    for a in range(4):
        if bin(adj[a]).count('1') == 1:
            n2, adj2 = add_vertex(4, adj, 1 << a)
            Ba = sum(hs[S] * St[15 & ~S][a] for S in range(16) if not (S >> a) & 1)
            print("strong 4-tournament, sole-exit cut at the out-degree-1 vertex a=%d: h=5, e_a=%d, B_a=%d (sigma=%d), H(T+v)=%d, increment %d; (h-e_a)/2 = %s" % (a, E[15][a], Ba, Ba - 5, H(n2, adj2), (H(n2, adj2) - 5) // 2, (5 - E[15][a]) / 2))


if __name__ == "__main__":
    check_a()
    check_d()
    check_f()
