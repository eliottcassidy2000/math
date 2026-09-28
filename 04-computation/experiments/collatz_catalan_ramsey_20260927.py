#!/usr/bin/env python3
"""collatz_catalan_ramsey_20260927.py -- Catalan numbers, the Ramsey rows and the Collatz spine
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, seventh note).

 (1) Catalan numbers are the spine-block counts of the critical member q = 4 of the family T_q (multiplier 4^o/2^k > 1 iff
     the +-1 walk stays positive): positive min-ending words of length 2n+1 are counted by C_n (a +1 step followed by a Dyck
     path), positive words by the central binomials C(k-1, floor((k-1)/2)), W = 1/(1 - P) as in THM-4495's ladder
     decomposition; the certification density at q = 4 decays like k^(-1/2) (the boundary of the fourth note's dichotomy),
     and C_n ~ 4^n n^(-3/2)/sqrt(pi) carries the k^(-3/2) exponent of THM-4495's W_k ~ 2^(h* k) k^(-3/2).
 (2) Catalan parity: C_n is odd iff n = 2^k - 1 (the Mersenne orders of the tournament tower); C_n is coprime to 6 iff
     n = 2^k - 1 and 2^k - 1 has no ternary digit 2 (Kummer): n = 1, 3, 31, 255 below 4000, and k = 1, 2, 5, 8 below 200 --
     the ternary-digit condition of Erdos's conjecture on 2^k, shifted by one.
 (3) Paley Ramsey rows: clique = independence number of the Paley graph P_q (q = 1 mod 4 prime, q <= 101) gives R(w+1, w+1) > q;
     the largest transitive subtournament of the Paley tournament P_q (q = 3 mod 4, q <= 31) gives R(t+1) > q for the
     tournament Ramsey numbers R(2..6) = 2, 4, 8, 14, 28 (Erdos-Moser, Reid-Parker, Sanchez-Flores; THM-455).
 (4) A graph doubling in the spirit of THM-447: for a graph G with Seidel matrix S (S_ij = -1 if adjacent, +1 if not, 0 on the
     diagonal), D(G) has Seidel matrix [[S, S + I], [S + I, -S]]: G, G and its complement on the copies, twins non-adjacent.
     The law for omega(D(G)) - omega(G) and alpha(D(G)) - alpha(G) over all labeled graphs on 5 and 6 vertices, and the values
     on P_5, P_13, P_17 (the Ramsey-extremal Paley graphs).
 (5) The value-time poset of an orbit as a permutation tournament: its largest transitive subtournament is the longer of the
     longest increasing and decreasing subsequences of the orbit's values (Erdos-Szekeres): 27, 703, 6171, 77031.
 (6) The exact law of the Seidel doubling (the graph analogue of THM-483's zigzag law): omega(D(G)) is the largest induced
     complete split subgraph (a clique completely joined to an independent set) and alpha(D(G)) = max(alpha(G) + 1, largest
     induced clique plus independent set with no edges between); verified exhaustively for n <= 6 and on random graphs for
     n = 7, 8; the split numbers of the Paley graphs P_5..P_37 and the Ramsey bounds their doubles give.
Usage: python3 collatz_catalan_ramsey_20260927.py
"""
import math, itertools, sys
from functools import lru_cache

sys.setrecursionlimit(10000)


def catalan(n):
    return math.comb(2 * n, n) // (n + 1)


def part1():
    print("== (1) Catalan numbers are the critical spine-block counts (q = 4, slope 1) ==")
    # positive words: prefix sums S_j = 2 o_j - j > 0 for all j >= 1 (odd = +1, even = -1)
    # positive min-ending words: positive and S_k = min_(1<=j<=k) S_j (= 1)
    P = {}; W = {}
    for k in range(1, 16):
        cntW = 0; cntP = 0
        for word in itertools.product((1, -1), repeat=k):
            s = 0; ok = True; mn = None
            for x in word:
                s += x
                if s <= 0:
                    ok = False; break
                mn = s if mn is None else min(mn, s)
            if ok:
                cntW += 1
                if s == mn:
                    cntP += 1
        W[k] = cntW; P[k] = cntP
    print(" positive min-ending words P_k, k = 1..15:", [P[k] for k in range(1, 16)])
    print(" Catalan C_n at k = 2n+1 (odd k), 0 at even k:", [catalan((k - 1) // 2) if k % 2 == 1 else 0 for k in range(1, 16)])
    print(" positive words W_k:", [W[k] for k in range(1, 16)], " central binomials C(k-1, floor((k-1)/2)):", [math.comb(k - 1, (k - 1) // 2) for k in range(1, 16)])
    # ladder identity W = 1/(1 - P) as power series
    ok = True
    for k in range(1, 16):
        conv = sum(P[j] * (W[k - j] if k - j > 0 else 1) for j in range(1, k + 1))
        ok &= (conv == W[k])
    print(" W = 1/(1 - P) (W_k = sum_j P_j W_(k-j)): %s" % ok)
    print(" certification density at q = 4: fraction of length-k words that are NOT positive = 1 - C(k-1,floor((k-1)/2))/2^k; the")
    print(" non-certified fraction ~ sqrt(2/(pi k)):", [round(math.comb(k - 1, (k - 1) // 2) / 2 ** k, 4) for k in (10, 100, 1000)], "vs", [round(math.sqrt(2 / (math.pi * k)), 4) for k in (10, 100, 1000)])
    print(" C_n ~ 4^n n^(-3/2)/sqrt(pi): ratio C_n n^(3/2) sqrt(pi)/4^n at n = 10, 100, 500:", [round(catalan(n) * n ** 1.5 * math.sqrt(math.pi) / 4 ** n, 4) for n in (10, 100, 500)])


def ternary_has_two(n):
    while n:
        if n % 3 == 2:
            return True
        n //= 3
    return False


def part2():
    print("== (2) Catalan parity, Mersenne indices and the ternary condition ==")
    N = 4000
    C = [1]
    for n in range(N):
        C.append(C[-1] * 2 * (2 * n + 1) // (n + 2))
    odd = [n for n in range(N + 1) if C[n] % 2 == 1]
    print(" n <= %d with C_n odd:" % N, odd[:12], "... = 2^k - 1: %s" % all(n + 1 & n == 0 for n in odd))
    cop6 = [n for n in range(N + 1) if C[n] % 2 == 1 and C[n] % 3 != 0]
    print(" n <= %d with C_n coprime to 6:" % N, cop6)
    ks = [k for k in range(1, 201) if not ternary_has_two(2 ** k - 1)]
    print(" k <= 200 with 2^k - 1 free of the ternary digit 2 (a sum of distinct powers of 3):", ks, "-> n = 2^k - 1 =", [2 ** k - 1 for k in ks])
    print(" (Kummer: v_3(C_n) = carries of n + n in base 3 minus v_3(n+1); no carries iff no ternary digit 2 in n)")
    print(" Erdos's conjecture concerns the ternary digits of 2^k (no digit 2 only for k = 0, 2, 8); here 2^k - 1, with k = 1, 2, 5, 8 below 200")


def paley_graph_omega(q):
    # q prime, q = 1 mod 4; clique number via Bron-Kerbosch with pivoting on bitsets
    sq = {(x * x) % q for x in range(1, q)}
    N = [0] * q
    for i in range(q):
        for j in range(q):
            if i != j and ((j - i) % q) in sq:
                N[i] |= 1 << j
    best = [0]

    def bk(R, P, X):
        if P == 0 and X == 0:
            best[0] = max(best[0], R)
            return
        if R + bin(P).count("1") <= best[0]:
            return
        # pivot
        PX = P | X
        u = max(range(q), key=lambda v: bin(P & N[v]).count("1") if (PX >> v) & 1 else -1)
        cand = P & ~N[u]
        while cand:
            v = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            bk(R + 1, P & N[v], X & N[v])
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << q) - 1, 0)
    return best[0]


def trans_tournament(n, out):
    # largest transitive subtournament: chains v1 -> v2 -> ... with each beating all later ones
    @lru_cache(maxsize=None)
    def f(S):
        best = 0; s = S
        while s:
            v = (s & -s).bit_length() - 1
            s &= s - 1
            best = max(best, 1 + f(S & out[v]))
        return best
    return f((1 << n) - 1)


def part3():
    print("== (3) Paley Ramsey rows ==")
    for q in (5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101):
        w = paley_graph_omega(q)
        print(" Paley graph P_%d: clique = independence number %d -> R(%d,%d) > %d" % (q, w, w + 1, w + 1, q))
    for q in (7, 11, 19, 23, 31):
        sq = {(x * x) % q for x in range(1, q)}
        out = [0] * q
        for i in range(q):
            for j in range(q):
                if i != j and ((j - i) % q) in sq:
                    out[i] |= 1 << j
        t = trans_tournament(q, out)
        print(" Paley tournament P_%d: largest transitive subtournament %d -> R(%d) > %d (tournament Ramsey)" % (q, t, t + 1, q))
    print(" known: R(3,3) = 6, R(4,4) = 18, 43 <= R(5,5) <= 46, 102 <= R(6,6); tournament R(2..6) = 2, 4, 8, 14, 28, 34 <= R(7) <= 47 (THM-455)")


def seidel_double(S):
    n = len(S)
    D = [[0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            D[i][j] = S[i][j]
            D[i][j + n] = S[i][j] + (1 if i == j else 0)
            D[i + n][j] = S[i][j] + (1 if i == j else 0)
            D[i + n][j + n] = -S[i][j]
    return D


def omega_alpha(D):
    n = len(D)
    adj = [0] * n; nadj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                if D[i][j] == -1:
                    adj[i] |= 1 << j
                else:
                    nadj[i] |= 1 << j

    def clique(N):
        best = [0]

        def bk(R, P, X):
            if P == 0 and X == 0:
                best[0] = max(best[0], R); return
            if R + bin(P).count("1") <= best[0]:
                return
            cand = P
            while cand:
                v = (cand & -cand).bit_length() - 1
                cand &= cand - 1
                bk(R + 1, P & N[v], X & N[v])
                P &= ~(1 << v); X |= 1 << v
        bk(0, (1 << n) - 1, 0)
        return best[0]
    return clique(adj), clique(nadj)


def part4():
    print("== (4) the Seidel doubling D(G): omega and alpha laws ==")
    for n in (5, 6):
        dist_w = {}; dist_a = {}; dist_m = {}
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            S = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                S[i][j] = S[j][i] = -1 if (mask >> k) & 1 else 1
            w, a = omega_alpha(S)
            W, A = omega_alpha(seidel_double(S))
            dist_w[W - w] = dist_w.get(W - w, 0) + 1
            dist_a[A - a] = dist_a.get(A - a, 0) + 1
            m = max(W, A) - max(w, a)
            dist_m[m] = dist_m.get(m, 0) + 1
        print(" all %d labeled graphs on %d vertices: omega(D) - omega: %s; alpha(D) - alpha: %s; max(omega,alpha) increment: %s" % (1 << len(pairs), n, dict(sorted(dist_w.items())), dict(sorted(dist_a.items())), dict(sorted(dist_m.items()))))
    for q in (5, 13, 17):
        sq = {(x * x) % q for x in range(1, q)}
        S = [[0] * q for _ in range(q)]
        for i in range(q):
            for j in range(q):
                if i != j:
                    S[i][j] = -1 if ((j - i) % q) in sq else 1
        w, a = omega_alpha(S)
        W, A = omega_alpha(seidel_double(S))
        print(" Paley P_%d: (omega, alpha) = (%d, %d); D(P_%d) on %d vertices: (omega, alpha) = (%d, %d) -> R(%d,%d) > %d from the double" % (q, w, a, q, 2 * q, W, A, max(W, A) + 1, max(W, A) + 1, 2 * q))


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def lis(seq):
    import bisect
    tails = []
    for x in seq:
        i = bisect.bisect_left(tails, x)
        if i == len(tails):
            tails.append(x)
        else:
            tails[i] = x
    return len(tails)


def part5():
    print("== (5) the orbit as a permutation tournament: Erdos-Szekeres ==")
    for n in (27, 703, 6171, 77031):
        vals = [n]; m = n
        while m != 1:
            m, _ = U(m); vals.append(m)
        L = len(vals); I = lis(vals); Dd = lis([-x for x in vals])
        r = math.isqrt(L - 1) + 1
        print(" n = %d: %d odd values; LIS = %d, LDS = %d (largest transitive subtournament of the permutation tournament = %d); Erdos-Szekeres guarantees a monotone subsequence of length %d" % (n, L, I, Dd, max(I, Dd), r))


def split_numbers(S):
    # G from the Seidel matrix; returns (largest |K|+|I| with K clique, I independent, disjoint, K x I all edges),
    #                                  (largest |I|+|K| with I independent, K clique, disjoint, no edges between)
    n = len(S)
    adj = [set(j for j in range(n) if j != i and S[i][j] == -1) for i in range(n)]
    def is_clique(T):
        return all(b in adj[a] for a, b in itertools.combinations(T, 2))
    def is_indep(T):
        return all(b not in adj[a] for a, b in itertools.combinations(T, 2))
    cliques = [T for r in range(0, n + 1) for T in itertools.combinations(range(n), r) if is_clique(T)]
    indeps = [T for r in range(0, n + 1) for T in itertools.combinations(range(n), r) if is_indep(T)]
    best_join = 0; best_sep = 0
    for K in cliques:
        Ks = set(K)
        for I in indeps:
            if Ks & set(I):
                continue
            if all(b in adj[a] for a in K for b in I):
                best_join = max(best_join, len(K) + len(I))
            if all(b not in adj[a] for a in K for b in I):
                best_sep = max(best_sep, len(K) + len(I))
    return best_join, best_sep


def part6():
    print("== (6) the exact law for the Seidel doubling: omega(D(G)) = largest induced complete split subgraph, alpha(D(G)) = max(alpha + 1, largest induced clique + independent set with no edges between) ==")
    import random
    ok = True; checked = 0
    for n in (3, 4, 5, 6):
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            S = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                S[i][j] = S[j][i] = -1 if (mask >> k) & 1 else 1
            W, A = omega_alpha(seidel_double(S))
            w, a = omega_alpha(S)
            J, Sep = split_numbers(S)
            ok &= (W == J) and (A == max(a + 1, Sep))
            checked += 1
    print(" law verified on all labeled graphs with 3 <= n <= 6 (%d graphs): %s" % (checked, ok))
    random.seed(9); ok2 = True
    for n, trials in ((7, 300), (8, 120)):
        for _ in range(trials):
            S = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    S[i][j] = S[j][i] = -1 if random.random() < 0.5 else 1
            W, A = omega_alpha(seidel_double(S)); w, a = omega_alpha(S); J, Sep = split_numbers(S)
            ok2 &= (W == J) and (A == max(a + 1, Sep))
    print(" law verified on 300 random graphs with n = 7 and 120 with n = 8: %s" % ok2)
    print(" proof sketch: a clique of D(G) is a clique K of G on the plain copy plus an independent set I of G on the prime copy (the prime copy carries the complement), twins never adjacent, cross pairs adjacent iff adjacent in G; an independent set of D(G) is an independent I on the plain copy plus a clique K on the prime copy, cross pairs non-adjacent iff non-adjacent in G, twins allowed, which forces |K| = 1 whenever K meets I")
    # Paley graphs: split numbers by clique enumeration in the common neighbourhood
    for q in (5, 13, 17, 29, 37):
        sq = {(x * x) % q for x in range(1, q)}
        adj = [set(j for j in range(q) if j != i and ((j - i) % q) in sq) for i in range(q)]
        # complete split: for each clique K (by size), maximum independent set inside the common neighbourhood
        def max_indep(vertices):
            vs = sorted(vertices); best = 0
            def rec(cands, size):
                nonlocal best
                if size + len(cands) <= best:
                    return
                if not cands:
                    best = max(best, size); return
                v = cands[0]
                rec([u for u in cands[1:] if u not in adj[v]], size + 1)
                rec(cands[1:], size)
            rec(vs, 0)
            return best
        def max_clique(vertices):
            vs = sorted(vertices); best = 0
            def rec(cands, size):
                nonlocal best
                if size + len(cands) <= best:
                    return
                if not cands:
                    best = max(best, size); return
                v = cands[0]
                rec([u for u in cands[1:] if u in adj[v]], size + 1)
                rec(cands[1:], size)
            rec(vs, 0)
            return best
        best_join = 0; best_sep = 0
        # enumerate cliques K up to the clique number (<= 5 here)
        def cliques_from(cands, cur):
            yield cur
            for idx, v in enumerate(cands):
                yield from cliques_from([u for u in cands[idx + 1:] if u in adj[v]], cur + [v])
        for K in cliques_from(list(range(q)), []):
            common = set(range(q)) - set(K)
            for v in K:
                common &= adj[v]
            best_join = max(best_join, len(K) + max_indep(common))
            far = set(range(q)) - set(K)
            for v in K:
                far -= adj[v]
            best_sep = max(best_sep, len(K) + max_indep(far))
        a = max_indep(range(q))
        print(" Paley P_%d: complete-split number %d = omega(D(P_%d)); separated number %d, alpha(D(P_%d)) = max(alpha + 1, sep) = %d -> the double on %d vertices gives R(%d,%d) > %d" % (q, best_join, q, best_sep, q, max(a + 1, best_sep), 2 * q, max(best_join, a + 1, best_sep) + 1, max(best_join, a + 1, best_sep) + 1, 2 * q))


def main():
    part1(); part2(); part3(); part4(); part5(); part6()


if __name__ == '__main__':
    main()
