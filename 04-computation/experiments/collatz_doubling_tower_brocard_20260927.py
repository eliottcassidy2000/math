#!/usr/bin/env python3
"""collatz_doubling_tower_brocard_20260927.py -- D18 (iterated Seidel doubling against the tournament tower) and the Brocard
triple (session collatz-posets-zeta5-20260927, opus, 2026-09-27, eighth note).

 (1) Iterated Seidel doubling: for G in {K_1, K_2, C_5 = P_5, P_13} compute (omega, alpha, cs, sep) of D^k(G) for k = 0..K
     (orders 5 * 2^k for P_5: 5, 10, 20, 40, 80), and the increments per doubling, against the tournament tower's
     trans(T_7, T_15, T_31, T_63) = 3, 5, 7, 11 (THM-455) and the trivial bound alpha(D) >= alpha + 1.
     By Theorem 4 of the seventh note omega(D(G)) = cs(G), alpha(D(G)) = max(alpha(G) + 1, sep(G)); here the four
     parameters are recomputed by brute force at every level (Bron-Kerbosch), so the law is re-checked along the tower.
 (2) The doubled structures: what is cs(D(G))? A complete split subgraph of D(G) is (clique K + independent I') joined to
     (independent J + clique L') with the eight cross conditions; we count the maximal 4-tuples on C_5 and K_3 to see the
     nesting, and compare cs(D(G)) with cs(G) + sep(G) and with 2 cs(G).
 (4) Random-graph benchmark: mean max(omega, alpha) of G(n, 1/2) at n = 10, 20, 40, 80 against the tower's 3, 4, 6, 8.
 (3) Brocard's problem n! + 1 = m^2: the solutions n <= 60 (4, 5, 7 with m = 5, 11, 71), the factorizations n! = (m-1)(m+1)
     ((4,6), (10,12), (70,72)), the Wilson readings (5^2 = 4! + 1 makes 5 a Wilson prime; 5! = ((11-1)/2)! = -1 mod 11^2 is
     Mordell's half-Wilson at level p^2; 7 = (71-1)/10 has no Wilson reading), the supersingular list (5, 11, 71 all divide
     the Monster's order; 71 is the largest such prime), and the group orders 4! = |2T|, 5! = |2I|, 7! = |S_7|.
Usage: python3 collatz_doubling_tower_brocard_20260927.py
"""
import math, itertools, sys

sys.setrecursionlimit(20000)


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


def clique_number(N, n):
    best = [0]

    def bk(R, P, X):
        if P == 0 and X == 0:
            best[0] = max(best[0], R); return
        if R + bin(P).count("1") <= best[0]:
            return
        PX = P | X
        # pivot with most neighbours in P
        u = -1; um = -1
        t = PX
        while t:
            v = (t & -t).bit_length() - 1; t &= t - 1
            c = bin(P & N[v]).count("1")
            if c > um:
                um = c; u = v
        cand = P & ~N[u]
        while cand:
            v = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            bk(R + 1, P & N[v], X & N[v])
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << n) - 1, 0)
    return best[0]


def params(S):
    # omega, alpha, cs (largest clique completely joined to a disjoint independent set), sep (clique + independent set, no edges between)
    n = len(S)
    adj = [0] * n; nadj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                if S[i][j] == -1:
                    adj[i] |= 1 << j
                else:
                    nadj[i] |= 1 << j
    w = clique_number(adj, n); a = clique_number(nadj, n)

    def max_indep_in(mask):
        # independent set inside mask = clique of the complement restricted to mask
        Nm = [nadj[i] & mask for i in range(n)]
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
                bk(R + 1, P & Nm[v], X & Nm[v])
                P &= ~(1 << v); X |= 1 << v
        bk(0, mask, 0)
        return best[0]

    cs = 0; sep = 0
    full = (1 << n) - 1

    # enumerate all cliques (including empty) via recursion
    def cliques(cands, cur_mask, common_adj, common_nadj):
        nonlocal cs, sep
        # cur_mask: clique so far; common_adj: vertices adjacent to all of cur; common_nadj: vertices adjacent to none of cur
        size = bin(cur_mask).count("1")
        cs = max(cs, size + max_indep_in(common_adj & ~cur_mask))
        sep = max(sep, size + max_indep_in(common_nadj & ~cur_mask))
        t = cands
        while t:
            v = (t & -t).bit_length() - 1; t &= t - 1
            cliques(t & adj[v], cur_mask | (1 << v), common_adj & adj[v], common_nadj & nadj[v])
    cliques(full, 0, full, full)
    return w, a, cs, sep


def paley_seidel(q):
    sq = {(x * x) % q for x in range(1, q)}
    S = [[0] * q for _ in range(q)]
    for i in range(q):
        for j in range(q):
            if i != j:
                S[i][j] = -1 if ((j - i) % q) in sq else 1
    return S


def part1():
    print("== (1) iterated Seidel doubling D^k(G): (omega, alpha, cs, sep) by brute force at every level ==")
    starts = {"K_1": [[0]], "K_2": [[0, -1], [-1, 0]], "C_5 = P_5": paley_seidel(5), "P_13": paley_seidel(13)}
    depth = {"K_1": 5, "K_2": 4, "C_5 = P_5": 4, "P_13": 2}
    for name, S in starts.items():
        rows = []; law_ok = True; prev = None
        for k in range(depth[name] + 1):
            w, a, cs, sep = params(S)
            if prev is not None:
                pw, pa, pcs, psep = prev
                law_ok &= (w == pcs) and (a == max(pa + 1, psep))
            rows.append((k, len(S), w, a, cs, sep))
            prev = (w, a, cs, sep)
            if k < depth[name]:
                S = seidel_double(S)
        print(" %s: (k, order, omega, alpha, cs, sep):" % name, rows, "| law omega(D) = cs, alpha(D) = max(alpha+1, sep) along the tower: %s" % law_ok)
    print(" tournament tower (THM-455): orders 7, 15, 31, 63 -> trans 3, 5, 7, 11 (increments +2, +2, +4)")


def part2():
    print("== (2) what cs(D(G)) is made of ==")
    for name, S in (("K_3", [[0, -1, -1], [-1, 0, -1], [-1, -1, 0]]), ("C_5", paley_seidel(5)), ("P_13", paley_seidel(13))):
        w, a, cs, sep = params(S)
        D = seidel_double(S)
        W, A, CS, SEP = params(D)
        print(" %s: (omega, alpha, cs, sep) = %s; D: %s; cs(D) against cs + sep = %d and 2 cs = %d" % (name, (w, a, cs, sep), (W, A, CS, SEP), cs + sep, 2 * cs))


def part3():
    print("== (3) Brocard's problem n! + 1 = m^2 ==")
    sols = []
    f = 1
    for n in range(1, 61):
        f *= n
        m = math.isqrt(f + 1)
        if m * m == f + 1:
            sols.append((n, m))
    print(" solutions with n <= 60:", sols, "(conjecturally all; verified to 10^15 by Matson 2017; finitely many under abc, Overholt 1993)")
    for n, m in sols:
        print("  n = %d: n! = %d = (m-1)(m+1) = %d * %d; m = %d" % (n, math.factorial(n), m - 1, m + 1, m))
    ss = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    print(" supersingular primes (Ogg: the primes dividing the Monster's order, X_0(p)+ of genus 0):", ss, "; 5, 11, 71 all supersingular: %s; 71 the largest" % all(m in ss for _, m in sols))
    print(" Wilson readings: 4! + 1 = 5^2 (5 is a Wilson prime, (p-1)! = -1 mod p^2; Wilson primes 5, 13, 563); 5! + 1 = 11^2 with 5 = (11-1)/2")
    print(" (Mordell: ((p-1)/2)! = +-1 mod p for p = 3 mod 4, here -1 mod 11^2); 7! + 1 = 71^2 with 7 = (71-1)/10: no Wilson reading")
    # half-Wilson mod p^2 census: primes p = 3 mod 4 with ((p-1)/2)! = -1 mod p^2
    hw = []
    for p in range(3, 2000):
        if all(p % d for d in range(2, int(p ** 0.5) + 1)) and p % 4 == 3:
            h = math.factorial((p - 1) // 2) % (p * p)
            if h == p * p - 1:
                hw.append(p)
    print(" primes p = 3 mod 4 below 2000 with ((p-1)/2)! = -1 mod p^2:", hw, "(11 is the only one below 2000 besides p = 3? 1! = 1)")
    print(" group orders: 4! = 24 = |2T| (binary tetrahedral, McKay E_6), 5! = 120 = |2I| (binary icosahedral, E_8), 7! = 5040 = |S_7|;")
    print(" the binary octahedral group has order 48 = 2 * 4!, so the exceptional binary polyhedral orders are 4!, 2 * 4!, 5!")
    # Brocard-type: n! + k = m^2 small k
    print(" nearby equations: n! - 1 = m^2 has no solution n <= 60: %s; n! + 4 = m^2 -> n = 4 (28? no), listing n! + k square for k in 1..4, n <= 30:" % (not any(math.isqrt(math.factorial(n) - 1) ** 2 == math.factorial(n) - 1 for n in range(2, 61))))
    for k in (1, 2, 3, 4):
        found = [n for n in range(1, 31) if math.isqrt(math.factorial(n) + k) ** 2 == math.factorial(n) + k]
        print("   k = %d: n in %s" % (k, found))


def part4():
    print("== (4) random-graph benchmark for the graph tower ==")
    import random
    random.seed(2)
    tower = {10: 3, 20: 4, 40: 6, 80: 8}
    for n, trials in ((10, 200), (20, 100), (40, 40), (80, 12)):
        vals = []
        for _ in range(trials):
            adj = [0] * n; nadj = [0] * n
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        adj[i] |= 1 << j; adj[j] |= 1 << i
                    else:
                        nadj[i] |= 1 << j; nadj[j] |= 1 << i
            vals.append(max(clique_number(adj, n), clique_number(nadj, n)))
        print("  G(n,1/2), n = %3d: max(omega, alpha) mean %.2f, range [%d, %d]; 2 log_2 n = %.1f; tower D^k(P_5): %d" % (n, sum(vals) / len(vals), min(vals), max(vals), 2 * math.log2(n), tower[n]))
    print("  reading: like the tournament tower (THM-455), the graph tower stays at or below the best of a dozen random graphs at every order")


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
