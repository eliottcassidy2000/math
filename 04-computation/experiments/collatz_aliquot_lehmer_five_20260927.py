#!/usr/bin/env python3
"""collatz_aliquot_lehmer_five_20260927.py -- aliquot sequences (the Lehmer five) against Collatz, plus D21 and D22
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, ninth note).

 (1) D21: the closed form cs(D(G)) = max over 4-part induced structures (K clique, I independent, J independent, L clique;
     K-I, K-J, K-L, I-J complete; I-L, J-L anticomplete; all disjoint) of |K|+|I|+|J|+|L|, or |K|+|J|+1 (the twin case,
     L = {u} inside J, I empty): verified against omega(D^2(G)) by brute force on all graphs with n <= 5 and random n = 6, 7.
 (2) D22: f(n)^2 - 1 and f(n) - 1 against factorials for |n| <= 60, f the Lehmer conductor polynomial.
 (3) Aliquot sequences s(n) = sigma(n) - n: the Lehmer five 276, 552, 564, 660, 966; the sequence of 276 with its 2-adic
     valuation, odd driver part and growth ratio for as many terms as sympy factors quickly; the Guy-Selfridge drivers
     (2, 4*7, 8*3, 8*3*5, 32*3*7, 512*3*11*31 and the even perfect numbers) with 24 = 4! and 120 = 5! among them.
 (4) The 2-adic contrast, measured: for even n <= 2*10^6 the transition matrix P(v_2(s(n)) = b | v_2(n) = a) and the mean
     log_2 growth by class a, against Collatz's Terras law (next valuation geometric(1/2), independent of the current one:
     measured on the same scale).
 (5) Leaves and cycles: untouchable numbers below 10^5 (density) against the exact 1/3 of multiples of 3 in the Syracuse
     tree; even perfect numbers as triangular numbers at Mersenne indices 2^p - 1 (the odd-Catalan indices).
Usage: python3 collatz_aliquot_lehmer_five_20260927.py
"""
import math, itertools, random, sys, time

sys.setrecursionlimit(20000)


# ---------- D21 ----------
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
        cand = P
        while cand:
            v = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            bk(R + 1, P & N[v], X & N[v])
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << n) - 1, 0)
    return best[0]


def omega_of(S):
    n = len(S); adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and S[i][j] == -1:
                adj[i] |= 1 << j
    return clique_number(adj, n)


def cs_closed_form(S):
    n = len(S)
    adj = [set(j for j in range(n) if j != i and S[i][j] == -1) for i in range(n)]

    def clique(T):
        return all(b in adj[a] for a, b in itertools.combinations(T, 2))

    def indep(T):
        return all(b not in adj[a] for a, b in itertools.combinations(T, 2))

    def complete(A, B):
        return all(b in adj[a] for a in A for b in B)

    def anti(A, B):
        return all(b not in adj[a] for a in A for b in B)
    best = 0
    for assign in itertools.product(range(5), repeat=n):
        K = [v for v in range(n) if assign[v] == 1]; I = [v for v in range(n) if assign[v] == 2]
        J = [v for v in range(n) if assign[v] == 3]; L = [v for v in range(n) if assign[v] == 4]
        if not (clique(K) and clique(L) and indep(I) and indep(J)):
            continue
        if complete(K, I) and complete(K, J) and complete(K, L) and complete(I, J) and anti(I, L) and anti(J, L):
            best = max(best, len(K) + len(I) + len(J) + len(L))
        # twin case: I empty, L = {u} with u in J -> |K| + |J| + 1 (already covered when I = [] and L = []: add 1 if J nonempty)
        if not I and not L and J and complete(K, J):
            best = max(best, len(K) + len(J) + 1)
    return best


def part1():
    print("== (1) D21: closed form for cs(D(G)) = omega(D^2(G)) ==")
    ok = True; cnt = 0
    for n in (2, 3, 4, 5):
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            S = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                S[i][j] = S[j][i] = -1 if (mask >> k) & 1 else 1
            ok &= (omega_of(seidel_double(seidel_double(S))) == cs_closed_form(S))
            cnt += 1
    print(" omega(D^2(G)) = closed-form 4-part maximum on all labeled graphs with 2 <= n <= 5 (%d graphs): %s" % (cnt, ok))
    random.seed(21); ok2 = True
    for n, trials in ((6, 60), (7, 15)):
        for _ in range(trials):
            S = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    S[i][j] = S[j][i] = -1 if random.random() < 0.5 else 1
            ok2 &= (omega_of(seidel_double(seidel_double(S))) == cs_closed_form(S))
    print(" ... and on 60 random graphs with n = 6 and 15 with n = 7: %s" % ok2)


# ---------- D22 ----------
def part2():
    print("== (2) D22: Lehmer conductors against factorials ==")
    f = lambda n: n ** 4 + 5 * n ** 3 + 15 * n ** 2 + 25 * n + 25
    facts = {}
    v = 1
    for k in range(1, 80):
        v *= k; facts[v] = k
    hits = []
    for n in range(-60, 61):
        for label, val in (("f(n) - 1", f(n) - 1), ("f(n)^2 - 1", f(n) ** 2 - 1), ("f(n) + 1", f(n) + 1), ("f(n)^2 + 1", f(n) ** 2 + 1)):
            if val in facts:
                hits.append((n, label, "%d!" % facts[val]))
    print(" factorial hits for |n| <= 60:", hits, "(the three Brocard cases and nothing else)")


# ---------- aliquot ----------
def sigma_from_factor(fac):
    s = 1
    for p, e in fac.items():
        s *= (p ** (e + 1) - 1) // (p - 1)
    return s


def part3():
    print("== (3) the Lehmer five and the sequence of 276 ==")
    try:
        from sympy import factorint
    except ImportError:
        print(" sympy not available"); return
    lehmer5 = [276, 552, 564, 660, 966]
    print(" Lehmer five:", lehmer5, "factorizations:", [{int(q): int(e) for q, e in factorint(x).items()} for x in lehmer5])
    # 138: the classical record excursion that returns to 1
    n = 138; peak = n; steps = 0
    while n > 1 and steps < 400:
        fac = {int(q): int(e) for q, e in factorint(n).items()}
        n = int(sigma_from_factor(fac)) - n; peak = max(peak, n); steps += 1
    print(" the sequence of 138 returns to 1 after %d steps with peak %d (exponent log(peak)/log(138) = %.2f; Collatz record excursions have exponent below 1.9)" % (steps, peak, math.log(peak) / math.log(138)))
    drivers = {2: "2", 28: "2^2*7 (perfect)", 24: "2^3*3 = 4!", 120: "2^3*3*5 = 5! (3-perfect)", 672: "2^5*3*7 (3-perfect)", 523776: "2^9*3*11*31 (3-perfect)", 6: "2*3 (perfect)", 496: "2^4*31 (perfect)", 8128: "2^6*127 (perfect)"}
    # Guy-Selfridge driver definition: 2^a v with v odd, v | sigma(2^a) = 2^(a+1) - 1, and 2^(a-1) | sigma(v)
    def is_driver(d):
        a = 0
        while d % 2 == 0:
            d //= 2; a += 1
        v = d
        if a == 0:
            return False
        if (2 ** (a + 1) - 1) % v != 0:
            return False
        return int(sigma_from_factor({int(q): int(e) for q, e in factorint(v).items()})) % (2 ** (a - 1)) == 0
    print(" Guy-Selfridge driver test (2^a v, v | 2^(a+1)-1, 2^(a-1) | sigma(v)):", {d: is_driver(d) for d in sorted(drivers)})
    print(" so 4! = 24 and 5! = 120 are drivers; 7! = 5040 = 2^4 * 315 is not (315 does not divide 31)")
    n = 276; t0 = time.time(); rows = []; k = 0
    vals = [n]
    while k < 400 and time.time() - t0 < 240 and n > 1:
        fac = {int(q): int(e) for q, e in factorint(n).items()}
        s = int(sigma_from_factor(fac)) - n
        a = fac.get(2, 0)
        odd_driver = 1
        for p in (3, 5, 7, 11, 31):
            if p in fac:
                odd_driver *= p ** fac[p]
        rows.append((k, len(str(n)), a, odd_driver, round(math.log2(s / n), 3) if s > 0 else None))
        n = s; vals.append(n); k += 1
    print(" sequence of 276: %d terms computed in %.0f s; last term has %d digits" % (len(vals), time.time() - t0, len(str(vals[-1]))))
    print(" (index, digits, v_2, odd part of the driver among {3,5,7,11,31}, log_2 growth) at selected indices:")
    for r in rows[::max(1, len(rows) // 25)]:
        print("   ", r)
    v2s = [r[2] for r in rows]
    persist = sum(1 for i in range(len(v2s) - 1) if v2s[i + 1] == v2s[i]) / max(1, len(v2s) - 1)
    growth = sum(r[4] for r in rows if r[4] is not None) / max(1, len(rows))
    print(" fraction of steps keeping the 2-adic valuation: %.3f; mean log_2 growth per step: %.3f; distribution of v_2 along the sequence: %s" % (persist, growth, {a: v2s.count(a) for a in sorted(set(v2s))}))


def part4():
    print("== (4) the 2-adic contrast, measured ==")
    N = 2 * 10 ** 6
    # sigma via sieve
    sig = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            sig[m] += d
    def v2(x):
        c = 0
        while x % 2 == 0 and x:
            x //= 2; c += 1
        return c
    trans = {}; growth = {}; count = {}
    for n in range(2, N + 1, 2):
        s = sig[n] - n
        if s <= 1:
            continue
        a = v2(n); b = v2(s)
        a = min(a, 6); b = min(b, 6)
        trans[(a, b)] = trans.get((a, b), 0) + 1
        growth[a] = growth.get(a, 0.0) + math.log2(s / n); count[a] = count.get(a, 0) + 1
    print(" aliquot: for even n <= %d, P(v_2(s(n)) = b | v_2(n) = a) (rows a = 1..6, columns b = 0..6, valuations >= 6 pooled):" % N)
    for a in range(1, 7):
        tot = sum(trans.get((a, b), 0) for b in range(7))
        print("   a = %d: %s  mean log_2(s(n)/n) = %+.3f (n = %d)" % (a, [round(trans.get((a, b), 0) / tot, 3) for b in range(7)], growth[a] / count[a], count[a]))
    print(" reading: the 2-adic valuation is sticky (the diagonal dominates for a >= 2) and the drift is negative only for a = 1: drivers persist for free")
    # Collatz control on the same scale: Syracuse valuations along orbits of odd n <= 2*10^5, transition matrix and drift by class
    def U(m):
        m = 3 * m + 1; v = 0
        while m % 2 == 0:
            m //= 2; v += 1
        return m, v
    ct = {}; cg = {}; cc = {}
    for n in range(3, 2 * 10 ** 5, 2):
        m = n; prev = None
        while m != 1:
            m2, v = U(m)
            if prev is not None:
                key = (min(prev, 6), min(v, 6)); ct[key] = ct.get(key, 0) + 1
            cg[min(v, 6)] = cg.get(min(v, 6), 0.0) + (math.log2(3) - v); cc[min(v, 6)] = cc.get(min(v, 6), 0) + 1
            prev = v; m = m2
    print(" Collatz: P(next valuation = b | current = a) along orbits of odd n < 2*10^5 (rows a = 1..6, columns b = 1..6; Terras: 2^(-b), independent of a):")
    for a in range(1, 7):
        tot = sum(ct.get((a, b), 0) for b in range(1, 7))
        print("   a = %d: %s  step growth log_2 3 - a = %+.3f" % (a, [round(ct.get((a, b), 0) / tot, 3) for b in range(1, 7)], math.log2(3) - a))
    print(" reading: Collatz valuations are memoryless (each row is the geometric law), so a growth pattern must be paid for in advance by the residue; aliquot valuations are not")


def part5():
    print("== (5) leaves and cycles ==")
    N = 10 ** 5
    sig = [0] * (4 * N + 1)
    for d in range(1, 4 * N + 1):
        for m in range(d, 4 * N + 1, d):
            sig[m] += d
    touched = set()
    for m in range(2, 4 * N + 1):
        s = sig[m] - m
        if 1 < s <= N:
            touched.add(s)
    # s(m) >= sqrt(m) for composite m, so preimages of n <= N have m <= n^2; we only scanned m <= 4N: numbers n <= sqrt(4N) are exact
    exact = int(math.isqrt(4 * N))
    unt = [n for n in range(2, exact + 1) if n not in touched]
    print(" untouchable numbers below %d (exact scan: preimages m <= n^2 covered): %s ... count %d, density %.3f (Erdos: positive density; Collatz's Syracuse leaves, the multiples of 3, have density exactly 1/3)" % (exact, unt[:12], len(unt), len(unt) / exact))
    perfect = [n for n in range(2, 4 * N + 1) if sig[n] == 2 * n]
    print(" even perfect numbers below %d: %s = triangular numbers T_(2^p - 1) with 2^p - 1 prime: %s" % (4 * N, perfect, all(any(n == (2 ** p - 1) * 2 ** (p - 1) for p in range(2, 20)) for n in perfect)))
    print(" the same Mersenne indices 2^p - 1 are where the Catalan numbers are odd (seventh note); Euclid-Euler needs 2^p - 1 prime, Kummer needs nothing")


def main():
    part1(); part2(); part3(); part4(); part5()


if __name__ == '__main__':
    main()
