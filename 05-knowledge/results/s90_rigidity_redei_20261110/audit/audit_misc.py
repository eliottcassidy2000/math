#!/usr/bin/env python3
"""Small independent checks of numbers and side claims in the S90 packet (pure Python, own code).

(a) An explicit ACTIVE additive exchange cube of dimension three (order-7 core), found by audit_exchange.py: the
    eight closure counts recomputed here by a subset dynamic programme (a different routine from audit_exchange's
    depth-first count).  Together with the exhaustive order-6 result this bounds what "additive cubes exist up to
    dimension n_-(W)" can mean.
(b) Statement 4 of THM-5124 as a special case of the inherited closure identity: S(A) is the protected closure of A
    with the full cut, so H(S(A)) = F(1) = h + 2 a.1 - (1/2) 1^T W 1 = H(A) + D(A); checked on every labelled A of
    order 1..5 by computing a and W from their combinatorial definitions (Session246 switches/PROOF.md section 2).
(c) Remark 4 constants: Busch floor f(m), f(n + 2) >= 5^((n+1)/3), 2^(floor(n/2)+1)/5^((n+1)/3) < 1.17 * 0.83^n,
    exponent (3/2) log_5 2 = 0.646; and the exact ratio 2^(floor(n/2)+1)/f(n + 2) for small n.
(d) README section 3 item 4: P(overlap >= 4) for two 9-subsets of the 50 odd numbers below 100 (hypergeometric);
    the nine members below 100 of each list; the shared members.
(e) Gauss: r_3(d) = 24 for d = 11, 19, 43, 67, 163 (brute force), r_3(3) = 8, and r_3(d) = 24 h(-d) for the
    d = 3 mod 8 squarefree below 400 against a brute-force class number (reduced forms).
(f) (q + 1)/4 prime for the Heegner q >= 7, and 1 for q = 3; x^2 + y^2 + 10 z^2: odd numbers below 3000 not
    represented (brute force).
"""
import itertools
import math
from math import comb


def hp_dp(out):
    n = len(out)
    cnt = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        cnt[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = cnt[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = out[v] & ~mask
                while nxt:
                    low = nxt & -nxt
                    nxt ^= low
                    cnt[mask | low][low.bit_length() - 1] += c
    return cnt


def H(out):
    n = len(out)
    if n == 0:
        return 1
    return sum(hp_dp(out)[(1 << n) - 1])


def closure(C, x):
    n = len(C)
    adj = [C[i] | (1 << n) | (0 if x >> i & 1 else 1 << (n + 1)) for i in range(n)]
    adj.append(1 << (n + 1))            # sink v -> r
    adj.append(x & ((1 << n) - 1))      # root r -> selected
    return adj


def is_strong(adj):
    n = len(adj)
    full = (1 << n) - 1

    def reach(a):
        seen = 1
        fr = 1
        while fr:
            low = fr & -fr
            fr ^= low
            new = a[low.bit_length() - 1] & ~seen
            seen |= new
            fr |= new
        return seen
    rev = [0] * n
    for i in range(n):
        for j in range(n):
            if adj[i] >> j & 1:
                rev[j] |= 1 << i
    return reach(adj) == full and reach(rev) == full


def main():
    ok = True
    # ---- (a)
    C = [68, 49, 18, 103, 105, 69, 6]
    n = 7
    assert all(((C[i] >> j) & 1) + ((C[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n))
    x0 = 0b0010111
    ex = [(0, 3), (1, 6), (2, 5)]
    vals = []
    strong = True
    for b in range(8):
        x = x0
        for e in range(3):
            if b >> e & 1:
                p, q = ex[e]
                assert x >> p & 1 and not x >> q & 1
                x ^= (1 << p) | (1 << q)
        adj = closure(C, x)
        vals.append(H(adj))
        strong &= is_strong(adj)
    g = [vals[1] - vals[0], vals[2] - vals[0], vals[4] - vals[0]]
    affine = all(vals[b] == vals[0] + sum(g[e] for e in range(3) if b >> e & 1) for b in range(8))
    print(f"(a) core {C} (H = {H(C)}), x0 = {x0:07b}, exchanges {ex}: eight counts {vals}; increments {g}; affine: {affine}; "
          f"all closures strong: {strong}; all increments nonzero: {all(g)}")
    ok &= affine and all(g)

    # ---- (b)
    bad = 0
    tot = 0
    for k in range(1, 6):
        prs = [(i, j) for j in range(1, k) for i in range(j)]
        for code in range(1 << len(prs)):
            A = [0] * k
            for t, (i, j) in enumerate(prs):
                if code >> t & 1:
                    A[i] |= 1 << j
                else:
                    A[j] |= 1 << i
            tot += 1
            cnt = hp_dp(A)
            rev = [0] * k
            for i in range(k):
                for j in range(k):
                    if A[i] >> j & 1:
                        rev[j] |= 1 << i
            cntr = hp_dp(rev)       # cntr[mask][v] = paths of A[mask] STARTING at v
            full = (1 << k) - 1
            Hs = [1] + [sum(cnt[m]) for m in range(1, 1 << k)]
            h = Hs[full]
            D = sum(Hs[U] * Hs[full ^ U] for U in range(1 << k))
            # a_i: ordered covers (P1, P2), P1 possibly empty, P2 starting at i
            a = [sum(Hs[U] * cntr[full ^ U][i] for U in range(1 << k) if not U >> i & 1) for i in range(k)]
            # W_ij: covers with P1 ending at i and P2 starting at j, plus the reversed endpoint version
            W = [[0] * k for _ in range(k)]
            for i in range(k):
                for j in range(k):
                    if i != j:
                        W[i][j] = (
                            sum(cnt[U][i] * cntr[full ^ U][j] for U in range(1 << k) if U >> i & 1 and not U >> j & 1)
                            + sum(cnt[U][j] * cntr[full ^ U][i] for U in range(1 << k) if U >> j & 1 and not U >> i & 1))
            F1 = h + 2 * sum(a) - sum(W[i][j] for i in range(k) for j in range(i + 1, k))
            S = closure(A, full)           # y0 = root r (beats all of A), x0 = sink v, x0 -> y0
            hS = H(S)
            if not (hS == F1 == h + D):
                bad += 1
            # the quadratic identity at every cut, as a check of a and W
            for x in range(1 << k):
                pred = h + 2 * sum(a[i] for i in range(k) if x >> i & 1) - sum(
                    W[i][j] for i in range(k) for j in range(i + 1, k) if x >> i & 1 and x >> j & 1)
                if H(closure(A, x)) != pred:
                    bad += 1
    print(f"(b) every labelled A of order 1..5 ({tot}): H(S(A)) = F(full cut) = h + 2 a.1 - sum_(i<j) W_ij = H(A) + D(A), "
          f"and F(x) = h + 2 a.x - (1/2) x^T W x at every cut: failures {bad}", "OK" if bad == 0 else "FAIL")
    ok &= bad == 0

    # ---- (c)
    def busch(m):
        if m % 3 == 0:
            return 3 * 5 ** ((m - 3) // 3)
        if m % 3 == 1:
            return 5 ** ((m - 1) // 3)
        return 9 * 5 ** ((m - 5) // 3)
    print("(c) Busch floor f(3..9):", [busch(m) for m in range(3, 10)])
    good = [busch(m) for m in range(3, 10)] == [3, 5, 9, 15, 25, 45, 75]
    good &= all(busch(m) ** 3 >= 5 ** (m - 1) for m in range(3, 200))          # exact: f(m) >= 5^((m-1)/3)
    cst = 2 / 5 ** (1 / 3)
    base = 2 ** 0.5 / 5 ** (1 / 3)
    print(f"    2/5^(1/3) = {cst:.4f} (< 1.17), sqrt2/5^(1/3) = {base:.4f} (< 0.83), (3/2) log_5 2 = {1.5 * math.log(2, 5):.4f}")
    good &= cst < 1.17 and base < 0.83 and abs(1.5 * math.log(2, 5) - 0.646) < 1e-3
    good &= all(2 ** (n // 2 + 1) / 5 ** ((n + 1) / 3) < 1.17 * 0.83 ** n for n in range(1, 400))
    print("    n, floor(n/2), f(n+2), 2^(floor(n/2)+1)/f(n+2):",
          [(n, n // 2, busch(n + 2), round(2 ** (n // 2 + 1) / busch(n + 2), 3)) for n in (4, 6, 8, 10, 12, 16, 20)])
    print("    constants of Remark 4:", "OK" if good else "FAIL")
    ok &= good

    # ---- (d)
    p = sum(comb(9, j) * comb(41, 9 - j) for j in range(4, 10)) / comb(50, 9)
    ram = [3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719]
    ozj = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637,
           821, 1067, 1339, 1861, 2647]
    print(f"(d) P(overlap >= 4) = {p:.4f}; members below 100: Ramanujan {len([x for x in ram if x < 100])}, O/Z/J "
          f"{len([x for x in ozj if x < 100])}; shared members: {sorted(set(ram) & set(ozj))}")
    ok &= abs(p - 0.043) < 0.002 and sorted(set(ram) & set(ozj)) == [7, 21, 43, 67]

    # ---- (e)
    def r3(d):
        c = 0
        m = int(d ** 0.5) + 1
        for x in range(-m, m + 1):
            for y in range(-m, m + 1):
                z2 = d - x * x - y * y
                if z2 < 0:
                    continue
                z = int(round(z2 ** 0.5))
                if z * z == z2:
                    c += 1 if z == 0 else 2
        return c

    def classno(D):          # D < 0 discriminant; number of reduced primitive forms
        h = 0
        b = D % 2
        while 3 * b * b <= -D:
            ac = (b * b - D) // 4
            a = max(b, 1)
            while a * a <= ac:
                if ac % a == 0:
                    c = ac // a
                    if math.gcd(math.gcd(a, b), c) == 1:
                        h += 1 if (b == 0 or a == b or a == c) else 2
                a += 1
            b += 2
        return h
    vals = {d: r3(d) for d in (3, 11, 19, 43, 67, 163)}
    bad = 0
    tested = 0
    for d in range(11, 400, 8):
        if all(d % (q * q) for q in range(2, 20)):
            tested += 1
            if r3(d) != 24 * classno(-d):
                bad += 1
    print(f"(e) r_3: {vals}; r_3(d) = 24 h(-d) for the {tested} squarefree d = 3 mod 8, 11 <= d < 400: failures {bad}",
          "OK" if bad == 0 and vals == {3: 8, 11: 24, 19: 24, 43: 24, 67: 24, 163: 24} else "FAIL")
    ok &= bad == 0 and vals == {3: 8, 11: 24, 19: 24, 43: 24, 67: 24, 163: 24}

    # ---- (f)
    def isprime(m):
        return m > 1 and all(m % q for q in range(2, int(m ** 0.5) + 1))
    print("(f) (q+1)/4 for q = 3, 7, 11, 19, 43, 67, 163:", [((q + 1) // 4, isprime((q + 1) // 4)) for q in (3, 7, 11, 19, 43, 67, 163)])
    lim = 3000
    rep = set()
    for x in range(0, 56):
        for y in range(x, 56):
            for z in range(0, 18):
                v = x * x + y * y + 10 * z * z
                if v < lim:
                    rep.add(v)
    miss = [m for m in range(1, lim, 2) if m not in rep]
    print("    odd numbers below 3000 not of the form x^2 + y^2 + 10 z^2:", miss)
    ok &= miss == ram
    print("AUDIT MISC:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
