#!/usr/bin/env python3
"""Auditor bridges-note, script A2 (independent code).

Section 6: signed-integer dictionary, favoured classes, head-series argument, fan primes, lines.
Section 5: packing polynomial of the F1 joins by brute force (directed odd cycles, disjoint packings),
           the digit claims for 17 and 45.
"""
import sys
from math import gcd
import numpy as np
from numba import njit


def psi_table(n):
    psi = [0] * n
    for h in range(1, n):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    return psi


def chi4(m):
    return 1 if m % 4 == 1 else -1


def n_of(m):
    return (chi4(m) * m - 1) // 4


def step(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def oddpart(n):
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    return n


def is_prime(p):
    if p < 2:
        return False
    i = 2
    while i * i <= p:
        if p % i == 0:
            return False
        i += 1
    return True


def ord2(p):
    k, v = 1, 2 % p
    while v != 1:
        v = (2 * v) % p
        k += 1
    return k


def section6():
    print("=== A2.1 dictionary (section 6) ===")
    R = range(-3000, 3001)
    ok_bij = sorted(abs(4 * n + 1) for n in R) == sorted(set(abs(4 * n + 1) for n in R))
    ok_inv = all(n_of(abs(4 * n + 1)) == n for n in R)
    print("n -> |4n+1| injective on [-3000,3000]:", ok_bij, "; inverse (chi_4(m) m - 1)/4:", ok_inv)
    ok_odd = all(4 * step(n) + 1 == 3 * (4 * n + 1) + 2 for n in R if n % 2)
    ok_even = all(2 * (4 * step(n) + 1) == (4 * n + 1) + 1 for n in R if n % 2 == 0)
    print("odd step m* -> 3 m* + 2:", ok_odd, "; even step m* -> (m* + 1)/2:", ok_even)
    ok_prod = all(abs(4 * n + 1) * abs(4 * k + 1) == abs(4 * (4 * n * k + n + k) + 1)
                  for n in range(-60, 61) for k in range(-60, 61))
    print("product of counts <-> n*k = 4nk+n+k:", ok_prod)
    ok_tri = all(3 * abs(4 * n + 1) == abs(4 * (-3 * n - 1) + 1) for n in R)
    print("tripling a count <-> n -> -3n-1:", ok_tri)
    ok_sheet = all((n >= 0) == (chi4(abs(4 * n + 1)) == 1) for n in R)
    print("sign of n = chi_4(m) (n >= 0 iff m = 1 mod 4):", ok_sheet)
    # note: -3n-1 = -(3n+1): tripling a count is the odd Collatz step followed by negation
    print("-3n-1 == -(3n+1) (odd step then negation), for all n: True by algebra")

    # favoured classes
    okA = okB = okC = True
    for p in [3, 5, 7, 11, 13, 17, 19, 23]:
        inv4 = pow(4, -1, p)
        inv2 = pow(2, -1, p)
        for n in R:
            N = abs(4 * n + 1)
            if (N % p == 0) != ((n - (-inv4)) % p == 0):
                okA = False
            if n >= 0:
                if ((N + 1) % p == 0) != ((n + inv2) % p == 0):
                    okB = False
            else:
                if ((N + 1) % p == 0) != (n % p == 0):
                    okC = False
    print("p | N iff n = -1/4 mod p:", okA, "; N = -1 iff n = -1/2 (n>=0):", okB, "; N = -1 iff n = 0 (n<0):", okC)
    okH = all(oddpart(abs(4 * n + 1) + 1) == (2 * n + 1 if n >= 0 else oddpart(n)) for n in R)
    print("odd part of N+1 = 2n+1 (n>=0), odd part of n (n<0):", okH)
    # tail sheet argument, for completeness: odd part of N = N itself (N odd)
    # fixed points: 3n+1 = n -> n = -1/2 ; on the negative side with x = -n: 3x-1 = x -> x = 1/2 -> n = -1/2 as well
    print("fixed point of n -> 3n+1 is n = -1/2 on both sides; the class n = 0 (negative side) is not that fixed point")

    print()
    print("=== A2.2 fan primes (section 6) ===")
    LIM = 1 << 22
    psi = psi_table(LIM)
    # balanced states x = y  <=>  psi(h) = h
    bal = [h for h in range(LIM) if psi[h] == h]
    print("balanced states psi(h) = h below 2^22:", bal)
    print("  all are 0 or powers of 2:", all(h == 0 or (h & (h - 1)) == 0 for h in bal))
    primes = [p for p in range(3, 200) if is_prime(p)]
    fan_first = []
    rows = []
    for p in primes:
        o = ord2(p)
        # first p-divisible state: p | gcd(2psi+1, 2h+1)
        first = None
        h = (p - 1) // 2
        while h < LIM:
            if (2 * psi[h] + 1) % p == 0:
                first = h
                break
            h += p
        if o % 2 == 0:
            t = o // 2
            assert pow(2, t, p) == p - 1
            # least t with p | 2^t + 1
            tmin = next(tt for tt in range(1, o + 1) if (pow(2, tt, p) + 1) % p == 0)
            assert tmin == t
            fan = 1 << (t - 1)
            is_first = (first == fan)
            if is_first:
                fan_first.append(p)
            rows.append((p, o, t, fan if t <= 40 else "2^%d" % (t - 1), first, is_first))
        else:
            assert all((pow(2, tt, p) + 1) % p for tt in range(1, o + 1))
            rows.append((p, o, None, None, first, None))
    print("p, ord_p(2), t = ord/2, fan h = 2^(t-1), first p-divisible state h, fan is first?")
    for r in rows:
        print("  ", r)
    print("primes p < 200 with ord_p(2) even:", [r[0] for r in rows if r[2] is not None])
    print("primes p < 200 with ord_p(2) odd (no balanced p-divisible state):", [r[0] for r in rows if r[2] is None])
    print("primes p < 200 whose first p-divisible state is the fan:", fan_first)
    print("first p-divisible state not found below 2^22 for:", [r[0] for r in rows if r[4] is None])

    print()
    print("=== A2.3 lines (section 6) ===")
    ok1 = all(oddpart(3 * n - 1) == oddpart(3 * n - 1) for n in range(1, 2001, 2))  # tautology: definition of T_-
    print("c_-(3n) = odd part of 3n-1 = T_-(n), c_+(3n) = odd part of 3n+1 = T_+(n): identities by definition")
    # the sentence: the join with a triangle moves a count along a Syracuse arc of ITS LINE LABEL ?
    tot = hit_m = hit_p = 0
    tot2 = hp = hm = 0
    for N in range(3, 20001, 2):
        c = oddpart(N - 1)
        c3 = oddpart(3 * N - 1)
        tot += 1
        if c3 == oddpart(3 * c + 1):
            hit_p += 1
        if c3 == oddpart(3 * c - 1):
            hit_m += 1
        d = oddpart(N + 1)
        d3 = oddpart(3 * N + 1)
        tot2 += 1
        if d3 == oddpart(3 * d + 1):
            hp += 1
        if d3 == oddpart(3 * d - 1):
            hm += 1
    print("odd N in [3,20000]: c_-(3N) = T_+(c_-(N)) for %d of %d; = T_-(c_-(N)) for %d" % (hit_p, tot, hit_m))
    print("                    c_+(3N) = T_+(c_+(N)) for %d of %d; = T_-(c_+(N)) for %d" % (hp, tot2, hm))
    # exact rule: N = 1 + c 2^s: 3N - 1 = 2 + 3c 2^s
    ex = [(N, oddpart(N - 1), oddpart(3 * N - 1), oddpart(3 * oddpart(N - 1) + 1)) for N in (3, 5, 7, 9, 11, 13, 15, 17, 19, 23)]
    print("examples (N, c_-(N), c_-(3N), T_+(c_-(N))):", ex)


# ------------------------------------------------------------------ packing polynomial by brute force
@njit(cache=True)
def packing_poly(adj, n):
    """coefficients alpha_k = number of sets of k pairwise vertex-disjoint directed odd cycles."""
    full = (1 << n) - 1
    hc = np.zeros(1 << n, np.int64)  # number of directed Hamiltonian cycles of T[S], S odd
    for s in range(n):
        # cycles whose lowest vertex is s; dp over masks of vertices > s
        m = n - s - 1
        if m < 2:
            continue
        dp = np.zeros((1 << m, m), np.int64)
        for v in range(m):
            if (adj[s] >> (s + 1 + v)) & 1:
                dp[1 << v, v] = 1
        for mask in range(1, 1 << m):
            pc = 0
            t = mask
            while t:
                pc += t & 1
                t >>= 1
            for v in range(m):
                c = dp[mask, v]
                if c == 0:
                    continue
                gv = s + 1 + v
                if pc % 2 == 0 and ((adj[gv] >> s) & 1):
                    hc[(mask << (s + 1)) | (1 << s)] += c
                for w in range(m):
                    if (mask >> w) & 1:
                        continue
                    if (adj[gv] >> (s + 1 + w)) & 1:
                        dp[mask | (1 << w), w] += c
    K = max(5, n // 3 + 2)
    P = np.zeros((1 << n, K), np.int64)
    P[0, 0] = 1
    for mask in range(1, full + 1):
        low = mask & (-mask)
        rest = mask ^ low
        for k in range(K):
            P[mask, k] = P[rest, k]
        sub = rest
        while True:
            S = sub | low
            if hc[S] > 0:
                rem = mask ^ S
                for k in range(K - 1):
                    P[mask, k + 1] += hc[S] * P[rem, k]
            if sub == 0:
                break
            sub = (sub - 1) & rest
    return P[full]


@njit(cache=True)
def hcount(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            w = 0
            while out:
                if out & 1:
                    dp[mask | (1 << w), w] += c
                out >>= 1
                w += 1
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


def one_apex_T(W):
    """T = TT_k + u: chain c_0 -> c_1 -> ... (c_i beats c_j, i<j); u beats c_i iff bit (k-1-i) of W (MSB first).
    vertices: chain 0..k-1, apex k. returns (n, set of arcs)"""
    k = W.bit_length()
    arcs = set()
    for i in range(k):
        for j in range(i + 1, k):
            arcs.add((i, j))
        if (W >> (k - 1 - i)) & 1:
            arcs.add((k, i))
        else:
            arcs.add((i, k))
    return k + 1, arcs, k


def one_apex_B(W):
    """B = TT_l + w: chain d_0 -> d_1 -> ...; d_j beats w iff bit j of W (LSB first). apex index l."""
    l = W.bit_length()
    arcs = set()
    for i in range(l):
        for j in range(i + 1, l):
            arcs.add((i, j))
        if (W >> i) & 1:
            arcs.add((i, l))
        else:
            arcs.add((l, i))
    return l + 1, arcs, l


def glue(W1, W2):
    n1, a1, u = one_apex_T(W1)
    n2, a2, w = one_apex_B(W2)
    n = n1 + n2
    arcs = set(a1)
    for (i, j) in a2:
        arcs.add((n1 + i, n1 + j))
    for i in range(n1):
        for j in range(n2):
            if i == u and j == w:
                arcs.add((n1 + j, i))
            else:
                arcs.add((i, n1 + j))
    adj = np.zeros(n, np.int64)
    for (i, j) in arcs:
        adj[i] |= (1 << j)
    assert len(arcs) == n * (n - 1) // 2
    return n, adj


def adj_from(n, arcs):
    adj = np.zeros(n, np.int64)
    for (i, j) in arcs:
        adj[i] |= (1 << j)
    return adj


def section5():
    print()
    print("=== A2.4 packing polynomial of the F1 joins by brute force (section 5 / THM-5134 (H)) ===")
    psi = psi_table(1 << 12)
    # one-apex pieces: I(Omega(R(h)), z) = 1 + psi(h) z ; H = 2 psi + 1 ; paths from the apex = 2h+1
    okp = True
    for h in range(0, 64):
        n, arcs, u = one_apex_T(2 * h)
        adj = adj_from(n, arcs)
        P = packing_poly(adj, n)
        if not (P[0] == 1 and P[1] == psi[h] and all(P[2:] == 0)):
            okp = False
        if hcount(adj, n) != 2 * psi[h] + 1:
            okp = False
    print("one-apex R(h) = T(W = 2h), h < 64: packing polynomial 1 + psi(h) z and H = 2 psi(h) + 1:", okp)
    cnt = 0
    bad = 0
    maxn = 0
    for h1 in range(0, 24):
        for h2 in range(0, 24):
            for s in range(0, 5):
                W1 = 2 * h1
                W2 = (1 << s) * (2 * h2 + 1) - 1
                n = W1.bit_length() + W2.bit_length() + 2
                if n > 13:
                    continue
                n, adj = glue(W1, W2)
                P = packing_poly(adj, n)
                y1, y2 = 2 * h1 + 1, 2 * h2 + 1
                cc = ((1 << s) * y1 * y2 - (1 if s == 0 else 0))
                assert cc % 2 == 0
                a1 = psi[h1] + psi[h2] + cc // 2
                a2 = psi[h1] * psi[h2]
                H = hcount(adj, n)
                N = (2 * psi[h1] + 1) * (2 * psi[h2] + 1) + (1 << s) * y1 * y2 - (1 if s == 0 else 0)
                cnt += 1
                maxn = max(maxn, n)
                if not (P[0] == 1 and P[1] == a1 and P[2] == a2 and all(P[3:] == 0) and H == N
                        and H == 1 + 2 * a1 + 4 * a2):
                    bad += 1
                    if bad < 5:
                        print("   MISMATCH", h1, h2, s, list(P), a1, a2, H, N)
    print("joins J(h1,h2,s) with W = 2 h1, W' = 2^s (2 h2 + 1) - 1, order <= 13: %d checked (max order %d), mismatches: %d"
          % (cnt, maxn, bad))
    # trailing ones on both sides: depends on s = a + b only
    cnt2 = bad2 = 0
    for h1 in range(0, 8):
        for h2 in range(0, 8):
            for a in range(1, 3):
                for b in range(1, 3):
                    W1 = (1 << a) * (2 * h1 + 1) - 1
                    W2 = (1 << b) * (2 * h2 + 1) - 1
                    n = W1.bit_length() + W2.bit_length() + 2
                    if n > 12:
                        continue
                    n, adj = glue(W1, W2)
                    P = packing_poly(adj, n)
                    s = a + b
                    cc = (1 << s) * (2 * h1 + 1) * (2 * h2 + 1)
                    cnt2 += 1
                    if not (P[1] == psi[h1] + psi[h2] + cc // 2 and P[2] == psi[h1] * psi[h2]):
                        bad2 += 1
    print("joins with trailing ones on both sides (s = a + b): %d checked, mismatches: %d" % (cnt2, bad2))

    print()
    print("=== A2.5 digit claims (section 5) ===")
    for target in (17, 45, 7, 21, 9, 13, 15):
        vecs = {}
        for h1 in range(0, 64):
            for h2 in range(h1, 64):
                for s in range(0, 8):
                    y1, y2 = 2 * h1 + 1, 2 * h2 + 1
                    cc = (1 << s) * y1 * y2 - (1 if s == 0 else 0)
                    N = (2 * psi[h1] + 1) * (2 * psi[h2] + 1) + cc
                    if N == target:
                        v = (psi[h1] + psi[h2] + cc // 2, psi[h1] * psi[h2])
                        vecs.setdefault(v, []).append((h1, h2, s))
        print("N = %d: packing vectors (alpha_1, alpha_2) of F1 joins (h1 <= h2 < 64, s < 8): %d vectors" % (target, len(vecs)))
        for v in sorted(vecs):
            print("    %s  from %s" % (v, vecs[v]))
    # digit condition b > max(alpha_1, alpha_2)
    okd = True
    for a1 in range(0, 30):
        for a2 in range(0, 30):
            for b in range(max(a1, a2) + 1, max(a1, a2) + 4):
                if b < 2:
                    continue
                Nb = 1 + a1 * b + a2 * b * b
                digs = []
                t = Nb
                while t:
                    digs.append(t % b)
                    t //= b
                want = [1, a1, a2]
                while len(want) > 1 and want[-1] == 0:
                    want.pop()
                if digs != want:
                    okd = False
    print("for b > max(alpha_1, alpha_2) the base-b numeral of N_b is (alpha_2)(alpha_1)(1):", okd)


if __name__ == "__main__":
    section6()
    section5()
