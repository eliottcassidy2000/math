#!/usr/bin/env python3
"""Skeptic check 1: the objects of THM-5141 (draft) and the inherited identification.

Independent of the coordinator's scripts.  Exact integer arithmetic.

 (a) q, x, g recursions; g(v) = v + 1 - 2 q(v) - (v mod 2); x(h) = g(2h); sizes.
 (b) brute force: Hamiltonian path counts of 'one-apex' tournaments (transitive chain + one extra vertex)
     by subset DP, compared with g under all four reading conventions; set equality with {x(h)}.
 (c) block law, exhaustive on small ranges; three-block formula (K independent of x), incl. s = 0, w = 1.
"""
import sys
from itertools import product

def q_rec(h, memo={0: 0}):
    if h in memo:
        return memo[h]
    if h % 2 == 0:
        r = q_rec(h // 2)
    else:
        k = h // 2
        r = q_rec(k) + (k + 1) // 2
    memo[h] = r
    return r

def g_rec(v, memo={0: 1}):
    if v in memo:
        return memo[v]
    if v % 2 == 1:
        r = g_rec(v // 2)
    else:
        u = v // 2
        r = g_rec(u) + u + (u % 2)
    memo[v] = r
    return r

def x_of(h):
    return 2 * h + 1 - 2 * q_rec(h)

def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r

def ham_paths(n, adj):
    """adj[i] = bitmask of out-neighbours of i; count Hamiltonian paths by subset DP."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c == 0:
                continue
            nxt = adj[last] & ~mask
            while nxt:
                low = nxt & (-nxt)
                j = low.bit_length() - 1
                dp[mask | low][j] += c
                nxt ^= low
    return sum(dp[full])

def main():
    out = []
    P = lambda *a: (print(*a), out.append(" ".join(str(t) for t in a)))
    N = 1 << 17
    bad = 0
    for v in range(N):
        if g_rec(v) != v + 1 - 2 * q_rec(v) - (v % 2):
            bad += 1
    P("(a) g(v) = v+1-2q(v)-(v mod 2) for v < 2^17: violations", bad)
    bad = sum(1 for h in range(N // 2) if x_of(h) != g_rec(2 * h))
    P("(a) x(h) = g(2h) for h < 2^16: violations", bad)
    bad = sum(1 for h in range(1, N) if not (0 <= q_rec(h) and 2 * q_rec(h) <= h - 1))
    P("(a) 0 <= q(h) <= (h-1)/2 for 1 <= h < 2^17: violations", bad)
    bad = sum(1 for h in range(1, N) if not (h + 2 <= x_of(h) <= 2 * h + 1))
    P("(a) h+2 <= x(h) <= 2h+1 for 1 <= h < 2^17: violations", bad)
    P("(a) x(0) =", x_of(0), " g(0) =", g_rec(0))
    # set {x(h)} versus set {g(v)}
    M = 1 << 14
    sx = set(x_of(h) for h in range(M))
    sg = set(g_rec(v) for v in range(2 * M))
    P("(a) {g(v): v < 2^15} == {x(h): h < 2^14}:", sx == sg)
    # unrolled formula for q
    bad = 0
    for h in range(1 << 14):
        s = 0
        i = 0
        while (h >> i):
            if (h >> i) & 1:
                s += ((h >> (i + 1)) + 1) // 2
            i += 1
        if s != q_rec(h):
            bad += 1
    P("(a) q(h) = sum_i h_i ceil(floor(h/2^(i+1))/2) for h < 2^14: violations", bad)

    # (b) brute force one-apex tournaments
    P("(b) one-apex tournaments: chain c_0 -> c_1 -> ... (transitive, c_i beats c_j for i < j) plus an apex;")
    P("    word bit = 1 iff the apex beats that chain vertex (convention A) or is beaten by it (convention B);")
    P("    the word is read with chain vertex 0 as most significant bit (MSB) or least significant bit (LSB).")
    for m in range(1, 10):
        counts = {}
        for word in range(1 << m):
            # chain vertices 0..m-1, apex = m
            adj = [0] * (m + 1)
            for i in range(m):
                for j in range(i + 1, m):
                    adj[i] |= 1 << j
            for i in range(m):
                if (word >> i) & 1:
                    adj[m] |= 1 << i     # apex beats chain vertex i
                else:
                    adj[i] |= 1 << m
            counts[word] = ham_paths(m + 1, adj)
        res = {}
        for name, f in (("A/LSB", lambda wd: wd),
                        ("A/MSB", lambda wd: rev(wd, m)),
                        ("B/LSB", lambda wd: wd ^ ((1 << m) - 1)),
                        ("B/MSB", lambda wd: rev(wd ^ ((1 << m) - 1), m))):
            res[name] = sum(1 for wd in range(1 << m) if counts[wd] != g_rec(f(wd)))
        setH = set(counts.values())
        setg = set(g_rec(v) for v in range(1 << m))
        P("    m =", m, " mismatches H != g(word) by convention:", res,
          " set{H} == {g(v): v < 2^m}:", setH == setg, " |set| =", len(setH))

    # (c) block law, exhaustive small
    bad = 0
    tot = 0
    for m in range(1, 9):
        for A in range(0, 1 << 7):
            for B in range(1 << m):
                tot += 1
                lhs = 2 * q_rec(A * (1 << m) + B)
                rhs = 2 * q_rec(A) + 2 * q_rec(B) + A * rev(B, m) + (A % 2) * ((B >> (m - 1)) & 1)
                if lhs != rhs:
                    bad += 1
    P("(c) block law 2q(A 2^m + B) = 2q(A)+2q(B)+A rev_m(B)+(A mod 2)B_(m-1): cases", tot, "violations", bad)
    # large random
    import random
    random.seed(20261111)
    bad = 0
    for _ in range(20000):
        m = random.randint(1, 60)
        A = random.getrandbits(random.randint(0, 60))
        B = random.getrandbits(m)
        lhs = 2 * q_rec(A * (1 << m) + B)
        rhs = 2 * q_rec(A) + 2 * q_rec(B) + A * rev(B, m) + (A % 2) * ((B >> (m - 1)) & 1)
        if lhs != rhs:
            bad += 1
    P("(c) block law, 20000 random cases with m <= 60, A < 2^60: violations", bad)

    # three-block formula
    bad = 0
    tot = 0
    for tp in range(1, 5):
        for w in range(1, 6):
            for s in range(0, 5):
                for Pp in range(1 << (tp - 1), 1 << tp):
                    for Z in range(1 << s):
                        c = (1 << (s + 1)) - rev(Z, s)
                        assert (1 << s) < c <= (1 << (s + 1))
                        if s >= 1:
                            assert c % 2 == (Z >> (s - 1)) & 1
                        else:
                            assert c == 2
                        Ks = set()
                        for xx in range(1 << w):
                            h = (Pp << (w + s)) | (xx << s) | Z
                            val = x_of(h) - (c * xx - Pp * rev(xx, w) - 2 * q_rec(xx)
                                             - (Pp % 2) * ((xx >> (w - 1)) & 1) - (c % 2) * (xx & 1))
                            Ks.add(val)
                        tot += 1
                        if len(Ks) != 1:
                            bad += 1
    P("(c) three-block formula: K(P,Z) independent of x; (P,Z,w) cases", tot, "violations", bad)
    # Z -> c bijection
    ok = all(sorted((1 << (s + 1)) - rev(Z, s) for Z in range(1 << s)) == list(range((1 << s) + 1, (1 << (s + 1)) + 1))
             for s in range(0, 12))
    P("(c) Z -> c = 2^(s+1) - rev_s(Z) is a bijection onto (2^s, 2^(s+1)] for s <= 11:", ok)
    with open(sys.argv[1] if len(sys.argv) > 1 else "s1_basics.out", "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
