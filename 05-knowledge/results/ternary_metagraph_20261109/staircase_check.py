"""Checks of the staircase constructions in THM-4951 (2a), (2b), (3b), (3d) and the tight-colouring counts (3c).

Convention: path 0 -> 1 -> ... -> n-1; tile (i, j), j >= i + 2, reversed in t iff j -> i.  u = t + 1 is the
complement tiling.  nu = maximum number of vertex-disjoint 3-cycles (= of disjoint cycles).

(a) Prefix staircase: reversed tiles = those inside [0, 3j + 1]; claim nu(t) = j, nu(u) = k - j (0 <= j <= k - 1).
(b) Multi-staircase at n = 3k + 2: for a block word beta in {0,1}^k, reversed tiles = those inside the vertex window
    [3s, 3s + 3r + 1] of each maximal run of 1-blocks [s, s + r - 1]; claim the triple word is the block expansion of
    beta, nu(t) = #1-blocks, nu(u) = #0-blocks, and {3b + 2 : beta_b = 0} is a feedback vertex set of u.
(c) Double staircase t = staircase(a) + transitive(c) + staircase(a), a = 3j + 2, n = 2a + c: blue (sigma-fixed),
    nu(t) = 2j, nu(u) = k - 2j (j >= 1, c >= 0).
(d) Tight colourings of [0, L - 1] (x_0 = 1) counted by a greedy-state DP, for L up to 40.
usage: python staircase_check.py
"""
import itertools


def tournament(n, rev):
    out = [0] * n
    for i in range(n - 1):
        out[i] |= 1 << (i + 1)
    for i in range(n):
        for j in range(i + 2, n):
            if (i, j) in rev:
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
    return out


def complement(n, rev):
    return {(i, j) for i in range(n) for j in range(i + 2, n)} - set(rev)


def nu(out, n, alive=None):
    alive = (1 << n) - 1 if alive is None else alive
    tri = []
    for a in range(n):
        for b in range(n):
            if out[a] >> b & 1:
                for c in range(n):
                    if out[b] >> c & 1 and out[c] >> a & 1 and a < b and a < c:
                        tri.append((1 << a) | (1 << b) | (1 << c))
    best = 0

    def rec(av, cur):
        nonlocal best
        best = max(best, cur)
        if not av or cur + bin(av).count("1") // 3 <= best:
            return
        lo = (av & -av)
        for T in tri:
            if T & lo and T & av == T:
                rec(av & ~T, cur + 1)
        rec(av & ~lo, cur)
    rec(alive, 0)
    return best


def acyclic(out, n, alive):
    """no directed cycle in the subtournament on alive (tournament: acyclic iff no 3-cycle)"""
    vs = [v for v in range(n) if alive >> v & 1]
    for a, b, c in itertools.combinations(vs, 3):
        if (out[a] >> b & 1 and out[b] >> c & 1 and out[c] >> a & 1) or (out[a] >> c & 1 and out[c] >> b & 1 and out[b] >> a & 1):
            return False
    return True


def inside(segments, n):
    return {(i, j) for (p, q) in segments for i in range(p, q + 1) for j in range(i + 2, q + 1)}


def sigma_fixed(n, rev):
    return all(((n - 1 - j, n - 1 - i) in rev) == ((i, j) in rev) for i in range(n) for j in range(i + 2, n))


def main():
    # (a)
    bad = 0
    for n in range(3, 19):
        k = n // 3
        for j in range(0, k):
            rev = inside([(0, 3 * j + 1)], n)
            t, u = tournament(n, rev), tournament(n, complement(n, rev))
            if (nu(t, n), nu(u, n)) != (j, k - j):
                bad += 1
                print("prefix staircase fails", n, j)
    print(f"(a) prefix staircase, n = 3..18, all j: failures {bad}")
    # (b)
    bad = 0
    cnt = 0
    for k in range(1, 6):
        n = 3 * k + 2
        for beta in itertools.product((0, 1), repeat=k):
            segs = []
            b = 0
            while b < k:
                if beta[b]:
                    s = b
                    while b < k and beta[b]:
                        b += 1
                    segs.append((3 * s, 3 * b + 1))
                else:
                    b += 1
            rev = inside(segs, n)
            word = tuple(1 if (a, a + 2) in rev else 0 for a in range(n - 2))
            expect = tuple(beta[a // 3] for a in range(n - 2))
            t, u = tournament(n, rev), tournament(n, complement(n, rev))
            F = [3 * bb + 2 for bb in range(k) if not beta[bb]]
            alive = ((1 << n) - 1) & ~sum(1 << f for f in F)
            ok = (word == expect and nu(t, n) == sum(beta) and nu(u, n) == k - sum(beta) and acyclic(u, n, alive)
                  and sigma_fixed(n, rev) == (tuple(beta) == tuple(reversed(beta))))
            cnt += 1
            if not ok:
                bad += 1
                print("multi-staircase fails", n, beta)
    print(f"(b) multi-staircase, n = 5, 8, 11, 14, 17, all {cnt} block words: failures {bad}")
    # (c)
    bad = 0
    for j in range(1, 3):
        a = 3 * j + 2
        for c in range(0, 7):
            n = 2 * a + c
            k = n // 3
            rev = inside([(0, a - 1), (a + c, n - 1)], n)
            t, u = tournament(n, rev), tournament(n, complement(n, rev))
            if not (sigma_fixed(n, rev) and nu(t, n) == 2 * j and nu(u, n) == k - 2 * j):
                bad += 1
                print("double staircase fails", n, j, c, nu(t, n), nu(u, n), k)
    print(f"(c) double staircase, j = 1, 2, c = 0..6: failures {bad}")
    # (d) tight colourings by DP: state (gap since last greedy pick in A, same in B, PA + PB)
    print("(d) tight colourings of [0, L-1] with x_0 = 1:")
    res = []
    for L in range(1, 41):
        k = -(-L // 3)
        st = {(3, 3, 0): 1}
        for pos in range(L):
            ns = {}
            for (ga, gb, s), c in st.items():
                for col in ((1,) if pos == 0 else (0, 1)):
                    if col:
                        nga, ngb, s2 = (1, min(gb + 1, 3), s + 1) if ga >= 3 else (min(ga + 1, 3), min(gb + 1, 3), s)
                    else:
                        nga, ngb, s2 = (min(ga + 1, 3), 1, s + 1) if gb >= 3 else (min(ga + 1, 3), min(gb + 1, 3), s)
                    key = (nga, ngb, s2)
                    ns[key] = ns.get(key, 0) + c
            st = ns
        res.append((L, sum(c for (ga, gb, s), c in st.items() if s == k)))
    for r in range(3):
        print("  L = 3k - " + str((3 - r) % 3) if r else "  L = 3k", [c for (L, c) in res if L % 3 == r][:13])


if __name__ == "__main__":
    main()
