"""Auditor's own check of the constructions in THM-4951 (2a), (2b), (3b), (3d), of the T5 first-hit witnesses, and a
brute-force count of tight colourings (3c).  Written independently of the packet's staircase_check.py.

Frame: path 0 -> 1 -> ... -> n-1, tile (i, j) (j >= i + 2) reversed in t means j -> i; t + 1 = complement of t;
sigma(i, j) = (n-1-j, n-1-i).  nu = max number of disjoint 3-cycles (memoised branch on the lowest free vertex).
Feedback claims are checked by testing that the tournament minus the given set has no 3-cycle (a tournament is
acyclic iff it has no 3-cycle), so nu <= tau <= |set| pins tau when nu = |set|.
usage: python staircase_audit.py PACKET_DIR
"""
import sys
import os
import itertools
from functools import lru_cache


def tourn(n, rev):
    """out-adjacency bitmasks of the tiling with reversed tile set rev"""
    out = [0] * n
    for v in range(n - 1):
        out[v] |= 1 << (v + 1)
    for i in range(n):
        for j in range(i + 2, n):
            if (i, j) in rev:
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
    return out


def all_tiles(n):
    return {(i, j) for i in range(n) for j in range(i + 2, n)}


def triangles(out, n):
    tri = []
    for a, b, c in itertools.combinations(range(n), 3):
        if (out[a] >> b & 1 and out[b] >> c & 1 and out[c] >> a & 1) or (out[a] >> c & 1 and out[c] >> b & 1 and out[b] >> a & 1):
            tri.append((1 << a) | (1 << b) | (1 << c))
    return tri


def nu(out, n):
    tri = triangles(out, n)
    by_low = {}
    for T in tri:
        lo = (T & -T).bit_length() - 1
        by_low.setdefault(lo, []).append(T)

    @lru_cache(maxsize=None)
    def best(free):
        if free.bit_count() < 3:
            return 0
        lo = (free & -free).bit_length() - 1
        r = best(free & ~(1 << lo))
        for T in by_low.get(lo, []):
            if T & free == T:
                r = max(r, 1 + best(free & ~T))
        return r
    return best((1 << n) - 1)


def acyclic_without(out, n, F):
    alive = [v for v in range(n) if v not in F]
    sub = [0] * n
    for v in alive:
        sub[v] = out[v] & sum(1 << w for w in alive)
    for a, b, c in itertools.combinations(alive, 3):
        if (sub[a] >> b & 1 and sub[b] >> c & 1 and sub[c] >> a & 1) or (sub[a] >> c & 1 and sub[c] >> b & 1 and sub[b] >> a & 1):
            return False
    return True


def inside(segs):
    return {(i, j) for (p, q) in segs for i in range(p, q + 1) for j in range(i + 2, q + 1)}


def blue(n, rev):
    return all(((n - 1 - j, n - 1 - i) in rev) == ((i, j) in rev) for (i, j) in all_tiles(n))


def P(X):
    c, last = 0, -10
    for a in sorted(X):
        if a - last >= 3:
            c, last = c + 1, a
    return c


def main():
    fails = []
    # (2a) prefix staircase, with the explicit feedback sets of the proof
    cnt = 0
    for n in range(3, 22):
        k = n // 3
        for j in range(0, k):
            a = 3 * j + 2
            rev = inside([(0, a - 1)])
            t, u = tourn(n, rev), tourn(n, all_tiles(n) - rev)
            Ft = {3 * i + 2 for i in range(j)}                        # every third vertex of L = [0, a-1]
            Fu = {a} | {a + 1 + 3 * i + 2 for i in range((n - a - 1) // 3)}  # a plus every third vertex of [a+1, n-1]
            ok = (nu(t, n), nu(u, n)) == (j, k - j) and len(Ft) == j and len(Fu) == k - j \
                and acyclic_without(t, n, Ft) and acyclic_without(u, n, Fu) and (blue(n, rev) == (j == 0))
            cnt += 1
            if not ok:
                fails.append(("2a", n, j))
    print(f"(2a) prefix staircase n = 3..21, all j: {cnt} cases, nu = (j, k-j), tau = nu via the proof's feedback sets, "
          f"black iff j >= 1: failures {[f for f in fails if f[0] == '2a']}")
    # (2b) double staircase
    cnt = 0
    for j in range(1, 4):
        a = 3 * j + 2
        for c in range(0, 10):
            n = 2 * a + c
            if n > 26:
                continue
            k = n // 3
            rev = inside([(0, a - 1), (a + c, n - 1)])
            t, u = tourn(n, rev), tourn(n, all_tiles(n) - rev)
            ok = blue(n, rev) and nu(t, n) == 2 * j and nu(u, n) == k - 2 * j
            cnt += 1
            if not ok:
                fails.append(("2b", j, c, nu(t, n), nu(u, n), k))
    print(f"(2b) double staircase j = 1..3, c = 0..9 (n <= 26): {cnt} cases: failures {[f for f in fails if f[0] == '2b']}")
    # (3b), (3d) multi-staircase at n = 3k + 2
    cnt = 0
    for k in range(1, 7):
        n = 3 * k + 2
        for beta in itertools.product((0, 1), repeat=k):
            segs, b = [], 0
            while b < k:
                if beta[b]:
                    s = b
                    while b < k and beta[b]:
                        b += 1
                    segs.append((3 * s, 3 * b + 1))
                else:
                    b += 1
            rev = inside(segs)
            word = [1 if (x, x + 2) in rev else 0 for x in range(n - 2)]
            F = {3 * bb + 2 for bb in range(k) if not beta[bb]}
            t, u = tourn(n, rev), tourn(n, all_tiles(n) - rev)
            ok = word == [beta[x // 3] for x in range(n - 2)] and nu(t, n) == sum(beta) and nu(u, n) == k - sum(beta) \
                and acyclic_without(u, n, F) and blue(n, rev) == (list(beta) == list(beta)[::-1])
            cnt += 1
            if not ok:
                fails.append(("3b", n, beta))
    print(f"(3b)/(3d) multi-staircase, all block words, k = 1..6 (n = 5..20): {cnt} words: failures "
          f"{[f for f in fails if f[0] == '3b']}")
    # (3c) tight colourings by brute force
    print("(3c) tight colourings of [0, L-1], x_0 = 1, brute force, against the stated closed forms:")
    for L in range(1, 20):
        k = -(-L // 3)
        c = 0
        for w in range(1 << (L - 1)):
            x = [1] + [(w >> i) & 1 for i in range(L - 1)]
            A = [a for a in range(L) if x[a]]
            B = [a for a in range(L) if not x[a]]
            if P(A) + P(B) == k:
                c += 1
        n = L + 2
        if n % 3 == 2:
            f = 2 ** (k - 1)
        elif n % 3 == 1:
            f = (k + 1) * 2 ** (k - 2) if k >= 2 else 1
        else:
            f = (k * k + 9 * k - 2) * 2 ** k // 16
        print(f"  n = {n}, k = {k}: {c} (formula {f}) {'ok' if c == f else 'MISMATCH'}")
    # T5 witnesses
    pk = sys.argv[1]
    bad, tot = 0, 0
    for line in open(os.path.join(pk, "t5_first_hits.out")):
        if not line.startswith("FOUND"):
            continue
        f = line.split()
        n = int(f[1].split("=")[1])
        w = f[3]
        nua, nub = int(f[5]), int(f[6])
        col = f[7]
        tiles = {tuple(int(z) for z in s.split(",")) for s in line.split("reversed):")[1].split()}
        t, u = tourn(n, tiles), tourn(n, all_tiles(n) - tiles)
        word = "".join("1" if (x, x + 2) in tiles else "0" for x in range(n - 2))
        a, b = nu(t, n), nu(u, n)
        ok = word == w and (a, b) == (nua, nub) and a + b == n // 3 and (col == "blue") == blue(n, tiles)
        tot += 1
        if not ok:
            bad += 1
            print("  T5 witness fails:", n, w, a, b, word)
    print(f"T5 first-hit witnesses (n = 13, 14): {tot} checked, failures {bad}")


if __name__ == "__main__":
    main()
