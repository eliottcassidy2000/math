#!/usr/bin/env python3
"""The square-sum graph Q_n (vertices 1..n, x ~ y iff x != y and x + y is a square): the owner's
observations of 2026-10-22, checked and explained (opus-2026-10-22-ssj, wave 2).

Inherited (not re-derived): Hamiltonian-path existence for n <= 40 is FINITE-EXACT and audited in
05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md (paths exactly at
n = 1, 15, 16, 17, 23, 25..40); Gerbicz (CITED there) for all n >= 25.

Parts
  1. components of Q_n for n <= 40 (counts, the n = 12 components and their sizes, the welding edges)
  2. new edges s(n) = #{a < n : a + n square} = floor(sqrt(2n-1)) - floor(sqrt(n)); run boundaries:
     gains at n = ceil((j^2+1)/2), losses at n = k^2; simultaneous gain and loss exactly at
     n = k^2 with k a Pell number (k^2 - 2t^2 = 1 or 2k^2 - m^2 = 1), i.e. n = 1, 9, 25, 289, 841, ...
  3. the onset lemma: if Q_n has a Hamiltonian path and Q_{n-1} has none, then s(n) >= 2
     (n is interior); the extension mechanism at 16, 17 and the loss at 18
  4. the owner's edge counts at n = 25 on {10,6,3,1,8,17}, {12,4,5,11}, {2,7,18,9,16}
  5. hostile probe for the 3-4-5 component sizes: other sum graphs (a + b in S)
  6. the Paley tournament QR_7: H = 189 and the free action of Aut = F_21 (9 orbits)
Run: python3 -B square_sum_structure.py > square_sum_structure.out
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import Counter


def require(cond, msg):
    if not cond:
        raise SystemExit("FAILED: " + msg)


def is_square(x: int) -> bool:
    return x >= 0 and math.isqrt(x) ** 2 == x


def edges(n, S=is_square):
    return [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1) if S(a + b)]


def components(n, S=is_square):
    parent = list(range(n + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges(n, S):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = {}
    for v in range(1, n + 1):
        comp.setdefault(find(v), []).append(v)
    return sorted(comp.values(), key=lambda c: (len(c), c))


def s_new(n):
    """number of squares k^2 in [n+1, 2n-1], i.e. of a < n with a + n square"""
    return sum(1 for k in range(math.isqrt(n) + 1, math.isqrt(max(2 * n - 1, 0)) + 1) if n + 1 <= k * k <= 2 * n - 1)


def s_brute(n):
    return sum(1 for a in range(1, n) if is_square(a + n))


# known Hamiltonian-path set for n <= 40 (inherited, FINITE-EXACT, audited)
HAM = {1, 15, 16, 17, 23} | set(range(25, 41))


def is_path(seq, n):
    return sorted(seq) == list(range(1, n + 1)) and all(is_square(seq[i] + seq[i + 1]) for i in range(n - 1))


def part1():
    print("== 1. components of Q_n ==")
    counts = []
    for n in range(1, 41):
        c = components(n)
        counts.append(len(c))
    print("  #components for n = 1..40:", counts)
    require(counts[0] == 1 and counts[1] == 2 and counts[2] == 2, "n = 1..3")
    require(all(counts[n - 1] == 3 for n in range(4, 13)), "3 components for 4..12")
    require(counts[12] == 2 and all(counts[n - 1] == 1 for n in range(14, 41)), "2 at 13, 1 from 14")
    c12 = components(12)
    print("  Q_12 components:", c12, " sizes", sorted(len(c) for c in c12))
    require(sorted(len(c) for c in c12) == [3, 4, 5], "sizes 3,4,5")
    require(max(a + b for a, b in edges(12)) == 16 and 11 + 12 < 25, "largest edge sum in Q_12 is 16; largest pair sum 23 < 25")
    print("  squares used inside Q_12:", sorted({a + b for a, b in edges(12)}), "(the square 25 = 3^2 + 4^2 is not yet available)")
    for n in (13, 14):
        new = [(a, n, a + n) for a in range(1, n) if is_square(a + n)]
        print(f"  welding edges of vertex {n}: {new}")
    print("  each Q_12 component is a path: ", [sorted(Counter(v for e in edges(12) for v in e if v in c).values()) for c in c12])


def part2():
    print("\n== 2. new edges s(n) and their runs ==")
    for n in range(1, 5001):
        require(s_brute(n) == s_new(n), f"brute force at {n}")
    for n in range(2, 1000001):
        require(s_new(n) == math.isqrt(2 * n - 1) - math.isqrt(n), f"formula at {n}")
    print("  s(n) = floor(sqrt(2n-1)) - floor(sqrt(n)): equals the brute-force count for n <= 5000 and the square count for 2 <= n <= 10^6")
    seq = [s_new(n) for n in range(1, 61)]
    runs, start = [], 1
    for n in range(2, 61):
        if seq[n - 1] != seq[n - 2]:
            runs.append((start, n - 1, seq[n - 2]))
            start = n
    runs.append((start, 60, seq[59]))
    print("  runs (first, last, s) for n <= 60:", runs)
    gains = sorted({(j * j + 2) // 2 if j % 2 else j * j // 2 + 1 for j in range(2, 80)})
    losses = [k * k for k in range(2, 60)]
    # boundary check: s(n) - s(n-1) = [n is a gain] - [n is a loss]
    G, L = set(gains), set(losses)
    for n in range(3, 3000):
        require(s_new(n) - s_new(n - 1) == (n in G) - (n in L), f"boundary law at {n}")
    print("  s(n) - s(n-1) = [n = ceil((j^2+1)/2)] - [n = k^2] verified for n < 3000")
    print("  first gains (square j^2 enters at n = ceil((j^2+1)/2)):", gains[:10])
    print("  first losses (square k^2 leaves at n = k^2):", losses[:8])
    both = sorted(G & L)
    print("  simultaneous gain and loss (invisible run boundaries) below 6400:", [b for b in both if b < 6400])
    # Pell characterisation, checked by brute force to 10^7
    coinc = []
    k = 1
    while k * k <= 10 ** 7:
        n = k * k
        j_odd = math.isqrt(2 * n - 1)
        hit = (j_odd * j_odd == 2 * n - 1) or (is_square(2 * n - 2) and math.isqrt(2 * n - 2) % 2 == 0)
        if hit:
            coinc.append(n)
        k += 1
    pell = []
    for a, b in [(1, 1)]:
        pass
    # k with k^2 - 2t^2 = 1 (k = 3, 17, 99, 577, 3363) or 2k^2 - m^2 = 1 (k = 1, 5, 29, 169, 985)
    ks = set()
    x, y = 3, 2
    while x * x <= 10 ** 7:
        ks.add(x)
        x, y = 3 * x + 4 * y, 2 * x + 3 * y
    m, k2 = 1, 1
    while k2 * k2 <= 10 ** 7:
        ks.add(k2)
        m, k2 = 3 * m + 4 * k2, 2 * m + 3 * k2
    require(coinc == sorted(kk * kk for kk in ks), f"Pell characterisation: {coinc} vs {sorted(kk*kk for kk in ks)}")
    print("  all coincidences n <= 10^7:", coinc, "= k^2 for k in", sorted(ks))
    print("  (at n = 25 the square 25 leaves the window [n+1, 2n-1] exactly as 49 enters: 7^2 - 2*5^2 = -1)")


def ham_path_witness(n):
    """Bitmask DFS with the sound pruning rules of the wave-six lane: start at a leaf if there is one;
    every unvisited vertex must keep a neighbour among the unvisited vertices or the current end;
    at most one unvisited vertex may have exactly one such neighbour; the unvisited set is connected."""
    adj = [0] * (n + 1)
    for a, b in edges(n):
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    full = sum(1 << v for v in range(1, n + 1))
    leaves = [v for v in range(1, n + 1) if bin(adj[v]).count("1") <= 1]
    if len(leaves) > 2 or any(adj[v] == 0 for v in range(1, n + 1)) and n > 1:
        return None
    starts = leaves if leaves else list(range(1, n + 1))

    def connected(mask):
        if mask == 0:
            return True
        start = mask & -mask
        seen, frontier = start, start
        while frontier:
            nxt = 0
            m = frontier
            while m:
                b = m & -m
                nxt |= adj[b.bit_length() - 1]
                m ^= b
            nxt &= mask & ~seen
            seen |= nxt
            frontier = nxt
        return seen == mask

    def dfs(v, visited, path):
        if visited == full:
            return list(path)
        rest = full & ~visited
        ones = 0
        m = rest
        while m:
            b = m & -m
            u = b.bit_length() - 1
            avail = adj[u] & (rest | (1 << v))
            c = bin(avail).count("1")
            if c == 0:
                return None
            if c == 1:
                ones += 1
                if ones > 1 and not (adj[u] & (1 << v)):
                    return None
            m ^= b
        if not connected(rest):
            return None
        m = adj[v] & rest
        while m:
            b = m & -m
            u = b.bit_length() - 1
            path.append(u)
            r = dfs(u, visited | b, path)
            if r:
                return r
            path.pop()
            m ^= b
        return None

    for st in starts:
        r = dfs(st, 1 << st, [st])
        if r:
            return r
    return None


def part3():
    print("\n== 3. the onset lemma and the extension mechanism ==")
    onsets = [n for n in range(2, 41) if n in HAM and (n - 1) not in HAM]
    losses = [n for n in range(2, 41) if n not in HAM and (n - 1) in HAM]
    print("  onsets (path at n, none at n-1):", onsets, "  s:", [s_new(n) for n in onsets])
    print("  losses (path at n-1, none at n):", losses, "  s:", [s_new(n) for n in losses])
    require(all(s_new(n) >= 2 for n in onsets), "onset lemma")
    dbl = [n for n in range(1, 26) if s_new(n) == 2]
    print("  vertices n <= 25 with s(n) = 2:", dbl, "; of these, onsets:", [n for n in dbl if n in onsets],
          "; non-onsets:", [n for n in dbl if n not in onsets])
    # the unique chains at 15, 16, 17 (inherited witnesses), and where 16, 17, 18 attach
    c15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    require(is_path(c15, 15), "15-chain")
    c16 = c15 + [16]
    c17 = [17] + c16
    require(is_path(c16, 16) and is_path(c17, 17), "16, 17 extensions")
    print("  unique 15-chain:", c15, " endpoints", (c15[0], c15[-1]))
    print("  16 + 9 = 25 extends at the end 9; 17 + 8 = 25 extends at the end 8:", c17)
    print("  18's only neighbour below 31:", [a for a in range(1, 31) if a != 18 and is_square(a + 18)],
          "- an interior vertex of the unique 17-chain, so 18 cannot be appended: the loss at 18")
    for n in (15, 16, 17, 23, 25):
        w = ham_path_witness(n)
        require(w is not None and is_path(w, n), f"witness at {n}")
    print("  independent DFS witnesses found for n = 15, 16, 17, 23, 25")
    for n in (18, 19, 24):
        require(ham_path_witness(n) is None, f"no path at {n}")
    print("  independent DFS confirms no Hamiltonian path at n = 18, 19, 24")


def part4():
    print("\n== 4. the owner's edge counts at n = 25 ==")
    E25 = set(edges(25))
    for S in ({10, 6, 3, 1, 8, 17}, {12, 4, 5, 11}, {2, 7, 18, 9, 16}):
        inner = sorted((a, b) for a, b in E25 if a in S and b in S)
        deg = Counter(v for e in inner for v in e)
        shape = "path" if max(deg.values()) <= 2 else "tree with a vertex of degree %d" % max(deg.values())
        print(f"  {sorted(S)}: {len(S)} vertices, {len(inner)} induced edges {inner}  ({shape})")
    print("  16, 17, 18 are the 25-partners of 9, 8, 7: the reflection x -> 25 - x on the n = 12 components")
    counts = sorted(len([e for e in E25 if e[0] in S and e[1] in S]) for S in ({10, 6, 3, 1, 8, 17}, {12, 4, 5, 11}, {2, 7, 18, 9, 16}))
    require(counts == [3, 4, 5], "edge counts 3,4,5")


def part5():
    print("\n== 5. hostile probe: component sizes of other sum graphs a + b in S ==")
    families = {
        "squares": is_square,
        "squares + 1": lambda x: is_square(x - 1),
        "squares - 1": lambda x: is_square(x + 1),
        "twice squares": lambda x: x % 2 == 0 and is_square(x // 2),
        "triangular": lambda x: is_square(8 * x + 1),
        "cubes": lambda x: round(x ** (1 / 3)) ** 3 == x,
        "pentagonal": lambda x: (lambda r: r * r == 24 * x + 1 and (r + 1) % 6 == 0)(math.isqrt(24 * x + 1)),
        "primes": lambda x: x > 1 and all(x % d for d in range(2, math.isqrt(x) + 1)),
    }
    for name, S in families.items():
        pyth = []
        last_multi = None
        for n in range(3, 120):
            c = components(n, S)
            if len(c) >= 3:
                last_multi = (n, sorted(len(x) for x in c))
                sz = sorted(len(x) for x in c)
                for a, b, cc in itertools.combinations(sz, 3):
                    if a * a + b * b == cc * cc:
                        pyth.append((n, sz))
                        break
        lm = last_multi if last_multi is None or len(last_multi[1]) <= 8 else (last_multi[0], f"{len(last_multi[1])} components")
        print(f"  {name:14s}: last n < 120 with >= 3 components: {lm}; stages whose sizes contain a Pythagorean triple: "
              f"{[(n, sz) for n, sz in pyth][:4]}{' ...' if len(pyth) > 4 else ''} ({len(pyth)} stages)")


def part6():
    print("\n== 6. the Paley tournament QR_7: H(QR_7) = 189 and the free action of F_21 ==")
    R = {1, 2, 4}
    V = range(7)
    arc = lambda u, v: (v - u) % 7 in R
    paths = [p for p in itertools.permutations(V) if all(arc(p[i], p[i + 1]) for i in range(6))]
    print("  Hamiltonian paths of QR_7:", len(paths))
    require(len(paths) == 189, "H(QR_7) = 189")
    aut = [(a, b) for a in (1, 2, 4) for b in range(7)]
    require(len(aut) == 21, "|F_21|")
    for a, b in aut:
        require(all(arc((a * u + b) % 7, (a * v + b) % 7) == arc(u, v) for u in V for v in V if u != v), "automorphism")
    orbits = set()
    for p in paths:
        orbits.add(min(tuple((a * x + b) % 7 for x in p) for a, b in aut))
    free = all(tuple((a * x + b) % 7 for x in p) != p for p in paths for a, b in aut if (a, b) != (1, 0))
    print(f"  Aut(QR_7) = x -> ax + b (a a square mod 7) has order 21; it acts freely on the 189 paths: {free}; orbits: {len(orbits)}")
    require(free and len(orbits) == 9, "free action, 9 orbits")


if __name__ == "__main__":
    print(f"python {sys.version.split()[0]}")
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    print("\nALL CHECKS PASSED")
