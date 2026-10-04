"""Exact finite audits for the owner's tiling, family, and atom proposal.

Run: python3 04-computation/experiments/tiling_modular_atoms_20261003.py
All arithmetic is integer; no external packages or inherited data filters.
"""
from itertools import combinations, permutations
from math import isqrt


def family(n):
    return n, (n + 1) // 2, 2 * n + 1


def fold(r, m):
    return (r, 1) if r <= m // 2 else (m - r, -1)


def tile(u, v, n):
    if u == v:
        return v, n
    if u == v + 1:
        return 1, v
    return v + 1, u - 1


def canonical(mask, n=5):
    pairs = list(combinations(range(n), 2))
    adj = [[0] * n for _ in range(n)]
    for bit, (u, v) in enumerate(pairs):
        adj[u][v] = (mask >> bit) & 1
        adj[v][u] = 1 - adj[u][v]
    return min(sum(adj[p[u]][p[v]] << bit
                   for bit, (u, v) in enumerate(pairs))
               for p in permutations(range(n)))


def main():
    cells = 0
    for m in range(2, 102):
        k = m // 2
        wedge = {(a, b): a * b % m
                 for a in range(1, k + 1) for b in range(a, k + 1)}
        for r in range(m):
            for s in range(m):
                if r == 0 or s == 0:
                    recovered = 0
                else:
                    a, er = fold(r, m)
                    b, es = fold(s, m)
                    recovered = er * es * wedge[tuple(sorted((a, b)))] % m
                assert recovered == r * s % m
                cells += 1
    print(f"Folding: every cell for 2 <= modulus <= 101: {cells} checks")
    boundary = [a * 5 % 10 for a in range(5, 0, -1)] + list(range(4, 0, -1))
    assert boundary == [5, 0, 5, 0, 5, 4, 3, 2, 1]
    print("mod 10 L boundary:", boundary)
    for n in range(2, 31):
        image = {tile(u, v, n) for u in range(1, n + 1)
                 for v in range(1, u + 1)}
        assert image == {(a, b) for a in range(1, n + 1)
                         for b in range(a, n + 1)}
    print("Pair/loop/path/wedge coordinate bijection: 2 <= n <= 30")
    for N in range(1, 501):
        assert family(2*N-1) == (2*N-1, N, 4*N-1)
        assert family(2*N) == (2*N, N, 4*N+1)
        for m, d, c in [(4*N-1, N, 3*N-1), (4*N+1, 3*N+1, N)]:
            h = (m-1)//2
            assert h*h % m == (h+1)*(h+1) % m == d
            assert h*(h+1) % m == c
            assert d+c == m
    print("Odd central blocks: both families, 1 <= atom N <= 500")
    for k in range(2, 501):
        m = 2*k
        center = k if k % 2 else 0
        assert k*k % m == center
        assert (k-1)**2 % m == (k+1)**2 % m == (center+1) % m
        assert (k-1)*(k+1) % m == (center-1) % m
    print("Even centers and four diagonal neighbors: 2 <= half-modulus <= 500")
    seen = set()
    for a in range(1, 101):
        for b in range(1, 101):
            s, t = a+b, (2*a+1)*b
            assert (s, t) not in seen
            seen.add((s, t))
            D = (2*s+1)**2 - 8*t
            assert D == (2*a-2*b+1)**2
            r = isqrt(D)
            candidates = [(2*s+1+sign*r)//4 for sign in (-1, 1)
                          if (2*s+1+sign*r) % 4 == 0]
            assert candidates == [b]
            assert 2*t == s*(s+1) - (a-b)*(a-b+1)
            assert 2*(a+t)+1 == (2*a+1)*(2*b+1)
            assert s % 2 == (a+b) % 2 and t % 2 == b % 2
    print("p: 10000 ordered input pairs, inverse, triangular identity, product law")
    for J in range(2, 101):
        actual = {(2*(I+J)-1, 4*J) for I in range(1, J)}
        expected = {(O, 4*J) for O in range(2*J+1, 4*J-2, 2)}
        assert actual == expected
        for O, E in actual:
            I = (2*O+2-E)//4
            assert 1 <= I < J and E//4 == J
            a, b = 2*I-1, 2*J
            assert (a+b, (2*a+1)*b) == (O, E*(2*O-E+1)//2)
    assert (5, 6) not in {(2*(I+J)-1, 4*J)
                         for J in range(2, 101) for I in range(1, J)}
    assert (2*7+2-8)//4 == 2  # (7,8) would require I=J, excluded.
    print("q: 4950 pairs, exact image/inverse; hostile (5,6) and boundary (7,8)")
    pairs = list(combinations(range(5), 2))
    free_bits = [i for i, (u, v) in enumerate(pairs) if v-u >= 2]
    tilings = [sum(((t >> j) & 1) << bit for j, bit in enumerate(free_bits))
               for t in range(64)]
    all_classes = {canonical(mask) for mask in range(1024)}
    tiled_classes = {canonical(mask) for mask in tilings}
    assert tiled_classes == all_classes and len(all_classes) == 12
    print("Order 5: all 64 fixed-path tilings cover all 12 classes of 1024 labeled tournaments")
    # Independent switching-gauge count: every labeled tournament has exactly two selfies.
    multiplicities = [0] * 1024
    for mask in tilings:
        for loops in range(32):
            cut = sum(((((loops >> u) ^ (loops >> v)) & 1) << i)
                      for i, (u, v) in enumerate(pairs))
            multiplicities[mask ^ cut] += 1
    assert set(multiplicities) == {2}
    print("Order 5 selfie gauge: 2048 inputs, every labeled tournament has exactly 2 preimages")
    print("ALL CHECKS PASSED")


if __name__ == '__main__':
    main()
