#!/usr/bin/env python3
"""Exact referee for THM-848's forward-polynomial/Mobius diagonalization.

The exhaustive lane checks every labelled tournament through n=5.  The
symbolic lane checks the generator eigenpolynomials through n=14.  This is a
tournament identity only; no runner, gap, or LRC metric data are present.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import comb
from pathlib import Path


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def phi(m, k):
    out = [1]
    for _ in range(k):
        out = poly_mul(out, [1, -1])
    for _ in range(m - k):
        out = poly_mul(out, [1, 1])
    return out


def generator(m, a):
    # L A=(1-x^2)A'-m(1-x)A.
    out = [0] * (m + 1)
    for j in range(1, m + 1):
        out[j - 1] += j * a[j]
        if j + 1 <= m:
            out[j + 1] -= j * a[j]
    for j in range(m + 1):
        out[j] -= m * a[j]
        if j + 1 <= m:
            out[j + 1] += m * a[j]
    return out


def krawtchouk(k, j, m):
    return sum(
        (-1) ** t * comb(j, t) * comb(m - j, k - t)
        for t in range(max(0, k - (m - j)), min(k, j) + 1)
    )


def edge_table(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return edges, {e: k for k, e in enumerate(edges)}


def forward(mask, u, v, index):
    if u < v:
        return (mask >> index[(u, v)]) & 1
    return 1 - ((mask >> index[(v, u)]) & 1)


def forward_counts(n, mask, perms, index):
    m = n - 1
    a = [0] * (m + 1)
    for p in perms:
        j = sum(forward(mask, p[t], p[t + 1], index) for t in range(m))
        a[j] += 1
    return a


def fwht(a):
    out = list(a)
    h = 1
    while h < len(out):
        for i in range(0, len(out), 2 * h):
            for j in range(i, i + h):
                x, y = out[j], out[j + h]
                out[j], out[j + h] = x + y, x - y
        h *= 2
    return out


def exhaustive(n):
    m = n - 1
    edges, index = edge_table(n)
    M = len(edges)
    N = 1 << M
    perms = list(permutations(range(n)))
    rows = [forward_counts(n, mask, perms, index) for mask in range(N)]
    h_values = [a[m] for a in rows]
    walsh_raw = fwht(h_values)

    # layer[k][T] is the degree-k homogeneous Walsh layer evaluated at T.
    layer = [[0] * N for _ in range(m + 1)]
    for k in range(m + 1):
        restricted = [v if s.bit_count() == k else 0 for s, v in enumerate(walsh_raw)]
        values = fwht(restricted)
        layer[k] = [Fraction(v, N) for v in values]

    distinct_B = set()
    for mask, a in enumerate(rows):
        c = [
            Fraction(sum(krawtchouk(k, j, m) * a[j] for j in range(m + 1)), 1 << m)
            for k in range(m + 1)
        ]
        assert all(c[k] == layer[k][mask] for k in range(m + 1))
        assert all(c[k] == 0 for k in range(1, m + 1, 2))
        rebuilt = [sum(c[k] * phi(m, k)[j] for k in range(m + 1)) for j in range(m + 1)]
        assert rebuilt == a
        distinct_B.add(tuple(c[::2]))

    return {
        "n": n,
        "tournaments": N,
        "modes": m // 2 + 1,
        "stalk_cells": len(distinct_B),
        "walsh_raw_sha256": sha256(
            ",".join(map(str, walsh_raw)).encode("ascii")
        ).hexdigest(),
    }


def main():
    print("THM-848 FORWARD-POLYNOMIAL / MOBIUS H-DRIFT REFEREE")
    print("=" * 78)
    print("A_T(x)=sum_r H_(2r)(T)(1-x)^(2r)(1+x)^(n-1-2r)")
    print("z=(1-x)/(1+x): A_T=(1+x)^(n-1) B_T(z)")
    print("S A_T=(1+x)^(n-1)(-2 z d/dz)B_T")
    print("e^(tS): B_T(z) -> B_T(exp(-2t)z)")
    print()

    symbolic = []
    for n in range(2, 15):
        m = n - 1
        for k in range(m + 1):
            p = phi(m, k)
            assert generator(m, p) == [-2 * k * x for x in p]
        symbolic.append((n, m + 1, tuple(-2 * k for k in range(m + 1))))
    print("symbolic_generator_sizes=2..14 all eigenpolynomial assertions passed")
    print("n=14 eigenvalues=" + ",".join(map(str, symbolic[-1][2])))

    print("exhaustive labelled tournament checks:")
    for n in range(2, 6):
        row = exhaustive(n)
        print(
            "  n={n} tournaments={tournaments} even_modes={modes} "
            "stalk_cells={stalk_cells} walsh_raw_sha256={walsh_raw_sha256}".format(**row)
        )

    print()
    print("Tournament Analysis:")
    print("  vertices=even Mobius/Walsh modes H_0,H_2,...")
    print("  pairwise observable=generator decay rate difference")
    print("  switch/gauge=orient slower decay toward faster decay (increasing degree)")
    print("  score histogram at n=14=0:1,1:1,2:1,3:1,4:1,5:1,6:1")
    print("  directed cycles=0; SCCs=7 singletons; Hamiltonian paths=1")
    print("  tie path=degree order (no ties because eigenvalues -4r are distinct)")
    print("  preserves=full averaged flip semigroup and every forward-count coefficient")
    print("  destroys=arc labels,directional currents,runner gaps,scale,owners,loneliness")
    print("  challenged assumption=tournament vertices are spectral proof obligations, not arcs")
    print("source_sha256=" + sha256(Path(__file__).read_bytes()).hexdigest())
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
