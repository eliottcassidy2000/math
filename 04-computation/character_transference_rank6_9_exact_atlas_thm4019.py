#!/usr/bin/env python3
"""Finite-exact hostile search at the arbitrary-lattice boundary of THM-4015.

For a positive-definite integral Gram matrix G and nonzero u in F_2^d, put

  A(G,u)=min_{x=u (mod 2)} x^T G x,
  B(G,u)=min_{u.z=1 (mod 2)} z^T G^{-1} z.

The sharper arbitrary-lattice candidate in THM-4015 is A(G,u)B(G,u)<=d.
Every minimum printed here is certified by rational LDL sphere enumeration.
Instances exceeding the declared node budget are skipped, never approximated.
The named A/D/E root lattices are required to finish without a node cutoff.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import isqrt
import random


SEED = 0x4015A11
NODE_LIMIT = 2_000_000
RANDOM_TRIALS = {6: 96, 7: 80, 8: 64, 9: 48}


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt]
            for row in a]


def inverse(a):
    n = len(a)
    aug = [[Fraction(a[i][j]) for j in range(n)]
           + [Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if aug[i][j]), None)
        if pivot is None:
            return None
        aug[j], aug[pivot] = aug[pivot], aug[j]
        q = aug[j][j]
        aug[j] = [x / q for x in aug[j]]
        for i in range(n):
            if i == j or not aug[i][j]:
                continue
            q = aug[i][j]
            aug[i] = [x - q * y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def quadratic(g, x):
    return sum(x[i] * sum(g[i][j] * x[j] for j in range(len(x)))
               for i in range(len(x)))


def ldl(g):
    """Return exact L,D with G=L diag(D) L^T, or None if not SPD."""
    n = len(g)
    ll = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    dd = [Fraction(0) for _ in range(n)]
    for j in range(n):
        dd[j] = Fraction(g[j][j]) - sum(
            ll[j][k] * ll[j][k] * dd[k] for k in range(j)
        )
        if dd[j] <= 0:
            return None
        for i in range(j + 1, n):
            ll[i][j] = (Fraction(g[i][j]) - sum(
                ll[i][k] * ll[j][k] * dd[k] for k in range(j)
            )) / dd[j]
    return ll, dd


def floor_fraction(q):
    return q.numerator // q.denominator


def floor_sqrt_fraction(q):
    assert q >= 0
    k = isqrt(q.numerator // q.denominator)
    while Fraction((k + 1) ** 2) <= q:
        k += 1
    while Fraction(k * k) > q:
        k -= 1
    return k


class NodeLimit(Exception):
    pass


def exact_minimum(g, initial, coordinate_parity=None, odd_u=None,
                  node_limit=NODE_LIMIT):
    """Certify a constrained minimum by exact LDL branch and bound."""
    decomp = ldl(g)
    assert decomp is not None
    ll, dd = decomp
    n = len(g)
    best = Fraction(quadratic(g, initial))
    best_x = tuple(initial)
    x = [0] * n
    nodes = 0

    def rec(i, spent):
        nonlocal best, best_x, nodes
        nodes += 1
        if nodes > node_limit:
            raise NodeLimit
        if i < 0:
            if odd_u is not None and not (sum(a * b for a, b in zip(odd_u, x)) & 1):
                return
            q = Fraction(quadratic(g, x))
            if q < best:
                best, best_x = q, tuple(x)
            return

        remaining = best - spent
        if remaining <= 0:
            return
        centre = sum(ll[j][i] * x[j] for j in range(i + 1, n))
        radius2 = remaining / dd[i]
        # Over-enumerate by two integers at both ends, then filter exactly.
        anchor = floor_fraction(-centre)
        radius = floor_sqrt_fraction(radius2) + 2
        candidates = []
        parity = None if coordinate_parity is None else coordinate_parity[i]
        for z in range(anchor - radius, anchor + radius + 2):
            if parity is not None and ((z - parity) & 1):
                continue
            cost = dd[i] * (Fraction(z) + centre) ** 2
            if cost < remaining:
                candidates.append((cost, z))
        candidates.sort(key=lambda item: (item[0], abs(item[1]), item[1]))
        for cost, z in candidates:
            x[i] = z
            rec(i - 1, spent + cost)

    rec(n - 1, Fraction(0))
    return best, best_x, nodes


def exact_pair(g, u, node_limit=NODE_LIMIT):
    gi = inverse(g)
    assert gi is not None
    primal_initial = tuple(int(a) for a in u)
    dual_i = next(i for i, a in enumerate(u) if a)
    dual_initial = tuple(int(i == dual_i) for i in range(len(u)))
    ap, xp, np = exact_minimum(
        g, primal_initial, coordinate_parity=u, node_limit=node_limit
    )
    bd, zd, nd = exact_minimum(
        gi, dual_initial, odd_u=u, node_limit=node_limit
    )
    return ap, bd, xp, zd, np + nd


def cartan_a(d):
    g = [[0] * d for _ in range(d)]
    for i in range(d):
        g[i][i] = 2
        if i:
            g[i][i - 1] = g[i - 1][i] = -1
    return g


def cartan_d(d):
    assert d >= 4
    g = [[0] * d for _ in range(d)]
    for i in range(d):
        g[i][i] = 2
    for i in range(d - 3):
        g[i][i + 1] = g[i + 1][i] = -1
    branch = d - 3
    g[branch][d - 2] = g[d - 2][branch] = -1
    g[branch][d - 1] = g[d - 1][branch] = -1
    return g


def cartan_e(d):
    assert d in (6, 7, 8)
    # Arms of lengths (2,d-4,1), hence E6/E7/E8.
    g = [[0] * d for _ in range(d)]
    for i in range(d):
        g[i][i] = 2
    for i in range(d - 2):
        g[i][i + 1] = g[i + 1][i] = -1
    g[2][d - 1] = g[d - 1][2] = -1
    return g


def named_scan(name, g):
    d = len(g)
    best = None
    nodes = 0
    for bits in product((0, 1), repeat=d):
        if not any(bits):
            continue
        ap, bd, xp, zd, used = exact_pair(g, bits)
        nodes += used
        score = ap * bd
        row = (score, bits, ap, bd, xp, zd)
        if best is None or row[:2] > best[:2]:
            best = row
    assert best is not None
    score, bits, ap, bd, xp, zd = best
    print(
        f"named={name};d={d};characters={2**d-1};nodes={nodes};"
        f"best={score};ratio={float(score/d):.12f};u={bits};"
        f"primal={ap}@{xp};dual={bd}@{zd}", flush=True
    )
    if score > d:
        print("HOSTILE_COUNTEREXAMPLE=YES", flush=True)
        print("G=" + repr(g), flush=True)
    return best


def diagonal_dominant_gram(rng, d):
    g = [[0] * d for _ in range(d)]
    for i in range(d):
        for j in range(i):
            g[i][j] = g[j][i] = rng.randint(-8, 8)
    for i in range(d):
        g[i][i] = rng.randint(2, 15) + sum(
            abs(g[i][j]) for j in range(d) if j != i
        )
    return g


def wishart_gram(rng, d):
    a = [[rng.randint(-3, 3) for _ in range(d)] for __ in range(d + 2)]
    g = matmul(transpose(a), a)
    for i in range(d):
        g[i][i] += rng.randint(1, 5)
    return g


def rank_one_gram(rng, d):
    v = [rng.randint(-5, 5) for _ in range(d)]
    scale = rng.randint(1, 8)
    return [[scale * int(i == j) + v[i] * v[j] for j in range(d)]
            for i in range(d)]


def random_scan(rng, d):
    identity = [[int(i == j) for j in range(d)] for i in range(d)]
    control = exact_pair(identity, [1] * d)
    best = (control[0] * control[1], identity, tuple([1] * d), *control[:4])
    checked = skipped = nodes = 0
    makers = (diagonal_dominant_gram, wishart_gram, rank_one_gram)
    for trial in range(RANDOM_TRIALS[d]):
        g = makers[trial % len(makers)](rng, d)
        supports = [[1] * d]
        for _ in range(3):
            u = [rng.randint(0, 1) for _ in range(d)]
            if not any(u):
                u[rng.randrange(d)] = 1
            supports.append(u)
        for u in supports:
            try:
                ap, bd, xp, zd, used = exact_pair(g, u)
            except NodeLimit:
                skipped += 1
                continue
            checked += 1
            nodes += used
            score = ap * bd
            row = (score, g, tuple(u), ap, bd, xp, zd)
            if row[0] > best[0]:
                best = row
    score, g, u, ap, bd, xp, zd = best
    print(
        f"random_d={d};checked={checked};skipped={skipped};nodes={nodes};"
        f"best={score};ratio={float(score/d):.12f};u={u};"
        f"primal={ap}@{xp};dual={bd}@{zd}", flush=True
    )
    if score > d:
        print("HOSTILE_COUNTEREXAMPLE=YES", flush=True)
        print("G=" + repr(g), flush=True)
    return checked, skipped, best


def invariance_controls():
    """Exercise scaling and a nontrivial unimodular coordinate shear."""
    g = cartan_e(6)
    u = (1, 0, 1, 1, 0, 1)
    base = exact_pair(g, u)
    scaled = exact_pair([[7 * a for a in row] for row in g], u)
    assert base[0] * base[1] == scaled[0] * scaled[1]

    q = [[int(i == j) for j in range(6)] for i in range(6)]
    q[0][1] = 1
    transformed = matmul(transpose(q), matmul(g, q))
    # x=Q^{-1}x' gives u'=Q^{-1}u modulo two; Q^{-1}=I-E_01.
    up = list(u)
    up[0] ^= up[1]
    sheared = exact_pair(transformed, up)
    assert base[0] * base[1] == sheared[0] * sheared[1]
    print(
        f"invariance_control=PASS;product={base[0]*base[1]};"
        f"base_nodes={base[4]};scaled_nodes={scaled[4]};sheared_nodes={sheared[4]}",
        flush=True,
    )


def main():
    rng = random.Random(SEED)
    print(
        f"SEED={SEED};status=FINITE-EXACT;node_limit={NODE_LIMIT};"
        "enumerator=rational_LDL_sphere", flush=True
    )
    invariance_controls()
    named_hostiles = []
    for d in range(6, 10):
        for name, g in ((f"A{d}", cartan_a(d)),
                        (f"D{d}", cartan_d(d))):
            best = named_scan(name, g)
            if best[0] > d:
                named_hostiles.append((name, best[0]))
        if d <= 8:
            name, g = f"E{d}", cartan_e(d)
            best = named_scan(name, g)
            if best[0] > d:
                named_hostiles.append((name, best[0]))
    total_checked = total_skipped = 0
    global_best = None
    for d in range(6, 10):
        checked, skipped, best = random_scan(rng, d)
        total_checked += checked
        total_skipped += skipped
        ratio = best[0] / d
        if global_best is None or ratio > global_best[0]:
            global_best = ratio, d, best[0]
    print(
        f"random_total_checked={total_checked};random_total_skipped={total_skipped};"
        f"best_ratio={global_best[0]};best_dimension={global_best[1]};"
        f"best_product={global_best[2]}", flush=True
    )
    print(f"named_hostiles={named_hostiles}", flush=True)
    if named_hostiles:
        print("RESULT=COUNTEREXAMPLE_FOUND_FINITE_EXACT", flush=True)
    else:
        print("RESULT=NO_COUNTEREXAMPLE_FOUND_FINITE_EXACT", flush=True)


if __name__ == "__main__":
    main()
