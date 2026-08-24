#!/usr/bin/env python3
"""Deterministic hostile search for parity-coset transference in ranks 4 and 5.

All reported minima are exact.  A trial Gram matrix is Q^T D Q with integral
Q and positive integral diagonal D.  Candidate upper bounds plus

    x_i^2 <= (x^T G x) (G^{-1})_{ii}

give a certified finite box for each primal parity minimum; the analogous
bound for G^{-1} uses G_{ii}.  Trials whose certified box is too large are
skipped, never approximated.
"""

from fractions import Fraction
from itertools import product
from math import isqrt
import random


SEED = 0xC0FFEE14
TRIALS = {4: 100, 5: 150}
NEAR_TRIALS = {4: 750, 5: 1000}
MAX_BOX = 20_000


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def inverse(a):
    n = len(a)
    aug = [[Fraction(a[i][j]) for j in range(n)]
           + [Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if aug[i][j]), None)
        if pivot is None:
            return None
        aug[j], aug[pivot] = aug[pivot], aug[j]
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(n):
            if i == j or not aug[i][j]:
                continue
            scale = aug[i][j]
            aug[i] = [x - scale * y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def quadratic(g, x):
    return sum(x[i] * sum(g[i][j] * x[j] for j in range(len(x)))
               for i in range(len(x)))


def floor_sqrt_fraction(q):
    assert q >= 0
    k = isqrt(q.numerator // q.denominator)
    while Fraction((k + 1) ** 2) <= q:
        k += 1
    while Fraction(k ** 2) > q:
        k -= 1
    return k


def box_size(ranges):
    ans = 1
    for r in ranges:
        ans *= len(r)
    return ans


def exact_pair(g, gi, u):
    d = len(u)

    # The 0/1 representative is a valid primal candidate.
    primal_best = quadratic(g, u)
    primal_caps = [floor_sqrt_fraction(primal_best * gi[i][i]) for i in range(d)]
    primal_ranges = [range(-cap + ((-cap - u[i]) & 1), cap + 1, 2)
                     for i, cap in enumerate(primal_caps)]

    # Any unit vector supported on u is a valid odd dual candidate.
    dual_best = min(gi[i][i] for i in range(d) if u[i])
    dual_caps = [floor_sqrt_fraction(dual_best * g[i][i]) for i in range(d)]
    dual_ranges = [range(-cap, cap + 1) for cap in dual_caps]

    if box_size(primal_ranges) + box_size(dual_ranges) > MAX_BOX:
        return None

    primal_arg = None
    for x in product(*primal_ranges):
        q = quadratic(g, x)
        if q < primal_best:
            primal_best, primal_arg = q, x
    if primal_arg is None:
        primal_arg = tuple(u)

    dual_arg = None
    for z in product(*dual_ranges):
        if sum(a * b for a, b in zip(u, z)) & 1:
            q = quadratic(gi, z)
            if q < dual_best:
                dual_best, dual_arg = q, z
    if dual_arg is None:
        i = next(i for i in range(d) if u[i] and gi[i][i] == dual_best)
        dual_arg = tuple(int(j == i) for j in range(d))

    return primal_best, dual_best, primal_arg, dual_arg,


def random_gram(rng, d):
    while True:
        q = [[rng.randint(-3, 3) for _ in range(d)] for _ in range(d)]
        # Bias away from singularity without forcing unimodularity.
        for i in range(d):
            q[i][i] += rng.choice((-4, 4))
        weights = [rng.randint(1, 8) for _ in range(d)]
        dq = [[weights[i] * q[i][j] for j in range(d)] for i in range(d)]
        g = matmul(transpose(q), dq)
        gi = inverse(g)
        if gi is not None:
            return g, gi


def near_orthogonal_gram(rng, d):
    g = [[0] * d for _ in range(d)]
    for i in range(d):
        for j in range(i):
            g[i][j] = g[j][i] = rng.randint(-6, 6)
    for i in range(d):
        # Strict diagonal dominance is an exact positive-definiteness proof.
        g[i][i] = 40 + rng.randint(0, 12) + sum(abs(g[i][j]) for j in range(d) if j != i)
    return g, inverse(g)


def main():
    rng = random.Random(SEED)
    print(f"SEED={SEED}")
    named = {
        "D4": [[2, -1, 0, 0], [-1, 2, -1, -1],
               [0, -1, 2, 0], [0, -1, 0, 2]],
        "A4": [[2, -1, 0, 0], [-1, 2, -1, 0],
               [0, -1, 2, -1], [0, 0, -1, 2]],
        "D5": [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0],
               [0, -1, 2, -1, -1], [0, 0, -1, 2, 0],
               [0, 0, -1, 0, 2]],
    }
    for name, g in named.items():
        gi = inverse(g)
        best = None
        for bits in product((0, 1), repeat=len(g)):
            if not any(bits):
                continue
            result = exact_pair(g, gi, list(bits))
            mp, md, xp, zd = result
            score = mp * md
            if best is None or score > best[0]:
                best = score, bits, mp, md, xp, zd
        print(f"named={name} all_characteristics={2**len(g)-1} "
              f"best_product={best[0]} u={best[1]} primal={best[2]} "
              f"dual_odd={best[3]}")
    for d in (4, 5):
        identity = [[int(i == j) for j in range(d)] for i in range(d)]
        control = exact_pair(identity, identity, [1] * d)
        best = (control[0] * control[1], identity, tuple([1] * d),
                control[0], control[1], control[2], control[3])
        checked = skipped = 0
        families = ((TRIALS[d], random_gram),
                    (NEAR_TRIALS[d], near_orthogonal_gram))
        for count, maker in families:
          for _ in range(count):
            g, gi = maker(rng, d)
            # Mix a random characteristic with the high-support one.
            us = [[1] * d]
            for __ in range(2):
                us.append([rng.randint(0, 1) for _ in range(d)])
                if not any(us[-1]):
                    us[-1][rng.randrange(d)] = 1
            for u in us:
                result = exact_pair(g, gi, u)
                if result is None:
                    skipped += 1
                    continue
                checked += 1
                mp, md, xp, zd = result
                score = mp * md
                if best is None or score > best[0]:
                    best = (score, g, tuple(u), mp, md, xp, zd)
        score, g, u, mp, md, xp, zd = best
        print(f"d={d} checked={checked} skipped={skipped} best_product={score} ratio_to_d={float(score/d):.12f}")
        print(f"u={u} primal={mp} at {xp} dual_odd={md} at {zd}")
        print("G=" + repr([[int(x) if Fraction(x).denominator == 1 else str(x)
                            for x in row] for row in g]))
        if score > d:
            print("HOSTILE_COUNTEREXAMPLE=YES")
        else:
            print("HOSTILE_COUNTEREXAMPLE=NO")


if __name__ == "__main__":
    main()
