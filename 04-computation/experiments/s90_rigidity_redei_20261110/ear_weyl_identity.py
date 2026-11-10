#!/usr/bin/env python3
"""The ear quadratic of a core and its Weyl-differenced exponential sum (README section 5, program item).

For a tournament C on n vertices and a cut x in {0,1}^n let E(x) = H(C + v), where the new vertex v beats exactly
the vertices with x_i = 1.  KNOWN (THM-4097 / THM-4114 section 6; THM-5122 eq. (1) for the protected closure):
E is a quadratic polynomial of x,

    E(x) = h + sum_i x_i s_i + sum_i (1 - x_i) e_i - x^T P (1 - x),   more precisely
    E(x) = sum_{i != j} (1 - x_i) x_j P_ij + sum_j x_j s_j + sum_i (1 - x_i) e_i,

with P_ij = #{ordered pairs of disjoint nonempty paths covering V(C), the first ending at i, the second starting
at j}, s_j = #{Hamiltonian paths of C starting at j}, e_i = #{... ending at i}.  Writing E(x) = c0 + L.x - (1/2) x^T W x
with W_ij = P_ij + P_ji (W_ii = 0), the squared modulus of S(theta) = sum_x exp(2 pi i theta E(x)) is

    |S(theta)|^2 = sum over d in {-1,0,1}^n of
                   exp(2 pi i theta (L.d - (1/2) sum_{i in supp d} (W d)_i)) * prod_{i not in supp d} (1 + exp(-2 pi i theta (W d)_i)).

This script checks (1) the quadratic formula against direct path counts for every cut of every labelled tournament
of orders 2-5, and (2) the differencing identity numerically for random cores of orders 3-6 and random theta.
Run:  python3 ear_weyl_identity.py > ear_weyl_identity.out
"""
import cmath
import itertools
import math
import random
import sys


def tables(n, beats):
    size = 1 << n
    end = [[0] * n for _ in range(size)]          # end[S][v] = paths of C[S] ending at v
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, size):
        row = end[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    end[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
    rb = [0] * n
    for i in range(n):
        for j in range(n):
            if beats[i] >> j & 1:
                rb[j] |= 1 << i
    start = [[0] * n for _ in range(size)]        # start[S][v] = paths of C[S] starting at v
    for v in range(n):
        start[1 << v][v] = 1
    for mask in range(1, size):
        row = start[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = rb[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    start[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
    return end, start


def direct(n, beats, x):
    """H(C + v) by subset DP on n + 1 vertices."""
    b = list(beats) + [0]
    for i in range(n):
        if x >> i & 1:
            b[n] |= 1 << i
        else:
            b[i] |= 1 << n
    m = n + 1
    size = 1 << m
    cnt = [[0] * m for _ in range(size)]
    for v in range(m):
        cnt[1 << v][v] = 1
    for mask in range(1, size):
        row = cnt[mask]
        for v in range(m):
            c = row[v]
            if c:
                nxt = b[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
    return sum(cnt[size - 1])


def quadratic_data(n, beats):
    end, start = tables(n, beats)
    full = (1 << n) - 1
    P = [[0] * n for _ in range(n)]
    for S in range(1, full):
        T = full ^ S
        for i in range(n):
            if end[S][i]:
                for j in range(n):
                    if start[T][j]:
                        P[i][j] += end[S][i] * start[T][j]
    s = [start[full][j] for j in range(n)]
    e = [end[full][i] for i in range(n)]
    return P, s, e


def ear(n, P, s, e, x):
    val = 0
    for i in range(n):
        xi = x >> i & 1
        val += xi * s[i] + (1 - xi) * e[i]
        if xi == 0:
            for j in range(n):
                if j != i and (x >> j & 1):
                    val += P[i][j]
    return val


def main():
    ok = True
    # (1) quadratic formula, exhaustive orders 2..5
    for n in range(2, 6):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        bad = 0
        for bits in range(1 << len(pairs)):
            beats = [0] * n
            for k, (i, j) in enumerate(pairs):
                if bits >> k & 1:
                    beats[i] |= 1 << j
                else:
                    beats[j] |= 1 << i
            P, s, e = quadratic_data(n, beats)
            for x in range(1 << n):
                if ear(n, P, s, e, x) != direct(n, beats, x):
                    bad += 1
        good = (bad == 0)
        ok &= good
        print(f"(1) ear quadratic = direct count, order {n}, all tournaments and cuts: mismatches {bad}", "OK" if good else "FAIL")
    # (2) differencing identity
    random.seed(20261110)
    worst = 0.0
    for trial in range(40):
        n = random.randint(3, 6)
        beats = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    beats[i] |= 1 << j
                else:
                    beats[j] |= 1 << i
        P, s, e = quadratic_data(n, beats)
        W = [[(P[i][j] + P[j][i]) if i != j else 0 for j in range(n)] for i in range(n)]
        # E(x) = c0 + L.x - (1/2) x^T W x :  c0 = sum e_i ;  L_j = s_j - e_j + sum_i P_ij   (from the formula above)
        L = [s[j] - e[j] + sum(P[i][j] for i in range(n) if i != j) for j in range(n)]
        c0 = sum(e)
        for x in range(1 << n):
            xv = [x >> i & 1 for i in range(n)]
            q = c0 + sum(L[i] * xv[i] for i in range(n)) - sum(W[i][j] * xv[i] * xv[j] for i in range(n) for j in range(n)) / 2
            if abs(q - ear(n, P, s, e, x)) > 1e-9:
                ok = False
        theta = random.random()
        S = sum(cmath.exp(2j * math.pi * theta * ear(n, P, s, e, x)) for x in range(1 << n))
        lhs = abs(S) ** 2
        rhs = 0
        for d in itertools.product((-1, 0, 1), repeat=n):
            Wd = [sum(W[i][j] * d[j] for j in range(n)) for i in range(n)]
            phase = sum(L[i] * d[i] for i in range(n)) - 0.5 * sum(Wd[i] for i in range(n) if d[i] != 0)
            term = cmath.exp(2j * math.pi * theta * phase)
            for i in range(n):
                if d[i] == 0:
                    term *= (1 + cmath.exp(-2j * math.pi * theta * Wd[i]))
            rhs += term
        worst = max(worst, abs(lhs - rhs.real) + abs(rhs.imag))
    good = worst < 1e-6
    ok &= good
    print(f"(2) |S(theta)|^2 = differenced sum on 40 random (core, theta), orders 3-6: max discrepancy {worst:.2e}", "OK" if good else "FAIL")
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
