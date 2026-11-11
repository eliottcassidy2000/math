#!/usr/bin/env python3
"""Scratch probe (not canon): the two-apex family G(2,0).  Spine s_1 -> ... -> s_n transitive, apexes u, v; letter
(a_i, b_i): a_i = 1 iff u -> s_i, b_i = 1 iff v -> s_i; arc u -> v (the other orientation is the swap of the roles).
H = number of Hamiltonian paths = sum over arrangements [R1] x [R2] y [R3], {x, y} = {u, v}, R_k increasing runs
partitioning the spine, with last(R1) -> x -> first(R2), last(R2) -> y -> first(R3); R2 empty needs x -> y.
Counts, for all words of length exactly n, the multiplicity of every odd value below a bound, and checks the
transfer count against brute force for small n.
"""
import itertools
import sys

import numpy as np
from numba import njit


def brute(n, a, b):
    V = n + 2
    u, v = n, n + 1
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if a[i]:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
        if b[i]:
            adj[v][i] = 1
        else:
            adj[i][v] = 1
    adj[u][v] = 1
    dp = [[0] * V for _ in range(1 << V)]
    for x in range(V):
        dp[1 << x][x] = 1
    for S in range(1, 1 << V):
        for x in range(V):
            c = dp[S][x]
            if c:
                for y in range(V):
                    if adj[x][y] and not (S >> y) & 1:
                        dp[S | (1 << y)][y] += c
    return sum(dp[(1 << V) - 1])


# state index: r1 in {0 empty, 1 last ok, 2 last bad}, r2 in {0, 1, 2}, r3 in {0, 1}; arrangement 0: x = u, y = v;
# arrangement 1: x = v, y = u.  36 states.
@njit(cache=True)
def step(vec, a, b, out):
    for k in range(36):
        out[k] = 0
    for arr in range(2):
        bx = a if arr == 0 else b            # bit of x at this vertex (1: x -> s)
        by = b if arr == 0 else a
        for r1 in range(3):
            for r2 in range(3):
                for r3 in range(2):
                    c = vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
                    if c == 0:
                        continue
                    # put the vertex in R1: allowed only while R2 and R3 are ... (no order constraint between runs)
                    n1 = 1 if bx == 0 else 2
                    out[arr * 18 + n1 * 6 + r2 * 2 + r3] += c
                    # put in R2: if R2 empty the vertex is first(R2): needs x -> s (bx = 1)
                    if r2 != 0 or bx == 1:
                        n2 = 1 if by == 0 else 2
                        out[arr * 18 + r1 * 6 + n2 * 2 + r3] += c
                    # put in R3: if R3 empty the vertex is first(R3): needs y -> s (by = 1)
                    if r3 != 0 or by == 1:
                        out[arr * 18 + r1 * 6 + r2 * 2 + 1] += c


@njit(cache=True)
def final(vec):
    tot = 0
    for arr in range(2):
        for r1 in range(3):
            if r1 == 2:
                continue
            for r2 in range(3):
                if r2 == 2:
                    continue
                if r2 == 0 and arr == 1:      # R2 empty needs x -> y; the arc is u -> v, so x = u only
                    continue
                for r3 in range(2):
                    tot += vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
    return tot


@njit(cache=True)
def count_word(n, a, b):
    vec = np.zeros(36, np.int64)
    out = np.zeros(36, np.int64)
    vec[0] = 1
    vec[18] = 1
    for i in range(n):
        step(vec, a[i], b[i], out)
        for k in range(36):
            vec[k] = out[k]
    return final(vec)


@njit(cache=True)
def census(n, bound, mult):
    """all 4^n words of length n: multiplicity of each odd value < bound."""
    stack = np.zeros((n + 1, 36), np.int64)
    stack[0, 0] = 1
    stack[0, 18] = 1
    letters = np.zeros(n, np.int64)
    depth = 0
    letters[0] = -1
    big = 0
    out = np.zeros(36, np.int64)
    while depth >= 0:
        letters[depth] += 1
        if letters[depth] == 4:
            depth -= 1
            continue
        l = letters[depth]
        step(stack[depth], l & 1, (l >> 1) & 1, out)
        for k in range(36):
            stack[depth + 1, k] = out[k]
        if depth + 1 == n:
            h = final(stack[n])
            if h < bound:
                mult[h >> 1] += 1
            else:
                big += 1
        else:
            depth += 1
            letters[depth] = -1
    return big


def main():
    ok = True
    for n in (1, 2, 3, 4, 5):
        for a in itertools.product((0, 1), repeat=n):
            for b in itertools.product((0, 1), repeat=n):
                h1 = brute(n, a, b)
                h2 = count_word(n, np.array(a, np.int64), np.array(b, np.int64))
                if h1 != h2:
                    ok = False
                    print("MISMATCH", n, a, b, h1, h2)
                    break
            if not ok:
                break
    print("transfer count = brute force for n <= 5:", ok)
    if not ok:
        sys.exit(1)
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    bound = 1 << 21
    total = np.zeros(bound >> 1, np.int64)
    for n in range(1, N + 1):
        mult = np.zeros(bound >> 1, np.int64)
        big = census(n, bound, mult)
        total += mult
        vals = np.nonzero(mult)[0]
        mx = 2 * vals.max() + 1 if len(vals) else 0
        print(f"n = {n:2d}: {4 ** n} words, {big} with H >= 2^21, max H below bound {mx}, distinct odd values {len(vals)}")
    # coverage and multiplicity by dyadic window, words of length <= N
    miss_all = [2 * i + 1 for i in range(min(len(total), 1 << 16)) if total[i] == 0]
    print(f"odd N < 2^17 not attained by words of length <= {N}: {miss_all[:40]}{' ...' if len(miss_all) > 40 else ''} ({len(miss_all)} in all)")
    for k in range(3, 21):
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg = total[lo:hi]
        print(f"  odd N in [2^{k}, 2^{k+1}): min multiplicity {int(seg.min())}, median {int(np.median(seg))}, mean {seg.mean():.1f}, zeros {int((seg == 0).sum())}")


if __name__ == "__main__":
    main()
