#!/usr/bin/env python3
"""HYP-9273 (B) at q = 31: the Paley orbit count H(P_q)/|Aff| modulo 256, by a subset DP from the arc 0 -> 1.

The stabiliser of the vertex 0 in Aff = {x -> a x + b : a a nonzero square} is {x -> a x}, which acts freely on the
Hamiltonian paths starting at 0 and transitively on their possible second vertices (the nonzero squares).  Hence

    H(P_q) / |Aff| = N_1 := #{Hamiltonian paths of P_q starting with the arc 0 -> 1}.

N_1 is counted by a pull DP over subsets S of the remaining q - 2 vertices (dp[S][j] = paths 0 -> 1 -> ... covering S,
ending at the j-th vertex), in uint8 arithmetic: wrap-around makes every count exact modulo 256.  The layers |S| = k
are processed in order and each layer in parallel (every row S is written by one thread only).  Memory 2^(q-2) (q-2)
bytes: 15.6 GB at q = 31.

Validation (same code): q = 7, 11, 19, 23 against the exact orbit counts 9, 1729, 6857869865, 62293308207033, and the
Paley tournament on GF(27) (modulus x^3 - x - 1) against the earlier result H(P_27)/351 = 1 mod 8.

Session opus-2026-10-14-S33.  Note: 05-knowledge/results/number_tournament_bridges_20261014.md (section 5c)
Reproduce: python3 04-computation/experiments/paley_orbit_count_mod256_20261014.py [qmax=31]   (q = 31: about 16 GB)
"""
import sys, time
import numpy as np
import numba as nb

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


@nb.njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(parallel=True, cache=True)
def count_from_arc(nw, inmask, first):
    """nw = number of remaining vertices W; inmask[j] = bitmask of i in W with W[i] -> W[j];
    first[j] = 1 if 1 -> W[j].  Returns N_1 mod 256 and the dp table is built layer by layer."""
    full = (1 << nw) - 1
    dp = np.zeros((1 << nw, nw), dtype=np.uint8)
    for j in range(nw):
        if first[j]:
            dp[1 << j, j] = 1
    for k in range(2, nw + 1):
        for SS in nb.prange(1 << nw):
            S = np.int64(SS)
            pc = 0
            x = S
            while x:
                x &= x - 1
                pc += 1
            if pc != k:
                continue
            m = S
            while m:
                j = 0
                low = m & (-m)
                while (low >> j) != 1:
                    j += 1
                m &= m - 1
                T = S ^ (1 << j)
                c = T & inmask[j]
                acc = 0
                while c:
                    lowi = c & (-c)
                    i = 0
                    while (lowi >> i) != 1:
                        i += 1
                    acc += dp[T, i]
                    c &= c - 1
                dp[S, j] = acc & 255
    tot = 0
    for j in range(nw):
        tot += dp[full, j]
    return tot & 255


def paley_prime(q):
    sq = {x * x % q for x in range(1, q)}
    verts = list(range(q))
    arc = lambda u, v: (v - u) % q in sq
    return verts, arc, 0, 1


def f27_elements():
    return [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]


def f27_mul(x, y):
    # x = a0 + a1 t + a2 t^2 with t^3 = t + 1
    p = [0] * 5
    for i in range(3):
        for j in range(3):
            p[i + j] = (p[i + j] + x[i] * y[j]) % 3
    for k in (4, 3):
        cf = p[k]
        if cf:
            p[k] = 0
            p[k - 3] = (p[k - 3] + cf) % 3
            p[k - 2] = (p[k - 2] + cf) % 3
    return tuple(p[:3])


def paley_27():
    els = f27_elements()
    zero = (0, 0, 0)
    sq = {f27_mul(x, x) for x in els if x != zero}
    assert len(sq) == 13
    sub = lambda x, y: tuple((x[i] - y[i]) % 3 for i in range(3))
    arc = lambda u, v: sub(v, u) in sq
    return els, arc, zero, (1, 0, 0)


def orbit_count_mod256(verts, arc, v0, v1):
    W = [v for v in verts if v != v0 and v != v1]
    nw = len(W)
    inmask = np.zeros(nw, dtype=np.int64)
    first = np.zeros(nw, dtype=np.int64)
    for j, w in enumerate(W):
        first[j] = 1 if arc(v1, w) else 0
        for i, u in enumerate(W):
            if i != j and arc(u, w):
                inmask[j] |= (1 << i)
    assert arc(v0, v1)
    return int(count_from_arc(nw, inmask, first))


EXACT = {7: 9, 11: 1729, 19: 6857869865, 23: 62293308207033}


def main():
    print('Paley orbit counts N_1 = H(P_q)/|Aff| modulo 256 (paths starting with the arc 0 -> 1)')
    for q in (3, 7, 11, 19, 23):
        t0 = time.time()
        r = orbit_count_mod256(*paley_prime(q))
        ex = EXACT.get(q)
        tag = f'exact {ex} = {ex % 256} mod 256: {check(f"q{q}", r == ex % 256)}' if ex else ''
        print(f'   q = {q:2d}: N_1 = {r:3d} (mod 256), {r % 8} mod 8;  {tag}  [{time.time() - t0:.1f} s]')
    t0 = time.time()
    r27 = orbit_count_mod256(*paley_27())
    print(f'   q = 27 (GF(27)): N_1 = H/351 = {r27:3d} (mod 256), {r27 % 8} mod 8 (earlier: 1 mod 8): '
          f'{check("q27", r27 % 8 == 1)}  [{time.time() - t0:.1f} s]')
    if QMAX >= 31:
        t0 = time.time()
        r31 = orbit_count_mod256(*paley_prime(31))
        print(f'   q = 31: N_1 = H(P_31)/465 = {r31:3d} (mod 256), {r31 % 8} mod 8;  HYP-9273 (B) at q = 31: '
              f'{"HOLDS" if r31 % 8 == 1 else "REFUTED"}  [{time.time() - t0:.1f} s]')
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
