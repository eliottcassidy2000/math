#!/usr/bin/env python3
"""Audit (a) of opus-2026-10-14-S35: independent recomputation of the Paley signed self-avoiding run sums.

A_L(q) = sum over walks v_0..v_L of pairwise distinct points of F_q of prod_i chi(v_{i+1} - v_i), chi quadratic.

Method (independent of the session's partition-lattice compiler): a signed subset DP over self-avoiding walks.
For even L the summand is invariant under AGL(1,q) (chi(a)^L = 1), which acts simply transitively on the
ordered pairs (v_0, v_1), v_0 != v_1.  Hence A_L = q(q-1) B_L with B_L the sum over walks 0 -> 1 -> v_2 -> ... -> v_L.
The DP runs over layers k = |{v_2..v_{k+1}}|; a state is (k-subset S of F_q minus {0,1}, last vertex), indexed by the
colex rank of S (Gosper enumeration) times n plus the last vertex.  Exact int64 arithmetic.
Works for prime fields and for F_27 (explicit field tables), i.e. also at a prime power outside the session's
interpolation set.  Validation: a numba DFS brute force (no symmetry reduction) at small q and L.

Comparisons: the closed forms of HYP-9274 (A_2..A_8), and S35's A_10, A_12, with a_q = coefficients of
eta(2z)^4 eta(4z)^4 computed here (Frobenius trace alpha^r + beta^r at q = p^r).
Audit of S35, 2026-10-14.  Independent code: nothing is imported from the session's scripts.
"""
import sys, time
import numpy as np
import numba as nb


# ---------------- finite fields ----------------
def prime_field_tables(p):
    sq = set((x * x) % p for x in range(1, p))
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, p)], dtype=np.int64)
    W = np.zeros((p, p), dtype=np.int64)
    for a in range(p):
        for b in range(p):
            W[a, b] = chi[(b - a) % p]
    return W


def f27_tables():
    # F_27 = F_3[t]/(t^3 - t + 1); element c0 + c1 t + c2 t^2 <-> c0 + 3 c1 + 9 c2.  t^3 = t - 1.
    def dig(a):
        return [a % 3, (a // 3) % 3, (a // 9) % 3]

    def undig(d):
        return d[0] + 3 * d[1] + 9 * d[2]

    def sub(a, b):
        da, db = dig(a), dig(b)
        return undig([(x - y) % 3 for x, y in zip(da, db)])

    def mul(a, b):
        da, db = dig(a), dig(b)
        prod = [0] * 5
        for i in range(3):
            for j in range(3):
                prod[i + j] += da[i] * db[j]
        # reduce t^4 = t^2 - t, t^3 = t - 1
        for deg in (4, 3):
            c = prod[deg]
            prod[deg] = 0
            prod[deg - 2] += c
            prod[deg - 3] -= c
        return undig([x % 3 for x in prod[:3]])

    # irreducibility check: no root in F_3
    for t in range(3):
        assert (t ** 3 - t + 1) % 3 != 0
    sq = set(mul(x, x) for x in range(1, 27))
    assert len(sq) == 13
    chi = [0] + [1 if x in sq else -1 for x in range(1, 27)]
    assert chi[sub(0, 1)] == -1  # -1 is a non-square (27 = 3 mod 4)
    W = np.zeros((27, 27), dtype=np.int64)
    for a in range(27):
        for b in range(27):
            W[a, b] = chi[sub(b, a)]
    return W


# ---------------- DP over self-avoiding walks ----------------
@nb.njit(cache=True)
def saw_B(q, W, L):
    """B_L = sum over walks 0 -> 1 -> v_2 .. v_L (v_2..v_L distinct, not in {0,1}) of prod_{i>=1} W[v_i, v_{i+1}]."""
    n = q - 2
    K = L - 1
    if K == 1:
        s = 0
        for x in range(n):
            s += W[1, x + 2]
        return s
    C = np.zeros((n + 3, n + 3), dtype=np.int64)
    for i in range(n + 3):
        C[i, 0] = 1
        for j in range(1, i + 1):
            C[i, j] = C[i - 1, j - 1] + C[i - 1, j]
    cur = np.zeros(n * n, dtype=np.int64)
    for x in range(n):
        cur[x * n + x] = W[1, x + 2]          # rank({x}) = C(x, 1) = x
    elems = np.zeros(n, dtype=np.int64)
    plow = np.zeros(n + 2, dtype=np.int64)
    phigh = np.zeros(n + 2, dtype=np.int64)
    k = 1
    while True:
        ncur = C[n, k]
        if k == K - 1:
            total = 0
            mask = (1 << k) - 1
            for r in range(ncur):
                for last in range(n):
                    if (mask >> last) & 1:
                        val = cur[r * n + last]
                        if val != 0:
                            s = 0
                            for x in range(n):
                                if not (mask >> x) & 1:
                                    s += W[last + 2, x + 2]
                            total += val * s
                c = mask & -mask
                rr = mask + c
                mask = (((rr ^ mask) >> 2) // c) | rr
            return total
        nxt = np.zeros(C[n, k + 1] * n, dtype=np.int64)
        mask = (1 << k) - 1
        for r in range(ncur):
            idx = 0
            for b in range(n):
                if (mask >> b) & 1:
                    elems[idx] = b
                    idx += 1
            plow[0] = 0
            for i in range(k):
                plow[i + 1] = plow[i] + C[elems[i], i + 1]
            phigh[k] = 0
            for i in range(k - 1, -1, -1):
                phigh[i] = phigh[i + 1] + C[elems[i], i + 2]
            for li in range(k):
                last = elems[li]
                val = cur[r * n + last]
                if val == 0:
                    continue
                j = 0
                for x in range(n):
                    if (mask >> x) & 1:
                        j += 1
                        continue
                    rn = plow[j] + C[x, j + 1] + phigh[j]
                    nxt[rn * n + x] += val * W[last + 2, x + 2]
            c = mask & -mask
            rr = mask + c
            mask = (((rr ^ mask) >> 2) // c) | rr
        cur = nxt
        k += 1


@nb.njit(cache=True)
def brute_A(q, W, L):
    """direct DFS over all injective walks v_0..v_L (no symmetry used)."""
    total = 0
    path = np.zeros(L + 1, dtype=np.int64)
    sign = np.zeros(L + 1, dtype=np.int64)
    choice = np.zeros(L + 1, dtype=np.int64)
    used = np.zeros(q, dtype=np.int64)
    for v0 in range(q):
        path[0] = v0
        sign[0] = 1
        used[v0] = 1
        depth = 1
        choice[1] = -1
        while depth >= 1:
            choice[depth] += 1
            if choice[depth] >= q:
                depth -= 1
                if depth >= 1:
                    used[path[depth]] = 0
                continue
            x = choice[depth]
            if used[x]:
                continue
            w = W[path[depth - 1], x]
            if w == 0:
                continue
            if depth == L:
                total += sign[depth - 1] * w
                continue
            path[depth] = x
            sign[depth] = sign[depth - 1] * w
            used[x] = 1
            depth += 1
            choice[depth] = -1
        used[v0] = 0
    return total


# ---------------- modular form a_q ----------------
def beukers_coeffs(N):
    """a_n for eta(2z)^4 eta(4z)^4 = q prod (1-q^{2m})^4 (1-q^{4m})^4, n < N."""
    P = [0] * N
    P[0] = 1
    for m in range(1, N):
        for step in (2 * m, 4 * m):
            if step >= N:
                continue
            for _ in range(4):
                for i in range(N - 1, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P[:N - 1]


def frob_trace(ap, p, r, k=4):
    """alpha^r + beta^r with alpha + beta = a_p, alpha beta = p^(k-1)."""
    s0, s1 = 2, ap
    if r == 0:
        return 2
    for _ in range(r - 1):
        s0, s1 = s1, ap * s1 - p ** (k - 1) * s0
    return s1


def closed_forms(q, aq):
    f = {}
    f[2] = q * (q - 1)
    f[4] = 2 * q * (q - 1) * (q - 3)
    f[6] = q * (q - 1) * (q - 3) * (5 * q - 29)
    f[8] = 2 * q * (q - 1) * (q - 3) * (q - 7) * (7 * q - 55)
    f[10] = 2 * q * (q - 1) * (21 * q ** 4 - 639 * q ** 3 + 7028 * q ** 2 - 32098 * q + 48662) + 34 * q * (q - 1) * aq
    P12 = 4 * q * (q - 1) * (33 * q ** 5 - 1540 * q ** 4 + 28243 * q ** 3 - 249550 * q ** 2 + 1027794 * q - 1488872)
    f[12] = P12 + 28 * q * (q - 1) * (17 * q - 251) * aq
    return f


if __name__ == '__main__':
    a = beukers_coeffs(200)
    print('eta(2z)^4 eta(4z)^4: a_3, a_5, a_7, a_11, a_19, a_23 =', a[3], a[5], a[7], a[11], a[19], a[23])
    # 1. validation of DP + normalisation against brute force
    t0 = time.time()
    for q in (7, 11, 13):
        W = prime_field_tables(q)
        for L in (2, 4, 6, 8):
            if L + 1 > q:
                continue
            bA = int(brute_A(q, W, L))
            dA = q * (q - 1) * int(saw_B(q, W, L))
            print(f'  validation q = {q:2d}, L = {L}: brute force {bA}, DP {dA}, equal: {bA == dA}', flush=True)
            assert bA == dA
    W11 = prime_field_tables(11)
    b10 = int(brute_A(11, W11, 10))
    print(f'  brute force A_10(11) = {b10} (S35 claims 478720; closed form {closed_forms(11, a[11])[10]})  [{time.time() - t0:.1f} s]', flush=True)
    # 2. the targets
    targets = [('19', 19, prime_field_tables(19), a[19]), ('23', 23, prime_field_tables(23), a[23]),
               ('27 = 3^3', 27, f27_tables(), frob_trace(a[3], 3, 3))]
    if len(sys.argv) > 1 and sys.argv[1] == 'quick':
        targets = targets[:2]
    for name, q, W, aq in targets:
        f = closed_forms(q, aq)
        for L in (2, 4, 6, 8, 10, 12):
            t0 = time.time()
            A = q * (q - 1) * int(saw_B(q, W, L))
            print(f'q = {name}: A_{L} = {A}  closed form = {f[L]}  agree: {A == f[L]}  (a_q = {aq})  [{time.time() - t0:.1f} s]', flush=True)
