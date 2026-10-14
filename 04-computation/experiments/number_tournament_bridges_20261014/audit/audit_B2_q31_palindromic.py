#!/usr/bin/env python3
"""Referee check of the late additions (sections 5b, 5c; HYP-9273 B at q = 31 and B').

(1) N_1 = number of Hamiltonian paths of P_q that start with the arc 0 -> 1 (= H/|Aff|), computed modulo 256 by a
    LAYERED pull DP over the subsets of W = F_q \\ {0, 1} indexed by their colex rank: only two adjacent layers are
    in memory (max about 2.3 GB at q = 31; the author's full table needs 15.6 GB).  Validated at q = 7, 11, 19, 23 against
    the exact orbit counts and at q = 27 against the referee's own full DP (H(P_27) = 95 mod 256, |Aff| = 351).
(2) pal(q) = number of (x -> -x)-palindromic Hamiltonian paths: own DP over the first half (one element from each
    pair {x, -x}, ending at a non-residue); brute-force enumeration of all Hamiltonian paths for q = 7, 11 as a check;
    pal mod 2^13 (q-1)/2 for q = 43, 47; the split H/|Aff| = 2F + P with P = 2 pal/(q-1).
"""
import sys, time
import numpy as np
import numba as nb

M = 29
BINOM = np.zeros((M + 2, M + 2), dtype=np.int64)
for n in range(M + 2):
    BINOM[n, 0] = 1
    for k in range(1, n + 1):
        BINOM[n, k] = BINOM[n - 1, k - 1] + (BINOM[n - 1, k] if k <= n - 1 else 0)


@nb.njit(cache=True)
def unrank_colex(r, size, binom):
    mask = np.int64(0)
    for j in range(size, 0, -1):
        s = j - 1
        while binom[s + 1, j] <= r:
            s += 1
        r -= binom[s, j]
        mask |= np.int64(1) << s
    return mask


@nb.njit(parallel=True, cache=True)
def layer_step(f_old, i, m, inmask, binom, nchunks):
    """from layer i (subsets of size i, shape (C(m,i), i)) to layer i+1"""
    size = i + 1
    total = binom[m, size]
    f_new = np.zeros((total, size), dtype=np.uint8)
    chunk = (total + nchunks - 1) // nchunks
    for c in nb.prange(nchunks):
        r0 = c * chunk
        r1 = min(total, r0 + chunk)
        if r0 >= r1:
            continue
        mask = unrank_colex(r0, size, binom)
        el = np.zeros(size, dtype=np.int64)
        pre = np.zeros(size + 1, dtype=np.int64)
        suf = np.zeros(size + 1, dtype=np.int64)
        for r in range(r0, r1):
            # elements in increasing order
            t = 0
            x = mask
            b = 0
            while x:
                if x & 1:
                    el[t] = b
                    t += 1
                x >>= 1
                b += 1
            # pre[k] = sum_{j<k} C(el_j, j+1);  suf[k] = sum_{j>k} C(el_j, j)
            pre[0] = 0
            for k in range(size):
                pre[k + 1] = pre[k] + binom[el[k], k + 1]
            suf[size - 1] = 0
            for k in range(size - 2, -1, -1):
                suf[k] = suf[k + 1] + binom[el[k + 1], k + 1]
            for k in range(size):
                v = el[k]
                rk = pre[k] + suf[k]
                acc = np.int64(0)
                im = inmask[v]
                for j in range(size):
                    if j == k:
                        continue
                    u = el[j]
                    if (im >> u) & 1:
                        pos = j if j < k else j - 1
                        acc += np.int64(f_old[rk, pos])
                f_new[r, k] = np.uint8(acc & 255)
            # Gosper: next mask with the same popcount
            if r + 1 < r1:
                cbit = mask & (-mask)
                rr = mask + cbit
                mask = (((rr ^ mask) >> 2) // cbit) | rr
    return f_new


def n1_mod256(elements, arc):
    """elements: list of field elements; arc(x, y) True iff x -> y; vertices '0' and '1' are elements[0], elements[1]"""
    z, one = elements[0], elements[1]
    W = elements[2:]
    m = len(W)
    inmask = np.zeros(m, dtype=np.int64)
    for a, v in enumerate(W):
        for b, u in enumerate(W):
            if a != b and arc(u, v):
                inmask[a] |= 1 << b
    binom = np.zeros((m + 2, m + 2), dtype=np.int64)
    for n in range(m + 2):
        binom[n, 0] = 1
        for k in range(1, n + 1):
            binom[n, k] = binom[n - 1, k - 1] + (binom[n - 1, k] if k <= n - 1 else 0)
    f = np.zeros((m, 1), dtype=np.uint8)
    for a, v in enumerate(W):
        f[a, 0] = 1 if arc(one, v) else 0
    for i in range(1, m):
        f = layer_step(f, i, m, inmask, binom, 192)
    return int(f[0].astype(np.int64).sum() % 256)


def prime_paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    return list(range(q)), (lambda x, y: (y - x) % q in sq)


def gf27_paley():
    def mul(a, b):
        A = [(a // 3 ** i) % 3 for i in range(3)]
        B = [(b // 3 ** i) % 3 for i in range(3)]
        pr = [0] * 5
        for i in range(3):
            for j in range(3):
                pr[i + j] += A[i] * B[j]
        for k in (4, 3):                    # y^3 = y + 2 for y^3 + 2y + 1
            c = pr[k]
            pr[k] = 0
            pr[k - 3] += 2 * c
            pr[k - 2] += c
        return sum((pr[i] % 3) * 3 ** i for i in range(3))

    def sub(a, b):
        return sum(((((a // 3 ** i) % 3) - ((b // 3 ** i) % 3)) % 3) * 3 ** i for i in range(3))
    sq = {mul(x, x) for x in range(1, 27)}
    return list(range(27)), (lambda x, y: sub(y, x) in sq)     # element 0 = zero, element 1 = one


@nb.njit(cache=True)
def pal_dp(q, isq, mod):
    """number of directed paths v_1 -> ... -> v_m (m = (q-1)/2) using one element of each pair {x, -x},
    ending at a non-residue; arithmetic modulo mod (mod <= 2^31)"""
    m = (q - 1) // 2
    pair = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        pair[x] = min(x, q - x) - 1
    N = 1 << m
    dp = np.zeros((N, q), dtype=np.int64)
    for v in range(1, q):
        dp[1 << pair[v], v] = 1
    for mask in range(1, N):
        for v in range(1, q):
            c = dp[mask, v]
            if c == 0:
                continue
            for w in range(1, q):
                if (mask >> pair[w]) & 1:
                    continue
                if isq[(w - v) % q]:
                    nm = mask | (1 << pair[w])
                    dp[nm, w] = (dp[nm, w] + c) % mod
    tot = 0
    for v in range(1, q):
        if not isq[v]:
            tot = (tot + dp[N - 1, v]) % mod
    return tot


def brute_pal(q):
    """enumerate all Hamiltonian paths of P_q by DFS, count H and the (x -> -x)-palindromic ones"""
    sq = {(x * x) % q for x in range(1, q)}
    H = 0
    pal = 0
    path = []
    used = [False] * q

    def dfs(v):
        nonlocal H, pal
        if len(path) == q:
            H += 1
            if all((path[q - 1 - i] + path[i]) % q == 0 for i in range(q)):
                pal += 1
            return
        for w in range(q):
            if not used[w] and (w - v) % q in sq:
                used[w] = True
                path.append(w)
                dfs(w)
                path.pop()
                used[w] = False
    for s in range(q):
        used[s] = True
        path.append(s)
        dfs(s)
        path.pop()
        used[s] = False
    return H, pal


def main():
    out = []
    exact = {7: 9, 11: 1729, 19: 6857869865, 23: 62293308207033}
    for q in (7, 11, 19, 23):
        els, arc = prime_paley(q)
        t0 = time.time()
        r = n1_mod256(els, arc)
        out.append(f'N_1(P_{q}) mod 256 = {r}; exact orbit count mod 256 = {exact[q] % 256}; match {r == exact[q] % 256}  ({time.time() - t0:.1f} s)')
    els, arc = gf27_paley()
    t0 = time.time()
    r27 = n1_mod256(els, arc)
    ref27 = (95 * pow(351, -1, 256)) % 256
    out.append(f'N_1(P_27) mod 256 = {r27}; referee full DP: H mod 256 = 95, H/351 mod 256 = {ref27}; match {r27 == ref27}  ({time.time() - t0:.1f} s)')
    print('\n'.join(out), flush=True)
    out = []
    if len(sys.argv) > 1 and sys.argv[1] == '31':
        els, arc = prime_paley(31)
        t0 = time.time()
        r31 = n1_mod256(els, arc)
        out.append(f'N_1(P_31) = H(P_31)/465 mod 256 = {r31}  (mod 8: {r31 % 8})  ({time.time() - t0:.0f} s)')
    # palindromic counts
    for q in (3, 7, 11, 19, 23, 31):
        isq = np.zeros(q, dtype=np.bool_)
        for x in range(1, q):
            isq[(x * x) % q] = True
        pal = int(pal_dp(q, isq, 2 ** 62))
        P = 2 * pal // (q - 1)
        line = f'q = {q}: pal = {pal} (mod 4: {pal % 4}), P = 2 pal/(q-1) = {P} (exact division {2 * pal % (q - 1) == 0}), P mod 8 = {P % 8}'
        if q in exact or q == 3:
            orb = exact.get(q, 1)
            F2 = orb - P
            line += f'; F = (H/|Aff| - P)/2 = {F2 // 2} (integer {F2 % 2 == 0 and F2 >= 0})'
        out.append(line)
    for q in (7, 11):
        H, palb = brute_pal(q)
        out.append(f'brute force q = {q}: H = {H}, palindromic = {palb}')
    for q, modulus in ((43, 2 ** 13 * 21), (47, 2 ** 13 * 23)):
        isq = np.zeros(q, dtype=np.bool_)
        for x in range(1, q):
            isq[(x * x) % q] = True
        t0 = time.time()
        r = int(pal_dp(q, isq, modulus))
        P8 = (r // ((q - 1) // 2)) % 8 if r % ((q - 1) // 2) == 0 else None
        out.append(f'q = {q}: pal mod {modulus} = {r} (mod 4: {r % 4}); P = pal/((q-1)/2) mod 8 = {P8}  ({time.time() - t0:.0f} s)')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
