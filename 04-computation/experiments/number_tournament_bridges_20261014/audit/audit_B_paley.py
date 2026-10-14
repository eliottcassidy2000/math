#!/usr/bin/env python3
"""Independent audit of HYP-9273 (B) and the Paley / circulant / doubly-regular controls (referee script).

Method differs from the author's: Hamiltonian paths are counted by a full subset DP over ALL start vertices
(dp[mask][v] = number of directed paths with vertex set mask ending at v), without using vertex-transitivity.
Automorphism groups are counted by an explicit backtracking search (no Burnside/Alspach shortcut).
P_27 is built on GF(27) = F_3[y]/(y^3 + 2y + 1) (a different modulus from the author's x^3 - x - 1) and H(P_27) is
computed modulo 256 in a uint8 table (|Aut| is odd, so H mod 8 already determines the orbit count mod 8).
"""
import sys, time
import numpy as np
import numba as nb


LOG37 = np.zeros(37, dtype=np.int64)
for _k in range(36):
    LOG37[(1 << _k) % 37] = _k


@nb.njit(cache=True)
def ham_all(outmask, n):
    """total number of Hamiltonian paths (all starts), int64"""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, N):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            nxt = outmask[v] & ~mask
            while nxt:
                b = nxt & (-nxt)
                w = LOG37[b % 37]
                dp[mask | b, w] += c
                nxt ^= b
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


@nb.njit(cache=True)
def ham_all_u8(outmask, n):
    """total number of Hamiltonian paths modulo 256 (uint8 wrap-around), all starts"""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.uint8)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, N):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            nxt = outmask[v] & ~mask
            while nxt:
                b = nxt & (-nxt)
                w = LOG37[b % 37]
                dp[mask | b, w] += c
                nxt ^= b
    tot = 0
    for v in range(n):
        tot += np.int64(dp[N - 1, v])
    return tot % 256


@nb.njit(cache=True)
def count_aut(adj, n):
    perm = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=np.uint8)
    outdeg = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if adj[i, j]:
                outdeg[i] += 1
    count = 0
    cand = np.zeros(n + 1, dtype=np.int64)
    depth = 0
    while depth >= 0:
        if depth == n:
            count += 1
            depth -= 1
            used[perm[depth]] = 0
            perm[depth] = -1
            continue
        found = False
        c = cand[depth]
        while c < n:
            if used[c] == 0 and outdeg[c] == outdeg[depth]:
                ok = True
                for k in range(depth):
                    if adj[k, depth] != adj[perm[k], c]:
                        ok = False
                        break
                if ok:
                    found = True
                    break
            c += 1
        if found:
            perm[depth] = c
            used[c] = 1
            cand[depth] = c + 1
            depth += 1
            cand[depth] = 0
        else:
            cand[depth] = 0
            depth -= 1
            if depth >= 0:
                used[perm[depth]] = 0
                perm[depth] = -1
    return count


def masks_from_adj(adj):
    n = adj.shape[0]
    return np.array([sum(1 << j for j in range(n) if adj[i, j]) for i in range(n)], dtype=np.int64)


def paley_prime(p):
    sq = {(x * x) % p for x in range(1, p)}
    adj = np.zeros((p, p), dtype=np.bool_)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in sq:
                adj[i, j] = True
    return adj


def is_tournament(adj):
    n = adj.shape[0]
    return all((adj[i, j] != adj[j, i]) for i in range(n) for j in range(n) if i != j) and not adj.diagonal().any()


# ---- GF(27) as F_3[y]/(y^3 + 2y + 1): elements are integers c0 + 3 c1 + 9 c2
def gf27_mul(a, b):
    A = [(a // 3 ** i) % 3 for i in range(3)]
    B = [(b // 3 ** i) % 3 for i in range(3)]
    pr = [0] * 5
    for i in range(3):
        for j in range(3):
            pr[i + j] += A[i] * B[j]
    # y^3 = -2y - 1 = y + 2  (mod 3)
    for k in (4, 3):
        c = pr[k]
        pr[k] = 0
        pr[k - 3] += 2 * c      # y^k = y^(k-3) * (y + 2)
        pr[k - 2] += c
    return sum((pr[i] % 3) * 3 ** i for i in range(3))


def gf27_sub(a, b):
    return sum(((((a // 3 ** i) % 3) - ((b // 3 ** i) % 3)) % 3) * 3 ** i for i in range(3))


def main():
    out = []
    t00 = time.time()
    # ---------- Paley primes
    out.append('== B1. H(P_q) by full subset DP over all start vertices; |Aut| by backtracking')
    res = {}
    for p in (3, 7, 11, 19, 23):
        adj = paley_prime(p)
        assert is_tournament(adj)
        t0 = time.time()
        H = int(ham_all(masks_from_adj(adj), p))
        aut = int(count_aut(adj, p)) if p <= 23 else None
        aff = p * (p - 1) // 2
        res[p] = H
        out.append(f'   P_{p}: H = {H}; |Aut| (backtracking) = {aut} (affine order {aff}); H/aff = {H // aff}, '
                   f'divisible {H % aff == 0}, orbit count mod 8 = {(H // aff) % 8};  {time.time() - t0:.1f} s')
    claimed = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}
    out.append(f'   matches the note: {all(res[p] == claimed[p] for p in claimed)}')

    # ---------- P_27 mod 256
    out.append('\n== B2. P_27 on F_3[y]/(y^3 + 2y + 1), H mod 256 via uint8 full DP')
    els = list(range(27))
    # irreducibility check of y^3 + 2y + 1: no root in F_3
    assert all((t ** 3 + 2 * t + 1) % 3 != 0 for t in range(3))
    sq = {gf27_mul(x, x) for x in range(1, 27)}
    assert len(sq) == 13 and gf27_sub(0, 1) not in sq       # -1 is a non-square
    adj = np.zeros((27, 27), dtype=np.bool_)
    for i in els:
        for j in els:
            if i != j and gf27_sub(j, i) in sq:
                adj[i, j] = True
    assert is_tournament(adj)
    t0 = time.time()
    aut27 = int(count_aut(adj, 27))
    out.append(f'   |Aut(P_27)| by backtracking = {aut27}  ({time.time() - t0:.1f} s)')
    t0 = time.time()
    h256 = int(ham_all_u8(masks_from_adj(adj), 27))
    inv1053 = pow(1053, -1, 256)
    inv351 = pow(351, -1, 256)
    out.append(f'   H(P_27) mod 256 = {h256};  H/1053 mod 8 = {(h256 * inv1053) % 8};  H/351 mod 8 = {(h256 * inv351) % 8};'
               f'  ({time.time() - t0:.0f} s)')
    out.append(f'   author: H mod 8424 = 3159 -> H mod 8 = {3159 % 8}; consistent with ours: {3159 % 8 == h256 % 8}')

    # ---------- all circulant tournaments on Z_p
    out.append('\n== B3. all circulant tournaments on Z_p: residues of H/|Aut| mod 8, |Aut| by backtracking (no Alspach)')
    for p in (3, 5, 7, 11, 13, 17, 19):
        pairs = [(x, p - x) for x in range(1, (p - 1) // 2 + 1)]
        resid = {}
        auts = {}
        t0 = time.time()
        for bits in range(1 << len(pairs)):
            S = {pairs[i][(bits >> i) & 1] for i in range(len(pairs))}
            adj = np.zeros((p, p), dtype=np.bool_)
            for i in range(p):
                for j in range(p):
                    if i != j and (j - i) % p in S:
                        adj[i, j] = True
            H = int(ham_all(masks_from_adj(adj), p))
            a = int(count_aut(adj, p))
            # cross-check with the multiplier stabilizer
            M = sum(1 for m in range(1, p) if {(m * x) % p for x in S} == S)
            assert a == p * M, (p, S, a, M)
            assert H % a == 0
            r = (H // a) % 8
            resid[r] = resid.get(r, 0) + 1
            auts[a] = auts.get(a, 0) + 1
        out.append(f'   p = {p:2d}: residues {dict(sorted(resid.items()))};  |Aut| distribution {dict(sorted(auts.items()))};'
                   f'  backtracking |Aut| = p |M| for all: True;  {time.time() - t0:.0f} s')

    # ---------- doubly regular tournament of order 15 (normalized skew-Hadamard core of the doubling tower)
    out.append('\n== B4. doubly regular tournament of order 15')
    Hm = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(3):
        Hm = np.block([[Hm, Hm], [-Hm.T, Hm.T]])
    assert (Hm @ Hm.T == 16 * np.eye(16, dtype=np.int64)).all()
    assert ((Hm - np.eye(16, dtype=np.int64)) == -(Hm - np.eye(16, dtype=np.int64)).T).all()
    A = (Hm[1:, 1:] == 1)
    np.fill_diagonal(A, False)
    assert is_tournament(A)
    n = 15
    scores = A.sum(axis=1)
    common = {int((A[i] & A[j]).sum()) for i in range(n) for j in range(n) if i != j}
    t0 = time.time()
    H15 = int(ham_all(masks_from_adj(A), 15))
    a15 = int(count_aut(A, 15))
    out.append(f'   scores {sorted(set(scores.tolist()))}, common out-neighbours of distinct pairs {common} (doubly regular iff {{3}})')
    out.append(f'   H = {H15}; |Aut| (backtracking) = {a15}; H/|Aut| = {H15 // a15 if H15 % a15 == 0 else "not integral"},'
               f' mod 8 = {(H15 // a15) % 8};  {time.time() - t0:.1f} s')
    # vertex-transitivity check: orbits of Aut on vertices (via count of automorphisms mapping 0 to v)
    out.append(f'\n(total {time.time() - t00:.0f} s)')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
