#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (D).
# (1) Pi_2(X) = sum over vertex subsets S of H(X[S]) H(X - S) (H(empty) = 1); check Pi_2 = (-1)^n H (mod 3) for
#     EVERY labelled tournament of order n <= 7 (brute force, bitmask DP).
# (2) The Gaussian kernel of THM-5133 (B): find every tournament K of order 4 with ordered ports (x, y) such that
#     H(K[A, B]) = Pi_2(A) Pi_2(B) + H(A) H(B) for all tournaments A, B of orders <= 3 (K[A, B]: x replaced by A,
#     y by B), by brute force on the substituted tournaments; then check the consequence modulo 3:
#     H(K[A, B]) = ((-1)^(|A|+|B|) + 1) H(A) H(B) (mod 3), in particular 3 | H(K[A, B]) when |A| + |B| is odd,
#     on all A, B of orders <= 4 (orders of K[A, B] up to 10).
# (3) The same for the hyperbolic kernel of order 5, H = 3 Pi_2 Pi_2' - H H'  (then H = -H H' mod 3: no constraint
#     from the parities; recorded for comparison).
import itertools, sys
import numpy as np
from numba import njit


@njit(cache=True)
def h_all_subsets(adj, n, Hs):
    # adj[v] = bitmask of out-neighbours. Hs[S] = number of Hamiltonian paths of the induced subtournament on S.
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full):
        for v in range(n):
            c = dp[S, v]
            if c:
                out = adj[v] & ~S
                w = 0
                while out:
                    if out & 1:
                        dp[S | (1 << w), w] += c
                    out >>= 1
                    w += 1
    Hs[0] = 1
    for S in range(1, full):
        t = 0
        for v in range(n):
            t += dp[S, v]
        Hs[S] = t


@njit(cache=True)
def check_all(n):
    m = n * (n - 1) // 2
    full = 1 << n
    Hs = np.zeros(full, dtype=np.int64)
    adj = np.zeros(n, dtype=np.int64)
    bad = 0
    cnt = 0
    hist = np.zeros((3, 3), dtype=np.int64)   # (H mod 3, Pi_2 mod 3)
    for code in range(1 << m):
        for v in range(n):
            adj[v] = 0
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> k) & 1:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                k += 1
        h_all_subsets(adj, n, Hs)
        pi2 = 0
        for S in range(full):
            pi2 += Hs[S] * Hs[(full - 1) ^ S]
        H = Hs[full - 1]
        sign = 1 if n % 2 == 0 else -1
        if (pi2 - sign * H) % 3 != 0:
            bad += 1
        hist[H % 3, pi2 % 3] += 1
        cnt += 1
    return cnt, bad, hist


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (code >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj


def H_and_Pi2(adj):
    n = len(adj)
    masks = np.array([sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)], dtype=np.int64)
    Hs = np.zeros(1 << n, dtype=np.int64)
    h_all_subsets(masks, n, Hs)
    full = (1 << n) - 1
    pi2 = sum(int(Hs[S]) * int(Hs[full ^ S]) for S in range(1 << n))
    return int(Hs[full]), pi2


def substitute(K, x, y, A, B):
    # K[A, B]: vertex x of K replaced by A, y by B; arcs between blocks as in K
    m = len(K)
    blocks = []
    for v in range(m):
        if v == x:
            blocks.append(("A", len(A)))
        elif v == y:
            blocks.append(("B", len(B)))
        else:
            blocks.append(("v", 1))
    verts = []
    for v, (kind, size) in enumerate(blocks):
        for i in range(size):
            verts.append((v, kind, i))
    N = len(verts)
    adj = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            va, ka, ia = verts[a]
            vb, kb, ib = verts[b]
            if va == vb:
                M = A if ka == "A" else B
                adj[a][b] = M[ia][ib]
            else:
                adj[a][b] = K[va][vb]
    return adj


def main():
    print("# pi2_mod3.py (auditor hypotheses; own code)")
    print("(1) Pi_2 = (-1)^n H (mod 3), every labelled tournament:")
    for n in range(1, 8):
        cnt, bad, hist = check_all(n)
        print("    n = %d: %8d tournaments, violations %d;  joint counts of (H mod 3, Pi_2 mod 3): %s" % (
            n, cnt, bad, {(a, b): int(hist[a, b]) for a in range(3) for b in range(3) if hist[a, b]}))
    # (2) the Gaussian kernel at order 4
    small = {k: list(tournaments(k)) for k in (1, 2, 3)}
    inv = {}
    for k in (1, 2, 3, 4):
        inv[k] = [(T, H_and_Pi2(T)) for T in tournaments(k)]
    print("(2) order-4 kernels with H(K[A, B]) = Pi_2(A) Pi_2(B) + H(A) H(B) for all |A|, |B| <= 3:")
    gauss = []
    for K in tournaments(4):
        for x in range(4):
            for y in range(4):
                if x == y:
                    continue
                ok = True
                for ka in (1, 2, 3):
                    for (A, (hA, pA)) in inv[ka]:
                        for kb in (1, 2, 3):
                            for (B, (hB, pB)) in inv[kb]:
                                if ka + kb > 5:
                                    continue
                                T = substitute(K, x, y, A, B)
                                hT, _ = H_and_Pi2(T)
                                if hT != pA * pB + hA * hB:
                                    ok = False
                                    break
                            if not ok:
                                break
                        if not ok:
                            break
                    if not ok:
                        break
                if ok:
                    gauss.append((K, x, y))
    print("    number of (labelled K, ordered ports) with the Gaussian law on |A| + |B| <= 5: %d" % len(gauss))
    if gauss:
        K, x, y = gauss[0]
        print("    example: K arcs %s, ports (x, y) = (%d, %d), H(K) = %d" % (
            [(i, j) for i in range(4) for j in range(4) if K[i][j]], x, y, H_and_Pi2(K)[0]))
        # full check and the consequence modulo 3 on |A|, |B| <= 4 (one representative kernel; random sample of pairs at 4+4)
        import random
        random.seed(1)
        tot = 0
        badlaw = 0
        badmod = 0
        odd_div = 0
        odd_tot = 0
        even_nonzero = 0
        even_tot = 0
        for ka in (1, 2, 3, 4):
            for kb in (1, 2, 3, 4):
                pairs = [(a, b) for a in inv[ka] for b in inv[kb]]
                if len(pairs) > 600:
                    pairs = random.sample(pairs, 600)
                for (A, (hA, pA)), (B, (hB, pB)) in pairs:
                    T = substitute(K, x, y, A, B)
                    hT, _ = H_and_Pi2(T)
                    tot += 1
                    if hT != pA * pB + hA * hB:
                        badlaw += 1
                    pred = (((-1) ** (ka + kb) + 1) * hA * hB) % 3
                    if hT % 3 != pred:
                        badmod += 1
                    if (ka + kb) % 2 == 1:
                        odd_tot += 1
                        odd_div += (hT % 3 == 0)
                    else:
                        even_tot += 1
                        even_nonzero += (hT % 3 != 0)
        print("    checked %d pairs (A, B) with |A|, |B| <= 4: law violated %d times; congruence mod 3 violated %d times" % (tot, badlaw, badmod))
        print("    |A| + |B| odd: %d pairs, 3 | H(K[A, B]) in %d of them;  |A| + |B| even: %d pairs, H != 0 mod 3 in %d" % (
            odd_tot, odd_div, even_tot, even_nonzero))
    # (3) hyperbolic kernel of order 5
    print("(3) order-5 kernels with H(K[A, B]) = 3 Pi_2(A) Pi_2(B) - H(A) H(B) for all |A| + |B| <= 4:")
    hyp = []
    for K in tournaments(5):
        for x in range(5):
            for y in range(5):
                if x == y:
                    continue
                ok = True
                for ka in (1, 2, 3):
                    for (A, (hA, pA)) in inv[ka]:
                        for kb in (1, 2, 3):
                            if ka + kb > 4:
                                continue
                            for (B, (hB, pB)) in inv[kb]:
                                T = substitute(K, x, y, A, B)
                                hT, _ = H_and_Pi2(T)
                                if hT != 3 * pA * pB - hA * hB:
                                    ok = False
                                    break
                            if not ok:
                                break
                        if not ok:
                            break
                    if not ok:
                        break
                if ok:
                    hyp.append((K, x, y))
    print("    number of (labelled K, ordered ports): %d" % len(hyp))


if __name__ == "__main__":
    main()
