#!/usr/bin/env python3
"""Auditor proof-5135, script a6: THM-5135 (6) and Remark 2.

  [B] Brahmagupta identity on one-apex states (x, y) = (2 psi(h) + 1, 2h + 1), psi from the recursion of the
      task statement; the F1 value N = M_s - [s = 0]; real and imaginary parts of z_A conj(z_B),
      z = x + sqrt(-2^s) y (exact arithmetic in Z[w], w^2 = -2^s).
  [V] Remark 2 says 'rank 3 occurs in the join calculus (kernels with three visits, THM-5133 (A))'.
      THM-5133 (A): H(K[A, B]) = sum N_ij P_i(A) P_j(B), N_ij = number of walks of K visiting every vertex other
      than the ports exactly once, the port x exactly i times and y exactly j times.
      Here: the visit matrices of all labelled tournaments of order 4 and 5 (and 6) with all ordered port pairs,
      validated against brute-force path counts of K[A, B]; among the kernels with at most three visits, how
      many have a 3 x 3 visit matrix of non-zero determinant (a seam of rank 3 at the primes not dividing it).
"""
import sys
import itertools
from math import gcd

import numpy as np


# ------------------------------------------------------------------ one-apex states
def psi(h, memo={0: 0}):
    if h in memo:
        return memo[h]
    if h % 2 == 0:
        r = psi(h // 2) + h // 2
    else:
        k = h // 2
        r = psi(k) + k // 2 + 1
    memo[h] = r
    return r


def ham_paths(adj, n):
    """number of Hamiltonian paths; adj[i] = bitmask of out-neighbours."""
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c == 0:
                continue
            nxt = adj[last] & ~mask
            while nxt:
                b = nxt & -nxt
                j = b.bit_length() - 1
                dp[mask | b][j] += c
                nxt ^= b
    return sum(dp[(1 << n) - 1])


def sub_adj(adj, verts):
    idx = {v: i for i, v in enumerate(verts)}
    out = []
    for v in verts:
        m = 0
        for w in verts:
            if (adj[v] >> w) & 1:
                m |= 1 << idx[w]
        out.append(m)
    return out


def path_cover_numbers(adj, n, imax):
    """P_i = number of ordered i-tuples of disjoint non-empty paths covering the tournament, i = 1..imax."""
    full = (1 << n) - 1
    Hsub = {}
    for mask in range(1, 1 << n):
        verts = [v for v in range(n) if (mask >> v) & 1]
        Hsub[mask] = ham_paths(sub_adj(adj, verts), len(verts))
    P = [0] * (imax + 1)
    # ordered set partitions into i non-empty blocks
    cur = {0: 1}
    for i in range(1, imax + 1):
        nxt = {}
        for mask, c in cur.items():
            rest = full & ~mask
            sub = rest
            while sub:
                nxt[mask | sub] = nxt.get(mask | sub, 0) + c * Hsub[sub]
                sub = (sub - 1) & rest
        cur = nxt
        P[i] = cur.get(full, 0)
    return P


def substitute(adjK, m, x, y, adjA, a, adjB, b):
    """K[A, B]: vertex x replaced by A, y by B."""
    others = [v for v in range(m) if v not in (x, y)]
    # new labels: others first, then A, then B
    lab = {}
    for i, v in enumerate(others):
        lab[("o", v)] = i
    for i in range(a):
        lab[("a", i)] = len(others) + i
    for i in range(b):
        lab[("b", i)] = len(others) + a + i
    n = len(others) + a + b
    adj = [0] * n

    def groups(v):
        if v == x:
            return [("a", i) for i in range(a)]
        if v == y:
            return [("b", i) for i in range(b)]
        return [("o", v)]
    for v in range(m):
        for w in range(m):
            if v != w and (adjK[v] >> w) & 1:
                for gv in groups(v):
                    for gw in groups(w):
                        adj[lab[gv]] |= 1 << lab[gw]
    for i in range(a):
        for j in range(a):
            if (adjA[i] >> j) & 1:
                adj[lab[("a", i)]] |= 1 << lab[("a", j)]
    for i in range(b):
        for j in range(b):
            if (adjB[i] >> j) & 1:
                adj[lab[("b", i)]] |= 1 << lab[("b", j)]
    return adj, n


def visit_matrix(adjK, m, x, y):
    """N[i][j], 0 <= i, j <= m - 1."""
    N = [[0] * m for _ in range(m)]
    others_mask = ((1 << m) - 1) & ~(1 << x) & ~(1 << y)

    def rec(v, mask, cx, cy):
        if mask == others_mask:
            N[cx][cy] += 1
        nxt = adjK[v]
        w = 0
        while nxt:
            if nxt & 1:
                if w == x:
                    if cx < m - 1:
                        rec(w, mask, cx + 1, cy)
                elif w == y:
                    if cy < m - 1:
                        rec(w, mask, cx, cy + 1)
                elif not (mask >> w) & 1:
                    rec(w, mask | (1 << w), cx, cy)
            nxt >>= 1
            w += 1

    for v in range(m):
        if v == x:
            rec(v, 0, 1, 0)
        elif v == y:
            rec(v, 0, 0, 1)
        else:
            rec(v, 1 << v, 0, 0)
    return N


def tournaments(m):
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    for bits in range(1 << len(pairs)):
        adj = [0] * m
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        yield bits, adj


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def main():
    out = []
    P = out.append

    # ------------------------------------------------------------------ [B]
    P("[B] Brahmagupta identity on one-apex states")
    # psi against the closed description q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2)
    okq = True
    for h in range(0, 3000):
        q2 = h - psi(h)
        if (2 * h - psi(2 * h)) != q2 or ((2 * h + 1) - psi(2 * h + 1)) != q2 + (h + 1) // 2:
            okq = False
    P("  recursion check q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), h < 3000: %s; first states (x,y): %s"
      % (okq, [(2 * psi(h) + 1, 2 * h + 1) for h in range(8)]))
    bad = 0
    bad_re = 0
    n = 0
    for hA in range(0, 400):
        xA, yA = 2 * psi(hA) + 1, 2 * hA + 1
        for hB in range(0, 400, 3):
            xB, yB = 2 * psi(hB) + 1, 2 * hB + 1
            for s in range(0, 7):
                t = 2 ** s
                MA = xA * xA + t * yA * yA
                MB = xB * xB + t * yB * yB
                MAB = xA * xB + t * yA * yB
                n += 1
                if MA * MB != MAB ** 2 + t * (xA * yB - xB * yA) ** 2:
                    bad += 1
                # z_A conj(z_B) in Z[w], w^2 = -t: (xA + w yA)(xB - w yB) = (xA xB + t yA yB) + w (xB yA - xA yB)
                re = xA * xB - (yA * (-yB)) * (-t) * (-1)   # written out below instead
                re = xA * xB + t * yA * yB
                im = xB * yA - xA * yB
                if re != MAB or im != -(xA * yB - xB * yA):
                    bad_re += 1
    P("  M_s(A,A) M_s(B,B) = M_s(A,B)^2 + 2^s (x_A y_B - x_B y_A)^2: %d cases, failures %d" % (n, bad))
    P("  z_A conj(z_B) = M_s(A,B) + w (x_B y_A - x_A y_B), w = sqrt(-2^s): the coefficient of w is MINUS the cross-determinant")
    P("  x_A y_B - x_B y_A, and the imaginary part proper is 2^(s/2) times it (failures of this bookkeeping: %d)" % bad_re)

    # ------------------------------------------------------------------ [V]
    P("")
    P("[V] visit matrices (THM-5133 (A)) and rank 3")
    # validation of the formula on random kernels of order 4, 5 with small modules
    rng = np.random.default_rng(5133)
    mods = []
    for a in (1, 2, 3, 4):
        seen = set()
        for bits, adj in tournaments(a):
            mods.append((a, adj))
    val_ok = True
    nval = 0
    for m in (4, 5):
        allK = list(tournaments(m))
        for _ in range(12):
            bits, adjK = allK[int(rng.integers(len(allK)))]
            x, y = [int(v) for v in rng.choice(m, size=2, replace=False)]
            N = visit_matrix(adjK, m, x, y)
            for _ in range(6):
                a, adjA = mods[int(rng.integers(len(mods)))]
                b, adjB = mods[int(rng.integers(len(mods)))]
                PA = path_cover_numbers(adjA, a, m - 1)
                PB = path_cover_numbers(adjB, b, m - 1)
                pred = sum(N[i][j] * PA[i] * PB[j] for i in range(1, m) for j in range(1, m))
                adj, n = substitute(adjK, m, x, y, adjA, a, adjB, b)
                true = ham_paths(adj, n)
                nval += 1
                if pred != true:
                    val_ok = False
    P("  H(K[A,B]) = sum N_ij P_i(A) P_j(B): %d random (K, ports, A, B) with |K| in {4,5}, |A|,|B| <= 4: all equal: %s" % (nval, val_ok))

    # the two kernels named in THM-5133 (B) / THM-5135 (4)
    found52 = None
    found116 = None
    for m in (4, 5):
        for bits, adjK in tournaments(m):
            for x in range(m):
                for y in range(m):
                    if x == y:
                        continue
                    N = visit_matrix(adjK, m, x, y)
                    mx = max([i for i in range(m) for j in range(m) if N[i][j]] + [0])
                    my = max([j for i in range(m) for j in range(m) if N[i][j]] + [0])
                    if mx <= 2 and my <= 2:
                        G = [[N[1][1], N[1][2]], [N[2][1], N[2][2]]]
                        if m == 4 and G == [[5, 2], [2, 1]] and found52 is None:
                            found52 = (bits, x, y)
                        if m == 5 and G == [[11, 6], [6, 3]] and found116 is None:
                            found116 = (bits, x, y)
    P("  a kernel of order 4 with visit matrix [[5,2],[2,1]] exists: %s; of order 5 with [[11,6],[6,3]]: %s"
      % (found52 is not None, found116 is not None))

    # census of three-visit kernels
    for m in (4, 5):
        n_all = n_le2 = n_3 = n_3x3 = n_det = 0
        dets = {}
        example = None
        for bits, adjK in tournaments(m):
            for x in range(m):
                for y in range(m):
                    if x == y:
                        continue
                    n_all += 1
                    N = visit_matrix(adjK, m, x, y)
                    mx = max([i for i in range(m) for j in range(m) if N[i][j]] + [0])
                    my = max([j for i in range(m) for j in range(m) if N[i][j]] + [0])
                    if mx <= 2 and my <= 2:
                        n_le2 += 1
                    elif mx <= 3 and my <= 3:
                        n_3 += 1
                        G = [[N[i][j] for j in (1, 2, 3)] for i in (1, 2, 3)]
                        if mx == 3 and my == 3:
                            n_3x3 += 1
                        d = det3(G)
                        if d != 0:
                            n_det += 1
                            dets[abs(d)] = dets.get(abs(d), 0) + 1
                            if example is None or abs(d) < abs(example[3]):
                                example = (bits, x, y, d, G)
        P("  order %d (labelled tournaments x ordered port pairs): %d kernels; at most two visits %d; at most three visits (and some port three times) %d,"
          % (m, n_all, n_le2, n_3))
        P("     of which both ports three times %d; 3 x 3 visit matrix of non-zero determinant: %d; |det| values (count): %s"
          % (n_3x3, n_det, sorted(dets.items())[:12]))
        if example:
            P("     example (least |det|): tournament bits %d, ports (%d, %d), visit matrix %s, det %d"
              % (example[0], example[1], example[2], example[4], example[3]))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
