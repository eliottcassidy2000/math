#!/usr/bin/env python3
"""gilbreath_fermat_tower_20260926_tournaments.py -- the doubling tower behind the Sierpinski sea, read as
tournaments (session gilbreath-fermat-platonic-20260926, opus, 2026-09-26).

 * H_2 = [[1,1],[-1,1]]; H_(2n) = [[H, H], [-H^T, H^T]]. Claim: H_(2^k) is a normalized skew Hadamard matrix
   (H + H^T = 2I, H H^T = 2^k I, first row all +1, first column -1 off the diagonal). Checked for k <= 8.
 * Deleting row/column 0 gives a tournament T_k on 2^k - 1 vertices (i -> j iff H[i][j] = +1). Claim: T_k is
   doubly regular (every pair has (n-3)/4 common out-neighbours). Checked for k <= 8, i.e. n = 3..255.
 * Recursion: T_(k+1) = T_k + {0'} + T_k' (arc-reversed copy), with i -> i', i -> j' iff i -> j, i' -> j iff i -> j,
   0' -> T_k, T_k' -> 0'. Checked.
 * Identification: automorphism-group orders of T_k for small k by individualization-refinement backtracking
   (node- and time-capped), and isomorphism tests against the Paley tournaments P_7 and P_31.
 * 4-vertex census (counts of the four tournament types on 4 vertices) as a cheap invariant.
Usage: python3 gilbreath_fermat_tower_20260926_tournaments.py
"""
import sys, time, itertools
import numpy as np

sys.setrecursionlimit(10000)


def skew_double(H):
    return np.vstack([np.hstack([H, H]), np.hstack([-H.T, H.T])])


def tower(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1):
        H = skew_double(H)
    return H


def tournament_from_skew(H):
    A = (H[1:, 1:] == 1).copy()
    np.fill_diagonal(A, False)
    return A


def check_drt(A):
    n = A.shape[0]
    J = np.ones((n, n), dtype=int); I = np.eye(n, dtype=int)
    Ai = A.astype(int)
    assert np.array_equal(Ai + Ai.T, J - I), "not a tournament"
    assert np.all(Ai.sum(1) == (n - 1) // 2), "not regular"
    C = Ai @ Ai.T
    off = C[~np.eye(n, dtype=bool)]
    assert np.all(off == (n - 3) // 4), "not doubly regular"
    return True


def paley(q):
    qr = set((x * x) % q for x in range(1, q))
    A = np.zeros((q, q), dtype=bool)
    for i in range(q):
        for j in range(q):
            if i != j and ((j - i) % q) in qr:
                A[i, j] = True
    return A


def masks(A):
    n = A.shape[0]
    out = [0] * n
    for i in range(n):
        m = 0
        for j in range(n):
            if A[i, j]:
                m |= 1 << j
        out[i] = m
    return out


def refine(out, n, colors):
    """equitable refinement: colour by (old colour, counts of out-neighbours in each class)"""
    while True:
        classes = {}
        for v in range(n):
            classes[colors[v]] = classes.get(colors[v], 0) | (1 << v)
        keys = sorted(classes)
        sig = [(colors[v],) + tuple(bin(out[v] & classes[c]).count('1') for c in keys) for v in range(n)]
        uniq = sorted(set(sig))
        index = {s: i for i, s in enumerate(uniq)}
        new = [index[s] for s in sig]
        if len(uniq) == len(classes):
            return new
        colors = new


class Cap(Exception):
    pass


class Found(Exception):
    pass


def iso_search(outA, outB, n, count_all, node_cap, seconds):
    """bijections phi with i -> j in A iff phi(i) -> phi(j) in B; count them (count_all) or return the first"""
    found = [0]; nodes = [0]; t0 = time.time()

    def rec(colA, colB):
        nodes[0] += 1
        if nodes[0] > node_cap or time.time() - t0 > seconds:
            raise Cap()
        colA = refine(outA, n, colA); colB = refine(outB, n, colB)
        if sorted(colA) != sorted(colB):
            return
        if len(set(colA)) == n:
            posB = {colB[v]: v for v in range(n)}
            phi = [posB[colA[v]] for v in range(n)]
            for i in range(n):
                oi = outA[i]; pi = phi[i]
                for j in range(n):
                    if i != j and (((oi >> j) & 1) != ((outB[pi] >> phi[j]) & 1)):
                        return
            found[0] += 1
            if not count_all:
                raise Found(phi)
            return
        cnt = {}
        for v in range(n):
            cnt[colA[v]] = cnt.get(colA[v], 0) + 1
        c = min(k for k in cnt if cnt[k] > 1)
        vA = next(v for v in range(n) if colA[v] == c)
        nc = max(colA) + 1
        for vB in range(n):
            if colB[vB] == c:
                cA = list(colA); cB = list(colB); cA[vA] = nc; cB[vB] = nc
                rec(cA, cB)

    try:
        rec([0] * n, [0] * n)
    except Found as e:
        return ('iso', e.args[0], nodes[0])
    except Cap:
        return ('cap', found[0], nodes[0])
    return ('done', found[0], nodes[0])


def census4(A):
    Ai = A.astype(int)
    n = Ai.shape[0]
    types = {}
    for S in itertools.combinations(range(n), 4):
        sub = Ai[np.ix_(S, S)]
        scores = tuple(sorted(sub.sum(1).tolist()))
        types[scores] = types.get(scores, 0) + 1
    return types


def check_recursion(Ak, Ak1):
    n = Ak.shape[0]; N = Ak1.shape[0]; assert N == 2 * n + 1
    for i in range(n):
        for j in range(n):
            if i != j:
                assert Ak1[i, j] == Ak[i, j]
                assert Ak1[n + 1 + i, n + 1 + j] == Ak[j, i]
                assert Ak1[i, n + 1 + j] == Ak[i, j]
                assert Ak1[n + 1 + i, j] == Ak[i, j]
        assert Ak1[i, n + 1 + i] and not Ak1[n + 1 + i, i]
        assert Ak1[n, i] and not Ak1[i, n]
        assert Ak1[n + 1 + i, n] and not Ak1[n, n + 1 + i]
    return True


def main():
    print("== skew Hadamard doubling tower ==")
    A_prev = None
    Ts = {}
    for k in range(1, 9):
        H = tower(k); n = H.shape[0]
        assert np.array_equal(H + H.T, 2 * np.eye(n, dtype=np.int64)), k
        assert np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64)), k
        assert np.all(H[0, :] == 1) and np.all(H[1:, 0] == -1)
        A = tournament_from_skew(H)
        if k >= 2:
            check_drt(A)
        if A_prev is not None and k >= 2:
            check_recursion(A_prev, A)
        Ts[k] = A; A_prev = A
        print(" k=%d: H_%d skew Hadamard, normalized; T_%d on %d vertices %s" % (k, n, k, n - 1, "doubly regular, recursion T_(k+1) = D(T_k) verified" if k >= 2 else "(single vertex)"))
    print("== identification ==")
    for k in (2, 3, 4, 5, 6, 7, 8):
        A = Ts[k]; n = A.shape[0]; out = masks(A)
        t0 = time.time()
        res = iso_search(out, out, n, True, 300000, 900)
        print(" T_%d (n=%d): |Aut| = %s (%s, %d nodes, %.1fs)" % (k, n, res[1] if res[0] == 'done' else '>=%d' % res[1], res[0], res[2], time.time() - t0))
    # orbit structure of Aut(T_5): enumerate the automorphisms (order 21) and their orbits
    A = Ts[5]; n = A.shape[0]; out = masks(A)
    auts = []

    def collect(colA, colB):
        colA = refine(out, n, colA); colB = refine(out, n, colB)
        if sorted(colA) != sorted(colB):
            return
        if len(set(colA)) == n:
            posB = {colB[v]: v for v in range(n)}
            phi = [posB[colA[v]] for v in range(n)]
            if all(((out[i] >> j) & 1) == ((out[phi[i]] >> phi[j]) & 1) for i in range(n) for j in range(n) if i != j):
                auts.append(tuple(phi))
            return
        cnt = {}
        for v in range(n):
            cnt[colA[v]] = cnt.get(colA[v], 0) + 1
        c = min(k for k in cnt if cnt[k] > 1)
        vA = next(v for v in range(n) if colA[v] == c)
        nc = max(colA) + 1
        for vB in range(n):
            if colB[vB] == c:
                cA = list(colA); cB = list(colB); cA[vA] = nc; cB[vB] = nc
                collect(cA, cB)

    collect([0] * n, [0] * n)
    orbits = []; seen = set()
    for v in range(n):
        if v not in seen:
            orb = sorted(set(phi[v] for phi in auts)); orbits.append(orb); seen |= set(orb)
    orders = {}
    for phi in auts:
        o = 1; x = phi
        while x != tuple(range(n)):
            x = tuple(phi[x[i]] for i in range(n)); o += 1
        orders[o] = orders.get(o, 0) + 1
    print(" Aut(T_5): %d automorphisms, element orders %s, orbit sizes %s (orbits: %s)" % (len(auts), dict(sorted(orders.items())), [len(o) for o in orbits], orbits))
    for k, q in ((3, 7), (5, 31)):
        A = Ts[k]; P = paley(q); check_drt(P)
        t0 = time.time()
        res = iso_search(masks(A), masks(P), q, False, 300000, 240)
        outP = masks(P)
        aut = iso_search(outP, outP, q, True, 300000, 240)
        print(" T_%d vs Paley P_%d: %s (%d nodes, %.1fs); |Aut(P_%d)| = %s" % (k, q, 'ISOMORPHIC' if res[0] == 'iso' else ('NOT isomorphic' if res[0] == 'done' else 'undecided (cap)'), res[2], time.time() - t0, q, aut[1] if aut[0] == 'done' else '>=%d (cap)' % aut[1]))
    print("== 4-vertex census: score multisets (0,1,2,3) transitive, (1,1,1,3) 3-cycle+source, (0,2,2,2) 3-cycle+sink, (1,1,2,2) strong ==")
    for k in (3, 4, 5):
        print(" T_%d:" % k, census4(Ts[k]))
    for q in (7, 31):
        print(" P_%d:" % q, census4(paley(q)))


if __name__ == '__main__':
    main()
