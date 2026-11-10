#!/usr/bin/env python3
"""Audit of the claims THM-5124 makes ABOUT THM-5122 (exchange cubes):
   status field : "THM-5122 constructs additive cubes of arbitrary dimension up to the negative inertia"
   Remark 2     : "additive cubes exist up to dimension n_-(W)"
   README sec 1 : "in the exchange cube of one core additive cubes exist"
THM-5122 itself only proves the bound k <= min(n_-(W), floor(n/2)) and exhibits ONE two-dimensional cube
(core [10,28,49,36,41,3], counts 345, 347, 347, 349).  This script asks what actually exists at small order.

Setting (Session246 switches/PROOF.md section 2).  Core C on n vertices, cut x in {0,1}^n.  Protected closure:
append a sink v (C -> v), then a root r with r -> c iff x_c = 1, c -> r otherwise, and v -> r.  F(x) = H(closure).
A pair exchange (p, q) moves the selection from p to q.  k pairwise disjoint pairs span a cube of 2^k cuts
(the remaining vertices keep a fixed background).  The cube is affine (additive) iff d_i^T W d_j = 0 for i != j,
W_ij = -(F(e_i + e_j) - F(e_i) - F(e_j) + F(0)).

Checks, every labelled core of order n = 4, 5, 6 (own depth-first path counting on the closure, all 2^n cuts):
  (Q) F is quadratic on the whole cube (every identity F(x) = F(0) + sum L_i x_i - sum_{i<j} W_ij x_i x_j);
      W_ij >= 2 H(C - {i, j}) > 0  (THM-5122 (5)).
  (E2) number of (core, two disjoint pairs, background) whose 4 cuts are affine with both increments nonzero.
  (E3) n = 6: number of (core, perfect matching) whose 8 cuts are affine; how many have all three increments
       nonzero; how many of those have all closures strong; examples.
  (W)  the THM-5122 witness: W matrix, 345, 347, 347, 349.
Then (E3) on random cores of orders 7 and 8 (three disjoint pairs, every background of the other vertices).
"""
import itertools
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def count_hp(n, adj, lowidx, st_mask, st_cand):
    """number of Hamiltonian paths, by depth-first enumeration."""
    full = (1 << n) - 1
    total = 0
    for s in range(n):
        m0 = 1 << s
        if m0 == full:
            total += 1
        depth = 0
        st_mask[0] = m0
        st_cand[0] = adj[s] & ~m0
        while depth >= 0:
            c = st_cand[depth]
            if c == 0:
                depth -= 1
            else:
                low = c & (-c)
                st_cand[depth] = c ^ low
                w = lowidx[low]
                nm = st_mask[depth] | low
                if nm == full:
                    total += 1
                else:
                    depth += 1
                    st_mask[depth] = nm
                    st_cand[depth] = adj[w] & ~nm
    return total


@njit(cache=True)
def closure_count(n, adjC, x, adj, lowidx, st_mask, st_cand):
    """C on 0..n-1, sink v = n, root r = n+1."""
    fullC = (1 << n) - 1
    for i in range(n):
        adj[i] = adjC[i] | (1 << n)
        if not (x >> i) & 1:
            adj[i] |= 1 << (n + 1)
    adj[n] = 1 << (n + 1)
    adj[n + 1] = x & fullC
    return count_hp(n + 2, adj, lowidx, st_mask, st_cand)


@njit(cache=True)
def is_strong(n, adj, lowidx):
    full = (1 << n) - 1
    seen = 1
    frontier = 1
    while frontier:
        low = frontier & (-frontier)
        frontier ^= low
        w = lowidx[low]
        new = adj[w] & ~seen
        seen |= new
        frontier |= new
    if seen != full:
        return False
    # co-reachability
    seen = 1
    changed = True
    while changed:
        changed = False
        for u in range(n):
            if not (seen >> u) & 1 and (adj[u] & seen) != 0:
                seen |= 1 << u
                changed = True
    return seen == full


@njit(cache=True)
def set_adj(n, code, adj):
    for i in range(n):
        adj[i] = 0
    k = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> k) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            k += 1


@njit(cache=True)
def all_F(n, adjC, F, adj, lowidx, st_mask, st_cand):
    for x in range(1 << n):
        F[x] = closure_count(n, adjC, x, adj, lowidx, st_mask, st_cand)


@njit(cache=True)
def exhaustive(n, lowidx, pairs2, match3):
    npairs = n * (n - 1) // 2
    adjC = np.zeros(n, np.int64)
    adj = np.zeros(n + 2, np.int64)
    F = np.zeros(1 << n, np.int64)
    W = np.zeros((n, n), np.int64)
    st_mask = np.zeros(n + 3, np.int64)
    st_cand = np.zeros(n + 3, np.int64)
    res = np.zeros(12, np.int64)
    ex = np.zeros((4, 12), np.int64)
    nex = 0
    sub = np.zeros(max(n - 2, 1), np.int64)
    for code in range(1 << npairs):
        set_adj(n, code, adjC)
        all_F(n, adjC, F, adj, lowidx, st_mask, st_cand)
        # quadratic fit
        for i in range(n):
            for j in range(n):
                W[i, j] = 0
        for i in range(n):
            for j in range(i + 1, n):
                w = -(F[(1 << i) | (1 << j)] - F[1 << i] - F[1 << j] + F[0])
                W[i, j] = w
                W[j, i] = w
        for x in range(1 << n):
            val = F[0]
            for i in range(n):
                if (x >> i) & 1:
                    val += F[1 << i] - F[0]
                    for j in range(i + 1, n):
                        if (x >> j) & 1:
                            val -= W[i, j]
            if val != F[x]:
                res[0] += 1
        # W_ij >= 2 H(C - {i, j})
        for i in range(n):
            for j in range(i + 1, n):
                m = 0
                idx = 0
                # induced subtournament on the other n - 2 vertices
                for u in range(n):
                    if u != i and u != j:
                        row = 0
                        idx2 = 0
                        for t in range(n):
                            if t != i and t != j:
                                if (adjC[u] >> t) & 1:
                                    row |= 1 << idx2
                                idx2 += 1
                        sub[idx] = row
                        idx += 1
                h2 = 1
                if n - 2 >= 1:
                    h2 = count_hp(n - 2, sub, lowidx, st_mask, st_cand)
                if W[i, j] < 2 * h2:
                    res[1] += 1
        # two disjoint pairs
        for t in range(pairs2.shape[0]):
            p1 = pairs2[t, 0]
            q1 = pairs2[t, 1]
            p2 = pairs2[t, 2]
            q2 = pairs2[t, 3]
            if W[q1, q2] - W[q1, p2] - W[p1, q2] + W[p1, p2] != 0:
                continue
            res[2] += 1        # orthogonal pair of exchanges (background-free condition)
            used = (1 << p1) | (1 << q1) | (1 << p2) | (1 << q2)
            rest = ((1 << n) - 1) ^ used
            bg = rest
            while True:
                x0 = bg | (1 << p1) | (1 << p2)
                g1 = F[x0 ^ (1 << p1) ^ (1 << q1)] - F[x0]
                g2 = F[x0 ^ (1 << p2) ^ (1 << q2)] - F[x0]
                if g1 != 0 and g2 != 0:
                    res[3] += 1
                if bg == 0:
                    break
                bg = (bg - 1) & rest
        # three disjoint pairs (n = 6: perfect matchings)
        for t in range(match3.shape[0]):
            p1 = match3[t, 0]
            q1 = match3[t, 1]
            p2 = match3[t, 2]
            q2 = match3[t, 3]
            p3 = match3[t, 4]
            q3 = match3[t, 5]
            if W[q1, q2] - W[q1, p2] - W[p1, q2] + W[p1, p2] != 0:
                continue
            if W[q1, q3] - W[q1, p3] - W[p1, q3] + W[p1, p3] != 0:
                continue
            if W[q2, q3] - W[q2, p3] - W[p2, q3] + W[p2, p3] != 0:
                continue
            used = (1 << p1) | (1 << q1) | (1 << p2) | (1 << q2) | (1 << p3) | (1 << q3)
            rest = ((1 << n) - 1) ^ used
            bg = rest
            while True:
                res[4] += 1      # affine 3-cubes (core, matching, background)
                x0 = bg | (1 << p1) | (1 << p2) | (1 << p3)
                g1 = F[x0 ^ (1 << p1) ^ (1 << q1)] - F[x0]
                g2 = F[x0 ^ (1 << p2) ^ (1 << q2)] - F[x0]
                g3 = F[x0 ^ (1 << p3) ^ (1 << q3)] - F[x0]
                if g1 != 0 and g2 != 0 and g3 != 0:
                    res[5] += 1
                    # all eight closures strong?
                    allstrong = True
                    for b in range(8):
                        x = x0
                        if b & 1:
                            x ^= (1 << p1) ^ (1 << q1)
                        if b & 2:
                            x ^= (1 << p2) ^ (1 << q2)
                        if b & 4:
                            x ^= (1 << p3) ^ (1 << q3)
                        closure_count(n, adjC, x, adj, lowidx, st_mask, st_cand)
                        if not is_strong(n + 2, adj, lowidx):
                            allstrong = False
                    if allstrong:
                        res[6] += 1
                    a1 = abs(g1) // 2
                    a2 = abs(g2) // 2
                    a3 = abs(g3) // 2
                    lo = min(a1, min(a2, a3))
                    hi = max(a1, max(a2, a3))
                    mid = a1 + a2 + a3 - lo - hi
                    complete = lo == 1 and mid <= 1 + lo and hi <= 1 + lo + mid
                    if complete:
                        res[7] += 1
                        if allstrong:
                            res[8] += 1
                    if nex < 4 and (allstrong or nex < 2):
                        ex[nex, 0] = code
                        ex[nex, 1] = x0
                        ex[nex, 2] = p1
                        ex[nex, 3] = q1
                        ex[nex, 4] = p2
                        ex[nex, 5] = q2
                        ex[nex, 6] = p3
                        ex[nex, 7] = q3
                        ex[nex, 8] = g1
                        ex[nex, 9] = g2
                        ex[nex, 10] = g3
                        ex[nex, 11] = 1 if allstrong else 0
                        nex += 1
                if bg == 0:
                    break
                bg = (bg - 1) & rest
    return res, ex, nex


@njit(cache=True)
def sampled3(n, trials, seed, lowidx, match3):
    np.random.seed(seed)
    adjC = np.zeros(n, np.int64)
    adj = np.zeros(n + 2, np.int64)
    W = np.zeros((n, n), np.int64)
    L = np.zeros(n, np.int64)
    st_mask = np.zeros(n + 3, np.int64)
    st_cand = np.zeros(n + 3, np.int64)
    res = np.zeros(6, np.int64)
    ex = np.zeros((4, 13 + n), np.int64)
    nex = 0
    for tr in range(trials):
        for i in range(n):
            adjC[i] = 0
        for j in range(1, n):
            for i in range(j):
                if np.random.random() < 0.5:
                    adjC[i] |= 1 << j
                else:
                    adjC[j] |= 1 << i
        F0 = closure_count(n, adjC, 0, adj, lowidx, st_mask, st_cand)
        for i in range(n):
            L[i] = closure_count(n, adjC, 1 << i, adj, lowidx, st_mask, st_cand) - F0
        for i in range(n):
            for j in range(i + 1, n):
                w = -(closure_count(n, adjC, (1 << i) | (1 << j), adj, lowidx, st_mask, st_cand) - L[i] - L[j] - F0)
                W[i, j] = w
                W[j, i] = w
        for t in range(match3.shape[0]):
            p1 = match3[t, 0]
            q1 = match3[t, 1]
            p2 = match3[t, 2]
            q2 = match3[t, 3]
            p3 = match3[t, 4]
            q3 = match3[t, 5]
            if W[q1, q2] - W[q1, p2] - W[p1, q2] + W[p1, p2] != 0:
                continue
            if W[q1, q3] - W[q1, p3] - W[p1, q3] + W[p1, p3] != 0:
                continue
            if W[q2, q3] - W[q2, p3] - W[p2, q3] + W[p2, p3] != 0:
                continue
            used = (1 << p1) | (1 << q1) | (1 << p2) | (1 << q2) | (1 << p3) | (1 << q3)
            rest = ((1 << n) - 1) ^ used
            bg = rest
            while True:
                res[0] += 1
                x0 = bg | (1 << p1) | (1 << p2) | (1 << p3)
                # increments from the quadratic form
                g = np.zeros(3, np.int64)
                ps = (p1, p2, p3)
                qs = (q1, q2, q3)
                for e in range(3):
                    p = ps[e]
                    q = qs[e]
                    # F(x0 - e_p + e_q) - F(x0) = L_q - L_p - sum_{s in x0, s != p} (W_qs - W_ps)   [W_pq term: +W_pq]
                    val = L[q] - L[p]
                    for s in range(n):
                        if (x0 >> s) & 1 and s != p:
                            val -= W[q, s] - W[p, s]
                    g[e] = val
                if g[0] != 0 and g[1] != 0 and g[2] != 0:
                    res[1] += 1
                    # verify the eight values directly and test strongness
                    allstrong = True
                    affine = True
                    base = closure_count(n, adjC, x0, adj, lowidx, st_mask, st_cand)
                    for b in range(8):
                        x = x0
                        pred = base
                        for e in range(3):
                            if (b >> e) & 1:
                                x ^= (1 << ps[e]) ^ (1 << qs[e])
                                pred += g[e]
                        val = closure_count(n, adjC, x, adj, lowidx, st_mask, st_cand)
                        if val != pred:
                            affine = False
                        if not is_strong(n + 2, adj, lowidx):
                            allstrong = False
                    if not affine:
                        res[2] += 1
                    if allstrong:
                        res[3] += 1
                    if nex < 4 and allstrong:
                        for i in range(n):
                            ex[nex, i] = adjC[i]
                        ex[nex, n] = x0
                        ex[nex, n + 1] = p1
                        ex[nex, n + 2] = q1
                        ex[nex, n + 3] = p2
                        ex[nex, n + 4] = q2
                        ex[nex, n + 5] = p3
                        ex[nex, n + 6] = q3
                        ex[nex, n + 7] = g[0]
                        ex[nex, n + 8] = g[1]
                        ex[nex, n + 9] = g[2]
                        ex[nex, n + 10] = base
                        nex += 1
                if bg == 0:
                    break
                bg = (bg - 1) & rest
    return res, ex, nex


def disjoint_pair_pairs(n):
    out = []
    for p1, q1 in itertools.combinations(range(n), 2):
        for p2, q2 in itertools.combinations(range(n), 2):
            if (p1, q1) < (p2, q2) and len({p1, q1, p2, q2}) == 4:
                out.append((p1, q1, p2, q2))
    return np.array(out, dtype=np.int64).reshape(-1, 4)


def disjoint_pair_triples(n):
    out = []
    prs = list(itertools.combinations(range(n), 2))
    for a, b, c in itertools.combinations(prs, 3):
        if len(set(a) | set(b) | set(c)) == 6:
            out.append(a + b + c)
    return np.array(out, dtype=np.int64).reshape(-1, 6)


def main():
    t0 = time.time()
    lowidx = np.zeros(1 << 11, np.int64)
    for i in range(11):
        lowidx[1 << i] = i
    ok = True

    # (W) the THM-5122 witness
    C = np.array([10, 28, 49, 36, 41, 3], dtype=np.int64)
    n = 6
    adj = np.zeros(n + 2, np.int64)
    sm = np.zeros(n + 3, np.int64)
    sc = np.zeros(n + 3, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            assert ((C[i] >> j) & 1) + ((C[j] >> i) & 1) == 1
    F = np.zeros(64, np.int64)
    all_F(n, C, F, adj, lowidx, sm, sc)
    x0 = 0b100101
    vals = [int(F[x0]), int(F[x0 ^ (1 << 2) ^ (1 << 4)]), int(F[x0 ^ (1 << 5) ^ (1 << 3)]),
            int(F[x0 ^ (1 << 2) ^ (1 << 4) ^ (1 << 5) ^ (1 << 3)])]
    Wm = [[0 if i == j else -int(F[(1 << i) | (1 << j)] - F[1 << i] - F[1 << j] + F[0]) for j in range(n)] for i in range(n)]
    good = vals == [345, 347, 347, 349] and Wm[0] == [0, 26, 14, 20, 26, 28]
    print("(W) THM-5122 witness: counts at bits 00, 10, 01, 11 =", vals, "; W row 0 =", Wm[0], "OK" if good else "FAIL")
    ok &= good
    ev = np.linalg.eigvalsh(np.array(Wm, dtype=float))
    print("    eigenvalues of W:", [round(float(e), 3) for e in ev], "; n_-(W) =", int((ev < -1e-9).sum()))

    for n in (4, 5, 6):
        p2 = disjoint_pair_pairs(n)
        m3 = disjoint_pair_triples(n) if n >= 6 else np.zeros((0, 6), dtype=np.int64)
        res, ex, nex = exhaustive(n, lowidx, p2, m3)
        good = res[0] == 0 and res[1] == 0
        print(f"(Q) n={n}: all {1 << (n*(n-1)//2)} labelled cores, all {1 << n} cuts: quadratic-identity failures {res[0]}; "
              f"violations of W_ij >= 2 H(C - i - j): {res[1]}", "OK" if good else "FAIL", f"[{time.time()-t0:.0f}s]")
        ok &= bool(good)
        print(f"(E2) n={n}: (core, two disjoint pairs) with zero interaction: {res[2]}; (core, pairs, background) "
              f"with both increments nonzero: {res[3]}")
        if n >= 6:
            print(f"(E3) n={n}: (core, perfect matching, background) affine: {res[4]}; with all three increments nonzero: "
                  f"{res[5]}; of those with all 8 closures strong: {res[6]}; interval-complete "
                  f"(weights 1, <=2, <=1+sum): {res[7]} (strong: {res[8]})")
            for i in range(nex):
                code = int(ex[i, 0])
                adjC = np.zeros(n, np.int64)
                set_adj(n, code, adjC)
                FF = np.zeros(1 << n, np.int64)
                a2 = np.zeros(n + 2, np.int64)
                all_F(n, adjC, FF, a2, lowidx, sm, sc)
                x0 = int(ex[i, 1])
                prs = [(int(ex[i, 2]), int(ex[i, 3])), (int(ex[i, 4]), int(ex[i, 5])), (int(ex[i, 6]), int(ex[i, 7]))]
                cube = []
                for b in range(8):
                    x = x0
                    for e in range(3):
                        if b >> e & 1:
                            x ^= (1 << prs[e][0]) ^ (1 << prs[e][1])
                    cube.append(int(FF[x]))
                print(f"     example: core out-masks {[int(v) for v in adjC]}, H(core) = "
                      f"{count_hp(n, adjC, lowidx, sm, sc)}, x0 = {x0:0{n}b}, exchanges {prs}, increments "
                      f"{[int(ex[i, 8]), int(ex[i, 9]), int(ex[i, 10])]}, eight counts {cube}, all closures strong: "
                      f"{bool(ex[i, 11])}")

    for n, trials in ((7, 20000), (8, 3000)):
        m3 = disjoint_pair_triples(n)
        res, ex, nex = sampled3(n, trials, 99 + n, lowidx, m3)
        print(f"(E3) n={n}: {trials} random cores, {len(m3)} triples of disjoint pairs each: affine (core, triple, "
              f"background): {res[0]}; all increments nonzero: {res[1]} (direct recount mismatches {res[2]}); "
              f"all 8 closures strong: {res[3]}", f"[{time.time()-t0:.0f}s]")
        ok &= res[2] == 0
        for i in range(nex):
            print(f"     example: core out-masks {[int(v) for v in ex[i, :n]]}, x0 = {int(ex[i, n]):0{n}b}, exchanges "
                  f"{[(int(ex[i, n+1]), int(ex[i, n+2])), (int(ex[i, n+3]), int(ex[i, n+4])), (int(ex[i, n+5]), int(ex[i, n+6]))]}, "
                  f"increments {[int(ex[i, n+7]), int(ex[i, n+8]), int(ex[i, n+9])]}, base count {int(ex[i, n+10])}")
    print("AUDIT EXCHANGE:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
