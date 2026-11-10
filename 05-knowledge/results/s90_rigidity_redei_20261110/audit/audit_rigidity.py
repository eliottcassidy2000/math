#!/usr/bin/env python3
"""Independent audit of THM-5124 (additive rigidity on presence cubes).  Own code, nothing imported from the packet.

Algorithms (deliberately different from the packet's subset DP and its quantified structure test):
  * H tables: one depth-first enumeration of ALL directed simple paths of T; Hs[mask] = number of directed paths
    whose vertex set is mask (each such path is a Hamiltonian path of T[mask]).  Cross-checked against brute-force
    permutation counting for every labelled tournament of order <= 5.
  * Structure (ii) of the pair lemma: tested through strong components.  Z = Y => Q => X with Q strong on four
    vertices  <=>  the strong component of Z containing a has exactly four vertices and contains b; then a, b are
    the two one-triangle vertices  <=>  each of a, b lies on exactly one cyclic triangle of that component.
  * "z lies on no cycle of T": tested directly (strong component of z is {z}), not through H.

Checks
  (P)  pair lemma, every labelled tournament of orders 3..7, every unordered pair {a, b}: hypothesis (i) <=> (ii),
       and the values c, 3c, 3c, 5c.  Counts expected 24, 240, 2160, 23520 at orders 4..7.
  (R)  statement 2, every labelled tournament of orders 3..7 and every ordered partition (B, W1, W2), W1, W2 nonempty:
       values 3c, 3c, 5c; exactly one vertex of each module active over B; structure (ii) on T[B + a + b] with hubs
       in B; every other module vertex is a singleton strong component of T; plus the proof step
       f(S1, S2) = f(S1, 0) + f(0, S2) - c and f(S1, 0) in {c, 3c} for all sub-selections (orders <= 6).
  (K)  statement 3, every labelled tournament of orders 3..7 and every unordered module triple over every base B:
       number of additive cubes with all three increments positive (expected 0).
  (T)  three pairwise additive active vertices over B = V - {a1, a2, a3} (orders 5..7): eight counts
       c (1; 3,3,3; 5,5,5; 9 or 15), 9 iff the triple is transitive.
  (S4) statement 4: H(S(A)) versus H(A) + D(A) for every labelled A of order 0..5 (order 0 is the degenerate case).
  (X)  random tournaments of orders 8, 9: (P) on all pairs, (R) on random partitions.
"""
import itertools
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def set_adj(n, code, adj, radj):
    for i in range(n):
        adj[i] = 0
        radj[i] = 0
    k = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> k) & 1:
                adj[i] |= 1 << j
                radj[j] |= 1 << i
            else:
                adj[j] |= 1 << i
                radj[i] |= 1 << j
            k += 1


@njit(cache=True)
def htable(n, adj, Hs, lowidx, st_mask, st_cand):
    size = 1 << n
    for m in range(size):
        Hs[m] = 0
    Hs[0] = 1
    for s in range(n):
        m0 = 1 << s
        Hs[m0] += 1
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
                Hs[nm] += 1
                depth += 1
                st_mask[depth] = nm
                st_cand[depth] = adj[w] & ~nm


@njit(cache=True)
def closure(adj, U, v, lowidx):
    seen = 1 << v
    frontier = seen
    while frontier:
        low = frontier & (-frontier)
        frontier ^= low
        w = lowidx[low]
        new = adj[w] & U & ~seen
        seen |= new
        frontier |= new
    return seen


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def component(adj, radj, U, v, lowidx):
    return closure(adj, U, v, lowidx) & closure(radj, U, v, lowidx)


@njit(cache=True)
def triangles_through(adj, radj, C, v, lowidx):
    t = 0
    outs = adj[v] & C
    ins = radj[v] & C
    while outs:
        low = outs & (-outs)
        outs ^= low
        x = lowidx[low]
        t += popcount(adj[x] & ins)
    return t


@njit(cache=True)
def structure_ii(adj, radj, U, a, b, lowidx):
    """T[U] = Y => Q => X, Q strong of order 4 containing a, b as its two one-triangle vertices.
    Returns the mask of the two hubs, or 0."""
    C = component(adj, radj, U, a, lowidx)
    if popcount(C) != 4 or not (C >> b) & 1:
        return 0
    if triangles_through(adj, radj, C, a, lowidx) != 1:
        return 0
    if triangles_through(adj, radj, C, b, lowidx) != 1:
        return 0
    return C & ~(1 << a) & ~(1 << b)


@njit(cache=True)
def run_pairs(n, lowidx):
    npairs = n * (n - 1) // 2
    size = 1 << n
    full = size - 1
    adj = np.zeros(n, np.int64)
    radj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    res = np.zeros(4, np.int64)     # hyp, struct, mismatch, value failures
    for code in range(1 << npairs):
        set_adj(n, code, adj, radj)
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for a in range(n):
            for b in range(a + 1, n):
                R = full ^ (1 << a) ^ (1 << b)
                c = Hs[R]
                ha = Hs[R | (1 << a)]
                hb = Hs[R | (1 << b)]
                hyp = ha > c and hb > c and (Hs[full] - ha - hb + c == 0)
                st = structure_ii(adj, radj, full, a, b, lowidx) != 0
                if hyp:
                    res[0] += 1
                if st:
                    res[1] += 1
                    if not (ha == 3 * c and hb == 3 * c and Hs[full] == 5 * c):
                        res[3] += 1
                if hyp != st:
                    res[2] += 1
    return res


@njit(cache=True)
def check_two_modules(n, adj, radj, Hs, B, W1, W2, lowidx, deep, res):
    """res[0] instances; res[1] value failures; res[2] active-vertex-count failures; res[3] structure failures;
    res[4] hub-not-in-B failures; res[5] inert-vertex failures; res[6] sub-selection additivity failures;
    res[7] sub-selection value failures (f(S,0) not in {c, 3c})."""
    c = Hs[B]
    f1 = Hs[B | W1]
    f2 = Hs[B | W2]
    f12 = Hs[B | W1 | W2]
    if not (f1 > c and f2 > c and f12 - f1 - f2 + c == 0):
        return
    res[0] += 1
    if not (f1 == 3 * c and f2 == 3 * c and f12 == 5 * c):
        res[1] += 1
    na = 0
    nb = 0
    a = -1
    b = -1
    for v in range(n):
        if (W1 >> v) & 1 and Hs[B | (1 << v)] > c:
            na += 1
            a = v
        if (W2 >> v) & 1 and Hs[B | (1 << v)] > c:
            nb += 1
            b = v
    if na != 1 or nb != 1:
        res[2] += 1
        return
    U = B | (1 << a) | (1 << b)
    hubs = structure_ii(adj, radj, U, a, b, lowidx)
    if hubs == 0:
        res[3] += 1
    elif (hubs & B) != hubs:
        res[4] += 1
    full = B | W1 | W2
    rest = (W1 | W2) & ~(1 << a) & ~(1 << b)
    for v in range(n):
        if (rest >> v) & 1:
            if component(adj, radj, full, v, lowidx) != (1 << v):
                res[5] += 1
    if deep:
        S1 = W1
        while True:
            fa = Hs[B | S1]
            if fa != c and fa != 3 * c:
                res[7] += 1
            S2 = W2
            while True:
                if Hs[B | S1 | S2] != fa + Hs[B | S2] - c:
                    res[6] += 1
                if S2 == 0:
                    break
                S2 = (S2 - 1) & W2
            if S1 == 0:
                break
            S1 = (S1 - 1) & W1
        S2 = W2
        while True:
            fb = Hs[B | S2]
            if fb != c and fb != 3 * c:
                res[7] += 1
            if S2 == 0:
                break
            S2 = (S2 - 1) & W2


@njit(cache=True)
def run_modules(n, lowidx, P2, P3, deep):
    npairs = n * (n - 1) // 2
    size = 1 << n
    adj = np.zeros(n, np.int64)
    radj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    res2 = np.zeros(8, np.int64)
    res3 = np.zeros(3, np.int64)    # triples examined, all three active, additive active cubes
    for code in range(1 << npairs):
        set_adj(n, code, adj, radj)
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for i in range(P2.shape[0]):
            check_two_modules(n, adj, radj, Hs, P2[i, 0], P2[i, 1], P2[i, 2], lowidx, deep, res2)
        for i in range(P3.shape[0]):
            B = P3[i, 0]
            W1 = P3[i, 1]
            W2 = P3[i, 2]
            W3 = P3[i, 3]
            c = Hs[B]
            d1 = Hs[B | W1] - c
            d2 = Hs[B | W2] - c
            d3 = Hs[B | W3] - c
            res3[0] += 1
            if d1 > 0 and d2 > 0 and d3 > 0:
                res3[1] += 1
                if (Hs[B | W1 | W2] == c + d1 + d2 and Hs[B | W1 | W3] == c + d1 + d3
                        and Hs[B | W2 | W3] == c + d2 + d3 and Hs[B | W1 | W2 | W3] == c + d1 + d2 + d3):
                    res3[2] += 1
    return res2, res3


@njit(cache=True)
def run_three_vertices(n, lowidx):
    """(T): triples pairwise in the situation of statement 1 over B = V - triple."""
    npairs = n * (n - 1) // 2
    size = 1 << n
    full = size - 1
    adj = np.zeros(n, np.int64)
    radj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    res = np.zeros(5, np.int64)     # instances, transitive (9c), cyclic (15c), failures, joint value 7c seen
    for code in range(1 << npairs):
        set_adj(n, code, adj, radj)
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for a in range(n):
            for b in range(a + 1, n):
                for d in range(b + 1, n):
                    ma = 1 << a
                    mb = 1 << b
                    md = 1 << d
                    B = full ^ ma ^ mb ^ md
                    c = Hs[B]
                    ha = Hs[B | ma]
                    hb = Hs[B | mb]
                    hd = Hs[B | md]
                    if not (ha > c and hb > c and hd > c):
                        continue
                    if Hs[B | ma | mb] - ha - hb + c != 0:
                        continue
                    if Hs[B | ma | md] - ha - hd + c != 0:
                        continue
                    if Hs[B | mb | md] - hb - hd + c != 0:
                        continue
                    res[0] += 1
                    cyc = Hs[ma | mb | md] == 3
                    top = Hs[full]
                    ok = (ha == 3 * c and hb == 3 * c and hd == 3 * c and Hs[B | ma | mb] == 5 * c
                          and Hs[B | ma | md] == 5 * c and Hs[B | mb | md] == 5 * c)
                    if cyc:
                        res[2] += 1
                        ok = ok and top == 15 * c
                    else:
                        res[1] += 1
                        ok = ok and top == 9 * c
                    if not ok:
                        res[3] += 1
                    if top == 7 * c:
                        res[4] += 1
    return res


@njit(cache=True)
def run_random(n, trials, seed, nparts, lowidx):
    np.random.seed(seed)
    size = 1 << n
    full = size - 1
    adj = np.zeros(n, np.int64)
    radj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    resp = np.zeros(4, np.int64)
    res2 = np.zeros(8, np.int64)
    for t in range(trials):
        for i in range(n):
            adj[i] = 0
            radj[i] = 0
        pfwd = 0.5 if t % 3 == 0 else (0.8 if t % 3 == 1 else 0.92)
        for j in range(1, n):
            for i in range(j):
                if np.random.random() < pfwd:
                    adj[i] |= 1 << j
                    radj[j] |= 1 << i
                else:
                    adj[j] |= 1 << i
                    radj[i] |= 1 << j
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for a in range(n):
            for b in range(a + 1, n):
                R = full ^ (1 << a) ^ (1 << b)
                c = Hs[R]
                ha = Hs[R | (1 << a)]
                hb = Hs[R | (1 << b)]
                hyp = ha > c and hb > c and (Hs[full] - ha - hb + c == 0)
                st = structure_ii(adj, radj, full, a, b, lowidx) != 0
                if hyp:
                    resp[0] += 1
                if st:
                    resp[1] += 1
                    if not (ha == 3 * c and hb == 3 * c and Hs[full] == 5 * c):
                        resp[3] += 1
                if hyp != st:
                    resp[2] += 1
        for q in range(nparts):
            B = 0
            W1 = 0
            W2 = 0
            for v in range(n):
                r = np.random.randint(0, 5)
                if r <= 2:
                    B |= 1 << v
                elif r == 3:
                    W1 |= 1 << v
                else:
                    W2 |= 1 << v
            if W1 == 0 or W2 == 0:
                continue
            check_two_modules(n, adj, radj, Hs, B, W1, W2, lowidx, True, res2)
    return resp, res2


def perm_htable(n, adj):
    """brute force: Hs[mask] = number of orderings of mask that are directed paths."""
    Hs = [0] * (1 << n)
    Hs[0] = 1
    for mask in range(1, 1 << n):
        vs = [v for v in range(n) if mask >> v & 1]
        cnt = 0
        for p in itertools.permutations(vs):
            good = True
            for i in range(len(p) - 1):
                if not (adj[p[i]] >> p[i + 1]) & 1:
                    good = False
                    break
            if good:
                cnt += 1
        Hs[mask] = cnt
    return Hs


def partitions2(n):
    out = []
    for code in itertools.product(range(3), repeat=n):
        B = sum(1 << v for v in range(n) if code[v] == 0)
        W1 = sum(1 << v for v in range(n) if code[v] == 1)
        W2 = sum(1 << v for v in range(n) if code[v] == 2)
        if W1 and W2:
            out.append((B, W1, W2))
    return np.array(out, dtype=np.int64).reshape(-1, 3)


def partitions3(n):
    out = []
    for code in itertools.product(range(4), repeat=n):
        W = [sum(1 << v for v in range(n) if code[v] == k) for k in range(4)]
        if W[1] and W[2] and W[3]:
            if (W[1] & -W[1]) < (W[2] & -W[2]) < (W[3] & -W[3]):
                # the universe is all of V; proper sub-universes are the same check at a smaller order
                out.append((W[0], W[1], W[2], W[3]))
    return np.array(out, dtype=np.int64).reshape(-1, 4)


def main():
    t0 = time.time()
    ok = True
    NMAX = 7
    lowidx = np.zeros(1 << 10, np.int64)
    for i in range(10):
        lowidx[1 << i] = i

    # ---- cross-check of the DFS table against permutation counting, all labelled tournaments of order <= 5
    bad = 0
    tot = 0
    for n in range(1, 6):
        adj = np.zeros(n, np.int64)
        radj = np.zeros(n, np.int64)
        Hs = np.zeros(1 << n, np.int64)
        sm = np.zeros(n + 1, np.int64)
        sc = np.zeros(n + 1, np.int64)
        for code in range(1 << (n * (n - 1) // 2)):
            set_adj(n, code, adj, radj)
            htable(n, adj, Hs, lowidx, sm, sc)
            ref = perm_htable(n, [int(x) for x in adj])
            tot += 1
            if list(Hs) != ref:
                bad += 1
    print(f"(0) DFS path tables vs brute-force permutations, all {tot} labelled tournaments of orders 1..5: "
          f"mismatches {bad}", "OK" if bad == 0 else "FAIL", flush=True)
    ok &= bad == 0

    expected = {3: 0, 4: 24, 5: 240, 6: 2160, 7: 23520}
    for n in range(3, NMAX + 1):
        res = run_pairs(n, lowidx)
        good = res[2] == 0 and res[3] == 0 and res[0] == res[1] == expected[n]
        print(f"(P) n={n}: pairs with hypothesis (i) {res[0]}, with structure (ii) {res[1]}, mismatches {res[2]}, "
              f"value failures {res[3]}", "OK" if good else "FAIL", f"[{time.time()-t0:.0f}s]", flush=True)
        ok &= bool(good)

    for n in range(3, NMAX + 1):
        P2 = partitions2(n)
        P3 = partitions3(n)
        res2, res3 = run_modules(n, lowidx, P2, P3, n <= 6)
        good = int(res2[1:].sum()) == 0 and res3[2] == 0
        print(f"(R) n={n}: ordered partitions per tournament {len(P2)}; active zero-interaction module pairs {res2[0]}; "
              f"failures [values, #active, structure, hubs in B, inert vertices, sub-additivity, sub-values] = "
              f"{[int(x) for x in res2[1:]]}", "OK" if int(res2[1:].sum()) == 0 else "FAIL", flush=True)
        print(f"(K) n={n}: unordered module triples per tournament {len(P3)}; examined {res3[0]}; all three active "
              f"{res3[1]}; additive active 3-cubes {res3[2]}", "OK" if res3[2] == 0 else "FAIL",
              f"[{time.time()-t0:.0f}s]", flush=True)
        ok &= bool(good)

    for n in range(5, NMAX + 1):
        res = run_three_vertices(n, lowidx)
        good = res[3] == 0 and res[4] == 0
        print(f"(T) n={n}: triples pairwise active with zero interaction over B = V - triple: {res[0]} "
              f"(transitive, 9c: {res[1]}; cyclic, 15c: {res[2]}); failures {res[3]}; joint value 7c: {res[4]}",
              "OK" if good else "FAIL", f"[{time.time()-t0:.0f}s]", flush=True)
        ok &= bool(good)

    # ---- statement 4: two-hub suspension
    for k in range(0, 6):
        n = k + 2
        adjA = np.zeros(max(k, 1), np.int64)
        radjA = np.zeros(max(k, 1), np.int64)
        HsA = np.zeros(1 << k, np.int64)
        adjS = np.zeros(n, np.int64)
        HsS = np.zeros(1 << n, np.int64)
        sm = np.zeros(n + 1, np.int64)
        sc = np.zeros(n + 1, np.int64)
        bad_plain = 0
        bad_fixed = 0
        cnt = 0
        fullA = (1 << k) - 1
        for code in range(1 << (k * (k - 1) // 2)):
            cnt += 1
            if k:
                set_adj(k, code, adjA, radjA)
                htable(k, adjA, HsA, lowidx, sm, sc)
            else:
                HsA[0] = 1
            D = sum(int(HsA[U]) * int(HsA[fullA ^ U]) for U in range(1 << k))
            # S(A): vertices 0..k-1 = A, k = x0, k+1 = y0; x0 -> y0, y0 -> A, A -> x0
            for i in range(k):
                adjS[i] = int(adjA[i]) | (1 << k)
            adjS[k] = 1 << (k + 1)
            adjS[k + 1] = fullA
            htable(n, adjS, HsS, lowidx, sm, sc)
            hS = int(HsS[(1 << n) - 1])
            if n <= 6:
                assert hS == perm_htable(n, [int(x) for x in adjS])[(1 << n) - 1]
            if hS != int(HsA[fullA]) + D:
                bad_plain += 1
            if hS != (int(HsA[fullA]) if k > 0 else 0) + D:
                bad_fixed += 1
        print(f"(S4) order of A = {k}: {cnt} labelled A; failures of H(S(A)) = H(A) + D(A) as stated: {bad_plain}; "
              f"failures with the first term dropped for empty A: {bad_fixed}",
              "OK" if bad_fixed == 0 else "FAIL", flush=True)
        ok &= bad_fixed == 0

    for n, trials, parts in ((8, 30000, 80), (9, 8000, 80)):
        resp, res2 = run_random(n, trials, 4242 + n, parts, lowidx)
        good = resp[2] == 0 and resp[3] == 0 and int(res2[1:].sum()) == 0
        print(f"(X) n={n}: {trials} random tournaments (three arc biases): pairs with (i) {resp[0]}, with (ii) {resp[1]}, "
              f"mismatches {resp[2]}; statement-2 instances {res2[0]}, failures {[int(x) for x in res2[1:]]}",
              "OK" if good else "FAIL", f"[{time.time()-t0:.0f}s]", flush=True)
        ok &= bool(good)

    print("AUDIT RIGIDITY:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
