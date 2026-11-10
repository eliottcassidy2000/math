#!/usr/bin/env python3
"""Audit of the items added to HYP-9657 / README after the pre-audit checkpoint (working-tree state of 2026-11-10,
about 12:05): (C) 3 "second-order gap" and (B') the O/Z/J mass profile.  Own code.

(G)  For two vertices a, b active over R = V - {a, b}: c = H(R), d_a, d_b the increments, mu the interaction.
     Candidate: mu = 0 or 2 mu >= max(d_a, d_b).
     * exhaustive, every labelled tournament of orders 4..7, EVERY pair (the packet's committed script computes
       orders 4-6 for the pair {0, 1} and prints the order-7 line as a constant): minimum of mu / max(d) over
       mu > 0, all equality cases 2 mu = max(d), number of violations;
     * the order-8 witness [246, 252, 248, 177, 160, 64, 152, 32], a = 1, b = 2: counts by depth-first enumeration;
     * random search at orders 8, 9, 10 (biased tournaments): least mu / max(d) found.
(O)  The S, B, F system: states (h, R, e, u) from (1, 1, 1, 1); m(N) = number of distinct states with h = N reached
     right after a B or F letter, N <= 12000; class means on odd N in [3000, 12000]; misses; primes; rough; Heegner.
"""
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit


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
def gap_exhaustive(n, lowidx):
    npairs = n * (n - 1) // 2
    size = 1 << n
    full = size - 1
    adj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    best_mu = 1
    best_d = 1
    best = np.zeros(4, np.int64)
    viol = 0
    inter = 0
    eq = np.zeros((64, 4), np.int64)
    neq = 0
    for code in range(1 << npairs):
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
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for a in range(n):
            for b in range(a + 1, n):
                R = full ^ (1 << a) ^ (1 << b)
                c = Hs[R]
                da = Hs[R | (1 << a)] - c
                db = Hs[R | (1 << b)] - c
                mu = Hs[full] - c - da - db
                if da > 0 and db > 0 and mu > 0:
                    inter += 1
                    d = max(da, db)
                    if 2 * mu < d:
                        viol += 1
                    if mu * best_d < best_mu * d:
                        best_mu = mu
                        best_d = d
                        best[0] = c
                        best[1] = da
                        best[2] = db
                        best[3] = mu
                    if 2 * mu == d:
                        lo = min(da, db)
                        found = False
                        for t in range(neq):
                            if eq[t, 0] == c and eq[t, 1] == d and eq[t, 2] == lo and eq[t, 3] == mu:
                                found = True
                        if not found and neq < 64:
                            eq[neq, 0] = c
                            eq[neq, 1] = d
                            eq[neq, 2] = lo
                            eq[neq, 3] = mu
                            neq += 1
    return inter, viol, best, eq[:neq]


@njit(cache=True)
def gap_random(n, trials, seed, lowidx):
    np.random.seed(seed)
    size = 1 << n
    full = size - 1
    adj = np.zeros(n, np.int64)
    Hs = np.zeros(size, np.int64)
    st_mask = np.zeros(n + 1, np.int64)
    st_cand = np.zeros(n + 1, np.int64)
    best_mu = 1
    best_d = 1
    best = np.zeros(4 + n + 2, np.int64)
    for t in range(trials):
        pf = 0.5 + 0.45 * np.random.random()
        for i in range(n):
            adj[i] = 0
        for j in range(1, n):
            for i in range(j):
                if np.random.random() < pf:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        htable(n, adj, Hs, lowidx, st_mask, st_cand)
        for a in range(n):
            for b in range(a + 1, n):
                R = full ^ (1 << a) ^ (1 << b)
                c = Hs[R]
                da = Hs[R | (1 << a)] - c
                db = Hs[R | (1 << b)] - c
                mu = Hs[full] - c - da - db
                if da > 0 and db > 0 and mu > 0:
                    d = max(da, db)
                    if mu * best_d < best_mu * d:
                        best_mu = mu
                        best_d = d
                        best[0] = c
                        best[1] = da
                        best[2] = db
                        best[3] = mu
                        for i in range(n):
                            best[4 + i] = adj[i]
                        best[4 + n] = a
                        best[5 + n] = b
    return best


def hp_count(out, mask):
    total = 0
    stack = []
    m = mask
    while m:
        low = m & -m
        m ^= low
        stack.append((low.bit_length() - 1, low))
    if mask == 0:
        return 1
    while stack:
        v, used = stack.pop()
        if used == mask:
            total += 1
            continue
        c = out[v] & mask & ~used
        while c:
            low = c & -c
            c ^= low
            stack.append((low.bit_length() - 1, used | low))
    return total


def main():
    t0 = time.time()
    ok = True
    lowidx = np.zeros(1 << 11, np.int64)
    for i in range(11):
        lowidx[1 << i] = i
    for n in range(4, 8):
        inter, viol, best, eq = gap_exhaustive(n, lowidx)
        if inter == 0:
            print(f"(G) n={n}: no pair with both vertices active and mu > 0")
            continue
        print(f"(G) n={n}: every labelled tournament, every pair: active pairs with mu > 0: {inter}; violations of "
              f"2 mu >= max(d): {viol}; least mu/max(d) = {best[3]}/{max(best[1], best[2])} at (c, d_a, d_b, mu) = "
              f"{tuple(int(x) for x in best)}; equality cases (c, max d, min d, mu): "
              f"{sorted(tuple(int(x) for x in r) for r in eq)}", f"[{time.time()-t0:.0f}s]", flush=True)
        ok &= viol == 0
    w = [246, 252, 248, 177, 160, 64, 152, 32]
    n = 8
    assert all(((w[i] >> j) & 1) + ((w[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n))
    a, b = 1, 2
    full = (1 << n) - 1
    R = full ^ (1 << a) ^ (1 << b)
    vals = [hp_count(w, R), hp_count(w, R | (1 << a)), hp_count(w, R | (1 << b)), hp_count(w, full)]
    c = vals[0]
    da, db = vals[1] - c, vals[2] - c
    mu = vals[3] - c - da - db
    good = vals == [15, 29, 29, 49] and 2 * mu < max(da, db)
    print(f"(G) order-8 witness {w}, a = {a}, b = {b}: counts {vals}, (c, d_a, d_b, mu) = {(c, da, db, mu)}; "
          f"2 mu < max(d): {2 * mu < max(da, db)}", "OK" if good else "FAIL")
    ok &= good
    for n, trials in ((8, 400000), (9, 150000), (10, 40000)):
        best = gap_random(n, trials, 777 + n, lowidx)
        c, da, db, mu = (int(x) for x in best[:4])
        print(f"(G) n={n}: {trials} random biased tournaments, all pairs: least mu/max(d) found = {mu}/{max(da, db)} = "
              f"{Fraction(mu, max(da, db))} ({mu / max(da, db):.4f}) at (c, d_a, d_b, mu) = {(c, da, db, mu)}, masks "
              f"{[int(x) for x in best[4:4 + n]]}, a = {int(best[4 + n])}, b = {int(best[5 + n])}",
              f"[{time.time()-t0:.0f}s]", flush=True)

    # ---- (O)
    X = 12000

    def step(st, ch):
        h, Rr, e, u = st
        if ch == 0:
            return (h, 2 * Rr, 0, u)
        if ch == 1:
            return (h + Rr - e, 2 * Rr - u, u, u)
        return (h + Rr - e, h + Rr, h - e, h)

    start = (1, 1, 1, 1)
    seen = {start}
    stack = [start]
    bf = set()
    while stack:
        st = stack.pop()
        for ch in (0, 1, 2):
            t = step(st, ch)
            if t[0] > X or t[1] > 8 * X:          # a wider cap than the packet's 4X
                continue
            if ch:
                bf.add(t)
            if t not in seen:
                seen.add(t)
                stack.append(t)
    m = np.zeros(X + 1, dtype=np.int64)
    for st in bf:
        m[st[0]] += 1
    odd = np.arange(1, X + 1, 2)
    miss = [int(N) for N in odd if m[N] == 0]
    s197 = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637,
            821, 1067, 1339, 1861, 2647]
    win = odd[odd >= X // 4]
    mw = m[win].astype(float)
    print(f"(O) S/B/F system: {len(bf)} states after B or F with h <= {X}; odd N <= {X} with m = 0: {len(miss)}, equal to "
          f"the 29 values: {miss == s197}; even h among B/F states: {sum(1 for s in bf if s[0] % 2 == 0)}")
    ok &= miss == s197
    print(f"    window [{X // 4}, {X}]: mean m = {mw.mean():.2f}")
    for p in (3, 5, 7):
        print(f"    p = {p}: class means", " ".join(f"{mw[win % p == c].mean():.1f}" for c in range(p)))

    def isprime(v):
        return v > 1 and all(v % q for q in range(2, int(v ** 0.5) + 1))
    pr = np.array([isprime(int(N)) for N in win])
    rough = np.array([isprime(int(N)) and N % 4 == 3 and isprime(int((N + 1) // 4)) for N in win])
    print(f"    primes: ratio {mw[pr].mean() / mw.mean():.3f}; N prime and (N+1)/4 prime: n = {rough.sum()}, ratio "
          f"{mw[rough].mean() / mw.mean():.3f}; m at Heegner numbers: {dict((d, int(m[d])) for d in (3, 7, 11, 19, 43, 67, 163))}")
    good = (abs(mw.mean() - 72.06) < 0.01 and abs(mw[win % 3 == 0].mean() - 101.3) < 0.06
            and abs(mw[pr].mean() / mw.mean() - 0.677) < 0.002 and abs(mw[rough].mean() / mw.mean() - 0.649) < 0.002)
    print("    packet numbers (72.06; 101.3; 0.677; 0.649) reproduced:", "OK" if good else "FAIL")
    ok &= good
    print("AUDIT NEW ITEMS:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
