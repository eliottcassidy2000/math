#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (c): growth exponent of the tower slopes (sign-flip towers of THM-5130 (C)).

Level map (THM-5130 (1)), parameters (delta, i, m), delta in {1, 2}, i >= 1, m >= 0, a = i + delta + m,
(alpha, gamma) = (3, 5) or (5, 13):
    r+ = -r,  h+ = r + 2^a M,  M+ = alpha 2^m (h + r) + 2 gamma 2^(i-1) 3^m M - [i = 1] E,
    E+ = alpha 2^(i-1) M + [m = 0] (r - h).
Bases: (TT_n, {p}) with state (1, 2^g, 1, [g = 0]), g = n - p >= 0;  (TT_a, {i, i + delta}) with state
(-1, alpha 2^m, 1, [m = 0]).
T0  the level map against a direct count (subset recursion on the explicit tournament), small towers.
T1  sigma_0: the root of  rho(s) = sum over (delta, i, m) of lambda^(-s) = 1,  lambda = gamma 2^i 3^m - [i = 1]
    (PROOFS.txt, Theorem T: M+ >= lambda M, so the number of tower slopes <= X is O(X^s) for every s > sigma_0).
T2  exact enumeration of all tower states with slope M <= X = 2^K: number of states, of distinct slopes of each
    sign, of distinct values 2^t M <= X (t >= 0), local exponents; odd parts c <= 4097 reached.
Usage: python3 c_towers.py K      (default 30)
"""
import math
import sys
import time

import numpy as np
from numba import njit

from seedlib import dp_tables, pair_mats, closed_out

ALPHA = (0, 3, 5)
GAMMA = (0, 5, 13)


def level(state, delta, i, m):
    r, M, h, E = state
    al, ga = ALPHA[delta], GAMMA[delta]
    a = i + delta + m
    return (-r, al * 2 ** m * (h + r) + 2 * ga * 2 ** (i - 1) * 3 ** m * M - (E if i == 1 else 0),
            r + 2 ** a * M, al * 2 ** (i - 1) * M + ((r - h) if m == 0 else 0))


def direct_state(n, out, O):
    fe = np.zeros((1 << n, n), np.int64)
    gs = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    Q = np.zeros((n, n), np.int64)
    U = np.zeros(n, np.int64)
    dp_tables(n, out, fe, gs, hS)
    pair_mats(n, fe, gs, hS, Q, U)
    vs = [v for v in range(n) if (O >> v) & 1]
    full = (1 << n) - 1
    h = int(hS[full])
    W = sum(int(Q[e, f]) for e in vs for f in vs if e != f)
    return (h - W, sum(int(U[v]) for v in vs), h, sum(int(fe[full, v]) for v in vs))


def check_level_map():
    cnt = 0
    ok = True
    for n in (1, 2, 3):
        for p in range(1, n + 1):
            out = np.zeros(n, np.int64)
            for x in range(n):
                for y in range(x + 1, n):
                    out[x] |= 1 << y
            O = 1 << (p - 1)
            st = (1, 2 ** (n - p), 1, 1 if p == n else 0)
            if direct_state(n, out, O) != st:
                ok = False
            todo = [(n, out, O, st, 0)]
            while todo:
                nn, oo, OO, s, depth = todo.pop()
                if depth == 2:
                    continue
                for delta in (1, 2):
                    for i in (1, 2, 3):
                        for m in (0, 1, 2):
                            a = i + delta + m
                            if nn + a + 1 > 13:
                                continue
                            res = np.zeros(nn + a + 1, np.int64)
                            N = closed_out(nn, oo, OO, a, res)
                            P = (1 << (nn + i - 1)) | (1 << (nn + i + delta - 1))
                            s2 = level(s, delta, i, m)
                            if direct_state(N, res, P) != s2:
                                ok = False
                            cnt += 1
                            if depth == 0 and N <= 8:
                                todo.append((N, res, P, s2, depth + 1))
    return cnt, ok


def rho(s, terms=200):
    tot = 0.0
    for delta in (1, 2):
        ga = GAMMA[delta]
        for m in range(terms):
            tot += (2 * ga * 3.0 ** m - 1) ** (-s)
        tot += ga ** (-s) * 4.0 ** (-s) / (1 - 2.0 ** (-s)) / (1 - 3.0 ** (-s))
    return tot


def rho_plain(s):
    return sum((2 * GAMMA[d]) ** (-s) for d in (1, 2)) / (1 - 2.0 ** (-s)) / (1 - 3.0 ** (-s))


def root(f):
    lo, hi = 0.3, 1.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


@njit(cache=True)
def enumerate_towers(X, outM, outr, counts_by_depth):
    """all tower states with M <= X; returns the number of states (records the first len(outM))."""
    al = np.array([0, 3, 5], np.int64)
    ga = np.array([0, 5, 13], np.int64)
    MAXD = 64
    sr = np.zeros(MAXD, np.int64)
    sM = np.zeros(MAXD, np.int64)
    sh = np.zeros(MAXD, np.int64)
    sE = np.zeros(MAXD, np.int64)
    sd = np.zeros(MAXD, np.int64)
    si = np.zeros(MAXD, np.int64)
    sm = np.zeros(MAXD, np.int64)
    cap = outM.shape[0]
    cnt = 0
    # base states
    nb = 0
    g = 0
    while (1 << g) <= X:
        nb += 1
        g += 1
    for b in range(nb + 2 * 56):
        if b < nb:
            r0, M0, h0, E0 = 1, 1 << b, 1, (1 if b == 0 else 0)
        else:
            bb = b - nb
            delta0 = 1 + (bb & 1)
            m0 = bb >> 1
            M0 = al[delta0] * (1 << m0)
            if M0 > X:
                continue
            r0, h0, E0 = -1, 1, (1 if m0 == 0 else 0)
        if cnt < cap:
            outM[cnt] = M0
            outr[cnt] = r0
        cnt += 1
        counts_by_depth[0] += 1
        depth = 0
        sr[0], sM[0], sh[0], sE[0] = r0, M0, h0, E0
        sd[0], si[0], sm[0] = 1, 1, 0
        while depth >= 0:
            delta = sd[depth]
            if delta > 2:
                depth -= 1
                continue
            i = si[depth]
            m = sm[depth]
            M = sM[depth]
            p3 = 1
            for _ in range(m):
                p3 *= 3
            lam = ga[delta] * (1 << i) * p3 - (1 if i == 1 else 0)
            if lam > X // M + 1 or lam * M > X:
                if m == 0:
                    sd[depth] = delta + 1
                    si[depth] = 1
                else:
                    si[depth] = i + 1
                sm[depth] = 0
                continue
            sm[depth] = m + 1
            r = sr[depth]
            h = sh[depth]
            E = sE[depth]
            a = i + delta + m
            M2 = al[delta] * (1 << m) * (h + r) + 2 * ga[delta] * (1 << (i - 1)) * p3 * M - (E if i == 1 else 0)
            if M2 > X:
                continue
            h2 = r + (1 << a) * M
            E2 = al[delta] * (1 << (i - 1)) * M + ((r - h) if m == 0 else 0)
            if cnt < cap:
                outM[cnt] = M2
                outr[cnt] = -r
            cnt += 1
            depth += 1
            counts_by_depth[depth] += 1
            sr[depth], sM[depth], sh[depth], sE[depth] = -r, M2, h2, E2
            sd[depth], si[depth], sm[depth] = 1, 1, 0
    return cnt


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    cnt, ok = check_level_map()
    P(f"T0  level map = direct count (subset recursion) on {cnt} one- and two-level towers of order <= 13: {ok}")
    s0 = root(rho)
    s1 = root(rho_plain)
    P(f"T1  sigma_0 (lambda = gamma 2^i 3^m - [i = 1]): {s0:.6f};  without the correction (lambda = gamma 2^i 3^m): {s1:.6f};  rho(1) = {rho(1.0):.6f}")
    X = 1 << K
    t0 = time.time()
    cbd = np.zeros(64, np.int64)
    cap = enumerate_towers(X, np.zeros(0, np.int64), np.zeros(0, np.int8), cbd)      # counting pass
    outM = np.zeros(cap, np.int64)
    outr = np.zeros(cap, np.int8)
    cbd = np.zeros(64, np.int64)
    n = enumerate_towers(X, outM, outr, cbd)
    P(f"T2  X = 2^{K}: {n} tower states with slope <= X  ({time.time() - t0:.1f} s); by depth: {[int(c) for c in cbd if c]}")
    if n > cap:
        P("    buffer overflow, distinct counts not computed")
    else:
        Mp = np.unique(outM[:n][outr[:n] == 1])
        Mm = np.unique(outM[:n][outr[:n] == -1])
        P(f"    states with r = +1: {int((outr[:n] == 1).sum())}, distinct slopes {Mp.shape[0]};  r = -1: {int((outr[:n] == -1).sum())}, distinct slopes {Mm.shape[0]}")
        P("    k, #states with M <= 2^k, #distinct slopes <= 2^k (sign +, sign -), #distinct values 2^t M <= 2^k (t >= 0; sign +, sign -), log2 of the last two")
        allM = np.sort(outM[:n])
        prev = None
        for k in range(8, K + 1, 2):
            b = 1 << k
            ns = int(np.searchsorted(allM, b, side="right"))
            dp = int(np.searchsorted(Mp, b, side="right"))
            dm = int(np.searchsorted(Mm, b, side="right"))
            vals = []
            for arr in (Mp, Mm):
                sub = arr[arr <= b]
                acc = [sub]
                cur = sub
                while True:
                    cur = cur[cur <= b // 2] * 2
                    if cur.shape[0] == 0:
                        break
                    acc.append(cur)
                vals.append(int(np.unique(np.concatenate(acc)).shape[0]))
            line = f"    k = {k:2d}: states {ns:10d}  distinct M {dp:9d} {dm:9d}  values {vals[0]:9d} {vals[1]:9d}  log2 {math.log2(max(vals[0], 1)):.3f} {math.log2(max(vals[1], 1)):.3f}"
            if prev is not None:
                line += f"  local exponents (states, +, -): {(math.log2(ns) - math.log2(prev[0])) / 2:.3f} {(math.log2(vals[0]) - math.log2(prev[1])) / 2:.3f} {(math.log2(vals[1]) - math.log2(prev[2])) / 2:.3f}"
            prev = (ns, vals[0], vals[1])
            P(line)
        for name, arr in (("+", Mp), ("-", Mm)):
            odd = arr.copy()
            while True:
                ev = (odd & 1) == 0
                if not ev.any():
                    break
                odd[ev] >>= 1
            odd = np.unique(odd)
            for bound in (129, 1025, 4097, 16385):
                c = int((odd <= bound).sum())
                P(f"    sign {name}: odd c <= {bound} that are the odd part of a tower slope M <= 2^{K}: {c} of {(bound + 1) // 2}")
            miss = sorted(set(range(1, 130, 2)) - set(int(x) for x in odd[odd <= 129]))
            P(f"    sign {name}: odd c <= 129 not reached: {miss}")
    with open(f"c_towers_{K}.out", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
