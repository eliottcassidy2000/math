#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): the twisted transfer operators T_k of the two-apex count modulo m, exactly.

Structure used (proved in PROOFS.txt, section C): in the coordinates x_j(w) = H(w s_j) the vector x_0 = (1, ..., 1)
is fixed by the four letters, so with y_j = x_j - x_1 (j = 2..17, the BASE state, 16 coordinates)

    y(w l) = y(w) B_l    (linear),          H(w l) = H(w) + delta_l(y(w))    (delta_l linear),

and  sum_{|w| = n} e(k H(w) / m) = e(k / m) * sum over the base orbit of (T_k^n nu_0),
    (T_k nu)(y') = sum over (y, l) with y B_l = y' of e(k delta_l(y) / m) nu(y).

usage: python c4_twisted_base.py m [NMAX] [maxstates]
Enumerates the base orbit modulo m, then iterates T_k / 4 for every k = 1 .. m - 1 (up to conjugation) and prints
|c_n(k)| = |sum nu_n| and the l2 norm of nu_n; the decay rate is fitted on the second half of the range.
"""
import json
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def twisted_step(nu, nxt, ph, table, out):
    n = nu.shape[0]
    for i in range(n):
        out[i] = 0.0
    for l in range(4):
        for i in range(n):
            v = nu[i]
            if v != 0.0:
                out[nxt[l, i]] += v * table[ph[l, i]]
    for i in range(n):
        out[i] *= 0.25


def main():
    m = int(sys.argv[1])
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 160
    maxstates = int(sys.argv[3]) if len(sys.argv) > 3 else 60_000_000
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    A = [np.array(M, dtype=np.int64) for M in d["A"]]
    x0 = np.array(d["x0"], dtype=np.int64)
    # check the fixed vector
    one = np.ones(r, np.int64)
    assert all((one @ A[l] == one).all() for l in range(4)), "x_0 = (1,...,1) is not fixed"
    C = [A[l][1:, :] for l in range(4)]                 # y C_l = z in Z^17 : delta = z_1, y' = z[1:] - z_1
    y0 = (x0[1:] - x0[0]) % m
    pw = m ** np.arange(r - 1, dtype=np.int64)

    def enc(Y):
        return Y @ pw

    def dec(codes):
        Y = np.empty((len(codes), r - 1), np.int64)
        c = codes.copy()
        for j in range(r - 1):
            Y[:, j] = c % m
            c //= m
        return Y

    def step(Y, l):
        Z = Y @ C[l]
        delta = Z[:, 0] % m
        Yn = (Z[:, 1:] - Z[:, :1]) % m
        return Yn, delta

    t0 = time.time()
    seen = np.array([int(enc(y0[None, :])[0])], dtype=np.int64)
    frontier = seen.copy()
    depth = 0
    closed = False
    while True:
        Y = dec(frontier)
        new = np.unique(np.concatenate([enc(step(Y, l)[0]) for l in range(4)]))
        new = new[~np.isin(new, seen, assume_unique=True)]
        depth += 1
        if len(new) == 0:
            closed = True
            break
        seen = np.union1d(seen, new)
        frontier = new
        if len(seen) > maxstates:
            break
    n = len(seen)
    print(f"m = {m}: base orbit {n} states, depth {depth}, {'closed' if closed else 'NOT closed'}, "
          f"{time.time() - t0:.1f} s", flush=True)
    if not closed:
        return
    nxt = np.empty((4, n), np.int32)
    ph = np.empty((4, n), np.uint8)
    CH = 2_000_000
    for lo in range(0, n, CH):
        Y = dec(seen[lo:lo + CH])
        for l in range(4):
            Yn, delta = step(Y, l)
            codes = enc(Yn)
            pos = np.searchsorted(seen, codes)
            assert (seen[pos] == codes).all()
            nxt[l, lo:lo + CH] = pos
            ph[l, lo:lo + CH] = delta
    start = int(np.searchsorted(seen, enc(y0[None, :])[0]))
    print("k, n: |c_n(k)| (twisted sum over all 4^n words, normalised) and the l2 norm of the twisted measure")
    results = {}
    ks = [k for k in range(1, m) if k <= m - k]
    for k in ks:
        if m % 2 == 0 and (2 * k) % m == 0:
            continue
        table = np.exp(2j * np.pi * k * np.arange(m) / m)
        nu = np.zeros(n, np.complex128)
        nu[start] = 1.0
        out = np.zeros(n, np.complex128)
        cs, l2 = [], []
        for stepn in range(1, nmax + 1):
            twisted_step(nu, nxt, ph, table, out)
            nu, out = out, nu
            cs.append(abs(nu.sum()))
            l2.append(float(np.sqrt((np.abs(nu) ** 2).sum())))
        cs = np.array(cs)
        l2 = np.array(l2)
        ns = np.arange(1, nmax + 1)
        half = ns >= nmax // 2
        ok = half & (cs > 1e-290)
        rate_c = float(np.exp(np.polyfit(ns[ok], np.log(cs[ok]), 1)[0])) if ok.sum() > 4 else float("nan")
        # the l2 norm of nu_n over a state space of size n: rate of ||T_k^n nu_0||
        rate_l2 = float(np.exp(np.polyfit(ns[half], np.log(l2[half]), 1)[0]))
        results[k] = (rate_c, rate_l2)
        pts = [t for t in (8, 12, 16, 17, 20, 24, 32, 48, 64, 96, 128, 160, 200) if t <= nmax]
        print(f"k = {k} (gcd {np.gcd(k, m)}): |c_n| at n = {pts}: " + " ".join(f"{cs[t - 1]:.3e}" for t in pts))
        print(f"   fitted on n >= {nmax // 2}: rate of |c_n| = {rate_c:.5f}, rate of the l2 norm = {rate_l2:.5f} "
              f"-> rho(T_k) = {4 * rate_c:.4f} (from c_n), {4 * rate_l2:.4f} (from the norm)", flush=True)
    print("summary for modulus", m, ":", {k: (round(v[0], 5), round(v[1], 5)) for k, v in results.items()})


if __name__ == "__main__":
    main()
