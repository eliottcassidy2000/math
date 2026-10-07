#!/usr/bin/env python3
"""Companion to twoorbit_exponents_20261007.py (mac-mini-2026-10-07-oaimath3).
  A. How fast the alignment coupling becomes Haar as the gap L grows: the excess
       eps(l) = E[2^(1-delta) | L_s = l] - 2/3   (aligned covariance = -3 eps)  and  P(delta >= 4 | L_s = l) / (1/8),
     for l = 1..40, from the same exact integer orbits (S19's lag-1 Mersenne pairs).
  B. The archimedean side: the normalized debt rho_t = D_t/2^(L_t) against the perpetuity
       Y = sum_(j>=0) 3^j / 2^(a_1+...+a_j),  a_i iid Geom(1/2)   (rho ~ -Y when L >> 0),
     E[(3/2^a)^theta] = 3^theta/(2^(theta+1) - 1) (= THM-4554's Moran function at theta + 1); Kesten index 1, so
     P(Y > u) is of order 1/u (lattice case: u P(Y > u) bounded above and below, log-periodic).
Usage: python3 twoorbit_decay_perpetuity_20261007.py [NPAIRS BITS SEED NPROC]   (defaults 3000 20000 2029 10)
"""
import sys, random, math, time
from fractions import Fraction
from multiprocessing import Pool
import numpy as np

args = [int(a) for a in sys.argv[1:5]]
NPAIRS, BITS, SEED, NPROC = (args + [3000, 20000, 2029, 10][len(args):])[:4]
MARGIN = 96
LMAX = 40


def v2(n):
    return (n & -n).bit_length() - 1


def run(idx):
    rng = random.Random(SEED * 1000003 + idx)
    y0 = rng.getrandbits(BITS) | 1
    v = 2
    while rng.random() < 0.5:
        v += 1
    x = 3 * (1 << v) * y0 + 1
    t3 = 3 * y0 + 1
    a0 = v2(t3)
    y = t3 >> a0
    L0 = v + a0
    L = L0
    ys, a_l, L_l = [], [], []
    A = B = 0
    Apos = {}
    consumed = a0
    sums = np.zeros((LMAX + 2, 3))          # per l: count, sum 2^(1-delta), count(delta >= 4)
    rhos = []
    for t in range(int((BITS - MARGIN) / 2.3)):
        if x == y:
            break
        Apos[A] = t
        ys.append(y)
        a = v2(3 * y + 1)
        b = v2(3 * x + 1)
        consumed += a
        if consumed > BITS - MARGIN:
            break
        s = Apos.get(B - L0)
        if s is not None and s < t:
            k = t - s
            E = 3 * (x - 3 ** k * ys[s]) + 1 - 3 ** k
            d = v2(E) if E != 0 else 60
            l = min(max(L_l[s], 0), LMAX + 1)
            sums[l] += (1, 2.0 ** (1 - min(d, 60)), d >= 4)
        # normalized debt every 7 steps when L >= 24 (float is enough: rho = D/2^L)
        if t % 7 == 0 and L >= 24:
            D = 3 * (x - (y << L)) + 1 - (1 << L)
            rhos.append(D / 2.0 ** L)
        a_l.append(a); L_l.append(L)
        A += a; B += b
        x = (3 * x + 1) >> b
        y = (3 * y + 1) >> a
        L = L + a - b
    return sums, rhos


def perpetuity(n, rng):
    """n samples of Y = sum_j 3^j / 2^(a_1+..+a_j), truncated when the term < 1e-18 (a.s. convergent)."""
    out = np.empty(n)
    for i in range(n):
        Y, term = 1.0, 1.0
        while term > 1e-18:
            a = 1
            while rng.random() < 0.5:
                a += 1
            term *= 3.0 / 2 ** a
            Y += term
            if Y > 1e300:
                break
        out[i] = Y
    return out


if __name__ == '__main__':
    t0 = time.time()
    with Pool(NPROC) as pool:
        res = pool.map(run, range(NPAIRS))
    sums = sum(r[0] for r in res)
    print(f"A. alignment coupling excess by gap L_s = l  (NPAIRS={NPAIRS}, BITS={BITS}; Haar: E[2^(1-delta)] = 2/3, P(delta>=4) = 1/8)")
    print("   l     n      E[2^(1-d)]-2/3   (s.e.)    P(d>=4)/(1/8)")
    for l in list(range(0, 21)) + [24, 28, 32, 36, 40]:
        n, s2, s4 = sums[l]
        if n < 50:
            continue
        m = s2 / n
        se = math.sqrt(max(m * (2 - m) - m * m, 1e-12) / n)    # crude: Var(2^(1-d)) <= E[2^(2-2d)]; use a safe bound
        print(f"   {l:3d} {int(n):8d}   {m - 2/3:+.5f}   ({se:.5f})    {s4 / n / 0.125:.3f}")
    rhos = np.array([r for res_ in res for r in res_[1]])
    print(f"B. normalized debt rho = D/2^L at L >= 24 (every 7th step): n = {len(rhos)}; P(rho < 0) = {np.mean(rhos < 0):.4f}; "
          f"min |rho| = {np.min(np.abs(rhos)):.3f}")
    rng = random.Random(SEED)
    Y = perpetuity(400000, rng)
    q = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    print("   quantiles of -rho:     " + " ".join(f"{np.quantile(-rhos, p):9.3f}" for p in q))
    print("   quantiles of Y (iid):  " + " ".join(f"{np.quantile(Y, p):9.3f}" for p in q))
    for u in (10, 30, 100, 300, 1000, 3000, 10000):
        print(f"   u = {u:6d}:  u P(-rho > u) = {u * np.mean(-rhos > u):.3f}    u P(Y > u) = {u * np.mean(Y > u):.3f}")
    th = np.linspace(0, 1.5, 7)
    print("   E[(3/2^a)^theta] = 3^theta/(2^(theta+1)-1): " + ", ".join(f"theta={t:.2f}: {3**t/(2**(t+1)-1):.4f}" for t in th)
          + "   (= 1 at theta = 0 and theta = 1: Kesten index 1)")
    print(f"[{time.time() - t0:.0f}s]")
