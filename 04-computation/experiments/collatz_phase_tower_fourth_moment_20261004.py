#!/usr/bin/env python3
"""
E5p (opus, 2026-10-04): second and fourth moments of f_n(0) for the random-start Pascal tower (the Collatz tower with
a random 2-adic start) against the i.i.d.-digit model, exactly, at small n (Fourier note 4k).

Pascal tower: the path phase is e(R Phi(a)), Phi(a) = sum_i 3^(n-i) 2^-D_i (mod 1), R uniform mod 2^M.
  E|f|^2 = sum over pairs of 1[Phi(a) = Phi(b)]      -- and Phi is injective on paths (PROVED), so E|f|^2 = 3^-n exactly;
  E|f|^4 = sum over quadruples w_a w_b w_c w_d 1[Phi(a) + Phi(b) = Phi(c) + Phi(d) mod 1]  (weighted additive energy).
i.i.d. model: independent strings per level; the quadruple condition is per level
  {D_i(a), D_i(b)} = {D_i(c), D_i(d)} as multisets for every i (dyadic sums 2^-x + 2^-y with x, y >= 1 coincide
  only as multisets), which is stricter than the single global equation of the Pascal tower; so E|f|^4_Pascal >=
  E|f|^4_iid: the same second moment, a heavier tail.
Both fourth moments are computed by hashing pair sums / pair-depth-multisets over all A^n paths.
Usage: python collatz_phase_tower_fourth_moment_20261004.py [NMAX=4] [A=8]
"""
import sys, math, time
from itertools import product
import numpy as np

def paths_and_phi(n, A, q=3):
    paths = list(product(range(1, A + 1), repeat=n))
    Dmax = n * A
    keys = np.zeros(len(paths), dtype=np.int64)
    w = np.zeros(len(paths))
    Dv = np.zeros((len(paths), n), dtype=np.int64)
    for p, a in enumerate(paths):
        D = 0; key = 0
        for i in range(n, 0, -1):
            D += a[i - 1]
            Dv[p, i - 1] = D
            key = (key + q ** (n - i) * (1 << (Dmax - D))) % (1 << Dmax)
        keys[p] = key; w[p] = 2.0 ** (-sum(a))
    return keys, w, Dv, Dmax

def energy(keys, w):
    """sum over groups of equal key of (sum of weights)^2."""
    order = np.argsort(keys, kind="stable")
    k = keys[order]; ww = w[order]
    # group boundaries
    boundaries = np.concatenate(([0], np.flatnonzero(k[1:] != k[:-1]) + 1, [len(k)]))
    sums = np.add.reduceat(ww, boundaries[:-1])
    return float(np.sum(sums ** 2))

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    qs = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [3]
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"Pascal tower (multiplier q) versus i.i.d. model, valuations <= {A}: moments scaled by 3^n (second) and 9^n (fourth)")
    P("q | n | #paths | 3^n E|f|^2 tower | 3^n E|f|^2 iid | 9^n E|f|^4 tower | 9^n E|f|^4 iid | ratio | (complex Gaussian reference 2)")
    t0 = time.time()
    for q in qs:
      for n in range(1, NMAX + 1):
        keys, w, Dv, Dmax = paths_and_phi(n, A, q)
        m2_pascal = energy(keys, w)                       # pairs with equal Phi
        m2_iid = float(np.sum(w ** 2))                    # diagonal only (independent strings)
        # fourth moments over all ordered pairs (a,b): key of the pair
        P_ = len(keys)
        ia, ib = np.meshgrid(np.arange(P_), np.arange(P_), indexing="ij")
        ia = ia.ravel(); ib = ib.ravel()
        wpair = w[ia] * w[ib]
        kp = (keys[ia] + keys[ib]) % (1 << Dmax)          # Phi(a) + Phi(b) mod 1 (scaled)
        m4_pascal = energy(kp, wpair)
        # iid: per-level unordered depth pairs; encode (min, max) per level into one integer key
        key_iid = np.zeros(len(ia), dtype=np.int64)
        for i in range(n):
            lo = np.minimum(Dv[ia, i], Dv[ib, i]); hi = np.maximum(Dv[ia, i], Dv[ib, i])
            key_iid = key_iid * (A + 1) ** 2 + lo * (A + 1) + hi
        m4_iid = energy(key_iid, wpair)
        P(f"{q} | {n} | {P_} | {3**n * m2_pascal:.6f} | {3**n * m2_iid:.6f} | {9**n * m4_pascal:.5f} | {9**n * m4_iid:.5f} | {m4_pascal / m4_iid:.4f}   [{time.time()-t0:.0f}s]")
    P("reading: equal second moments (no path collisions), and the Pascal tower's fourth moment exceeds the i.i.d. model's: heavier tails from a single shared string.")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
