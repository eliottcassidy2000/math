#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): residue classes of H over ALL 4^n words of each length n (by count), from the
exact histograms of H modulo 720720 = 16 * 9 * 5 * 7 * 11 * 13 written by b1_census.py.

usage: python c3_local_from_census.py STATE_DIR NMAX
For every modulus m in (3, 5, 7, 9, 11, 13, 15, 21, 35, 8, 16):
  class factors f_n(c) = m * #{w : H(w) = c mod m} / 4^n  (for even m: over the odd classes, factor relative to the
  uniform law on the odd classes);
  twisted sums c_n(k) = 4^-n sum_w e(k H(w) / m), their moduli for the characters k, and the decay rate estimated
  by least squares of log |c_n(k)| over n in [NFIT0, NMAX] (the rate is rho(T_k) / 4 if a single mode dominates).
"""
import sys

import numpy as np

MOD = 720720


def main():
    sdir, nmax = sys.argv[1], int(sys.argv[2])
    hists = {}
    for n in range(1, nmax + 1):
        h = np.load(f"{sdir}/modhist_n{n:02d}.npy")
        assert int(h.sum()) == 4 ** n
        hists[n] = h.astype(np.float64)
    res = np.arange(MOD)
    summary = []
    for m in (3, 5, 7, 9, 11, 13, 15, 21, 35, 27 if MOD % 27 == 0 else 3, 8, 16):
        if MOD % m:
            continue
        print(f"===== modulus {m} =====")
        cls = {n: np.bincount(res % m, weights=hists[n], minlength=m) / 4.0 ** n for n in hists}
        odd_only = m % 2 == 0
        active = [c for c in range(m) if (c % 2 == 1 or not odd_only)]
        norm = len(active)
        print("class factors (share times the number of admissible classes):")
        for n in sorted(hists):
            if n >= nmax - 5 or n in (8, 10):
                f = cls[n][active] * norm
                print(f"  n = {n:2d}: " + " ".join(f"{c}:{v:.5f}" for c, v in zip(active, f))
                      + f"   max |f - 1| = {np.abs(f - 1).max():.5f}")
        # twisted sums
        ks = [k for k in range(1, m // 2 + 1)]
        print("twisted sums |c_n(k)| and fitted rates:")
        worst = (0.0, None)
        for k in ks:
            if odd_only and (2 * k) % m == 0:
                continue            # k = m/2: e(kH/m) = -1 identically on odd H
            ph = np.exp(2j * np.pi * k * np.arange(m) / m)
            a = np.array([abs(np.sum(cls[n] * ph)) for n in sorted(hists)])
            ns = np.array(sorted(hists))
            n0 = max(8, nmax - 8)
            sel = (ns >= n0) & (a > 0)
            rate = float(np.exp(np.polyfit(ns[sel], np.log(a[sel]), 1)[0])) if sel.sum() >= 3 else float("nan")
            last_ratio = a[-1] / a[-2] if a[-2] > 0 else float("nan")
            two_step = (a[-1] / a[-3]) ** 0.5 if a[-3] > 0 else float("nan")
            print(f"  k = {k:2d} (gcd with m: {np.gcd(k, m)}): |c_n| at n = {nmax - 3}..{nmax}: "
                  + " ".join(f"{v:.3e}" for v in a[-4:])
                  + f" | fit n >= {n0}: {rate:.4f}; last ratio {last_ratio:.4f}; two-step {two_step:.4f}")
            if rate > worst[0]:
                worst = (rate, k)
        summary.append((m, worst[1], worst[0]))
        print(f"  slowest character of modulus {m}: k = {worst[1]}, fitted rate {worst[0]:.4f} "
              f"(twisted spectral radius about {4 * worst[0]:.3f} of 4)")
    print()
    print("summary: modulus, slowest character, fitted rate per letter, 4 x rate")
    for m, k, rate in summary:
        print(f"  {m:3d}  k = {k}  {rate:.4f}  {4 * rate:.3f}")
    # ---- primitive characters only: root mean square of |c_n(k)| against the baseline 2^-n of 4^n independent
    #      uniform residues (for which E |c_n(k)|^2 = 4^-n)
    print()
    print("primitive characters: R_n(m) = rms over k prime to m of |c_n(k)|, shown as R_n * 2^n, and the rate fitted")
    print("on n = 10..NMAX (rate 0.5 = square-root cancellation = spectral radius 2 of 4)")
    for m in (5, 7, 11, 13, 35, 55, 65, 77, 91, 143, 9, 16):
        if MOD % m:
            continue
        ks = [k for k in range(1, m) if np.gcd(k, m) == 1]
        R = []
        for n in sorted(hists):
            cl = np.bincount(res % m, weights=hists[n], minlength=m) / 4.0 ** n
            cf = np.fft.fft(cl)
            R.append(np.sqrt(np.mean(np.abs(cf[ks]) ** 2)))
        R = np.array(R)
        ns = np.array(sorted(hists))
        sel = ns >= 10
        rate = float(np.exp(np.polyfit(ns[sel], np.log(R[sel]), 1)[0]))
        print(f"  m = {m:3d}: R_n 2^n for n = 10..{nmax}: " + " ".join(f"{R[i] * 2.0 ** ns[i]:6.2f}" for i in range(len(ns)) if ns[i] >= 10)
              + f" | fitted rate {rate:.4f} (rho about {4 * rate:.3f})")


if __name__ == "__main__":
    main()
