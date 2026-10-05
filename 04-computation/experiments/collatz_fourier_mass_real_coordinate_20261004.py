#!/usr/bin/env python3
"""
Experiment E3 (opus, 2026-10-04): where does the Fourier mass of the 3-adic Syracuse law sit, in the REAL
frequency coordinate theta = t/3^n in R/Z?

Motivation (E1): frequencies that are 3-adically fixed (t = u for a fixed unit u) decay at ~0.57 per level, below
the Parseval rate 3^(-1/2); S20/S23 found the maxima at t = +-2^s with 2^s/3^n ~ 2^-6, i.e. at frequencies that
converge ARCHIMEDEANLY (theta -> ~1/64).  Hypothesis: the Parseval mass at level n is carried by frequencies t
whose real position theta = t/3^n is close to a dyadic rational of small denominator, not by 3-adically fixed ones.

For n <= NMAX we compute the full law (dense vector, valuations <= A), its FFT (all 3^n coefficients), and report:
  * the mass |mu_hat|^2 summed over the units, and its distribution in 64 bins of theta;
  * the top-30 coefficients with t, theta (binary), and whether t is +-2^s;
  * the mass fraction carried by frequencies within 2^-m of a dyadic rational with denominator <= 2^m, for m = 4..8,
    against the Haar (uniform) expectation of that set;
  * the mass carried by the 3-adically small frequencies |t| <= T0 (fixed units), against the uniform expectation.
Usage: python collatz_fourier_mass_real_coordinate_20261004.py [NMAX=13] [A=40]
"""
import sys, math, time
import numpy as np

def law_vector(n, A, prev=None):
    mod = 3 ** n
    if prev is None:
        prev = np.array([1.0])
    z = (3 * np.arange(3 ** (n - 1), dtype=np.int64) + 1) % mod
    inv2 = pow(2, -1, mod)
    cur = np.zeros(mod)
    r = z.copy()
    for a in range(1, A + 1):
        r = (r * inv2) % mod
        cur += np.bincount(r, weights=prev * (2.0 ** (-a)), minlength=mod)
    return cur

def is_pm_power_of_two(t, mod):
    """return s if t = +-2^s mod 3^n for some 0 <= s < L, else None (brute force via a dict built once per n)."""
    return POW2.get(t % mod)

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    prev = None
    for n in range(1, NMAX + 1):
        mu = law_vector(n, A, prev); prev = mu
        if n < 8:
            continue
        mod = 3 ** n; L = 2 * 3 ** (n - 1)
        t0 = time.time()
        F = np.fft.fft(mu)                        # F[t] = sum_y mu(y) e(-t y / 3^n); |F[t]| = |mu_hat_n(t)|
        t = np.arange(mod)
        units = (t % 3) != 0
        E = np.abs(F) ** 2
        mass_units = float(E[units].sum())
        theta = t / mod
        # 64 bins of theta over the units
        bins = np.minimum((theta * 64).astype(int), 63)
        massbin = np.bincount(bins[units], weights=E[units], minlength=64) / mass_units
        P(f"n={n}: mass over units {mass_units:.4f}; top theta-bins (of 64, uniform = 0.0156 each): " +
          ", ".join(f"[{b/64:.3f},{(b+1)/64:.3f}):{massbin[b]:.4f}" for b in np.argsort(-massbin)[:8]))
        # powers of two table
        POW2 = {}
        r = 1
        for s in range(L):
            POW2.setdefault(r, s); POW2.setdefault((mod - r) % mod, -s if s else 0)  # sign-marked: negative s means -2^s (s=0: both 1 and -1)
            r = (2 * r) % mod
        # top coefficients
        idx = np.argsort(-E * units)[:30]
        lines = []
        for tt in idx:
            s = POW2.get(int(tt))
            th = theta[tt]
            thb = format(int(th * 2 ** 12), '012b')
            tag = f"{'+' if s is None or s >= 0 else '-'}2^{abs(s)}" if s is not None else "   --  "
            lines.append(f"t={int(tt):>10d} theta={th:.6f} (bin {thb[:8]}..) |mu_hat|={math.sqrt(E[tt]):.5f} {tag}")
        P("   top coefficients:")
        for l in lines[:16]:
            P("     " + l)
        # dyadic neighbourhoods
        for m in (3, 4, 5, 6, 7, 8):
            near = np.zeros(mod, dtype=bool)
            for j in range(0, 2 ** m):
                c = j / 2 ** m
                d = np.abs(((theta - c + 0.5) % 1.0) - 0.5)
                near |= d <= 2.0 ** (-m) / 4      # within a quarter spacing of the dyadic rational j/2^m
            frac_mass = float(E[units & near].sum()) / mass_units
            frac_haar = float((units & near).sum()) / float(units.sum())
            P(f"   m={m}: mass within 2^-{m}/4 of a dyadic rational j/2^{m}: {frac_mass:.4f} (uniform expectation {frac_haar:.4f}, ratio {frac_mass/frac_haar:.2f})")
        # 3-adically small frequencies
        for T0 in (100, 1000, 10000):
            sel = units & ((t <= T0) | (t >= mod - T0))
            fm = float(E[sel].sum()) / mass_units; fh = float(sel.sum()) / float(units.sum())
            P(f"   fixed units |t| <= {T0}: mass fraction {fm:.3e} (uniform {fh:.3e}, ratio {fm/fh:.3f})")
        # exponent coordinate: mass by exponent window (|s| small vs s near n log2 3)
        P(f"   [level {n} done in {time.time()-t0:.1f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
