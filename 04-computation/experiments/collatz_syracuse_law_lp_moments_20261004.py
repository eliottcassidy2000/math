#!/usr/bin/env python3
"""
Experiment E2 (opus, 2026-10-04): L^p criticality of the 3-adic Syracuse law (prediction P7 of the coalescence note).

The law mu_n of Y_n on Z/3^n (Y_0 = 0, Y_n = 2^-a (3 Y_{n-1} + 1) mod 3^n, a geometric of mean 2) is computed as a
dense vector by the exact recursion (valuations a <= A, dropped mass < n 2^-A).  The reference density on the units is
rho_n(y) = (2/3) 3^n mu_n(y) (S19), of mean 1 over the units.  We report, per level n:
    E[rho_n^p] over the units for p in {1.25, 1.5, 1.75, 1.9, 2, 2.1, 2.25, 2.5},
    the tail P(rho_n > t) for t = 2, 4, 8, 16, 32, 64, 128 and the local tail exponents,
    the maximal atom and its location (expected: the class of -1, S19), and the atoms at -1, -5, -7, -17, 1.
P7 predicts: E[rho^p] bounded in n for p < 2, linear growth at p = 2 (S19: 0.31 per level), exponential growth for
p > 2; equivalently P(rho > t) ~ c t^-2.
Usage: python collatz_syracuse_law_lp_moments_20261004.py [NMAX=14] [A=40]
"""
import sys, math, time
import numpy as np

def law_vector(n, A, prev=None):
    """mu_n as a float64 vector of length 3^n, from mu_{n-1} (vector of length 3^(n-1)) or from scratch."""
    mod = 3 ** n
    if prev is None:
        prev = np.array([1.0])
    # lift: a residue y mod 3^(n-1) splits into y, y + 3^(n-1), y + 2 3^(n-1); the recursion Y_n = 2^-a (3 Y_{n-1} + 1)
    # needs Y_{n-1} as a residue mod 3^n?  No: 3 Y_{n-1} + 1 mod 3^n depends only on Y_{n-1} mod 3^(n-1).  So:
    z = (3 * np.arange(3 ** (n - 1), dtype=np.int64) + 1) % mod      # z = 3y+1 for y mod 3^(n-1)
    inv2 = pow(2, -1, mod)
    cur = np.zeros(mod)
    r = z.copy()
    for a in range(1, A + 1):
        r = (r * inv2) % mod
        cur += np.bincount(r, weights=prev * (2.0 ** (-a)), minlength=mod)
    return cur

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    ps = [1.25, 1.5, 1.75, 1.9, 2.0, 2.1, 2.25, 2.5]
    ts = [2, 4, 8, 16, 32, 64, 128, 256]
    P(f"Syracuse law L^p moments, A = {A}, levels 1..{NMAX}")
    P("n | mass | " + " ".join(f"E[rho^{p}]" for p in ps) + " | max rho (at) | rho(-1) rho(-5) rho(-7) rho(-17) rho(1)")
    prev = None
    prevE2 = None
    tails = {}
    for n in range(1, NMAX + 1):
        t0 = time.time()
        mu = law_vector(n, A, prev)
        prev = mu
        mod = 3 ** n
        units = np.arange(mod) % 3 != 0
        rho = (2.0 / 3.0) * mod * mu[units]
        mass = mu.sum()
        moms = [float(np.mean(rho ** p)) for p in ps]
        imax = int(np.argmax(mu));
        atoms = [float((2.0 / 3.0) * mod * mu[x % mod]) for x in (-1, -5, -7, -17, 1)]
        tail = [float(np.mean(rho > t)) for t in ts]
        tails[n] = tail
        # THM-4263's uniform-integrability quantity (14): |Y_n|^-1 sum_{rho_n(y) > M} rho_n(y) = E[rho 1_{rho > M}]
        ui = [float(np.mean(rho * (rho > M))) for M in (4, 8, 16, 32, 64, 128)]
        P(f"     UI tails E[rho 1(rho > M)] for M = 4, 8, 16, 32, 64, 128: " + " ".join(f"{x:.4f}" for x in ui) + "   (t^-2 tail with c: 2c/M)")
        P(f"{n:2d} | {mass:.6f} | " + " ".join(f"{m:.4f}" for m in moms) + f" | {(2/3)*mod*mu[imax]:.3f} (y={imax}, = {imax - mod if imax > mod//2 else imax}) | " + " ".join(f"{a:.3f}" for a in atoms) + f"   [{time.time()-t0:.1f}s]")
        if prevE2 is not None:
            P(f"     increment of E[rho^2]: {moms[4] - prevE2:+.4f}")
        prevE2 = moms[4]
    P("tails P(rho_n > t):")
    P("n | " + " ".join(f"t={t}" for t in ts))
    for n in range(4, NMAX + 1):
        P(f"{n:2d} | " + " ".join(f"{x:.2e}" for x in tails[n]))
    P("local tail exponents -log2(P(rho>2t)/P(rho>t)) at the top level (prediction: 2):")
    tl = tails[NMAX]
    P("  " + " ".join(f"t={ts[i]}: {(-math.log2(tl[i+1]/tl[i]) if tl[i+1] > 0 and tl[i] > 0 else float('nan')):.2f}" for i in range(len(ts) - 1)))
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
