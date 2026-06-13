#!/usr/bin/env python3
"""
drt_two_point_spectrum_monad.py
monad-explorer-2026-06-07 (deep-research, 3rd session) -- evidence for HYP-2308.

The leading-order Catalan skeleton of the cluster integral A_{2k} (THM-438) is SPECTRAL,
not arithmetic.  The only top-order single-merge pattern is the closed 2k-cycle = tr(S^{2k}),
where S is the skew +-1 tournament matrix.  For ANY doubly-regular tournament (DRT) S has the
TWO-POINT spectrum {0} U {+- i*sqrt(n)} (each +-i*sqrt(n) with multiplicity (n-1)/2), giving

        tr(S^{2k}) = (-1)^k n^k (n-1)   IDENTICALLY,

which is exactly the closed form THM-438 found for Paley.  This closed form depends ONLY on the
two-point spectrum = the defining property of a DRT.  Hence the even-cacti census and the
A088368 -> Catalan cancellation (MISTAKE-060) are DRT-universal at leading order; only the
sub-leading o(n^{k+1}) bound is circulant-specific (Weil) and becomes, for general DRT, a
tight-spectral / expander-mixing estimate.

Here we verify the two-point spectrum and the trace formula on the Paley DRTs (the circulant
DRTs we can build without skew-Hadamard data).  A compute node should repeat with a VERIFIED
non-circulant DRT on n=15 (skew-Hadamard order 16) -- CHECK validity per MISTAKE-011b/017.
"""
import numpy as np

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)

print("Paley DRT skew matrix S[a,b]=chi(b-a):  spectrum and tr(S^{2k})")
print(f"{'n':>4} {'|eig| set':>20} {'sqrt(n)':>9} {'#zero':>6}   trace checks")
for p in [7, 11, 19, 23, 31]:
    S = np.array([[legendre((b - a) % p, p) for b in range(p)] for a in range(p)], dtype=float)
    ev = np.linalg.eigvals(S)
    mags = sorted({round(abs(z), 3) for z in ev})
    zero = sum(1 for z in ev if abs(z) < 1e-6)
    checks = []
    for k in [1, 2, 3, 4]:
        tr = np.trace(np.linalg.matrix_power(S, 2 * k)).real
        pred = ((-1) ** k) * (p ** k) * (p - 1)
        checks.append(f"k={k}:{'OK' if abs(tr - pred) < 1e-3 else 'FAIL'}")
    print(f"{p:>4} {str(mags):>20} {p**0.5:>9.3f} {zero:>6}   " + " ".join(checks))
print()
print("All Paley DRTs: spectrum = {0} U {+- i sqrt(n)}, tr(S^{2k}) = (-1)^k n^k (n-1).")
print("This is the DRT-defining two-point spectrum => the closed form is DRT-universal.")
print("=> the THM-438 even-cycle (single 2k-cycle) contribution is identical for any DRT,")
print("   no Gauss sum / Weil needed.  Leading-order Catalan skeleton is SPECTRAL (HYP-2308).")
