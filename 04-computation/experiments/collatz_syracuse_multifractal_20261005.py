#!/usr/bin/env python3
"""
Multifractal spectrum of the Syracuse law (opus, 2026-10-05).

The Syracuse law mu_n on Z/3^n is the law of the n-th accelerated odd iterate Y_n = 2^-a_n (3 Y_(n-1) + 1) mod 3^n
for a Haar-random odd 2-adic start (valuations a_i iid geometric(1/2)).  Exact recursion on residues:
    mu_n(c) = sum_{a >= 1, 2^a c = 1 mod 3} 2^-a mu_(n-1)( (2^a c - 1)/3 mod 3^(n-1) ),   mu_0 = delta.
Partition function Z_n(q) = sum_c mu_n(c)^q ~ 3^(-n tau(q)); generalized dimension D_q = tau(q)/(q-1); the density
rho_n = (2/3) 3^n mu_n on 3-units has E[rho_n^q] = 3^(n((q-1) - tau(q))) up to constants.

CONJECTURED SPECTRUM (derived in the note from the word-energy heuristic: in the concentrated regime each residue is
carried by few parity words, so Z_n(q) ~ sum_words 2^(-q S_n) = (2^q - 1)^-n; in the diffuse regime Z_n(q) ~ 3^(-n(q-1))):
    tau(q) = q - 1                 for q <= 2   (D_q = 1),
    tau(q) = log_3(2^q - 1)        for q >= 2   (D_q = log_3(2^q - 1)/(q - 1), D_inf = log_3 2),
with the transition exactly at q = 2 because sum_a 4^-a = 1/3 = 3^-(2-1).  Legendre: alpha(q) = tau'(q) =
2^q ln 2 / ((2^q - 1) ln 3), f = q alpha - tau; at q = 2: alpha = 4/(3 log_2 3) = 1/(p log_2 3) with p = 3/4, and
2 - H(3/4) = (3/4) log_2 3 exactly (the critical valuation-one frequency 3/4 = the energy-weighted step law 3 4^-a).
This script evaluates a valuation-truncated exact recursion in float64 to n = N.
The partition functions and slopes are numerical, not exact rational certificates.
Usage: python collatz_syracuse_multifractal_20261005.py [N=13] [A=60]
"""
import sys, math, time
import numpy as np

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    qs = [0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0, 4.0, 6.0, 8.0]
    def tau_conj(q):
        return (q - 1) if q <= 2 else math.log(2 ** q - 1) / math.log(3)
    P(f"Syracuse law mu_n on Z/3^n, exact recursion (valuations <= {A}), n <= {N}; conjecture tau(q) = q-1 (q<=2), log_3(2^q-1) (q>=2)")
    mu = np.array([1.0])            # n = 0
    Z_prev = None
    rows = []
    for n in range(1, N + 1):
        M = 3 ** n; Mp = 3 ** (n - 1)
        c = np.arange(M, dtype=np.int64)
        new = np.zeros(M)
        for a in range(1, A + 1):
            pa = pow(2, a, 3)
            # need 2^a c = 1 mod 3  <=>  c = pa^-1 mod 3
            cres = (pow(pa, -1, 3)) % 3
            sel = c[c % 3 == cres]
            y = ((pow(2, a, 3 * M) * sel - 1) // 3) % Mp          # exact: 2^a sel - 1 divisible by 3
            new[sel] += (2.0 ** -a) * mu[y]
        mu = new
        tot = mu.sum()
        Z = {q: float(np.sum(mu[mu > 0] ** q)) for q in qs}
        amax = mu.max(); cmax = int(np.argmax(mu))
        line = f"n={n:2d} mass {tot:.12f} max atom {amax:.6e} at c={cmax} (c=-1: {cmax == M - 1}), 2^-n ratio {amax * 2.0 ** n:.4f}"
        if Z_prev:
            slopes = {q: -math.log(Z[q] / Z_prev[q]) / math.log(3) for q in qs}
            line += " | tau slopes: " + " ".join(f"q{q}:{slopes[q]:.4f}" for q in qs)
            rows.append((n, slopes))
        P(line)
        Z_prev = Z
    P("")
    P("per-level tau estimates at the last three levels vs the conjecture (and D_q = tau/(q-1)):")
    P("  q      conj tau   " + " ".join(f"n={n:2d}      " for n, _ in rows[-3:]) + " conj D_q  ")
    for q in qs:
        vals = " ".join(f"{s[q]:.4f}    " for _, s in rows[-3:])
        Dq = tau_conj(q) / (q - 1) if q != 1 else 1.0
        P(f"  {q:<5} {tau_conj(q):.4f}     {vals} {Dq:.4f}")
    P("")
    P(f"max atom at c = -1: mu_n(-1) 2^n -> {rows and amax * 2.0 ** N:.4f} (the L^inf local dimension log_3 2 = {math.log(2)/math.log(3):.4f}); D_inf conjectured = log_3 2")
    H = lambda p: -p * math.log2(p) - (1 - p) * math.log2(1 - p)
    a = math.log2(3)
    P(f"exact identities: 2 - H(3/4) = {2 - H(0.75):.10f} = (3/4) log_2 3 = {0.75 * a:.10f}; alpha(2) = 4/(3 log_2 3) = {4 / (3 * a):.6f}; f(alpha(2)) = H(3/4)/((3/4) log_2 3) = {H(0.75) / (0.75 * a):.6f}; h* = H(1/log_2 3) = {H(1 / a):.6f}")
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
