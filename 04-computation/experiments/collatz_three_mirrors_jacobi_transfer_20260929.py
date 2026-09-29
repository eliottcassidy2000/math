#!/usr/bin/env python3
"""The Jacobi-sum transfer of the multiplicative spectrum of the Syracuse law (S23, obligation (a)).

Y_n = 2^-a (3 Y_(n-1) + 1) mod 3^n, so for a character psi of (Z/3^n)^x
    E[psi(Y_n)] = G_psi * E[psi(3 Y_(n-1) + 1)],   G_psi = sum_(a>=1) 2^-a psi(2)^-a = (psi(2)^-1/2)/(1 - psi(2)^-1/2),
and y -> psi(3y + 1) (a function of y mod 3^(n-1), y a unit since Y_(n-1) is) expands in the characters psi' of
(Z/3^(n-1))^x with the Jacobi-type coefficients
    c(psi, psi') = (1/L_(n-1)) sum_(y unit mod 3^(n-1)) psi(3y + 1) conj(psi'(y)),
so E[psi(Y_n)] = G_psi sum_(psi') c(psi, psi') E[psi'(Y_(n-1))]: an exact linear recursion on the spectrum.
This script computes the law by forward DP for n <= NMAX, the spectra S_n, checks the recursion exactly, and reports
the structure of |c(psi, psi')|: the row sums of |c|^2 (Parseval: = mean |psi(3y+1)|^2 = 1), the concentration
(largest |c| per row and how many psi' carry 90% of the row's mass), and whether the dominant psi' of a given psi
is related to psi by a simple rule (e.g. psi' = psi restricted / psi^2 / the 3-adic-log characters).
Run: python 04-computation/experiments/collatz_three_mirrors_jacobi_transfer_20260929.py [NMAX]
"""
from __future__ import annotations

import math
import sys

import numpy as np


def law(n: int) -> np.ndarray:
    mod = 3 ** n
    mu = np.zeros(mod)
    mu[0] = 1.0
    inv2 = [pow(2, -a, mod) for a in range(61)]
    ys = np.arange(mod)
    src = (3 * ys + 1) % mod
    for _ in range(n):
        new = np.zeros(mod)
        for a in range(1, 61):
            np.add.at(new, (src * inv2[a]) % mod, 2.0 ** (-a) * mu)
        mu = new
    return mu


def dlog_table(n: int) -> np.ndarray:
    """discrete log base 2 of the units mod 3^n: dlog[y] = k with 2^k = y, -1 for non-units."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    d = -np.ones(mod, dtype=np.int64)
    x = 1
    for k in range(L):
        d[x] = k
        x = (2 * x) % mod
    return d


if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for n in range(2, NMAX + 1):
        mod, modp = 3 ** n, 3 ** (n - 1)
        L, Lp = 2 * 3 ** (n - 1), 2 * 3 ** (n - 2)
        mu_n, mu_p = law(n), law(n - 1)
        dl, dlp = dlog_table(n), dlog_table(n - 1)
        units = np.nonzero(dl >= 0)[0]
        unitsp = np.nonzero(dlp >= 0)[0]
        # spectra: S_n(psi_j) = sum_y mu_n(y) e(j dlog(y)/L)
        S_n = np.array([np.sum(mu_n[units] * np.exp(2j * np.pi * j * dl[units] / L)) for j in range(L)])
        S_p = np.array([np.sum(mu_p[unitsp] * np.exp(2j * np.pi * j * dlp[unitsp] / Lp)) for j in range(Lp)])
        # transfer coefficients c(psi_j, psi'_j') = (1/Lp) sum_(y unit mod 3^(n-1)) psi_j(3y+1) conj psi'_j'(y)
        y = unitsp                                   # units mod 3^(n-1)
        w = (3 * y + 1) % mod                        # 1-units mod 3^n
        dl_w = dl[w]                                 # their discrete logs mod L
        C = np.zeros((L, Lp), dtype=complex)
        for j in range(L):
            phase = np.exp(2j * np.pi * j * dl_w / L)
            for jp in range(Lp):
                C[j, jp] = np.mean(phase * np.exp(-2j * np.pi * jp * dlp[y] / Lp))
        # G_psi
        G = np.array([(np.exp(-2j * np.pi * j / L) / 2) / (1 - np.exp(-2j * np.pi * j / L) / 2) for j in range(L)])
        pred = G * (C @ S_p)
        err = np.max(np.abs(pred - S_n))
        rowsum2 = np.sum(np.abs(C) ** 2, axis=1)
        conc = []
        for j in range(L):
            a = np.sort(np.abs(C[j]) ** 2)[::-1]
            cs = np.cumsum(a) / a.sum()
            conc.append(int(np.searchsorted(cs, 0.9) + 1))
        prim = np.array([j % 3 != 0 for j in range(L)])
        print(f"n={n}: L={L}, L'={Lp}; recursion check max |pred - S_n| = {err:.2e}; row sums of |c|^2: min {rowsum2.min():.4f} max {rowsum2.max():.4f} (Parseval 1)")
        print(f"   concentration: number of psi' carrying 90% of a row's mass: min {min(conc)}, median {int(np.median(conc))}, max {max(conc)} of {Lp}; largest |c| per row: mean {np.mean(np.max(np.abs(C), axis=1)):.3f}")
        # which psi' dominates for the 3-adic-log characters psi_(+-2) and for a generic primitive psi
        for j in [2 % L, (L - 2) % L, 1, L // 2 + 1] if L > 4 else [1]:
            top = np.argsort(-np.abs(C[j]))[:4]
            print(f"   psi_{j}: |S_n| 3^(n/2) = {abs(S_n[j]) * 3 ** (n / 2):.3f}; dominant psi'_j': " + ", ".join(f"j'={jp} (|c|={abs(C[j, jp]):.3f}, |S'| 3^((n-1)/2)={abs(S_p[jp]) * 3 ** ((n - 1) / 2):.3f})" for jp in top))
        # is the dominant j' a fixed function of j? test the rule j' = j mod L' (restriction) and j' = j/3 for 3 | j
        # exact structure: for primitive psi (3 does not divide j) and primitive psi' (3 does not divide j'), |c| = (#prim')^(-1/2); imprimitive columns vanish
        primp = np.array([jp % 3 != 0 for jp in range(Lp)]) if n >= 3 else np.array([jp != 0 for jp in range(Lp)])
        nprimp = int(primp.sum())
        sub_pp = np.abs(C[np.ix_(prim, primp)]) * math.sqrt(nprimp)
        sub_pi = np.abs(C[np.ix_(prim, ~primp)])
        Gabs = np.abs(G)
        print(f"   exact structure: over primitive (psi, psi') pairs, |c| sqrt(#prim') in [{sub_pp.min():.6f}, {sub_pp.max():.6f}] (1 if flat); max |c| on imprimitive psi' columns {sub_pi.max():.2e}; "
              f"|G_psi| in [{Gabs.min():.4f}, {Gabs.max():.4f}], rms over primitive psi {math.sqrt(np.mean(Gabs[prim] ** 2)):.4f} (1/sqrt 3 = 0.5774); |G| at j = +-2: {Gabs[2 % L]:.4f}")
        # prediction of the spectrum's size from the recursion: |S_n(psi)| ~ |G_psi| * (flat mix of S_(n-1)); check corr(|S_n| , |G|) over primitive psi
        corr = np.corrcoef(np.abs(S_n[prim]), Gabs[prim])[0, 1]
        print(f"   correlation of |S_n(psi)| with |G_psi| over primitive psi: {corr:.3f}; mean |S_n| 3^(n/2) for |G| > 0.9: {np.mean(np.abs(S_n[prim & (Gabs > 0.9)]) * 3 ** (n / 2)) if (prim & (Gabs > 0.9)).any() else float('nan'):.3f}, for |G| < 0.45: {np.mean(np.abs(S_n[prim & (Gabs < 0.45)]) * 3 ** (n / 2)) if (prim & (Gabs < 0.45)).any() else float('nan'):.3f}")
        rule_hits = sum(1 for j in range(L) if np.argmax(np.abs(C[j])) == (j % Lp))
        print(f"   rule j' = j mod L' is the argmax for {rule_hits} of {L} rows")
    print("DONE")
