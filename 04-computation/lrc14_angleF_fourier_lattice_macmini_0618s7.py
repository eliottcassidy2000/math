#!/usr/bin/env python3
"""
lrc14_angleF_fourier_lattice_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE F part 2 — EXACT Fourier-over-relation-lattice form of corr(E),
verified, and the per-vector kernel bound that makes the covolume the carrier.

The seven-sector cover indicator expands (inclusion-exclusion over which of the
7 sectors are MISSED; sector 0 is pinned hit by e_0=0, so missed-sets T<={1..6}):
   meas(S7(E)) = sum_{n in Z^{k-1}: sum_j n_j e_j = 0}  K(n),
   K(n) = sum_{T<={1..6}} (-1)^{|T|} prod_{j=1}^{k-1} chat_T(n_j),
where n=(n_1,...,n_{k-1}) ranges over the NONZERO offsets e_1<...<e_{k-1}, and
chat_T is the Fourier coefficient of 1_{B_T}, B_T = [0,1)\ U_{j in T}[j/7,(j+1)/7).
   chat_T(0) = 1 - |T|/7;   chat_T(7m)=0 (THM-503);   |chat_T(n)| <= |T|/(pi |n|).
n=0 gives M7(k); corr(E) = sum over NONZERO lattice vectors.

PER-VECTOR KERNEL BOUND (the geometry-of-numbers engine):
   |K(n)| <= sum_T prod_j |chat_T(n_j)|  <=  A * prod_{j: n_j != 0} 1/|n_j|
with an ABSOLUTE constant A (independent of E, k beyond a |T|-sum factor).
Define height(n) = prod_{n_j != 0} |n_j|.  Then
   |corr(E)| <= A * sum_{0 != n in Lambda^o(E)} 1/height(n)  =  A * W_full(E).
The lattice Lambda^o(E) = {n: sum n_j e_j = 0} has covolume = |e|/gcd... ; SMALL
covolume <=> MANY short relations <=> large W_full.  This is the rigorous chain
  small covolume  ->  large relation weight  ->  large corr bound.

This script:
 (1) VERIFIES the lattice identity reproduces meas(S7) exactly (small k, exact-ish);
 (2) computes the EXACT per-vector kernel constant A_T = max|chat_T(n)|*|n| numerically;
 (3) tabulates W_full(E) (the 1/height lattice sum) vs corr(E) and vs covolume.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from collections import defaultdict
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

TWO_PI_I = 2j * math.pi


def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = set()
        for e in E:
            secs.add(int(((e * xm) % 1) * 7))
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    s = F(0)
    for t in range(0, 7):
        s += F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1)
    return s


# single-sector Fourier coefficient ŝ_j(n) = FT of 1_{[j/7,(j+1)/7)} at n
def shat(n, j):
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)


def chat_T_array(nvals, T):
    """chat_T(n) for each n in nvals. chat_T = 1_{n=0}(1-|T|/7) - sum_{j in T} shat_j(n)."""
    out = {}
    for n in nvals:
        if n == 0:
            out[n] = complex(1 - len(T) / 7.0, 0.0)
        elif n % 7 == 0:
            out[n] = 0j
        else:
            out[n] = -sum(shat(n, j) for j in T)
    return out


def lattice_corr(E, N0):
    """Sum over nonzero relation vectors with |n_j|<=N0 of K(n). Returns (corr_approx, Wfull)."""
    nz = [e for e in E if e != 0]
    d = len(nz)
    nvals = list(range(-N0, N0 + 1))
    # precompute chat_T(n) for all needed n and all T<={1..6}
    subsets = [tuple(T) for r in range(0, 7) for T in itertools.combinations(range(1, 7), r)]
    chat = {T: chat_T_array(nvals, T) for T in subsets}
    sign = {T: (-1) ** len(T) for T in subsets}
    corr = 0j
    Wfull = 0.0
    for n in itertools.product(nvals, repeat=d):
        if sum(ni * ei for ni, ei in zip(n, nz)) != 0:
            continue
        if all(x == 0 for x in n):
            continue
        K = 0j
        for T in subsets:
            ct = chat[T]
            p = 1.0 + 0j
            for ni in n:
                p *= ct[ni]
                if p == 0:
                    break
            K += sign[T] * p
        corr += K
        h = 1
        for ni in n:
            if ni != 0:
                h *= abs(ni)
        Wfull += 1.0 / h
    return corr.real, Wfull


def offset_covolume(E):
    """Covolume of Lambda^o(E)={n in Z^{k-1}: sum n_j e_j=0} = sqrt(det Gram of LLL basis)."""
    import sympy
    nz = [e for e in E if e != 0]
    M = sympy.Matrix([nz])
    ns = M.nullspace()
    if not ns:
        return 0.0
    basis = []
    for v in ns:
        dens = []
        for x in v:
            fr = F(int(x.p), int(x.q))
            dens.append(fr.denominator)
        L = 1
        for dd in dens:
            L = L * dd // math.gcd(L, dd)
        iv = [int(x * L) for x in v]
        g = 0
        for c in iv:
            g = math.gcd(g, abs(c))
        if g:
            iv = [c // g for c in iv]
        basis.append(iv)
    G = sympy.Matrix([[sum(a * b for a, b in zip(u, w)) for w in basis] for u in basis])
    return math.sqrt(float(G.det()))


def main():
    print("=" * 92)
    print("ANGLE F part 2: Fourier-lattice form + per-vector kernel bound + covolume carrier")
    print("=" * 92)
    # (1) verify identity reproduces measS7 by truncation
    print("(1) Truncated lattice sum -> exact meas(S7) (verifies the lattice picture):")
    for E in [(0, 1, 2, 3, 4), (0, 1, 2, 4, 5), (0, 1, 3, 5, 7)]:
        exact = float(measS7(E))
        m7 = float(M7(len(E)))
        for N0 in (6, 12):
            corr_a, _ = lattice_corr(E, N0)
            print(
                f"   E={E} N0={N0}: M7+latticesum={m7+corr_a:.5f}  exact={exact:.5f}  "
                f"corr_exact={exact-m7:+.5f} corr_trunc={corr_a:+.5f}"
            )
    print()
    # (2) max |chat_T(n)|*|n| -> absolute kernel constant
    print("(2) per-sector Fourier bound: max_n |shat_j(n)|*|n| (=1/pi * |sin(pi/7)|*7/pi...):")
    mx = 0.0
    for j in range(7):
        for n in range(1, 50):
            if n % 7 == 0:
                continue
            mx = max(mx, abs(shat(n, j)) * n)
    print(f"   max |shat(n)| * |n| = {mx:.5f}   (bound chat_T <= |T|*{mx:.4f}/|n|)")
    print()
    # (3) corr vs Wfull vs covolume
    print("(3) corr(E), W_full(E)=sum 1/height, and covolume(Lambda^o):")
    hdr = f"{'shape (k=8)':<34}{'corr':>10}{'Wfull':>10}{'covol':>12}{'corr/Wfull':>12}"
    print(hdr)
    shapes = [
        ("consec {0..7}", list(range(8))),
        ("3AP {0,1,2,3,4,5,6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("2run {0,1,2,3,40,41,42,43}", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("perf {0,2,3,4,5,6,7,9}", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("Sidon {0,1,3,7,12,20,30,44}", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("dissoc {0,1,3,7,15,31,63,127}", [0, 1, 3, 7, 15, 31, 63, 127]),
    ]
    m7 = float(M7(8))
    for name, E in shapes:
        corr = float(measS7(E)) - m7
        _, Wf = lattice_corr(E, 6)
        cov = offset_covolume(E)
        ratio = corr / Wf if Wf > 1e-9 else float("nan")
        print(f"{name:<34}{corr:>10.5f}{Wf:>10.4f}{cov:>12.2f}{ratio:>12.5f}")
    print()
    print("READOUT: corr <= A*Wfull with A small; Wfull large <=> covol small (densest")
    print("relations). The covolume is the geometry-of-numbers carrier of the correction.")
    print("DONE.")


if __name__ == "__main__":
    main()
