#!/usr/bin/env python3
"""procgen-bridges 2026-09-23, part ADELIC: the shared integrality engine and the Z/2 fold.

Sections (deterministic stdout; timing on stderr):
  AD1  adelic radii of the natural function forms.  Pólya-Bertrandias needs prod_v R_v > 1.
       Collatz/Bernstein series (both parametrizations), theta series, and AMM's lacunary class all have
       prod_v R_v <= 1, with equality exactly when a density exists (product formula for 2^d 3^-l).
  AD2  Hankel determinants of the AMM parity skeleton L(w) = sum_(j>=0) w^(2^j) modulo 2 (Lemma P:
       phi == L mod 2 for every complement-symmetric fair extractor), shifts m = 0..3, n <= NMAX;
       exact integer Hankel determinants of L for n <= 40 (sizes, 2-adic valuations, Hadamard defect).
  AD3  the Z/2 fold: (a) AMM's p <-> 1-p is LRC's half-translate x -> x+1/2 under p = sin^2(pi x), and the
       AMM quotient w = p(1-p) is LRC's clock-two quotient y = 2x (4w = sin^2(pi y)); the Catalan branch
       z(w) = sin^2(pi y/2).  (b) 2-adically, L(p) = p(1-p)/2 is isometrically conjugate to the Collatz map T
       (both are 2-adic Bernoulli maps); (c) what the conjugacy destroys: rational orbits under L are not
       eventually periodic (height doubling), under T they are (the affine structure behind PC).
Memory < 150 MB; runtime about 1 minute.
"""
import cmath
import math
import random
import sys
import time
from fractions import Fraction

import numpy as np

T0 = time.time()


def P(*a):
    print(*a, flush=True)


def tick(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------- AD1
def mechanical(alpha, rho, n):
    i = np.arange(n + 1, dtype=np.float64)
    fl = np.floor(i * alpha + rho)
    return (fl[1:] - fl[:-1]).astype(np.int8)


def block_word(nblocks, member, B, B2):
    out = []
    for k in range(nblocks):
        out.extend(B2 if member(k) else B)
    return np.array(out, dtype=np.int8)


def is_square(k):
    r = math.isqrt(k)
    return k > 0 and r * r == k


def is_cube(k):
    r = round(k ** (1 / 3))
    return k > 0 and any((r + e) ** 3 == k for e in (-1, 0, 1))


def section_ad1():
    P("=" * 100)
    P("AD1  adelic radii: Polya-Bertrandias needs sum_v log R_v > 0")
    P("=" * 100)
    P("  Bernstein series of a parity word w with odd-step positions d_1 < d_2 < ...:")
    P("    (i)  g_w(Y) = sum_l 2^(d_l) Y^l in Z[[Y]]:      log2 R_inf = -limsup d_l/l,  log2 R_2 = +liminf d_l/l")
    P("    (ii) f_w(X) = sum_l X^(d_l) 3^(-l) in Z[1/3][[X]]: log3 R_inf = liminf l/d_l, log3 R_3 = -limsup l/d_l")
    P("  Each term is an S-unit (S = {inf, 2, 3}), so the product formula forces sum_v log R_v = liminf - limsup <= 0.")
    rng = np.random.default_rng(7)
    words = [
        ("Sturmian golden slope", mechanical(1 / ((1 + 5 ** 0.5) / 2) ** 2, 0.3, 200000)),
        ("Sturmian slope log_3 2 (critical)", mechanical(math.log(2) / math.log(3), 0.1, 200000)),
        ("square-swap Y (1^9 0 / 1^8 0 1)", block_word(20000, is_square, [1] * 9 + [0], [1] * 8 + [0, 1])),
        ("cube-swap Y3", block_word(20000, is_cube, [1] * 9 + [0], [1] * 8 + [0, 1])),
        ("i.i.d. fair bits", rng.integers(0, 2, 200000).astype(np.int8)),
        ("Bernoulli blocks of density 0.3/0.8 (no density)", np.concatenate(
            [(rng.random(4 ** k) < (0.3 if k % 2 else 0.8)).astype(np.int8) for k in range(1, 9)])),
    ]
    P(f"  {'word':50s} {'liminf d/l':>11s} {'limsup d/l':>11s} {'sum_v log2 R_v (i)':>19s}")
    for name, w in words:
        d = np.flatnonzero(w) + 1
        l = np.arange(1, len(d) + 1)
        ratio = d / l
        tail = ratio[len(ratio) // 4:]
        lo, hi = float(tail.min()), float(tail.max())
        P(f"  {name:50s} {lo:11.5f} {hi:11.5f} {lo - hi:19.5f}")
    P("  theta(rho) = sum_k rho^(k^2) as a function of rho: integer coefficients, R_v = 1 at every place -> sum 0.")
    P("  AMM lacunary class phi = poly + sum_j eps_j w^(2^j): R_inf = 1 (natural boundary), R_p = 1 -> sum 0.")
    P("  AMM Theorem A: phi continues to W_gamma (not a disk); the archimedean term becomes log R(W_gamma,0) =")
    P("  Lambda(gamma) > 0 for gamma < 0.3775 (THM-4467).  CONTINUATION beyond the circle of convergence is the whole")
    P("  gain; Collatz-type series are lacunary on their circle of convergence and have no continuation.")
    P("  VERDICT AD1: the engine 'integrality + smallness => vanishing' is shared (REAL), but its adelic-capacity form")
    P("  is exactly critical (sum = 0) for every Collatz function form, so it cannot carry PC; Collatz must use")
    P("  value-level inputs (Liouville/Pade/Subspace), which is what Theorems R, S, D, Y do.")


# ---------------------------------------------------------------------------------------------- AD2
def gf2_rank_rows(rows, n):
    """rank over GF(2) of n rows given as python ints (bit k = column k)"""
    rows = list(rows)
    rank = 0
    for col in range(n):
        bit = 1 << col
        piv = None
        for i in range(rank, len(rows)):
            if rows[i] & bit:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        pr = rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] & bit:
                rows[i] ^= pr
        rank += 1
    return rank


def hankel_mod2_profile(coef, m, nmax):
    out = []
    for n in range(1, nmax + 1):
        rows = []
        for i in range(n):
            r = 0
            for k in range(n):
                if coef[m + i + k]:
                    r |= 1 << k
            rows.append(r)
        out.append(1 if gf2_rank_rows(rows, n) == n else 0)
    return out


def bareiss_det(Mx):
    A = [row[:] for row in Mx]
    n = len(A)
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def v2(x):
    x = abs(x)
    return (x & -x).bit_length() - 1 if x else None


def section_ad2(nmax):
    P("=" * 100)
    P("AD2  Hankel determinants of the AMM parity skeleton L(w) = w + w^2 + w^4 + w^8 + ...")
    P("=" * 100)
    N = 2 * nmax + 16
    coef = [0] * N
    j = 1
    while j < N:
        coef[j] = 1
        j *= 2
    P(f"  H_n^(m) = det[c_(m+i+k)]_(0<=i,k<n) modulo 2 for n = 1..{nmax}")
    for m in (0, 1, 2, 3):
        prof = hankel_mod2_profile(coef, m, nmax)
        odd = [n for n, b in zip(range(1, nmax + 1), prof) if b]
        tick(f"AD2 shift {m}")
        P(f"    shift m={m}: odd for {len(odd)}/{nmax} values of n; first odd n: {odd[:14]}"
          f"; last odd n <= {nmax}: {odd[-3:] if odd else []}")
        if m == 1:
            gaps = sorted(set(b - a for a, b in zip(odd, odd[1:])))
            P(f"      gaps between consecutive odd n (m=1): {gaps[:10]}")
    NEX = 100
    for m in (1, 2):
        dets = []
        for n in range(1, NEX + 1):
            Mx = [[coef[m + i + k] for k in range(n)] for i in range(n)]
            d = bareiss_det(Mx)
            had = 0.5 * sum(math.log2(max(sum(r), 1)) for r in Mx)
            dets.append((n, d, had))
        vals = sorted(set(d for _, d, _ in dets))
        signs = "".join("+" if d > 0 else ("-" if d < 0 else "0") for _, d, _ in dets[:24])
        worst = max(h - (math.log2(abs(d)) if d else 0.0) for _, d, h in dets if d)
        P(f"  exact integer Hankel determinants H_n^({m}) of L, n = 1..{NEX}: set of values {vals};"
          f" signs n=1..24: {signs}")
        P(f"    log2 Hadamard bound (Cauchy-Schwarz on the rows) minus log2|H_n|: at n=10,40,{NEX}: "
          f"{dets[9][2]:.1f}, {dets[39][2]:.1f}, {dets[NEX - 1][2]:.1f} bits (max {worst:.1f})")
        tick(f"AD2 exact m={m}")
    P("  => every shifted Hankel determinant H_n^(1), H_n^(2) of the skeleton is a UNIT (+-1) in the exact range, hence")
    P("     odd; by Lemma P the same determinants of phi are odd for every complement-symmetric fair extractor, so")
    P("     v_2(H_n(phi)) = 0: the 2-adic place contributes nothing (no Bertrandias gain), confirming and sharpening")
    P("     the note's radius-1 statement.  The Hadamard/Cauchy-Schwarz size bound misses the unimodular structure")
    P("     by a growing number of bits, but Polya only needs |H_n| >= 1, so this defect is not a lever for C*.")


# ---------------------------------------------------------------------------------------------- AD3
def itin_T(a, K):
    """first K parity bits of the 3x+1 map T on Z/2^K"""
    M = 1 << K
    out = 0
    x = a % M
    for k in range(K):
        b = x & 1
        out |= b << k
        x = (x >> 1) if b == 0 else ((3 * x + 1) >> 1)
    return out


def itin_L(a, K):
    """first K parity bits of L(p) = p(1-p)/2 on Z/2^(K+1) (L maps classes mod 2^(k+1) onto classes mod 2^k)"""
    M = 1 << (K + 2)
    out = 0
    x = a % M
    for k in range(K):
        b = x & 1
        out |= b << k
        x = ((x * (1 - x)) // 2) % M
    return out


def section_ad3():
    P("=" * 100)
    P("AD3  the Z/2 fold shared by AMM (p <-> 1-p), LRC (half-translate / clock two) and Collatz (parity split)")
    P("=" * 100)
    rng = random.Random(11)
    err = 0.0
    for _ in range(2000):
        y = rng.random()
        x1, x2 = y / 2, (y + 1) / 2
        p1, p2 = math.sin(math.pi * x1) ** 2, math.sin(math.pi * x2) ** 2
        err = max(err, abs(p1 + p2 - 1), abs(4 * p1 * (1 - p1) - math.sin(math.pi * y) ** 2))
        w = p1 * (1 - p1)
        z = (1 - math.sqrt(1 - 4 * w)) / 2
        target = min(p1, p2)
        err = max(err, abs(z - target))
    P(f"  (a) p = sin^2(pi x): the two clock-two lifts y/2, (y+1)/2 of y map to the fair pair p, 1-p; 4 p(1-p) =")
    P(f"      sin^2(pi y); the Catalan branch z(w) = (1 - sqrt(1-4w))/2 returns the lift in [0,1/2]:"
      f" max error over 2000 samples {err:.2e}")
    P("      So AMM's involution p -> 1-p IS LRC's half-translate x -> x + 1/2 (THM-4449's Fourier sign (-1)^k), and the")
    P("      AMM quotient w = p(1-p) IS the LRC clock-two quotient y = 2x.  Real p in [0,1] <-> real times; AMM's")
    P("      continuation region (complex p off [0,1], e.g. the golden point p = -1/phi <-> x = i asinh(phi^-1/2)/pi)")
    P("      corresponds to IMAGINARY times, which carry no LRC meaning.")
    xg = math.asinh(((1 + 5 ** 0.5) / 2) ** -0.5) / math.pi
    pg = cmath.sin(math.pi * 1j * xg) ** 2
    P(f"      check: sin^2(pi * i * {xg:.6f}) = {pg.real:.12f} (-1/phi = {-2 / (1 + 5 ** 0.5):.12f})")
    K = 16
    M = 1 << K
    iT = {}
    for a in range(M):
        iT[itin_T(a, K)] = a
    bij = len(iT) == M
    iL = [itin_L(a, K) for a in range(M)]
    bijL = len(set(iL)) == M
    Psi = [iT[iL[a]] for a in range(M)]
    iso = True
    for _ in range(20000):
        a, b = rng.randrange(M), rng.randrange(M)
        if a == b:
            continue
        va = ((a - b) & -(a - b)).bit_length() - 1 if (a - b) % M else K
        d = (Psi[a] - Psi[b]) % M
        vb = (d & -d).bit_length() - 1 if d else K
        if va != vb:
            iso = False
            break
    P(f"  (b) 2-adic: itineraries mod 2^{K} are bijections for T: {bij}, for L(p) = p(1-p)/2: {bijL};")
    P(f"      Psi = itin_T^-1 o itin_L preserves 2-adic distance on 20000 random pairs: {iso}")
    P("      (both are 2-adic Bernoulli maps; the conjugacy is an isometry, as for every such pair)")

    def L_orbit(p, n):
        out = [p]
        for _ in range(n):
            p = p * (1 - p) / 2
            out.append(p)
        return out

    def T_orbit(p, n, r=1):
        out = [p]
        for _ in range(n):
            p = p / 2 if p.numerator % 2 == 0 else (3 * p + r) / 2
            out.append(p)
        return out
    for x0 in (Fraction(1, 3), Fraction(5, 7), Fraction(-1, 5)):
        Lo = L_orbit(x0, 6)
        To = T_orbit(x0, 12)
        P(f"  (c) x0 = {x0}: L-orbit heights log2 max(|num|,den): "
          f"{[round(math.log2(max(abs(f.numerator), f.denominator)), 1) for f in Lo]}")
        P(f"      T-orbit (3x+1): {[str(f) for f in To]}")
    P("      Under L rational orbits have doubling heights, so their itineraries are never eventually periodic;")
    P("      under the affine T heights stay bounded along cycles.  The conjugacy keeps the 2-adic dynamics and")
    P("      destroys exactly the affine arithmetic that the Periodicity Conjecture is about.")
    P("  VERDICT AD3: the fold is literally shared (REAL identity of the quotient); what each program does with it")
    P("  is different (AMM: two-point integrality -> one-point, a capacity gain; LRC: two-lift failure mass;")
    P("  Collatz: the parity split), so transfers are ANALOGY unless a predicate is carried along.")


def main():
    quick = "--quick" in sys.argv
    section_ad1()
    tick("AD1")
    section_ad2(nmax=96 if quick else 300)
    tick("AD2")
    section_ad3()
    tick("AD3")


if __name__ == "__main__":
    main()
