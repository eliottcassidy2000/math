#!/usr/bin/env python3
"""
lrc14_angleF_fourier_lattice_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE F part 2 — the EXACT Fourier-over-relation-lattice form of corr(E),
and the geometry-of-numbers bound test.

We verify the identity
   meas(S7(E)) = sum_{n in Z^{k} : sum_e n_e e = 0}  K(n),   K(n) = sum_T (-1)^|T| prod_e chat_T(n_e)
where chat_T is the Fourier coefficient of 1_{B_T}, B_T = [0,1)\ U_{j in T}[j/7,(j+1)/7).
The n=0 vector gives the main term M7(k); corr(E) = sum over NONZERO lattice vectors.

We then evaluate the lattice sum by enumerating all n with |n_e|<=N0 in the
offset relation lattice Lambda^o(E) = {m: sum_{e!=0} m_e e = 0} (n_0 free, but its
contribution to K factors as a scalar). We check truncated lattice sum -> meas(S7).

KEY GEOMETRY-OF-NUMBERS QUANTITY.  Each lattice vector n contributes |K(n)| <=
prod_e (bound on |chat_T(n_e)|) ~ prod_{n_e != 0} c/|n_e|.  So the dominant
contributions come from SHORT lattice vectors (small support, small entries) =
short relations.  The total |corr| <= sum_{0!=n in Lambda^o} |K(n)|, and the
number of lattice vectors in a box of radius R is ~ vol(box)/covol(Lambda^o).
Hence small covolume => more short vectors => larger corr bound.  We make the
per-vector kernel bound explicit and compute the partial sums by shell.

EXACT chat_T(n):
  1_{B_T}(y) = 1 - sum_{j in T} 1_{[j/7,(j+1)/7)}(y).
  Fourier of 1_{[a,a+L)}: at freq n!=0 is (e^{-2pi i n a} - e^{-2pi i n (a+L)})/(2 pi i n).
  With a=j/7, L=1/7:  chat_{single j}(n) = e^{-2pi i n j/7} (1 - e^{-2pi i n/7})/(2 pi i n).
  chat_T(0) = 1 - |T|/7.  chat_T(7m)=0 for m!=0 (since 1-e^{-2pi i (7m)/7}=0). [THM-503]
All exact via cmath at moderate precision; we only need to CONFIRM the lattice
picture and the short-vector dominance, not exact rationals here.
"""
from __future__ import annotations
import sys, itertools, math, cmath
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
            y = (e * xm) % 1
            secs.add(int(y * 7))
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    s = F(0)
    for t in range(0, 7):
        s += F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1)
    return s


def chat_T(n, T):
    """Fourier coeff of 1_{B_T} at frequency n. T subset of {0..6}."""
    if n == 0:
        return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0:
        return 0j  # THM-503 vanishing
    s = 0j
    for j in T:
        # Fourier of 1_{[j/7,(j+1)/7)} at n
        a = j / 7.0
        s += (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)
    return -s  # 1_{B_T} = 1 - sum_j 1_{sector j}; constant 1 only contributes at n=0


def K_kernel(n_vec):
    """K(n) = sum_T (-1)^|T| prod_e chat_T(n_e), T over subsets of {0..6}.
    Sector 0 is auto-hit by e=0 so we can restrict missed sectors to {1..6}? No:
    inclusion-exclusion is over ALL 7 sectors; but e=0 pins sector 0 so any T
    containing 0 gives chat_{T}(0-coordinate via e0)... handled by n_0 freedom.
    Here n_vec already EXCLUDES e0; we fold e0 by noting e0=0 => 1_{B_T}(0).
    1_{B_T}(0): is 0 in B_T? 0 in [0,1/7) which is sector 0, so 0 in B_T iff 0 not in T.
    So the e0 factor = [0 not in T]. Thus only T<={1..6} survive."""
    total = 0j
    for r in range(0, 7):
        for T in itertools.combinations(range(1, 7), r):  # 0 excluded (e0 pins it)
            prod = 1.0 + 0j
            for ne in n_vec:
                prod *= chat_T(ne, set(T))
            total += ((-1) ** len(T)) * prod
    return total


def offset_relation_vectors(E, N0):
    """All integer vectors m over nonzero offsets with |m_j|<=N0 and sum m_j e_j=0.
    Returns list including the zero vector."""
    nz = [e for e in E if e != 0]
    out = []
    ranges = [range(-N0, N0 + 1)] * len(nz)
    for m in itertools.product(*ranges):
        if sum(mi * ei for mi, ei in zip(m, nz)) == 0:
            out.append(m)
    return out, nz


def main():
    print("=" * 90)
    print("ANGLE F part 2: exact Fourier-over-relation-lattice form of corr(E)")
    print("=" * 90)
    # Verify lattice sum reproduces measS7 for a small shape, by truncation.
    print("Truncated lattice sum vs exact meas(S7) (k=5, small offsets):")
    for E in [(0, 1, 2, 3, 4), (0, 1, 2, 4, 5), (0, 1, 3, 5, 7)]:
        exact = float(measS7(E))
        m7 = float(M7(len(E)))
        for N0 in (3, 6, 10):
            vecs, nz = offset_relation_vectors(E, N0)
            s = 0j
            for m in vecs:
                s += K_kernel(m)
            print(
                f"  E={E} N0={N0:>2}: lattice_sum={s.real:.5f} (im {s.imag:+.1e}) "
                f"exact={exact:.5f} M7={m7:.5f} corr_exact={exact-m7:+.5f}"
            )
        print()

    # Per-shell decomposition: how much of corr lives in short relations.
    print("=" * 90)
    print("Shell decomposition: |corr| concentration by max-entry of relation vector")
    print("=" * 90)
    for E in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 3, 7, 12, 20, 30, 44)]:
        m7 = float(M7(len(E)))
        exact = float(measS7(E))
        print(f"E={E}  corr_exact = {exact-m7:+.5f}")
        vecs, nz = offset_relation_vectors(E, 8)
        from collections import defaultdict
        shell = defaultdict(float)
        for m in vecs:
            if all(x == 0 for x in m):
                continue
            sh = max(abs(x) for x in m)
            shell[sh] += K_kernel(m).real
        cum = 0.0
        for sh in sorted(shell):
            cum += shell[sh]
            print(f"    shell maxentry<= {sh}: shell_contrib={shell[sh]:+.5f}  cumulative={cum:+.5f}")
        print()
    print("DONE.")


if __name__ == "__main__":
    main()
