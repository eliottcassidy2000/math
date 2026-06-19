#!/usr/bin/env python3
"""
lrc14_angleF_verify_identity_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

Tight verification that the Fourier-over-relation-lattice identity
   meas(S7(E)) = M7(k) + corr(E),
   corr(E) = sum_{0 != n in Lambda^o(E)} K(n),  K(n)=sum_{T<={1..6}}(-1)^|T| prod_j chat_T(n_j)
reproduces the EXACT measS7 for genuine k>=7 shapes, by ENUMERATING ONLY the
lattice (not the ambient box) -> fast, high truncation possible.

We enumerate Lambda^o(E) = {n: sum n_j e_j = 0} up to a box radius using the
LLL-reduced basis and coefficient enumeration (so the cost is ~ (#lattice pts),
not (box)^d). This lets us push N0 high and confirm convergence to measS7.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from fractions import Fraction as F
import sympy
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
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    return sum(F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


def shat(n, j):
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)


SUBSETS = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SIGN = {T: (-1) ** len(T) for T in SUBSETS}


def chat_T(n, T):
    if n == 0:
        return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0:
        return 0j
    return -sum(shat(n, j) for j in T)


def lll_basis(nz):
    M = sympy.Matrix([nz])
    ns = M.nullspace()
    basis = []
    for v in ns:
        L = 1
        for x in v:
            fr = F(int(x.p), int(x.q))
            L = L * fr.denominator // math.gcd(L, fr.denominator)
        iv = [int(x * L) for x in v]
        g = 0
        for c in iv:
            g = math.gcd(g, abs(c))
        if g:
            iv = [c // g for c in iv]
        basis.append(iv)
    # crude LLL via sympy's not available; use Gram-Schmidt size reduction loop
    return basis


def enum_lattice(nz, N0, max_coeff=4):
    """Enumerate lattice vectors n with |n_j|<=N0 by integer combos of basis."""
    basis = lll_basis(nz)
    d = len(basis)
    seen = set()
    out = []
    for coeffs in itertools.product(range(-max_coeff, max_coeff + 1), repeat=d):
        v = tuple(sum(coeffs[i] * basis[i][t] for i in range(d)) for t in range(len(nz)))
        if all(abs(x) <= N0 for x in v) and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def kernel(n):
    K = 0j
    for T in SUBSETS:
        p = 1.0 + 0j
        for ni in n:
            p *= chat_T(ni, T)
            if p == 0:
                break
        K += SIGN[T] * p
    return K


def main():
    print("Verify  meas(S7) = M7(k) + sum_{0!=n in Lambda^o} K(n)  for genuine k>=7:")
    print("(enumerate ONLY lattice points => high truncation)")
    print("-" * 80)
    cases = [
        (0, 1, 2, 3, 4, 5, 6),          # k=7 AP
        (0, 1, 2, 3, 4, 5, 7),          # k=7 perturbed
        (0, 1, 2, 3, 4, 5, 6, 7),       # k=8 AP
        (0, 1, 2, 3, 5, 6, 7, 8),       # k=8 two-run
    ]
    for E in cases:
        k = len(E)
        nz = [e for e in E if e != 0]
        exact = float(measS7(E))
        m7 = float(M7(k))
        for (N0, mc) in [(20, 5), (40, 7)]:
            vecs = enum_lattice(nz, N0, mc)
            s = 0j
            for n in vecs:
                if all(x == 0 for x in n):
                    continue
                s += kernel(n)
            approx = m7 + s.real
            print(
                f"E={E} k={k} N0={N0}: M7+corr={approx:.6f} exact={exact:.6f} "
                f"|diff|={abs(approx-exact):.2e} (#latpts={len(vecs)})"
            )
        print()
    print("DONE.")


if __name__ == "__main__":
    main()
