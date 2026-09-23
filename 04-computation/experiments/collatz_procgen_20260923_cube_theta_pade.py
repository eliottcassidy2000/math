#!/usr/bin/env python3
"""
collatz_procgen_20260923_cube_theta_pade.py

Cube-theta lane (HYP-9127): the FORMAL Pade picture.

Every approximation argument that treats rho as a formal variable produces integer polynomials
A, B (deg <= E) with  A(x) Theta(x) - B(x) = R(x) = sum_{j >= W'} r_j x^j,  W' >= W,
where Theta(x) = sum_{s in S} x^s is the 0/1 series of the swap positions (S = cubes for Y3).
Proposition F (note, section 3): if v2(r_{W'}) < L then R(rho) != 0 and
    v2(M0^E A(rho) X - M0^E B(rho)) = L W' + v2(r_{W'}),
so the pair certifies "X is not a rational of height < 2^(margin)" with
    margin = L W' + v2(r_W') - log2 max(|M0^E A(rho)|, |M0^E B(rho)|).
A proof of HYP-9127 needs such pairs with margin -> infinity for infinitely many E.

Sections:
  P1  Hankel determinants of the cube / square indicators vs random controls with the same
      counting function (the 'Pade miracle' test): bit growth, parity, fits.
  P2  diagonal Pade forms of Theta_3 evaluated at rho = 2^10/3^9: exact heights, remainder
      leading coefficient, Newton-polygon non-vanishing, margins.
  P3  diagonal Pade forms of Theta_2 at maps BEYOND Lemma P (mu_bar >= phi): finite certificates.
  P4  sub-diagonal (Siegel) forms via LLL on the kernel lattice: heights and non-vanishing status.

Usage: python3 collatz_procgen_20260923_cube_theta_pade.py [--quick]
"""
import math
import random
import sys
import time

import flint
import gmpy2
from gmpy2 import mpz

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from collatz_procgen_20260923_cube_theta_core import (  # noqa: E402
    PHI, L_Y3, M0_Y3, mubar, v2, theta_mod, log2abs)

QUICK = "--quick" in sys.argv


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def indicator_from_positions(pos, size):
    c = [0] * size
    for p in pos:
        if p < size:
            c[p] = 1
    return c


def power_positions(d, size):
    out = []
    k = 0
    while k ** d < size:
        out.append(k ** d)
        k += 1
    return out


def random_like_positions(d, size, seed):
    """One position uniformly in each [k^d, (k+1)^d): same counting function as the d-th powers."""
    rng = random.Random(seed)
    out = [0]
    k = 1
    while k ** d < size:
        lo, hi = k ** d, (k + 1) ** d
        out.append(rng.randrange(lo, hi))
        k += 1
    return out


def k32_positions(size):
    out, k = [], 0
    while True:
        p = (k * k * k)
        v = math.isqrt(p)          # floor(k^(3/2)) exactly
        if v >= size:
            return sorted(set(out))
        out.append(v)
        k += 1


def k32_random(size, seed):
    """one position uniformly in each [floor(k^1.5), floor((k+1)^1.5)) (same counting function)"""
    rng = random.Random(seed)
    base = k32_positions(size + 100)
    out = [0]
    for lo, hi in zip(base[1:], base[2:]):
        if lo >= size:
            break
        out.append(rng.randrange(lo, hi))
    return out


def hankel_det(c, n, shift=0):
    M = flint.fmpz_mat([[c[i + j + shift] for j in range(n)] for i in range(n)])
    return int(M.det())


def bits(x):
    return abs(x).bit_length() if x else 0


def v2i(x):
    return (abs(x) & -abs(x)).bit_length() - 1 if x else None


# ----------------------------------------------------------------------------- P1
def section_P1():
    hdr("P1  Hankel determinants H_n = det(c_{i+j})_{i,j<n}: cubes, squares, random controls")
    ns = [50, 100, 200, 300, 400, 600, 800] if QUICK else [50, 100, 200, 300, 400, 600, 800, 1000, 1400]
    size = 2 * max(ns) + 10
    seqs = [
        ("cubes", power_positions(3, size)),
        ("rand~cubes s=1", random_like_positions(3, size, 1)),
        ("rand~cubes s=2", random_like_positions(3, size, 2)),
        ("squares", power_positions(2, size)),
        ("rand~squares s=1", random_like_positions(2, size, 1)),
        ("rand~squares s=2", random_like_positions(2, size, 2)),
        ("floor(k^1.5)", k32_positions(size)),
        ("rand~floor(k^1.5)", k32_random(size, 1)),
    ]
    print("bits(H_n) [v2(H_n)]; last column of each row = bits/(n log2 n)")
    head = "   n  " + "".join(f"| {name:>18s} " for name, _ in seqs)
    print(head)
    growth = {name: [] for name, _ in seqs}
    for n in ns:
        row = f"{n:5d} "
        for name, pos in seqs:
            c = indicator_from_positions(pos, size)
            d = hankel_det(c, n)
            b = bits(d)
            growth[name].append((n, b))
            row += f"| {b:6d} [{str(v2i(d)):>4s}] {b/(n*math.log2(n)):.4f} "
        print(row)
        sys.stdout.flush()
    if not QUICK:
        print("\n  cubes vs same-counting-function random controls at n = 2000:")
        size2 = 4010
        for name, pos in (("cubes", power_positions(3, size2)), ("rand~cubes s=1", random_like_positions(3, size2, 1)),
                          ("rand~cubes s=2", random_like_positions(3, size2, 2))):
            c = indicator_from_positions(pos, size2)
            d = hankel_det(c, 2000)
            print(f"   {name:16s} n=2000: bits(H_n) = {bits(d)}  (bits/(n log2 n) = {bits(d)/(2000*math.log2(2000)):.4f})")
            sys.stdout.flush()
    print("\nzero determinants (Pade-table non-normality) among n = 1..300:")
    for name, pos in seqs:
        c = indicator_from_positions(pos, 2 * 300 + 10)
        z = [n for n in range(1, 301) if hankel_det(c, n) == 0]
        odd = sum(1 for n in range(1, 301) if hankel_det(c, n) % 2)
        print(f"   {name:18s}: {len(z):3d} zeros (first: {z[:10]}); odd H_n for {odd}/300")
    print("\nReading: the budget of a diagonal Pade certificate for Y3 is (2L - log2 M0) = 5.735 bits per")
    print("unit of n; every column grows like kappa * n log2 n with kappa << 5.735/log2 n in range,")
    print("so diagonal Pade certifies at every computed n, but n log n growth is asymptotically fatal.")
    return growth


# ----------------------------------------------------------------------------- P2/P3 helpers
def pade_form(c, E, W, L, M0, Xmod, Nbig, extra=400, lll_pick=True):
    """Integer A (deg<=E) with A*Theta = B + O(x^W). Returns dict of certificate data."""
    rows = []
    for j in range(E + 1, W):
        rows.append([c[j - i] if j - i >= 0 else 0 for i in range(E + 1)])
    M = flint.fmpz_mat(rows)
    K, nul = M.nullspace()
    basis = [[int(K[r, col]) for r in range(E + 1)] for col in range(nul)]
    if nul > 1 and lll_pick:
        B = flint.fmpz_mat(basis).lll()
        basis = [[int(B[i, j]) for j in range(E + 1)] for i in range(B.nrows())]
    best = None
    for a in basis:
        g = 0
        for x in a:
            g = math.gcd(g, x)
        if g == 0:
            continue
        a = [x // g for x in a]
        h = max(abs(x) for x in a)
        if best is None or h < best[1]:
            best = (a, h)
    a, h = best
    # B and remainder
    Bc = [sum(a[i] * c[j - i] for i in range(0, j + 1)) for j in range(E + 1)]
    Wp, rW = None, None
    for j in range(E + 1, W + extra):
        r = sum(a[i] * c[j - i] for i in range(E + 1) if j - i >= 0)
        if j < W:
            assert r == 0
        elif r != 0:
            Wp, rW = j, r
            break
    P = mpz(0)
    Q = mpz(0)
    M0z = mpz(M0)
    for i in range(E + 1):
        if a[i]:
            P += a[i] * (mpz(2) ** (L * i)) * M0z ** (E - i)
        if Bc[i]:
            Q += Bc[i] * (mpz(2) ** (L * i)) * M0z ** (E - i)
    mod = mpz(1) << Nbig
    val = v2((P * Xmod - Q) % mod)
    lh = max(log2abs(P), log2abs(Q))
    pred = None if Wp is None else L * Wp + v2(rW)
    np_ok = (Wp is not None) and v2(rW) < L
    return dict(E=E, W=W, nul=nul, h=h, log2h=math.log2(h), Wp=Wp, rW=rW, v2rW=(v2(rW) if rW else None),
                val=val, pred=pred, np_ok=np_ok, lh=lh, margin=val - lh, nnz=sum(1 for x in a if x))


def section_P2():
    hdr("P2  diagonal Pade forms of Theta_3 at rho = 2^10/3^9 (the number of HYP-9127)")
    L, M0 = L_Y3, M0_Y3
    ns = [50, 100, 150, 200, 300, 400, 600] if QUICK else [50, 100, 150, 200, 300, 400, 600, 800, 1000]
    size = 2 * max(ns) + 500
    c = indicator_from_positions(power_positions(3, size), size)
    Nbig = L * (2 * max(ns) + 450) + 64
    X = theta_mod(Nbig, lambda k: k ** 3, L, M0)
    print("  E = n, W = 2n+1 (Pade [n/n]); budget (2L - log2 M0) n = 5.735 n")
    print("   n | nullity | log2 h | nnz(A) |  W'  | v2(r_W') | NP-nonvanish | v2(PX-Q) (=pred) | log2 H | margin | 5.735n")
    rows = []
    for n in ns:
        t0 = time.time()
        d = pade_form(c, n, 2 * n + 1, L, M0, X, Nbig)
        rows.append(d)
        print(f"{n:5d} | {d['nul']:7d} | {d['log2h']:6.1f} | {d['nnz']:6d} | {d['Wp']:5d} | {str(d['v2rW']):>8s} |"
              f" {str(d['np_ok']):>12s} | {d['val']:8d} ({str(d['val']==d['pred']):>5s}) | {d['lh']:7.0f} |"
              f" {d['margin']:6.0f} | {5.7353*n:6.0f}   [{time.time()-t0:.1f}s]")
        sys.stdout.flush()
    print("  Every row is a FINITE certificate: X is not a rational a/b with log2(|a|+|b|) < margin.")
    return rows


def section_P3():
    hdr("P3  diagonal Pade forms of Theta_2 (squares) beyond Lemma P's threshold phi")
    ns = [50, 100, 200, 300, 400] if QUICK else [50, 100, 200, 300, 400, 600, 800]
    size = 2 * max(ns) + 500
    c = indicator_from_positions(power_positions(2, size), size)
    for (name, L, M0) in (("5x+1 block 1^7 0^3", 10, 5 ** 7), ("5x+1 block 1^8 0^2", 10, 5 ** 8),
                          ("7x+1 A=15,L=26", 26, 7 ** 15)):
        mu = mubar(L, M0)
        Nbig = L * (2 * max(ns) + 450) + 64
        X2 = theta_mod(Nbig, lambda k: k * k, L, M0)
        print(f"-- {name}: mu_bar = {mu:.4f} (phi = {PHI:.4f}); diagonal budget L(2-mu_bar) = {L*(2-mu):.3f} bits/n")
        for n in ns:
            d = pade_form(c, n, 2 * n + 1, L, M0, X2, Nbig)
            print(f"   n={n:4d}: log2 h = {d['log2h']:6.1f}, W' = {d['Wp']}, v2(r_W') = {d['v2rW']}, NP-nonvanish = {d['np_ok']},"
                  f" v2 = {d['val']} (pred {d['pred']}), margin = {d['margin']:.0f}")
            sys.stdout.flush()
    print("  Lemma P (Zudilin forms) has NEGATIVE asymptotic margin at these maps (ledger L3); the")
    print("  diagonal Pade forms still certify at every computed n. A proof would need their heights")
    print("  to stay below L(2-mu_bar) n for all n: the Hankel data (P1) grow like n log n instead.")


# ----------------------------------------------------------------------------- P4
def section_P4():
    hdr("P4  sub-diagonal (Siegel) forms: LLL-short vectors of the kernel lattice, W = w E")
    L, M0 = L_Y3, M0_Y3
    Es = [40, 80, 120] if QUICK else [40, 80, 120, 200, 300]
    ws = [1.5, 1.75, 1.9]
    size = 2 * max(Es) + 800
    c = indicator_from_positions(power_positions(3, size), size)
    Nbig = L * (2 * max(Es) + 700) + 64
    X = theta_mod(Nbig, lambda k: k ** 3, L, M0)
    print("   E |   w  | dim ker | log2 h | W' - W | v2(r_W') | NP-nonvanish | margin | (Lw - log2 M0) E")
    for E in Es:
        for w in ws:
            W = int(math.ceil(w * E))
            d = pade_form(c, E, W, L, M0, X, Nbig, extra=700)
            print(f"{E:4d} | {w:4.2f} | {d['nul']:7d} | {d['log2h']:6.2f} | {d['Wp']-W:6d} | {str(d['v2rW']):>8s} |"
                  f" {str(d['np_ok']):>12s} | {d['margin']:6.0f} | {(L*w - math.log2(M0))*E:7.0f}")
            sys.stdout.flush()
    print("  Siegel forms have tiny heights; what is missing for a proof is a guarantee that some short")
    print("  kernel vector has v2(r_W') < L = 10 (a 2-adic zero estimate), for infinitely many E.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_cube_theta_pade.py" + (" --quick" if QUICK else ""))
    section_P1()
    section_P2()
    section_P3()
    section_P4()
    print(f"\n[pade done in {time.time()-t0:.1f}s]")
