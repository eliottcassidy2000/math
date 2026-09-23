#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta_zero.py

HYP-9130 (2-adic zero estimate for the cube Pade lattices), restricted forms.

Lattice:  Lambda_(E,W) = { A in Z^(E+1) : c_j . A = 0 for E < j < W },  c_j(i) = 1_cube(j - i),
so that A(x) Theta_3(x) = B(x) + O(x^W) with deg B <= E; r_W(A) = c_W . A is the x^W coefficient.

  Z1  F2-rank frontier: W_2(E) = max W such that c_(E+1), ..., c_W are independent over F_2.
      Lemma F2 (note, section 3.1): then some A in Lambda_(E,W) has r_W(A) odd (for every W <= W_2(E)).
  Z2  heights: LLL basis of Lambda_(E, ceil(wE)); the least-height vector whose first nonzero remainder
      coefficient has v2 < 10 (the Newton-polygon condition for rho = 2^10/3^9); also the largest basis
      vector (a proxy for the last successive minimum).
  Z3  the echelon case: (E, W] inside one cube interval  =>  independence (checked).

Usage: python3 collatz_procgen_20260923_theta_zero.py [--quick]
"""
import math
import sys
import time

import flint

QUICK = "--quick" in sys.argv
L_Y3 = 10


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def cube_set(limit):
    s = set()
    k = 0
    while k ** 3 <= limit:
        s.add(k ** 3)
        k += 1
    return s


def row_mask(j, E, cubes_sorted):
    m = 0
    for s in cubes_sorted:
        if s > j:
            break
        i = j - s
        if i <= E:
            m |= 1 << i
    return m


def f2_frontier(E, cubes_sorted, wmax=2.2):
    """Largest W such that rows E+1..W are F2-independent (scan up to wmax*E)."""
    piv = {}
    W = E
    for j in range(E + 1, int(wmax * E) + 3):
        r = row_mask(j, E, cubes_sorted)
        while r:
            hb = r.bit_length() - 1
            if hb in piv:
                r ^= piv[hb]
            else:
                piv[hb] = r
                break
        if r == 0:
            return W
        W = j
    return W


def v2(x):
    x = abs(int(x))
    return (x & -x).bit_length() - 1 if x else None


def lattice_basis(E, W, c):
    rows = [[c[j - i] if 0 <= j - i else 0 for i in range(E + 1)] for j in range(E + 1, W)]
    M = flint.fmpz_mat(rows)
    K, nul = M.nullspace()
    basis = flint.fmpz_mat([[int(K[r, col]) for r in range(E + 1)] for col in range(nul)])
    B = basis.lll()
    return [[int(B[i, j]) for j in range(E + 1)] for i in range(B.nrows())]


def lead(a, c, W, E, extra):
    for j in range(W, W + extra):
        r = 0
        for i, ai in enumerate(a):
            if ai and 0 <= j - i < len(c):
                r += ai * c[j - i]
        if r:
            return j, r
    return None, None


def section_Z1():
    hdr("Z1  F2-rank frontier W_2(E)/E for the cube indicator (Lemma F2: parity half of HC-CT1)")
    Es = list(range(20, 301, 20)) + [400, 600, 800] if QUICK else list(range(10, 401, 1)) + list(range(420, 1501, 20))
    cubes = sorted(cube_set(4 * max(Es) + 10))
    rat = []
    t0 = time.time()
    for E in Es:
        W2 = f2_frontier(E, cubes)
        rat.append((W2 / E, E, W2))
    rmin = min(rat)
    print(f"  tested {len(Es)} values of E in [{min(Es)}, {max(Es)}]  ({time.time()-t0:.1f}s)")
    print(f"  min W_2(E)/E = {rmin[0]:.4f} (E = {rmin[1]}, W_2 = {rmin[2]});  mean = {sum(r for r,_,_ in rat)/len(rat):.4f}")
    for thr in (1.4265, 1.5, 1.75, 1.9, 1.95):
        cnt = sum(1 for r, _, _ in rat if r >= thr)
        print(f"  #E with W_2(E) >= {thr} E: {cnt}/{len(rat)}")
    big = [x for x in rat if x[1] >= 100]
    print(f"  E >= 100: min ratio {min(big)[0]:.4f}, max {max(big)[0]:.4f}; 2E - W_2(E) ranges over "
          f"[{min(2*E - W for _, E, W in big)}, {max(2*E - W for _, E, W in big)}]")
    return rat


def section_Z2():
    hdr("Z2  heights of odd-leading (v2 < 10) forms in Lambda_(E, ceil(wE)) for the cube series")
    Es = [40, 80, 120, 160, 200] if QUICK else list(range(20, 301, 5)) + [350, 400, 450, 500, 600]
    ws = [1.5, 1.75, 1.9]
    size = 3 * max(Es) + 2000
    cs = cube_set(size)
    c = [1 if n in cs else 0 for n in range(size)]
    stats = {w: [] for w in ws}
    t0 = time.time()
    for E in Es:
        for w in ws:
            W = int(math.ceil(w * E))
            basis = lattice_basis(E, W, c)
            best = None
            hmax = 0
            for a in basis:
                h = max(abs(x) for x in a)
                hmax = max(hmax, h)
                Wp, r = lead(a, c, W, E, 3000)
                if Wp is None:
                    continue
                if v2(r) < L_Y3 and (best is None or h < best[0]):
                    marg = 10 * Wp + v2(r) - 14.265 * E - math.log2(h * (E + 1) * max(1, len(basis)))
                    best = (h, Wp - W, v2(r), marg)
            stats[w].append((E, best, hmax, len(basis)))
    print(f"  ({time.time()-t0:.1f}s)")
    for w in ws:
        rows = stats[w]
        found = sum(1 for _, b, _, _ in rows if b is not None)
        hs = [b[0] for _, b, _, _ in rows if b is not None]
        hm = [hm for _, _, hm, _ in rows]
        print(f"  w = {w}: HC-CT1 witness (v2(lead) < 10) found for {found}/{len(rows)} values of E;"
              f" max witness height {max(hs)} (log2 {math.log2(max(hs)):.1f}); max LLL-basis height log2 {math.log2(max(hm)):.1f}")
        worst = max(rows, key=lambda t: (t[1][0] if t[1] else float('inf')))
        mm = min(rows, key=lambda t: (t[1][3] if t[1] else -float('inf')))
        print(f"      worst height at E = {worst[0]}: (h, W'-W, v2, margin) = ({worst[1][0]}, {worst[1][1]}, {worst[1][2]},"
              f" {worst[1][3]:.0f}), dim {worst[3]}")
        print(f"      minimal margin 10 W' + v2 - log2(3^9) E - log2(h (E+1) dim) over all E: {mm[1][3]:.1f} at E = {mm[0]}")
    print("  every witness certifies 'X is not a rational of height < 2^margin' (FINITE certificates only).")
    return stats


def section_Z3():
    hdr("Z3  echelon case: (E, W] inside one cube interval [N^3, (N+1)^3) => rows independent")
    ok = True
    for N in range(4, 13):
        E = N ** 3 - 1
        W = (N + 1) ** 3 - 1          # <= 2E + 1 exactly when N >= 4
        cubes = sorted(cube_set(4 * W))
        W2 = f2_frontier(E, cubes, wmax=2.2)
        ok = ok and (W2 >= W)
    print(f"  for E = N^3 - 1, N = 4..12: W_2(E) >= (N+1)^3 - 1 (distinct minimal columns j - N^3): {ok}")
    print("  (ratio (N+1)^3/N^3 -> 1: a provable but weak restricted form; the data of Z1 show W_2(E) ~ 2E)")


def section_Z4():
    hdr("Z4  F2 linear complexity profile of the cube indicator (Berlekamp-Massey)")
    N = 6000 if QUICK else 40000
    cs = cube_set(N + 10)
    bits = [1 if n in cs else 0 for n in range(N)]
    C, B, Lc, m = 1, 1, 0, 1
    S = 0
    maxdev = 0.0
    jumps = []
    t0 = time.time()
    for n in range(N):
        S = (S << 1) | bits[n]
        d = bin(C & S).count("1") & 1
        if d:
            T = C
            C ^= B << m
            if 2 * Lc <= n:
                jumps.append(n + 1 - 2 * Lc)
                Lc = n + 1 - Lc
                B = T
                m = 1
            else:
                m += 1
        else:
            m += 1
        maxdev = max(maxdev, abs(Lc - (n + 1) / 2))
    print(f"  N = {N}: final linear complexity {Lc} (N/2 = {N/2:.0f}); max |L_n - n/2| = {maxdev:.1f}; "
          f"#jumps {len(jumps)}, max jump {max(jumps)}  ({time.time()-t0:.1f}s)")
    from collections import Counter
    cnt = Counter(jumps)
    tot = sum(cnt.values())
    print("  LCP increments (= degrees of the F2 partial quotients) distribution: " +
          ", ".join(f"{k}:{cnt[k]}" for k in sorted(cnt)[:10]))
    print(f"  random binary sequences: increment k with frequency 2^-k; here k=1: {cnt.get(1,0)/tot:.3f},"
          f" k=2: {cnt.get(2,0)/tot:.3f}, k=3: {cnt.get(3,0)/tot:.3f}, k=4: {cnt.get(4,0)/tot:.3f}")
    print("  A near-perfect profile (bounded deviation) is the F2 behaviour behind W_2(E) ~ 2E in Z1.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta_zero.py" + (" --quick" if QUICK else ""))
    section_Z1()
    section_Z2()
    section_Z3()
    section_Z4()
    print(f"\n[zero done in {time.time()-t0:.1f}s]")
