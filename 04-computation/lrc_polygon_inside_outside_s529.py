#!/usr/bin/env python3
"""
lrc_polygon_inside_outside_s529.py    oracle-2026-06-01-S528/S529

The regular polygon's OUTSIDE (sides) vs its HIDDEN INSIDE (diagonals), the
dihedral group, and how this organizes the LRC covering sum.

THESIS
------
A tournament is an orientation of K_n = the 1-skeleton of the (n-1)-simplex.
The maximally symmetric tournament is the REGULAR POLYGON one (rotational /
quadratic-residue), the LRC tight witness. Its chords split into:
  * OUTSIDE  = sides       = skip-1 chords = the Hamiltonian cycle = ranking
               (the cut / score space, the base path of the tiling model).
  * INSIDE   = diagonals   = skip>=2 chords = the genuinely cyclic content
               (the cycle / tile space; the 3-cycles live here).
The dihedral group D_n permutes sides as one orbit and diagonals as shells
(skip-j orbits). Two classical facts tie this to LRC's threshold 1/n:
  (A) CHORD PRODUCT:  prod_{k=1}^{n-1} |1 - w^k| = n     (w = e^{2pi i/n}).
      The single-sector / nearest-neighbour arc gap of the regular n-gon is
      exactly 1/n -- the LRC threshold IS the n-gon's own gap.
  (B) For prime n, the diagonal orientation is the Legendre symbol; the
      balance of the inside is the GAUSS SUM, |sum chi(k) w^k| = sqrt(n).

LRC COVERING SUM, organized by the polygon shells
--------------------------------------------------
Observer lonely at time t  <=>  every runner v_i t lies in the SAFE band
B = [1/n, 1-1/n]. With  f = 1_B,  f^(0) = 1-2/n,  f^(m) = -sin(2pi m/n)/(pi m),
    |LONELY(v)| = INT_0^1 prod_i f(v_i t) dt
                = SUM over (m_1..m_{n-1}) with  SUM m_i v_i = 0  of  prod f^(m_i).
The resonance condition SUM m_i v_i = 0 is the algebraic shadow of the polygon
diagonals. The MAIN term (all m_i = 0) = (1-2/n)^{n-1} = the "outside"
/independence value. Corrections are graded by RESONANCE ORDER r = #{i: m_i!=0}:
  r=2  pairwise  (the skip-1 / outer relations),
  r>=3 multi-way (the deep inside diagonals = the "higher-resonance debt"
       that blocks n>=4, S526/S527).

This script:
  (1) verifies the cyclotomic chord-product = n and prints the skip-shell sizes;
  (2) builds the rotational regular-polygon tournament, verifies it is regular
      for odd n, splits sides vs diagonals, and (prime n) matches the QR/Paley
      orientation; computes the Gauss-sum magnitude sqrt(n);
  (3) computes |LONELY(v)| by direct integration (ground truth) and by the
      resonance sum graded by order; shows the AP (regular polygon) -> 0 (tight),
      and that the n=3 order-2 term reproduces S526's Legendre closed form;
  (4) isolates the INSIDE DEBT = sum of order>=3 resonance terms, showing it is
      identically 0 for n=3 (two speeds admit no 3-term resonance) and first
      switches on at n=4 -- the geometric birth of the open case.
"""
from math import sin, pi, gcd, cos
from fractions import Fraction
from itertools import product as iproduct
import cmath

# ----------------------------------------------------------------------
# (1) The polygon: chord product = n, and the skip-shells
# ----------------------------------------------------------------------
def chord_product(n):
    w = cmath.exp(2j*pi/n)
    p = 1.0
    for k in range(1, n):
        p *= abs(1 - w**k)
    return p

def skip_shells(n):
    """skip j connects vertices i,i+j; j=1..floor(n/2). Returns {j: (count, chord_len, kind)}."""
    shells = {}
    for j in range(1, n//2 + 1):
        # number of distinct chords of skip j on n-gon
        cnt = n if (n % 2 == 1 or j != n//2) else n//2   # diameters at even n are shared
        length = 2*abs(sin(pi*j/n))
        kind = "SIDE (outside)" if j == 1 else ("DIAMETER (inside)" if (n%2==0 and j==n//2) else "DIAGONAL (inside)")
        shells[j] = (cnt, length, kind)
    return shells

# ----------------------------------------------------------------------
# (2) The regular-polygon tournament (rotational), sides vs diagonals, QR
# ----------------------------------------------------------------------
def rotational_tournament(n):
    """i -> j iff (j-i) mod n in {1,..,(n-1)//2} (ahead within leading semicircle).
    Defined for odd n (regular tournament). Returns adjacency dict of out-sets."""
    half = (n-1)//2
    out = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (j - i) % n in set(range(1, half+1)):
                out[i].add(j)
    return out

def is_regular(out, n):
    degs = {len(out[i]) for i in range(n)}
    return len(degs) == 1, degs

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def qr_tournament_matches_rotational(p):
    """For prime p: i->j iff (j-i) is a QR mod p. Compare arc structure to rotational
    by checking whether the 'forward skips' {1..(p-1)/2} are exactly the QRs (true iff
    -1 is a non-residue, i.e. p=3 mod 4); otherwise they are isomorphic rotational tours.
    We just report the Gauss sum magnitude (the inside balance)."""
    g = sum(legendre(k, p) * cmath.exp(2j*pi*k/p) for k in range(1, p))
    return abs(g)

# ----------------------------------------------------------------------
# (3) LRC covering sum, graded by resonance order = polygon shell depth
# ----------------------------------------------------------------------
def fhat(m, n):
    if m == 0:
        return 1.0 - 2.0/n
    return -sin(2*pi*m/n)/(pi*m)

def lonely_measure_direct(speeds, n, G=400000):
    """ground truth: fraction of t with all runners in OPEN safe band (1/n,1-1/n)."""
    lo, hi = 1.0/n, 1.0 - 1.0/n
    cnt = 0
    for i in range(G):
        t = (i + 0.5)/G
        if all(lo < ((s*t) % 1.0) < hi for s in speeds):
            cnt += 1
    return cnt/G

def lonely_measure_resonance(speeds, n, M=12):
    """Sum over (m_i) with sum m_i v_i = 0, |m_i|<=M, of prod fhat(m_i),
    graded by resonance order r = #nonzero m_i. Returns dict order->partial sum."""
    by_order = {}
    rng = range(-M, M+1)
    for ms in iproduct(rng, repeat=len(speeds)):
        if sum(m*v for m, v in zip(ms, speeds)) != 0:
            continue
        r = sum(1 for m in ms if m != 0)
        term = 1.0
        for m in ms:
            term *= fhat(m, n)
        by_order[r] = by_order.get(r, 0.0) + term
    return by_order

def n3_legendre_closed(a, b):
    """S526 closed form for n=3 pairwise overlap of SAFE: 1/9 + (2/9) chi(a)chi(b)/(ab)."""
    def chi3(k):
        k %= 3
        return 0 if k == 0 else (1 if k == 1 else -1)
    return Fraction(1,9) + Fraction(2,9)*chi3(a)*chi3(b)*Fraction(1, a*b)

# ----------------------------------------------------------------------
def main():
    print("="*72)
    print("(1) THE POLYGON: chord product = n, and the skip-shells (outside/inside)")
    print("="*72)
    for n in range(3, 11):
        p = chord_product(n)
        gap = 1.0/n
        print(f"  n={n:2d}: prod|1-w^k| = {p:.6f}  (= n? {abs(p-n)<1e-6})   nearest gap = 1/n = {gap:.4f}")
    print()
    for n in (3,4,5,6,7,8):
        sh = skip_shells(n)
        print(f"  n={n} shells:")
        for j,(cnt,length,kind) in sh.items():
            print(f"      skip {j}: {cnt:2d} chords, len {length:.4f}  {kind}")
    print()

    print("="*72)
    print("(2) REGULAR-POLYGON TOURNAMENT: regular?, sides vs diagonals, Gauss sum")
    print("="*72)
    for n in (3,5,7,9,11):
        out = rotational_tournament(n)
        reg, degs = is_regular(out, n)
        # classify each arc by skip
        side_arcs = sum(1 for i in out for j in out[i] if (j-i)%n==1 or (i-j)%n==1)
        diag_arcs = sum(len(out[i]) for i in out) - side_arcs
        print(f"  n={n:2d}: regular={reg} (out-degrees={degs}); "
              f"arcs on sides(skip1)={side_arcs}, on diagonals(inside)={diag_arcs}")
    print("  Gauss sum magnitude |sum chi(k) w^k| (inside balance), should be sqrt(p):")
    for p in (3,5,7,11,13):
        g = qr_tournament_matches_rotational(p)
        print(f"      p={p:2d}: |g| = {g:.6f}  vs sqrt(p) = {p**0.5:.6f}")
    print()

    print("="*72)
    print("(3) LRC COVERING SUM graded by resonance order (= shell depth)")
    print("="*72)
    tests = {
        3: [(1,2),(1,3),(2,3),(1,4)],
        4: [(1,2,3),(1,2,4),(1,3,5)],
        5: [(1,2,3,4),(1,2,3,5)],
    }
    for n, sets in tests.items():
        print(f"  --- n={n}  (main term (1-2/n)^(n-1) = {(1-2.0/n)**(n-1):.6f}) ---")
        for v in sets:
            direct = lonely_measure_direct(v, n)
            by_ord = lonely_measure_resonance(v, n, M=14)
            total = sum(by_ord.values())
            ap = " [AP=regular polygon]" if list(v)==list(range(1,n)) else ""
            ordstr = "  ".join(f"r{r}:{by_ord[r]:+.5f}" for r in sorted(by_ord))
            print(f"    v={v}{ap}")
            print(f"        direct |LONELY|={direct:.6f}   resonance-sum={total:.6f}")
            print(f"        by order: {ordstr}")
    print()

    print("="*72)
    print("(4) n=3: order-2 term reproduces S526 Legendre form; INSIDE DEBT birth")
    print("="*72)
    print("  n=3 pairwise (order-2) overlap vs S526 closed form 1/9+(2/9)chi(a)chi(b)/(ab):")
    for (a,b) in [(1,2),(1,3),(2,3),(1,4),(2,5)]:
        # SAFE pairwise measure = |LONELY| for two speeds (n=3): direct
        direct = lonely_measure_direct((a,b), 3)
        closed = float(n3_legendre_closed(a,b))
        print(f"    (a,b)=({a},{b}): direct|LONELY|={direct:.6f}  S526 closed={closed:.6f}  "
              f"match={abs(direct-closed)<2e-3}")
    print()
    print("  INSIDE DEBT = sum of resonance orders r>=3 (the deep diagonals):")
    for n, sets in tests.items():
        for v in sets:
            by_ord = lonely_measure_resonance(v, n, M=14)
            inside = sum(val for r,val in by_ord.items() if r>=3)
            print(f"    n={n} v={v}: inside-debt(r>=3) = {inside:+.6f}"
                  + ("   (=0: no 3-term resonance possible with 2 speeds)" if n==3 else ""))
    print()
    print("  => n=3 has NO inside debt (the triangle's only inside is its single")
    print("     diagonal-class = the 3-cycle, captured entirely by the order-2 Legendre")
    print("     term). The inside debt first switches on at n=4: that is the geometric")
    print("     birth of the open case -- the polygon acquires genuine interior diagonals.")

if __name__ == "__main__":
    main()
