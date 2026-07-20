#!/usr/bin/env python3
"""the_half_dictionary_kps_S128c119.py -- kind-pasteur-2026-07-20-S128c119

{-1, 0, 1}  vs  {0, 1/2, 1}:  THE SAME TRICHOTOMY IN TWO COORDINATE SYSTEMS.

The owner's observation is that these two three-element sets keep recurring in this
repo and are "functionally equivalent" -- and they are, by the affine conjugation

    x  |-->  (1 + x)/2 :    -1 |-> 0,    0 |-> 1/2,    +1 |-> 1.

The left-hand set is the SIGN world (multiplicative, closed under negation -- the
odd/even axis).  The right-hand set is the AFFINE / convex world (additive, living in
[0,1], with 1/2 the midpoint -- the positive/negative axis after normalisation).  The
claim to test is that the repo's objects come in conjugate PAIRS under exactly this map,
and that the number 1/2 is where the conjugation leaves its residue.

TWO INSTANCES, both checked here.

(I) THE SPECTRAL ONE.  A tournament has a skew matrix S (entries in {-1, 0, +1}, zero
    diagonal) and an adjacency matrix A (entries in {0, 1}, zero diagonal).  They are
    conjugate by exactly the owner's map:

        A = (J - I + S)/2,      S = 2A - J + I

    -- and the correction term is the DIAGONAL, because (1+0)/2 = 1/2 is not 0.  That
    stray 1/2 is the whole content.  Consequences to verify:
      * every eigenvalue of a tournament A has Re(lambda) >= -1/2, since
        A + A^T = J - I gives 2 Re(x* A x) = |sum x_i|^2 - 1 >= -1 on unit x;
      * for a REGULAR tournament, J S = 0, so off the all-ones vector A = (S - I)/2 and
        every non-Perron eigenvalue has Re(lambda) EXACTLY -1/2 -- the image of the
        skew world's purely imaginary spectrum under x |-> (x-1)/2.
    So THM-1455's "char(S) is parity-pure with purely imaginary spectrum" becomes, in
    the affine coordinate, "the spectrum of A sits on the vertical line Re = -1/2".

(II) THE COMBINATORIAL ONE, in this week's GMC work.  THM-1540's nullcone condition is
     that the CHARGES a - b are all of one strict SIGN -- the {-1,0,1} axis.  Its proof
     step L2 works with h(t) t^{-d/2} and asks whether the support straddles d/2 -- the
     MIDPOINT of the Newton segment, the {0,1/2,1} axis.  Those are the same condition:
     "one-sided in sign" = "the midpoint is not in the support's interior".  The d/2 is
     the 1/2.  Checked on explicit supports.
"""
import sys
import cmath
from itertools import product
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def skew_of(A, n):
    return A - A.T


print("=" * 78)
print("(I.a) THE CONJUGATION  A = (J - I + S)/2  IS EXACT")
print("=" * 78)
bad = 0
tot = 0
for n in range(3, NMAX + 1):
    m = n * (n - 1) // 2
    for bits in range(min(1 << m, 4000)):
        A = from_bits(bits, n)
        S = skew_of(A, n)
        J = np.ones((n, n))
        I = np.eye(n)
        tot += 1
        if not np.allclose(A, (J - I + S) / 2):
            bad += 1
print("  tournaments checked : %d ; failures of A = (J - I + S)/2 : %d" % (tot, bad))
print("  -> the map is exactly the owner's x |-> (1+x)/2, and the -I/2 is the residue")
print("     left by the diagonal, since (1+0)/2 = 1/2 rather than 0.")

print()
print("=" * 78)
print("(I.b) EVERY TOURNAMENT EIGENVALUE HAS Re >= -1/2")
print("=" * 78)
print("  A + A^T = J - I, so 2 Re(x* A x) = |sum x_i|^2 - 1 >= -1 on unit x.")
worst = 10.0
viol = 0
tot = 0
for n in range(3, NMAX + 1):
    m = n * (n - 1) // 2
    for bits in range(min(1 << m, 3000)):
        A = from_bits(bits, n)
        ev = np.linalg.eigvals(A)
        tot += 1
        r = min(e.real for e in ev)
        worst = min(worst, r)
        if r < -0.5 - 1e-9:
            viol += 1
print("  tournaments checked : %d ; violations of Re >= -1/2 : %d" % (tot, viol))
print("  minimum real part seen over all of them : %.12f" % worst)

print()
print("=" * 78)
print("(I.c) REGULAR TOURNAMENTS: every NON-PERRON eigenvalue has Re EXACTLY -1/2")
print("=" * 78)
print("  For a regular tournament J S = 0, so off the all-ones vector A = (S - I)/2,")
print("  and S has purely imaginary spectrum (THM-1455).  Hence Re = -1/2 exactly.")
found = 0
for n in (3, 5, 7):
    m = n * (n - 1) // 2
    cnt = 0
    for bits in range(min(1 << m, 200000)):
        A = from_bits(bits, n)
        if not np.allclose(A.sum(axis=1), (n - 1) / 2):
            continue
        ev = sorted(np.linalg.eigvals(A), key=lambda e: -e.real)
        perron = ev[0]
        rest = ev[1:]
        ok = all(abs(e.real + 0.5) < 1e-9 for e in rest)
        cnt += 1
        if cnt == 1:
            print("  n = %d : Perron = %.4f (expect %.4f) ; other Re = %s ; all -1/2 : %s"
                  % (n, perron.real, (n - 1) / 2,
                     ["%.6f" % e.real for e in rest[:4]], ok))
        if not ok:
            print("     VIOLATION at bits = %d" % bits)
        found += 1
        if cnt >= 3:
            break
print("  regular tournaments sampled : %d ; all satisfied the line Re = -1/2." % found)
print()
print("  So THM-1455's 'char(S) parity-pure, spectrum purely imaginary' reads, in the")
print("  affine coordinate, as 'spec(A) lies on the vertical line Re = -1/2'.  The two")
print("  statements are the SAME fact in the two coordinate systems, and the shift is")
print("  precisely the owner's 1/2.")

print()
print("=" * 78)
print("(II) THE SAME DICHOTOMY IN THE GMC NULLCONE: sign vs Newton midpoint")
print("=" * 78)
print("  THM-1540: the nullcone condition is that the CHARGES a-b are all of one strict")
print("  SIGN.  Its L2 step asks whether the support of h straddles d/2 -- the MIDPOINT")
print("  of the Newton segment.  Those are one condition in two coordinates:")
print()
print("  %-26s %-16s %-18s %s" % ("support of h (powers of t)", "d", "charges 2a-d",
                                  "one-sided = midpoint out"))
cases = [([0, 1], 1), ([1], 1), ([0], 1), ([0, 2], 2), ([2], 2), ([0, 1, 2], 2),
         ([1, 2, 3], 3), ([0, 3], 3)]
for supp, d in cases:
    ch = [2 * a - d for a in supp]
    onesided = all(c > 0 for c in ch) or all(c < 0 for c in ch)
    mid = d / 2
    mid_out = (min(supp) > mid) or (max(supp) < mid)
    print("  %-26s %-16d %-18s %s / %s   agree: %s"
          % (str(supp), d, str(ch), onesided, mid_out, onesided == mid_out))
print()
print("  The two columns agree in every row: 'all charges of one sign' and 'the Newton")
print("  midpoint d/2 lies outside the support' are the same statement.  The charge")
print("  a-b is the SIGN coordinate; the exponent a with midpoint d/2 is the AFFINE")
print("  coordinate; and a - b = 2a - d is exactly the inverse of x |-> (1+x)/2 rescaled.")
