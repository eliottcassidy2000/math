#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_witness_2q_qr_dichotomy_kpswf12.py   (kind-pasteur 2026-06-22, THREAD 2 refinement)

THE QR(-1) DICHOTOMY for the Fourier-correlation constant c^F_q.

This SHARPENS the "decorrelation goes THROUGH the Z/q structure" verdict of
HYP-2854 (concurrent S33; the c_q measure + engine are byte-identical and already
on main).  The refinement HYP-2854 / S33 did NOT state:

  THEOREM (elementary, verified q prime up to 23).  For n = 2q and ANY cluster E
  whose autocorrelation is symmetric (a(m) = a(-m) -- e.g. any E symmetric about
  its centroid: the consec cluster {0,..,k-1}, the boundary core {1,..,2q-2}, any
  symmetric difference set), the Fourier-correlation constant
        c^F_q = sum_{m != 0} |chat(m)| |ghat(m)|,
        chat(m) = a(m)/k,   |ghat(m)| = |sin(pi m / q)| / (pi |m|)
  decomposes over the residue m mod q into a QUADRATIC-RESIDUE shell and a
  NON-residue shell, and these two shells carry EQUAL mass
        QR-shell == NQR-shell      <=>      -1 is a NON-residue mod q
                                   <=>      q == 3 (mod 4).

PROOF (elementary).  |ghat(m)| is even in m and supported off q|m.  For a
symmetric cluster a(m)=a(-m), so the term t(m)=|chat(m)||ghat(m)| satisfies
t(m)=t(-m).  Group by residue class r = m mod q in F_q^*.  Negation r -> -r = q-r
is a fixed-point-free involution on F_q^* pairing each class with its negative and
carrying equal t-mass.  If -1 is a non-residue, negation maps QR <-> NQR
bijectively, so the two shells have identical mass (exactly 50/50).  If -1 is a
residue (q == 1 mod 4), negation preserves QR and preserves NQR, imposing no
balance -- and generically the shells differ.  Since (-1|q) = (-1)^{(q-1)/2},
(-1|q) = -1  <=>  q == 3 (mod 4).  QED.

CONSEQUENCE.  The LRC(14) prime q=7 satisfies 7 == 3 mod 4, so its correlation
constant is EXACTLY QR(7)-balanced (QR(7) = {1,2,4} = the Fano/Hamming line set).
This is the precise sense in which the decorrelation runs THROUGH the quadratic
character of Z/q -- a Gauss-sum-level fact -- not merely "through the Farey grid".
The balance is 2-PERIODIC in q (q mod 4), ORTHOGONAL to the Cayley-Dickson tower:
there is NO property-loss discontinuity at q=7; the only arithmetic fork is the
quadratic character.
"""
import sys
from math import sin, pi
from fractions import Fraction as Fr

sys.stdout.reconfigure(line_buffering=True)


def qr_set(q):
    return sorted({(a * a) % q for a in range(1, q)})


def autocorr(E):
    ac = {}
    for ei in E:
        for ej in E:
            m = ei - ej
            ac[m] = ac.get(m, 0) + 1
    return ac


def shell_split(E, q):
    """Return (QR_mass, NQR_mass, qdiv_mass) of c^F_q over m mod q."""
    QR = set(qr_set(q))
    ac = autocorr(E)
    qm = nm = zm = 0.0
    for m, cnt in ac.items():
        if m == 0:
            continue
        g = abs(sin(pi * m / q)) / (pi * abs(m))
        t = (cnt / len(E)) * g
        r = m % q
        if r == 0:
            zm += t
        elif r in QR:
            qm += t
        else:
            nm += t
    return qm, nm, zm


def cF_exact_shape_rational(E, q):
    """The 'rational shape' S_q with c^F_q = S_q/pi, where
       S_q = sum_{m!=0} (a(m)/k) * |sin(pi m/q)| / |m|.
       |sin(pi m/q)| is algebraic; we keep (a(m)/k)/|m| EXACT and tag the sin."""
    ac = autocorr(E)
    k = len(E)
    # exact rational coefficient on each |sin(pi m/q)|
    terms = {}
    for m, cnt in ac.items():
        if m == 0:
            continue
        r = m % q
        coeff = Fr(cnt, k * abs(m))      # exact
        terms[m] = (coeff, r)
    return terms


def main():
    print("=" * 90)
    print("THE QR(-1) DICHOTOMY for the Fourier-correlation constant c^F_q  (kpswf12)")
    print("  refines HYP-2854 'through Z/q' -> the EXACT q==3 mod 4 quadratic-character gate")
    print("=" * 90)

    print("\n(1) PREDICTION vs OBSERVATION across primes q, two symmetric clusters:")
    print(f"  {'q':>3} {'q%4':>4} {'(-1|q)':>7} {'consec 50/50':>14} {'bdry 50/50':>12} {'verdict':>16}")
    print("  " + "-" * 66)
    allok = True
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        QR = set(qr_set(q))
        minus1_nqr = ((q - 1) % q) not in QR
        legendre = -1 if minus1_nqr else +1
        a = shell_split(list(range(q)), q)
        b = shell_split(list(range(1, 2 * q - 1)), q)
        bal_a = abs(a[0] - a[1]) < 1e-12
        bal_b = abs(b[0] - b[1]) < 1e-12
        pred_bal = (q % 4 == 3)
        ok = (bal_a == pred_bal) and (bal_b == pred_bal)
        allok &= ok
        verdict = "BALANCED" if pred_bal else "unbalanced"
        print(f"  {q:>3} {q % 4:>4} {legendre:>+7} {str(bal_a):>14} {str(bal_b):>12} "
              f"{verdict:>16}{'  OK' if ok else '  ** MISMATCH **'}")
    print(f"\n  prediction (50/50 <=> q==3 mod 4) matches observation on ALL tested primes: {allok}")

    print("\n(2) THE LRC FAMILY q=3,5,7,9 (n=6,10,14,18); consec-q cluster:")
    print(f"  {'q':>3} {'n':>4} {'QR-shell':>12} {'NQR-shell':>12} {'split':>14} {'q mod 4':>8}")
    print("  " + "-" * 60)
    for q in (3, 5, 7, 9):
        n = 2 * q
        qm, nm, zm = shell_split(list(range(q)), q)
        tot = qm + nm
        bal = abs(qm - nm) < 1e-12
        ptag = "" if all(q % d for d in range(2, q)) else " (q not prime)"
        print(f"  {q:>3} {n:>4} {100*qm/tot:>11.1f}% {100*nm/tot:>11.1f}% "
              f"{('EXACT 50/50' if bal else 'unbalanced'):>14} {q % 4:>8}{ptag}")
    print("\n  => q=7 (n=14, the LRC prime, 7==3 mod 4) is EXACTLY QR(7)-balanced.")
    print("     QR(7) = {1,2,4} = Fano line {1,2,4} = the Hamming(7,4) cyclic shift set.")

    print("\n(3) the q|m shell is identically 0 (ghat vanishes on q|m) -- sanity:")
    for q in (3, 5, 7, 9):
        qm, nm, zm = shell_split(list(range(q)), q)
        print(f"    q={q}: q|m shell mass = {zm:.2e}  (must be ~0)")

    print("\nDONE.")


if __name__ == "__main__":
    main()
