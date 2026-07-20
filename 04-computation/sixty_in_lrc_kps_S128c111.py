#!/usr/bin/env python3
"""sixty_in_lrc_kps_S128c111.py -- kind-pasteur-2026-07-20-S128c111

DOES 60 APPEAR IN THE n=12 AP-UNIQUENESS PROBLEM, AND IS IT THE SAME 60?

THM-1415 deflated the "three sixties" -- ord_1001(2) = lcm(3,10,12), the Pisano
period pi(10) = lcm(3,20), and |A_5| = 5!/2 -- by showing the causes differ and that
60 is the second most frequent lcm of subsets of {1..12} (4.90%), because
60 = lcm(1..6).  The owner now asks how those relate to the LRC n=12 AP uniqueness.

THIS SCRIPT LOOKS RATHER THAN ASSERTS.  Four places a 60 could enter, all checked:

 (1) THE DISTANCE SPECTRUM.  For LRC(n) the extremal family is the AP {1,...,n-1}
     with witness t = 1/n, and the distances are ||k/n|| = min(k, n-k)/n.  So the
     NUMERATORS are exactly {1, ..., floor(n/2)} and the natural constant is
     lcm(1..floor(n/2)).  At n = 13 that is lcm(1..6) = 60.  If so, 60 is the n=13
     member of a family, not a bridge -- and the value AT THE ACTUAL TARGET n = 14
     is lcm(1..7), a different number.  Computed for a range of n.

 (2) THE PRIMITIVE-ROOT FACT.  ord_13(2) = 12 = the number of speeds, i.e. 2 is a
     primitive root mod 13.  That is NOT an lcm coincidence, so it is the one
     candidate with independent content.  Tested by asking whether the AP's WITNESS
     SET depends on it: if every t = j/n with gcd(j,n) = 1 is a witness, then no
     generator is distinguished and ord_n(2) can play no role.

 (3) THE FAREY GAP.  The crux interval is (1/13, 2/25).  Checked: is 2/25 the
     MEDIANT of 1/13 and 1/12?  If so the interval is pure Farey structure with no
     room for a 60.

 (4) THE COVERING MODULI.  The covering-system reduction tests families mod
     lcm(2..Q_0).  lcm(2..5) = 60 exactly -- so a 60 does appear, but as another
     initial-segment lcm, i.e. the SAME attractor mechanism.  Tabulated.
"""
import sys
from math import gcd
from fractions import Fraction as F

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20


def lcm(a, b):
    return a * b // gcd(a, b)


def lcm_range(lo, hi):
    L = 1
    for x in range(lo, hi + 1):
        L = lcm(L, x)
    return L


def multord(a, m):
    if gcd(a, m) != 1:
        return None
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def is_prime(p):
    return p > 1 and all(p % d for d in range(2, int(p ** .5) + 1))


print("=" * 78)
print("(1) THE DISTANCE SPECTRUM of the extremal AP -- where an lcm naturally sits")
print("=" * 78)
print("  LRC(n): speeds {1,...,n-1}, witness t = 1/n, distances ||k/n|| = min(k,n-k)/n.")
print("  So the numerators are exactly {1, ..., floor(n/2)}.")
print()
print("  %-5s %-14s %-10s %-14s %s" % ("n", "#speeds", "spectrum", "lcm(numerators)", ""))
for n in range(5, NMAX + 1):
    half = n // 2
    nums = sorted({min(k, n - k) for k in range(1, n)})
    L = lcm_range(1, half)
    mark = ""
    if n == 13:
        mark = "  <-- the n=12 AP-uniqueness problem"
    if n == 14:
        mark = "  <-- LRC(14), the actual target"
    print("  %-5d %-14d {1..%-6d} %-14d%s" % (n, n - 1, half, L, mark))
    assert nums == list(range(1, half + 1))
print()
print("  So the 60 attached to the 12-speed problem is lcm(1..6) -- the n = 13 member")
print("  of a family.  At the real target n = 14 the same construction gives")
print("  lcm(1..7) = %d, NOT 60." % lcm_range(1, 7))

print()
print("=" * 78)
print("(2) THE PRIMITIVE-ROOT FACT: ord_13(2) = 12.  Does the witness set see it?")
print("=" * 78)
print("  ord_13(2) = %s   (2 is a primitive root mod 13: %s)"
      % (multord(2, 13), multord(2, 13) == 12))
print("  powers of 2 mod 13 : %s"
      % sorted((pow(2, k, 13)) for k in range(12)))
print("  -> they are ALL of {1,...,12}, i.e. the AP itself.  But so is any full")
print("     residue system, so this alone distinguishes nothing.  The real test:")
print()
for n in (11, 13, 17):
    good = []
    for j in range(1, n):
        if gcd(j, n) != 1:
            continue
        M = min(min((k * j) % n, n - (k * j) % n) for k in range(1, n))
        if M == 1:
            good.append(j)
    print("  n = %-3d : witnesses t = j/n attaining the AP optimum 1/%d : j in %s"
          % (n, n, good))
    print("            that is ALL %d units mod %d ; ord_%d(2) = %s"
          % (len(good), n, n, multord(2, n)))
print()
print("  Every unit is a witness, at every n.  No generator is distinguished, so")
print("  ord_n(2) -- and hence the primitive-root status of 2 -- has NO role in the")
print("  AP's extremal structure.  The one non-lcm sixty does not enter LRC.")

print()
print("=" * 78)
print("(3) THE FAREY GAP (1/13, 2/25) -- is it pure mediant structure?")
print("=" * 78)
a, b = F(1, 13), F(1, 12)
med = F(a.numerator + b.numerator, a.denominator + b.denominator)
print("  mediant(1/13, 1/12) = %s   equals 2/25 : %s" % (med, med == F(2, 25)))
print("  |1/13 * 12 - 1/12 * 13| = %s  (Farey neighbours have determinant 1: %s)"
      % (abs(1 * 12 - 1 * 13), abs(1 * 12 - 1 * 13) == 1))
print("  runner-up gap quoted in canon: 1/156, and 156 = 12*13 = %d : %s"
      % (12 * 13, 12 * 13 == 156))
print("  -> the interval is Farey structure on 1/13 and its neighbour 1/12.")
print("     Denominators 13, 12, 25, 156 -- no 60, and no room for one.")

print()
print("=" * 78)
print("(4) THE COVERING MODULI lcm(2..Q) -- where a 60 DOES appear")
print("=" * 78)
for Q in range(2, 13):
    L = lcm_range(2, Q)
    mark = "   <-- 60" if L == 60 else ""
    print("  Q = %-3d : lcm(2..Q) = %-12d%s" % (Q, L, mark))
print()
print("  So a 60 does occur, at Q = 5 (and Q = 6), as another INITIAL-SEGMENT LCM.")
print("  That is the same attractor THM-1415 quantified, not a new phenomenon.")

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
print("  60 appears in the 12-speed problem TWICE, and both are initial-segment lcms:")
print("     * lcm(1..6) = 60, the distance-spectrum numerators of the extremal AP;")
print("     * lcm(2..5) = 60, a covering modulus.")
print("  Neither is a bridge to the other sixties -- they are further instances of")
print("  the SAME cheap mechanism, 60 = lcm(1..6) being the smallest highly-composite")
print("  target for lcms of small integers.")
print()
print("  The only sixty with independent content, ord_13(2) = 12 (2 primitive mod 13),")
print("  provably does NOT enter: every unit j/n is an equally good witness, so no")
print("  generator is distinguished.")
print()
print("  AND THE PARAMETRISED FORM IS THE USEFUL OUTPUT: the constant attached to")
print("  LRC(n) is lcm(1..floor(n/2)).  It is 60 at n = 13 and %d at n = 14."
      % lcm_range(1, 7))
print("  Anything built on '60' would not survive the move to the actual target.")
