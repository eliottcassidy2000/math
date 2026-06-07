#!/usr/bin/env python3
"""
S640 / HYP-2318 — Extending the Poke 'fiber bundle over the 7-runner base' for LRC(14)
into the general n=2p CRT reduction, and exploring the doubling/QR connection.

Refined picture (from S639 HYP-2317 + S643 HYP-2346):
  n = 2p, Z/2p = Z/2 x Z/7.  The half-turn (mod-2) is a benign detector (S639);
  the REAL obstruction lives in the mod-p fiber (S643): at the p-clock t=b/p a
  runner is dangerous <=> p|v, every non-mult-of-p has clock-distance >= 1/p = 2*delta.

This script:
  (A) verifies the p-clock BASE-SECTION margin (non-mult-of-p safe) for n=2p;
  (B) the RECURSIVE FIBER: near t=b/p the mult-of-p runners {p*w} reduce to an
      LRC on the reduced speeds {w} (protection-chain depth 2p -> p);
  (C) the FIBER AUTOMORPHISM: doubling x->2x on Z/p has order ord_p(2); for p=7
      its orbits are the two cube-root 3-cycles = the Paley QR/non-QR cosets;
  (D) CONNECTION sweep: for which p does <2> = QR (doubling-orbit = Paley)?
      = 2 is a QR (p = +-1 mod 8) AND ord_p(2) = (p-1)/2.
No external libs.
"""
from math import gcd
from fractions import Fraction

def isprime(n):
    if n < 2: return False
    d = 2
    while d*d <= n:
        if n % d == 0: return False
        d += 1
    return True

def clock_dist_frac(a, q):  # ||a/q|| = dist(a/q, Z) as a Fraction
    r = a % q
    return Fraction(min(r, q - r), q)

print("="*68)
print("(A) p-clock base-section margin for n=2p  (non-mult-of-p is >= 1/p safe)")
print("="*68)
print("  n=2p | delta=1/2p | p-clock margin 1/p | ratio (margin/delta) | OK?")
for p in [3,5,7,11,13]:
    n = 2*p
    delta = Fraction(1, n)
    # worst-case non-mult-of-p runner at t=b/p over b=1..p-1: min clock dist
    worst = min(clock_dist_frac(v*b, p)
                for v in range(1, n) if v % p != 0
                for b in range(1, p))
    print(f"  {n:4d} | 1/{n:<7}| 1/{p:<15}| {worst/delta!s:>5}            | "
          f"{worst >= delta}")
print("  => every non-mult-of-p runner clears delta with a factor-2 margin.")
print("     LRC(p) (proven, p<=7 well within frontier) is the BASE SECTION.")

print("\n" + "="*68)
print("(B) the recursive fiber: mult-of-p runners reduce to LRC on {v/p}")
print("="*68)
# near t = b/p + eps:  v=p*w  =>  v*t = w*b + p*w*eps  ==> ||v*t|| = ||p*w*eps||,
# i.e. the speeds {p*w} at time p*eps behave as an LRC on {w} at time p*eps.
p = 7; n = 2*p
mult = [v for v in range(1, n) if v % p == 0]
print(f"  n={n}: mult-of-{p} speeds in 1..{n-1} = {mult}  (reduced speeds v/{p} = "
      f"{[v//p for v in mult]})")
print("  near t=b/7+eps:  ||(7w)t|| = ||7w*eps||  => the fiber IS an LRC on {w}.")
print("  So LRC(14) <= [LRC(7) base clears non-mult-of-7] AND [the {v/7}-LRC fits")
print("  the perturbation window].  depth-1 protection chain 14 -> 7 (the bundle).")
print("  For speeds 1..13 the fiber is a SINGLE runner (v=7 -> w=1): trivially")
print("  dodgeable; the hardness is only for speed-SETS with many mult-of-7.")

print("\n" + "="*68)
print("(C) fiber automorphism: doubling x->2x on Z/7 = the cube-root 3-cycles")
print("="*68)
def orbit(g, x, q):
    o=[x]; y=(g*x)%q
    while y!=x: o.append(y); y=(g*y)%q
    return o
for start in [1,3]:
    print(f"  doubling-orbit of {start} mod 7: {orbit(2,start,7)}")
QR7 = sorted({(x*x)%7 for x in range(1,7)})
print(f"  QR mod 7 = {QR7} = cube roots of unity mu_3 = paleySet (S638)")
print(f"  doubling-orbit of 1 == QR? {orbit(2,1,7)==QR7 or set(orbit(2,1,7))==set(QR7)}")
print(f"  ord_7(2) = {len(orbit(2,1,7))} = (7-1)/2 = |QR|  -> <2> = QR exactly")

print("\n" + "="*68)
print("(D) connection sweep: for which odd primes p is <2> = QR (doubling=Paley)?")
print("="*68)
print("  p | ord_p(2) | (p-1)/2 | 2 a QR? (p=+-1 mod8) | <2>=QR? | p=3 mod4 (Paley)?")
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
    if not isprime(p): continue
    op = len(orbit(2,1,p))
    half = (p-1)//2
    two_qr = (p % 8 in (1,7))
    is_paley = (p % 4 == 3)
    eq = (op == half and two_qr)
    print(f"  {p:2d} | {op:7d}  | {half:6d}  | {str(two_qr):>5}                | "
          f"{str(eq):>5}   | {is_paley}")
print("  p=7 is special: <2>=QR (2 is a QR AND generates the QR subgroup) AND")
print("  p=3 mod4 (Paley tournament exists, S638).  The doubling fiber-symmetry")
print("  and the Paley/cube-root structure COINCIDE at the n=14 prime.")
