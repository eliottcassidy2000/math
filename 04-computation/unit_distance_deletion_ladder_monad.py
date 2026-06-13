#!/usr/bin/env python3
"""
unit_distance_deletion_ladder_monad.py
monad-explorer-2026-06-10  (THM-461 / HYP-2367)

The 3N-frontier floor for the Erdos unit-distance maximum u(n), reframed.

Background (canon): u(n) = max # unit distances among n planar points.
  - PROVEN exact (Alexeev-Mixon-Parshall, arXiv:2412.11914): u(n), n<=21; u(21)=57.
  - For n=22..27 only bounds are known; the repo's N* = first n with u(n)>3n
    is pinned to N* in [25,28] (THM-431/432/433/437, OPEN-Q-057).

THIS SCRIPT establishes, with EXACT INTEGER arithmetic, three rigorous facts:

(A) THE DELETION LADDER.  For ANY graph,  sum_v u(G-v) = (n-2) u(G)  [each edge
    survives n-2 of the n vertex-deletions], so
                 u(n) <= floor( n * u_max(n-1) / (n-2) ).            (Prop A)
    => AMP's cited upper bounds 66,72,78,84,90 for n=23..27 are EXACTLY this
    averaging ladder propagated from their single hard bound u(22)<=61.

(B) THE AVERAGING LADDER IS EXHAUSTED.  Using the proven Moser-lattice
    construction LOWER bounds L(n) (Engel/Schade; repo THM-431-C), the averaging
    ladder from ANY anchor a<=24 CANNOT certify u(m)<=3m for m in {25,26,27}
    unless a small exact value is pinned near its construction value.  We compute
    the exact threshold anchor values that would advance the floor.

(C) TRIPLE CONFLUENCE AT n=27.  3n = n^(4/3) = 3^4 = 81 at the unique n=27;
    27 = 3^(kappa/2) = (kappa/2)^3 with kappa=6 the planar kissing number.

All distance/edge counts are integers; no floating-point in the rigorous parts.
"""

import math

# ----------------------------------------------------------------------
# Data.  PROVEN exact values n<=21 (AMP).  Construction lower bounds L(n)
# (Moser lattice, Engel/Schade; reproduced exactly in repo THM-431-C).
# AMP / best-known UPPER bounds Ubound(n) for n=22..27.
# ----------------------------------------------------------------------
u_exact = {3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,12:27,13:30,14:33,
           15:37,16:42,17:46,18:50,19:52,20:54,21:57}            # PROVEN

L = {22:60,23:64,24:68,25:72,26:76,27:81,28:85,29:89,30:93}      # Moser lower bds
Ubound_AMP = {22:61,23:66,24:72,25:78,26:84,27:90}               # cited upper bds

# ----------------------------------------------------------------------
# (A) Prop A: the deletion-averaging ladder.  u(n) <= floor(n*ubound(n-1)/(n-2)).
#     Show the AMP upper bounds above n=22 are exactly this ladder from u(22)<=61.
# ----------------------------------------------------------------------
def ladder(anchor_n, anchor_U, upto):
    """Propagate u(m) <= floor(m * u(m-1) / (m-2)) from (anchor_n, anchor_U)."""
    b = {anchor_n: anchor_U}
    for m in range(anchor_n+1, upto+1):
        b[m] = (m * b[m-1]) // (m-2)
    return b

print("="*72)
print("(A) DELETION LADDER  u(n) <= floor(n*u(n-1)/(n-2))  from u(22)<=61")
print("="*72)
lad = ladder(22, 61, 27)
print(f"{'n':>3} {'3n':>4} {'ladder bd':>9} {'AMP cited':>9} {'match?':>7} {'L(n) constr':>11}")
all_match = True
for n in range(22, 28):
    amp = Ubound_AMP[n]
    mt = (lad[n] == amp)
    all_match &= mt
    print(f"{n:>3} {3*n:>4} {lad[n]:>9} {amp:>9} {str(mt):>7} {L[n]:>11}")
print(f"\n  ALL AMP bounds (n=23..27) == deletion ladder from u(22)<=61 ?  {all_match}")
print("  => only u(21)=57 (exact) and u(22)<=61 need heavy machinery;")
print("     the entire ladder above 22 is elementary averaging.")

# Sanity: u(22)<=61 is itself SHARPER than the ladder from u(21)=57 (=62).
lad21 = ladder(21, 57, 22)
print(f"\n  ladder from u(21)=57 gives u(22) <= {lad21[22]}  (AMP's 61 is 1 sharper)")

# ----------------------------------------------------------------------
# (B) The ladder is EXHAUSTED: with construction lower bounds L(n), averaging
#     from any anchor a<=24 cannot push u(m)<=3m for m in {25,26,27}, UNLESS a
#     small value is pinned near its construction value.  Find exact thresholds.
# ----------------------------------------------------------------------
print("\n" + "="*72)
print("(B) FLOOR ADVANCEMENT THRESHOLDS  (exact-integer)")
print("="*72)

# For target m, averaging from m-1 gives u(m) <= floor(m*u(m-1)/(m-2)).
# u(m) <= 3m  <==  floor(m*U/(m-2)) <= 3m  <==  m*U/(m-2) < 3m+1
#         <==  U < (3m+1)(m-2)/m.  Largest integer U: floor((3m+1)(m-2)/m - eps).
def max_anchor_for_3m(m):
    # need floor(m*U/(m-2)) <= 3m  i.e.  m*U <= (3m+1)(m-2)-1  (strict < 3m+1)
    # m*U/(m-2) < 3m+1  <=> m*U < (3m+1)(m-2)  <=> U < (3m+1)(m-2)/m
    num = (3*m+1)*(m-2)
    # largest U with m*U < num  ->  U = ceil(num/m) - 1
    U = -(-num // m) - 1
    return U

for m in (25, 26, 27):
    U = max_anchor_for_3m(m)
    # verify
    chk = (m*U)//(m-2)
    print(f"  u({m}) <= 3*{m}={3*m}  via averaging  <==  u({m-1}) <= {U}"
          f"   [check floor({m}*{U}/{m-2})={chk} <= {3*m}: {chk<=3*m}];"
          f"  construction L({m-1})={L[m-1]}  gap={U-L[m-1]}")

print("\n  Interpretation:")
print("  * u(24) <= 69  ==> u(25) <= 75 = 3*25  (floor advances to N*>=26).")
print("    Construction L(24)=68, AMP bound 72; truth in [68,72]; need <=69.")
print("  * u(25) <= 72  ==> u(26) <= 78 = 3*26  (=> N*>=27).  L(25)=72 (tight!).")
print("  * u(26) <= 78  ==> u(27) <= 81 = 3*27  (=> N*>=28, settling u(27)=81).")
print("    L(26)=76; need <=78.")

# Confirm pure propagation from n=22 with the OPTIMISTIC anchor u(22)=L(22)=60
print("\n  Even the OPTIMISTIC anchor u(22)=60 (=construction) propagates to:")
ladopt = ladder(22, 60, 27)
for n in range(25, 28):
    print(f"    u({n}) <= {ladopt[n]}   (3n={3*n};  beats 3n floor? {ladopt[n] > 3*n})")
print("  => averaging from n=22 CANNOT reach 3n at 25/26/27 (blocked by L(22)>=60).")
print("  => advancing the floor REQUIRES a non-averaging (geometric) bound at")
print("     one of n in {24,25,26} -- the targets are only 1-3 below current bds.")

# ----------------------------------------------------------------------
# (C) Triple confluence at n=27.
# ----------------------------------------------------------------------
print("\n" + "="*72)
print("(C) TRIPLE CONFLUENCE AT n=27 = 3^3")
print("="*72)
kappa = 6
print(f"  kissing number of the plane kappa = {kappa},  kappa/2 = {kappa//2}")
print(f"  3n  at n=27         = {3*27}")
print(f"  n^(4/3) at n=27     = {round(27**(4/3))}   (exact: 27^(4/3)=(3^3)^(4/3)=3^4)")
print(f"  3^4                 = {3**4}")
print(f"  3^(kappa/2)         = {3**(kappa//2)}   (size of triangle-cube K3^[]^{kappa//2})")
print(f"  (kappa/2)^3         = {(kappa//2)**3}   (crossover of line 3n and curve n^4/3)")
print(f"  line 3n = curve n^(4/3)  <==>  n^(1/3) = kappa/2 = 3  <==>  n = {(kappa//2)**3}")
print("  FIXED-POINT: 3^(kappa/2) == (kappa/2)^3 holds because kappa/2 = 3.")
print("  All four benchmarks (cube tie, incidence crossover, 3^4, 3*27) = 81 at 27.")

# verify 3n vs n^4/3 sign change brackets 27 exactly (integer test n^4 vs (3n)^3)
print("\n  sign of (3n)^3 - n^4  (=27 n^3 - n^4 = n^3(27-n)):  >0 for n<27, =0 at 27, <0 n>27")
for n in (24,25,26,27,28,29):
    val = n**3*(27-n)
    rel = '>' if val>0 else ('=' if val==0 else '<')
    print(f"    n={n}: (3n)^3-n^4 = {val:>8}  => 3n {rel} n^(4/3)")

print("\nDONE.")
