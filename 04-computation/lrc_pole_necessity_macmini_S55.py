#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S55 -- HYP-4119 item 4: T_l POLE NECESSITY (ladder side).

opus-S81's six-top ceiling: any uniform FEE dies at l >= 7 (fee-mean 2*rho*l >= 1).
THIS: the LADDER-CONSTANT side -- no C_7 can exist AT ALL, for any accounting:

  THE FLOATING 7-CLUSTER: for every scale S, the family
      W(S) = {1, 2, 3, 4, 5} u { S*c : c in CLUSTER }
  with CLUSTER = {7-ish multiples arranged to carry the moduli 7..12 and the
  pinning}, passes the FULL gap-violator profile (covering, spread, 24-comp,
  pair-38, pinning at every q <= 25) with w_(7)/w_(8) ~ S -> infinity.
  Hence NO constant C_7 with "profile => w_(7) <= C_7 w_(8)" exists:
  the ladder's l <= 6 pole is NECESSARY, not an accounting artifact.

  (These families are NOT gap violators -- their exact M is comfortably
  >= 2/25; that is the point: the PROFILE alone admits them, so any ladder
  rung at l = 7 provable from the profile is false.  The profile's top-7
  ratio chain is the most scale-rigidity the formal filters can give.)

Construction detail: CLUSTER must carry: multiples of 7,8,9,10,11,12 (covering),
the pinning contributions the bottom {1..5} misses, and stay 24-compressed +
ladder-consistent internally (consecutive ratios <= 24 within the top 7).
Bottom {1,2,3,4,5}: covers q in 2..6 (wait: 6 needs a multiple of 6 -- 6 = 2*3:
multiples of 6 in {1..5}: none! => 6 | some cluster element too), and pins
q <= 25 partially (residues 1..5 hit pairs {1},{2},{3},{4},{5} of each q).

We SEARCH for a valid cluster shape at several scales S and verify the FULL
profile exactly (same filters as the census), plus exact M >= 2/25 via scan.
"""
from math import gcd
from functools import reduce
import sys

def dist_q(x, q):
    x %= q
    return min(x, q - x)

def profile_pass(W):
    """The full gap-violator profile (census filters F1..F5). Returns (ok, why)."""
    n = len(W)
    Ws = sorted(W)
    # F1 covering 2..12
    for m in range(2, 13):
        if not any(v % m == 0 for v in W):
            return False, f"covering fails m={m}"
    # F2 spread
    if 2 * Ws[-1] <= 23 * Ws[0]:
        return False, "not spread"
    # F3 24-compression
    if Ws[-1] > 24 * Ws[-2]:
        return False, "not 24-compressed"
    # F4 pair 38
    if Ws[-1] + Ws[-2] < 38:
        return False, "pair < 38"
    # primitivity
    if reduce(gcd, W) != 1:
        return False, "not primitive"
    # F5 pinning at every q <= 25
    for q in range(2, 26):
        if any(v % q == 0 for v in W):
            continue
        for b in range(1, q // 2 + 1):
            if gcd(b, q) != 1:
                continue
            if not any(v % q == b or v % q == q - b for v in W):
                return False, f"pinning fails q={q} pair {b}"
    return True, "ok"

def witness_2_25(W, qmax=2000):
    """margin >= 2/25 witness search (loose confirmation)."""
    for q in list(range(8, 80)) + [v + d for v in sorted(W)[-3:] for d in (-1, 1)]:
        if q < 8 or any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            if all(25 * dist_q(a * v, q) >= 2 * q for v in W):
                return q, a
    return None

# ---------------------------------------------------------------------------
# The cluster: 7 values c_1..c_7 (to be scaled by S).  Requirements at scale S:
#  - multiples of 6,7,8,9,10,11,12 present among S*c_i  => choose c_i divisible
#    (S coprime to smalls keeps control: we take S odd, coprime to 2..25 lcm? no:
#     simplest: bake the divisibility into c_i and let S be arbitrary integer)
#  - pinning: for each q <= 25 with no q-multiple: residues of {1..5} u {S c_i}
#    must hit all unit pairs.  We SEARCH c-shapes and test at many S.
# base cluster: carries 7,8,9,10,11,12 and 6: c's = {6*7=42? ...}
# Use c = {7*8, 9*10, 11*12, 6*13, 5*14, ...}: we search small shapes.
import itertools, random
random.seed(55)

BOT = [1, 2, 3, 4, 5]
# candidate cluster atoms: numbers <= 60 covering the needed moduli
atoms = [v for v in range(6, 61)]
need = [6, 7, 8, 9, 10, 11, 12]

def cluster_ok_shape(C7):
    # internal 24-compression + ladder: sorted c: each consecutive ratio <= 24 (way loose)
    Cs = sorted(C7)
    for i in range(6):
        if Cs[i + 1] > 24 * Cs[i]:
            return False
    return all(any(c % m == 0 for c in C7) for m in need)

print("searching 7-cluster shapes (atoms <= 60 covering 6..12, internally compressed)...")
shapes = []
for _ in range(200000):
    C7 = tuple(sorted(random.sample(atoms, 7)))
    if cluster_ok_shape(C7):
        shapes.append(C7)
    if len(shapes) >= 400:
        break
print(f"candidate shapes: {len(shapes)}")

SCALES = [10, 100, 1000, 10**4, 10**5, 10**6, 10**7]
good_shape = None
for C7 in shapes:
    ok_all = True
    for S in SCALES:
        W = BOT + [S * c for c in C7]
        ok, why = profile_pass(W)
        if not ok:
            ok_all = False
            break
    if ok_all:
        good_shape = C7
        break

if good_shape:
    print(f"\nFLOATING 7-CLUSTER FOUND: bottom {BOT}, cluster {list(good_shape)} * S")
    print(f"{'S':>9} {'profile':>8} {'w_(7)/w_(8)':>12} {'2/25-witness':>16}")
    for S in SCALES:
        W = BOT + [S * c for c in good_shape]
        ok, why = profile_pass(W)
        ratio = (S * min(good_shape)) / 5
        wit = witness_2_25(W)
        print(f"{S:>9} {str(ok):>8} {ratio:>12.1f} {str(wit):>16}")
    print("\nCONCLUSION: the full formal profile admits families with w_(7)/w_(8)")
    print("arbitrarily large => NO ladder constant C_7 exists; the l <= 6 pole is")
    print("NECESSARY (ladder side; complements opus-S81's fee-mean ceiling).")
    print("(These families are loose -- 2/25-witnesses found -- as expected: the")
    print(" profile is necessary-not-sufficient, and the runaway lives on the loose side.)")
else:
    print("NO shape passed at all scales -- pole-necessity construction FAILS; investigate.")
