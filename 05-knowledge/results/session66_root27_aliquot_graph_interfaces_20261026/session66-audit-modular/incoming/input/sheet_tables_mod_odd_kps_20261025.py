#!/usr/bin/env python3
"""
sheet_tables_mod_odd_kps_20261025.py  (kind-pasteur-2026-10-25-S4)

Modular edge tables of the two Collatz sheets.  For an odd modulus m and k >= 1 halvings define the
edge table of the plus sheet  E_k^+(m) = {(x mod m, (3x+1) 2^{-k} mod m) : x odd}  and of the minus
sheet  E_k^-(m) = {(x mod m, (3x-1) 2^{-k} mod m) : x odd}  (2 is invertible mod odd m; the halving
count k is not determined by x mod m, so the table is indexed by k).

Exact statements checked:
  T1  (affine step only) x -> 3x+1 and x -> 3x-1 are conjugate by the unit translation x -> x-1 over
      Z and hence over every Z/m; the translation does not commute with halving, so this is not a
      conjugacy of the Syracuse maps.
  T2  (negation, both coordinates) E_k^-(m) = {(-x, -y) : (x,y) in E_k^+(m)} for every odd m, every k:
      the sheet involution of THM-4474 (E) / THM-4610 read modulo m.
  T3  (alternation by a halving shift) if -1 lies in the subgroup <2> of (Z/m)^*, say 2^t = -1 (mod m)
      with t = ord_m(2)/2, then E_{k+t}^-(m) = {(-x, y) : (x,y) in E_k^+(m)}: the minus table with t
      more halvings is the plus table read from the negated source, with the SAME target residue.
      For m = 3: t = 1 (the tables alternate with the parity of k); m = 5: t = 2; m = 9: t = 3; m = 27: t = 9.
  T4  (hostiles) no halving shift relates the sheets when -1 is not a power of 2 mod m: m = 7 (ord 3),
      m = 15 (the shifts t = 1 mod 3 and t = 2 mod 5 are incompatible), 21, 23, 31, 45, 63.
"""
import sys
from math import gcd

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def edge_table(m, k, sign):
    inv2k = pow(pow(2, k, m), -1, m)
    S = set()
    for x in range(1, 2 * m, 2):          # odd residues x mod 2m determine x mod m over all odd classes
        y = ((3 * x + sign) * inv2k) % m
        S.add((x % m, y))
    return S


def order2(m):
    t, v = 1, 2 % m
    while v != 1:
        v = (v * 2) % m
        t += 1
    return t


def neg_one_power(m):
    """Return t with 2^t = -1 mod m, or None."""
    o = order2(m)
    if o % 2 == 0 and pow(2, o // 2, m) == m - 1:
        return o // 2
    return None


# T1: translation conjugacy of the affine step, all m <= 60; and failure for the Syracuse step
ok = True
for m in range(2, 61):
    for x in range(m):
        ok &= ((3 * ((x - 1) % m) + 1 + 1) % m == (3 * x - 1) % m)   # phi^{-1} f phi with phi(x) = x - 1
check("T1 affine step: x->3x-1 equals (x->x+1) o (x->3x+1) o (x->x-1) over Z/m for all m <= 60", ok)
# the translation does not conjugate the Syracuse maps: S^-(x) vs S^+(x-1)+1 on odd x
def S(x, sign):
    y = 3 * x + sign
    while y % 2 == 0:
        y //= 2
    return y
mism = sum(1 for x in range(1, 2001, 2) if (x - 1) % 2 == 1 and S(x, -1) != S(x - 1, 1) + 1)
mism2 = sum(1 for x in range(3, 2001, 2) if S(x, -1) != S(x - 2, 1) + 2)
check("T1b the unit translation is not a conjugacy of the Syracuse maps (counterexamples exist among odd x < 2000)", mism2 > 0, mism2)

# T2: negation conjugacy of the full tables, odd m <= 101, k <= 6
ok = True
for m in range(3, 102, 2):
    for k in range(1, 7):
        Ep = edge_table(m, k, +1)
        Em = edge_table(m, k, -1)
        ok &= (Em == {((-x) % m, (-y) % m) for (x, y) in Ep})
check("T2 E_k^-(m) = negation of E_k^+(m) in both coordinates, all odd m <= 101, k = 1..6", ok)

# T3: alternation by the halving shift t when 2^t = -1
alts = {}
ok = True
for m in range(3, 102, 2):
    t = neg_one_power(m)
    if t is None:
        continue
    alts[m] = t
    for k in range(1, 7):
        Ep = edge_table(m, k, +1)
        Em = edge_table(m, k + t, -1)
        ok &= (Em == {((-x) % m, y) for (x, y) in Ep})
check("T3 E_{k+t}^-(m) = E_k^+(m) with the source negated and the target kept, for every odd m <= 101 with 2^t = -1", ok)
print("odd m <= 101 with -1 in <2> and the shift t:", alts)
check("T3b m = 3, 5, 9, 27 have t = 1, 2, 3, 9", alts.get(3) == 1 and alts.get(5) == 2 and alts.get(9) == 3 and alts.get(27) == 9)
check("T3c mod 3 the plus table with k halvings lands in residue 2 for odd k and 1 for even k, the minus table the reverse",
      all(edge_table(3, k, +1) == {(1, 2 if k % 2 else 1), (0, 2 if k % 2 else 1), (2, 2 if k % 2 else 1)} and
          edge_table(3, k, -1) == {(1, 1 if k % 2 else 2), (0, 1 if k % 2 else 2), (2, 1 if k % 2 else 2)} for k in range(1, 7)))

# T4: hostiles — no halving shift relates the sheets when -1 is not a power of 2
ok = True
for m in (7, 15, 21, 23, 31, 45, 63):
    assert neg_one_power(m) is None
    o = order2(m)
    found = False
    for k in range(1, 7):
        Ep = edge_table(m, k, +1)
        for sh in range(0, o):
            Em = edge_table(m, k + sh, -1)
            if Em == {((-x) % m, y) for (x, y) in Ep} or Em == Ep:
                found = True
    ok &= (not found)
check("T4 hostiles m = 7, 15, 21, 23, 31, 45, 63: no halving shift makes the minus table the (source-negated) plus table", ok)
print("ord_m(2) for m = 3,5,7,9,15,27:", {m: order2(m) for m in (3, 5, 7, 9, 15, 27)})

bad = [n for n, ok in CHECKS if not ok]
print("%d checks, %d failed" % (len(CHECKS), len(bad)))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")
