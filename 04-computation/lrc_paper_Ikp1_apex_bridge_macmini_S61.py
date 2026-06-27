#!/usr/bin/env python3
"""The paper's I(k,p,1) descent 7 -> 14 IS the project's covering condition (exact bridge).

Paper: Sungkawichai-Trakulthongchai, arXiv:2604.23906 ("Eleven, twelve, thirteen lonely runners",
proves LRC k<=12). A tuple v in Z^k has the LR property iff exists t with ||t v_i|| >= 1/(k+1) for all i.
I(k,p,l) = the (k,p,l)-IMPROPER tuples = no witness t in (1/(lp))Z (and gcd condition fails, vacuous at l=1).
Their composite-(k+1) fallback: lift at the prime factors c of k+1. For k=13, k+1=14=2*7 -> lifts c=7,2.

THIS SCRIPT verifies the exact dictionary (mac-mini-S61), extending kps-S31ag's mapping:
  (1) LEVEL 7  [ansatz (1/7)Z]:  no witness  <=>  some v_i == 0 (mod 7)   [= covering mod 7]
  (2) LIFT c=2 [ansatz (1/14)Z]: no witness  <=>  some v_i == 0 (mod 14)  [= project covering condition]
      (odd multiples of 7 are RESCUED by odd j/14; even multiples of 7 = multiples of 14 SURVIVE)
  (3) the (1/14)Z ansatz fails EXACTLY on covering tuples -> a real/finer witness is needed ->
      this is where the project's CONTINUOUS bound (p0<=cap, gK8) replaces the paper's enumeration.
  (4) small-k sanity: the improper FRACTION at ansatz (1/d)Z decreases as d refines (covering shrinks),
      the discrete shadow of the continuous lonely measure.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(11)

def ansatz_witness(v, d, kk1):
    """Exists t = j/d (j=0..d-1) with ||t v_i|| >= 1/kk1 for all i?  ||m/d||=min(m,d-m)/d."""
    for j in range(d):
        ok = True
        for vi in v:
            m = (j * vi) % d
            nd = min(m, d - m)
            if nd * kk1 < d:          # nd/d < 1/kk1
                ok = False; break
        if ok:
            return True
    return False

print("=" * 76)
print(" THE I(k,p,1) DESCENT 7 -> 14 = THE PROJECT'S COVERING CONDITION (k=13, k+1=14)")
print("=" * 76)
k, kk1 = 13, 14
N = 40000
bad7 = bad14 = 0
frac7_improper = frac14_improper = 0
for _ in range(N):
    v = [random.randrange(1, 14*6) for _ in range(k)]      # generic speeds, residues mod 14 generic
    # LEVEL 7
    w7 = ansatz_witness(v, 7, kk1)
    cov7 = any(x % 7 == 0 for x in v)
    if (not w7) != cov7: bad7 += 1
    if not w7: frac7_improper += 1
    # LIFT to 14
    w14 = ansatz_witness(v, 14, kk1)
    cov14 = any(x % 14 == 0 for x in v)
    if (not w14) != cov14: bad14 += 1
    if not w14: frac14_improper += 1
print(f"(1) LEVEL 7 :  [no (1/7)-witness  <=>  some v_i==0 mod 7 ]  counterexamples = {bad7}/{N}")
print(f"(2) LIFT 14 :  [no (1/14)-witness <=>  some v_i==0 mod 14]  counterexamples = {bad14}/{N}")
print(f"    improper fraction:  level-7 = {frac7_improper/N:.4f}   level-14 = {frac14_improper/N:.4f}")
print(f"    (theory: 1-(6/7)^13 = {1-(6/7)**13:.4f} ;  1-(13/14)^13 = {1-(13/14)**13:.4f})")

print("\n--- the RESCUE mechanism (odd mult of 7 vs mult of 14) ---")
for r in (7, 14):
    # a tuple whose ONLY non-generic coordinate is == r-residue; rest are units mod 14
    base = [1, 3, 5, 9, 11, 13, 1, 3, 5, 9, 11, 13]      # 12 units mod 14
    v = base + [r]                                         # 13th coord == 7 or 14 (mod 14)
    w14 = ansatz_witness(v, 14, kk1)
    print(f"   coord == {r} mod 14 (={'odd mult of 7' if r==7 else 'mult of 14'}): "
          f"has (1/14)-witness? {w14}  -> {'RESCUED (proper)' if w14 else 'SURVIVES improper = COVERING'}")

print("\n" + "=" * 76)
print(" (4) SMALL-k: improper fraction shrinks as the ansatz (1/d)Z refines (-> continuous)")
print("=" * 76)
for kk in (3, 4, 5):
    kk1s = kk + 1
    for d in (kk1s, 2*kk1s, 6*kk1s, 30*kk1s):
        imp = 0; M = 6000
        for _ in range(M):
            v = [random.randrange(1, 12*kk1s) for _ in range(kk)]
            if not ansatz_witness(v, d, kk1s): imp += 1
        print(f"   k={kk} (1/{kk1s} target), ansatz (1/{d})Z: improper fraction = {imp/M:.4f}")
    print()
print("INTERPRETATION: level-7 improper = covering mod 7 (coarse); one c=2 lift lands on covering mod 14")
print("= the project's exact covering condition. The (1/14)Z ansatz then FAILS on covering tuples, so a")
print("finer witness is needed -- where p0<=cap (gK8, continuous) REPLACES the paper's large-p enumeration.")
