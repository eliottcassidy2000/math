"""
lrc_exact_Aprime_exhaustive_opus_S136.py

EXHAUSTIVE EXACT verification of (A') per k, with certified tight loci and exact
isolation gaps — powered by the order-cell exact-mu engine (lrc_exact_mu_ordercells).

(A') per k: mu_{1/7}(E) >= mu_{1/7}(AP_k) for every k-element integer set E, with
equality iff E is affine-AP.  This is the load-bearing open lemma of the density floor
(monad HYP-4787).  Grids can verify ">=" to ~4 digits but can NEVER certify the equality
case or the isolation gap; exact rationals can.

For each k we sweep ALL normalized (min=1, gcd of differences=1) k-subsets of {1..H_k}
and report:  #below bar (expect 0) / equality classes (expect only AP_k) /
the exact runner-up value and ISOLATION GAP gamma_k = runnerup - bar.

Boxes chosen for ~minutes each: k=8:H=18, k=9:H=17, k=10:H=17, k=11:H=16, k=12:H=16,
k=13:H=17.  (Beyond the box, kps-S60's intersected diameter ledger + the diameter floor
cover the large-diameter regimes at the DAG bars; this run targets the exact-extremality
statement itself on the near-AP boxes where the tight locus lives.)
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import time, sys

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import mu_exact, normalize

BARS = {8: F(691,735), 9: F(247,294), 10: F(38,49),
        11: F(1381,2205), 12: F(13823,24255), 13: F(477,1078)}
BOXES = {8: 18, 9: 17, 10: 17, 11: 16, 12: 16, 13: 17}

def main():
    print("="*100)
    print("EXHAUSTIVE EXACT (A') per k: 0-below check + certified equality classes + isolation gaps")
    print("="*100)
    grand = True
    for k in range(8, 14):
        H = BOXES[k]
        bar = BARS[k]
        t0 = time.time()
        seen = set()
        below = []
        eqcls = []
        runner = None
        for comb in combinations(range(1, H+1), k):
            En = normalize(comb)
            if En in seen: continue
            seen.add(En)
            m = mu_exact(list(En))
            if m < bar:
                below.append((En, m))
            elif m == bar:
                eqcls.append(En)
            else:
                if runner is None or m < runner[1]:
                    runner = (En, m)
        ok = (len(below) == 0)
        grand = grand and ok
        gap = runner[1] - bar if runner else None
        print(f"\n k={k}  box {{1..{H}}}  classes={len(seen)}  bar={bar}")
        print(f"   below bar: {len(below)}   {'(A) VERIFIED EXACT on box' if ok else '*** VIOLATION ***'}")
        for En, m in below[:4]:
            print(f"      VIOLATOR {En}: {m}")
        print(f"   equality classes: {len(eqcls)} -> {eqcls[:4]}")
        if runner:
            print(f"   runner-up: {runner[0]}  mu = {runner[1]} = {float(runner[1]):.6f}")
            print(f"   ISOLATION GAP gamma_{k} = {gap} = {float(gap):.6f}")
        print(f"   [{time.time()-t0:.0f}s]")
    print("\n" + "="*100)
    print(f"GRAND VERDICT: (A') exact on all boxes: {grand}; tight locus certified per k above.")

if __name__ == "__main__":
    main()
