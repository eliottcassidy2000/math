#!/usr/bin/env python3
"""
lrc14_angleE_bonferroni_prefix_macmini_0620s8.py  (mac-mini-2026-06-20-S8)  ANGLE E part 3

BREAKTHROUGH SIGNAL (part 2/quick test): the CUMULATIVE Bonferroni prefix of the signed sum
  P_L(E) = Sum_{M subset {1..6}, |M|<=L} (-1)^|M| J(M,E)
is consec-EXTREMAL (max) at EVERY level L in {1,2,4,5,6} -- only L=3 fails (14/330 at k=8).

This is the irreducibly-aggregate C1 path made ALMOST term-group-wise: the prefix telescoping
works at 5 of 6 cut points.  measS7 = P_6, P_5, ... are nested.

This script:
  (A) reconfirms prefix-extremality per L for k=8,9,10 over the certified box.
  (B) isolates the L=3 violators and tests REPAIR: is the level-{3,4} JOINT step extremal?
      i.e. is P_4 = P_2 + (level3 + level4) extremal AND is (level3+level4) alone extremal?
      Also test pairing (1,2),(3,4),(5,6) -- the "even-closure" Bonferroni pairs.
  (C) the KEY structural test: are the EVEN prefixes P_0,P_2,P_4,P_6 a DECREASING chain of
      UPPER bounds, each consec-extremal, converging down to measS7?  And ODD prefixes
      P_1,P_3,P_5 increasing LOWER bounds?  (classical Bonferroni sandwich.)  If every EVEN
      prefix is a consec-extremal UPPER bound and they decrease to measS7=P_6, that ALONE
      doesn't prove P_6 extremal; but if P_4 (consec-extremal upper bd) >= measS7(adv) for the
      RIGHT reason... we need P_6 itself.  The clean statement: P_6 extremal needs all 6 levels.
      BUT: pairing levels into (2j-1,2j) "Kounias/Hunter" blocks may give monotone consec-extremal
      block sums.  Test block sums b_j = (level 2j-1)+(level 2j) signed contributions.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def J(E, M):
    Mset = set(M)
    if 0 in Mset: return F(0)
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; ok = True
        for e in E:
            if e == 0: continue
            if int(((e * xm) % 1) * 7) in Mset: ok = False; break
        if ok: tot += x1 - x0
    return tot

def levels(E):
    """signed contribution of each level L=0..6: c_L = sum_{|M|=L} (-1)^L J(M)."""
    c = [F(0)] * 7
    for L in range(7):
        for M in itertools.combinations(range(1, 7), L):
            c[L] += F((-1) ** L) * J(E, M)
    return c

def main():
    boxes = {8: 11, 9: 12, 10: 13}  # span box; certified region per crux assembly (max(E)<=13 overall)
    for k in [8, 9, 10]:
        N = boxes[k]
        consec = list(range(k))
        cc = levels(consec)
        Pc = [sum(cc[:L + 1]) for L in range(7)]
        print("=" * 96)
        print(f"k={k}, box 0..{N}, measS7(consec)=P6={Pc[6]}={float(Pc[6]):.6f}")
        print("  consec level contributions c_L: " + ", ".join(f"c{L}={float(cc[L]):+.4f}" for L in range(7)))
        print("  consec prefixes  P_L:           " + ", ".join(f"P{L}={float(Pc[L]):+.4f}" for L in range(7)))

        pool = [(0,) + c for c in itertools.combinations(range(1, N + 1), k - 1)]
        print(f"  pool size: {len(pool)}")

        # (A) prefix extremality per L
        prefix_fail = {L: 0 for L in range(1, 7)}
        l3_violators = []
        # (B) block sums b_j = c_{2j-1}+c_{2j}: b1=c1+c2, b2=c3+c4, b3=c5+c6
        bc = [cc[1] + cc[2], cc[3] + cc[4], cc[5] + cc[6]]
        block_fail = [0, 0, 0]
        # joint level (3,4)
        lvl34_c = cc[3] + cc[4]; lvl34_fail = 0
        for E in pool:
            ce = levels(E)
            Pe = [sum(ce[:L + 1]) for L in range(7)]
            for L in range(1, 7):
                if Pe[L] > Pc[L] + F(0):
                    prefix_fail[L] += 1
                    if L == 3: l3_violators.append((E, Pe[3] - Pc[3]))
            be = [ce[1] + ce[2], ce[3] + ce[4], ce[5] + ce[6]]
            for j in range(3):
                if be[j] > bc[j]: block_fail[j] += 1
            if ce[3] + ce[4] > lvl34_c: lvl34_fail += 1

        print("\n  (A) cumulative prefix extremality (adversaries beating consec):")
        for L in range(1, 7):
            tag = "  <-- ONLY failure" if prefix_fail[L] > 0 else ""
            print(f"      P_{L}: {prefix_fail[L]}/{len(pool)}{tag}")

        print("\n  (B) Kounias/Hunter BLOCK sums b_j=c_{2j-1}+c_{2j} extremality:")
        print(f"      consec blocks: b1=c1+c2={float(bc[0]):+.4f}, b2=c3+c4={float(bc[1]):+.4f}, b3=c5+c6={float(bc[2]):+.4f}")
        for j in range(3):
            tag = "  <-- FAIL" if block_fail[j] > 0 else "  OK"
            print(f"      b{j+1}: {block_fail[j]}/{len(pool)}{tag}")
        print(f"      joint level (3+4) c3+c4={float(lvl34_c):+.4f}: {lvl34_fail}/{len(pool)} beat  {'<-- FAIL' if lvl34_fail else 'OK'}")

        if l3_violators:
            l3_violators.sort(key=lambda t: -t[1])
            print(f"\n  L=3 prefix violators (top 5 of {len(l3_violators)}):")
            for E, d in l3_violators[:5]:
                print(f"      E={E}  P3 excess={float(d):+.6f}")

if __name__ == "__main__":
    main()
