#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_lattice_rearrangement_opus_0620.py   (opus-2026-06-20, THREAD D part iii)

After (D-i)(D-ii) refuted the gap-multiset Schur order and the per-S permanent linearity,
this is the CONSTRUCTIVE refinement: measS7(E) = iid_k + corr(E) and corr(E) lives ENTIRELY
in the relation lattice Lambda(E) = { n in Z^E : sum_i n_i e_i = 0 }.  The Fourier identity
(kps-S6/wide-bound):

   corr(E) = Sum_{S subset Z/7, 0 notin S}  (-1)^|S|  Sum_{0 != n in Lambda(E)}  prod_i chat_S(n_i)

where chat_S(0)=1-|S|/7 and chat_S(7m)=0 (apex-prime vanishing), so ONLY relation vectors n
with NO coordinate a nonzero multiple of 7 contribute (call these EFFECTIVE relations).

THE REARRANGEMENT QUESTION (the right one): rank shapes by a relation-lattice MONOTONE.
 (L1) Verify the Fourier corr identity reproduces measS7 exactly (closes the bridge in code).
 (L2) Define LATTICE-ENERGY functionals and test which one consec MAXIMIZES over the pool:
        - eff_count(E)  = # effective short relation vectors (||n||_1 <= T) -- proxy for
          "how many low-height resonances".  consec packs the MOST low-height relations.
        - the leading 1-norm shell: the n with ||n||_1 = 2 are the pair relations e_i = e_j
          (none, distinct) and e_i + e_j = e_l etc.  consec has e_a+e_b=e_{a+b}: a TON of
          additive relations.  Count additive triples a+b=c within E.
 (L3) THE CLEAN INEQUALITY CANDIDATE: corr(E) is dominated by the leading lattice shell.
      Test  measS7(E) <= iid_k + Phi(N_2(E))  where N_2(E)=# additive triples (Schur-like
      energy) and Phi is the per-triple contribution at consec.  Is N_2 maximized at consec?
      Is corr monotone in N_2 across the pool (a genuine 1-D rearrangement certificate)?
 (L4) HONEST refutation check: if corr is NOT monotone in any single lattice-energy
      coordinate, report the scatter (corr vs N_2) and the rank correlation, and the
      decisive non-monotone pair.  A clean "monotone in additive-triple count" would be a
      WORKING extremality inequality; a clean scatter is the refutation.

EXACT Fractions; stdlib only.  Restrict to a feasible pool (k=8, span<=11; k=9 span<=11).
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

def cell_decomp(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hits = {0}
        for e in E:
            if e == 0: continue
            hits.add(int(((e * xm) % 1) * 7))
        cells.append((x1 - x0, frozenset(hits)))
    return cells

def measS7(E):
    return sum(L for L, h in cell_decomp(E) if len(h) == 7)

def stirling2(n, k):
    return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)
def measS7_iid(k):
    return F(factorial(7)*stirling2(k,7), 7**k)

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def shapes(k, span):
    return [(0,)+combo for combo in itertools.combinations(range(1, span+1), k-1)]

# ----- lattice-energy proxies (no full lattice enumeration; use structural counts) -----
def additive_triples(E):
    """# ordered (a,b,c) with a,b,c in E, a+b=c, a<=b (unordered pairs). The dominant
    rank-2 relation shell e_a + e_b - e_c = 0.  consec is RICH in these."""
    Es = set(E); cnt = 0
    El = sorted(E)
    for i in range(len(El)):
        for j in range(i, len(El)):
            if El[i]+El[j] in Es:
                cnt += 1
    return cnt

def difference_collisions(E):
    """# distinct difference VALUES with multiplicity (Sidon-defect): sum over d of C(rep_d,2).
    consec has the MOST repeated differences (anti-Sidon)."""
    El = sorted(E); reps = defaultdict(int)
    for i in range(len(El)):
        for j in range(i+1, len(El)):
            reps[El[j]-El[i]] += 1
    return sum(r*(r-1)//2 for r in reps.values())

def main():
    print("="*92)
    print("THREAD D (part iii): relation-lattice REARRANGEMENT -- the right coordinate for corr(E)")
    print("="*92)

    CONFIG = {8: 11, 9: 11}  # feasible boxes

    for k in (8, 9):
        span = CONFIG[k]
        consec = tuple(range(k))
        iid = measS7_iid(k)
        cap = CAP[k]
        pool = shapes(k, span)
        print("\n" + "#"*92)
        print(f"### k={k}, span box [0..{span}], |pool|={len(pool)}, iid={float(iid):.5f}, cap={float(cap):.4f}")
        print("#"*92)

        data = []
        for E in pool:
            m = measS7(E)
            data.append((E, m, m-iid, additive_triples(E), difference_collisions(E)))

        cons = next(d for d in data if d[0]==consec)
        print(f"\n[L2] consec: measS7={float(cons[1]):.5f} corr=+{float(cons[2]):.5f} "
              f"add-triples={cons[3]} diff-collisions={cons[4]}")
        # is consec the max of each lattice-energy proxy?
        mx_at = max(data, key=lambda d:d[3]); mc_at = max(data, key=lambda d:d[4])
        mxc = max(data, key=lambda d:d[1])
        print(f"     max add-triples: {mx_at[3]} at {mx_at[0]} (consec? {mx_at[0]==consec})")
        print(f"     max diff-collisions: {mc_at[4]} at {mc_at[0]} (consec? {mc_at[0]==consec})")
        print(f"     max measS7: {float(mxc[1]):.5f} at {mxc[0]} (consec? {mxc[0]==consec})")

        # [L3]/[L4] monotonicity of corr in the lattice-energy coordinates
        # sort by add-triples, see if measS7 is (weakly) monotone increasing
        def mono_check(idx, label):
            srt = sorted(data, key=lambda d:(d[idx], d[1]))
            inv = 0; pairs = 0; worst = None
            # count inversions: higher energy but LOWER measS7
            for a in range(len(srt)):
                for b in range(a+1, len(srt)):
                    da, db = srt[a], srt[b]
                    if da[idx] < db[idx]:
                        pairs += 1
                        if da[1] > db[1]:  # lower energy but higher measS7 = inversion
                            inv += 1
                            d = da[1]-db[1]
                            if worst is None or d>worst[0]: worst=(d,da,db)
            rho_num = pairs - 2*inv
            print(f"\n[L3:{label}] strictly-comparable pairs={pairs}, MONOTONE-violations(inversions)={inv} "
                  f"({100*inv/max(pairs,1):.1f}%)")
            if worst:
                d,da,db = worst
                print(f"     worst inversion: {label}={da[idx]} measS7={float(da[1]):.4f}  >  "
                      f"{label}={db[idx]} measS7={float(db[1]):.4f}  (lower energy, higher measS7 by +{float(d):.4f})")
            return inv, pairs
        mono_check(3, "add-triples")
        mono_check(4, "diff-collisions")

    print("\n" + "="*92)
    print("READOUT (Thread D part iii)")
    print("="*92)
    print("If consec maximizes a lattice-energy proxy AND corr is monotone in it -> WORKING")
    print("rearrangement inequality. If monotone fails -> the relation lattice is also the wrong")
    print("single coordinate; the extremality is genuinely a MULTI-shell aggregate (matches the")
    print("CLAUDE.md 'consec-max is irreducibly aggregate' note).")

if __name__ == "__main__":
    main()
