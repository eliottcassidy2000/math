#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHARP FREIMAN DICHOTOMY for LRC(14) (opus-2026-06-19, second pass).

First pass (lrc14_freiman_partition_rigor) showed the naive GAP detector absorbs
everything into category (a). That is uninformative because it allows GAP size up
to 8k. Here we SHARPEN:

  * For each primitive non-AP E (|E|=k), compute the SMALLEST-size proper 2-or-3
    dim GAP that contains E, and the ratio rho_GAP = minGAPsize / k.
  * Compute doubling sigma = |E+E|/k.
  * A set is DISSOCIATED-type if it has a Sidon-like stranger: some w in E whose
    pairwise sums with E are all distinct from the rest, formalized via the
    "no bounded relation" test at height H.

KEY QUESTIONS:
  (Q1) Plemann/Freiman: does small sigma  <=> small min-GAP-size ratio?
       Find the empirical curve sigma -> min achievable GAP-size/k.
  (Q2) Is there ANY set with small sigma but NO GAP of size <= C*k for moderate C?
       (the genuine third pocket). Push C down and see where sets fall out.
  (Q3) For the LRC margin: ALL non-AP sets have L_y well below cap. Confirm the
       L_y vs (sigma, GAP-dim) relationship: higher dim / higher sigma => lower L_y.

The point for the proof program step 3(a)/3(b):
  - If small-doubling => bounded GAP (dim>=2) with explicit constant C, then 3(a)
    (dimension penalty) applies and L_y << cap.
  - If large-doubling => dissociated stranger => contraction (HYP-2610) -> 3(b).
  - The threshold sigma* is where GAP-size ratio explodes.
"""
import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            val = Fraction((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = Fraction(-(t-2)*(t-3)*(t-6), 36)
        else:
            val = Fraction((t-3)*(t-4), 12)
        g.append(val)
    return g

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7)), p

CAP = {8: 0.38153, 9: 0.49426, 10: 0.6044}

def doubling_size(E):
    E = sorted(set(E))
    return len({a + b for a in E for b in E})

def is_full_AP(E):
    E = sorted(set(E))
    if len(E) < 2: return True
    d = E[1] - E[0]
    return all(E[i+1] - E[i] == d for i in range(len(E)-1))

def primitive(E):
    E = sorted(set(E))
    return len(E) > 0 and E[0] == 0 and reduce(gcd, E) == 1

def min_gap2_size(E):
    """
    Smallest size of a PROPER 2-dim GAP {x*q1 + y*q2 : 0<=x<l1, 0<=y<l2} containing
    E (after translating min to 0). Proper = all l1*l2 lattice points distinct.
    Returns (size, q1,q2,l1,l2) or None if no proper 2-dim GAP of size <= span+1 found
    (full AP cover is the dim-1 fallback handled separately).
    """
    E = sorted(set(E))
    E = [e - E[0] for e in E]
    n = len(E)
    span = E[-1]
    best = None
    # candidate generators: all differences
    diffs = sorted({E[j]-E[i] for i in range(n) for j in range(i+1, n)})
    for q1 in diffs:
        if q1 <= 0: continue
        # reduce mod q1: residues = q2-direction candidates
        res = sorted({e % q1 for e in E})
        quo = sorted({e // q1 for e in E})
        # rows (quo) must be an AP; cols (res) must be an AP
        def ap_info(vals):
            vals = sorted(set(vals))
            if len(vals) == 1: return (1, 1)
            ds = [vals[i+1]-vals[i] for i in range(len(vals)-1)]
            g = reduce(gcd, ds)
            return (g, (vals[-1]-vals[0])//g + 1)
        rstep, rlen = ap_info(quo)
        cstep, clen = ap_info(res)
        size = rlen * clen
        # properness: q1 (the row jump * rstep) and cstep generate distinct points.
        # 2-dim genuine: both lens >= 2.
        if rlen >= 2 and clen >= 2:
            # verify it actually covers E and is proper
            q1eff = q1 * rstep
            q2eff = cstep
            pts = set()
            ok = True
            for x in range(rlen):
                for y in range(clen):
                    p = x * q1eff + y * q2eff
                    if p in pts:
                        ok = False; break
                    pts.add(p)
                if not ok: break
            if ok and all(e in pts for e in E):
                if best is None or size < best[0]:
                    best = (size, q1eff, q2eff, rlen, clen)
    return best

def has_bounded_relation_for_each(E, H, maxsupp=3):
    """
    A 'dissociated stranger' w in E\{0} has NO bounded integer relation
        c*w + sum_{v in S} m_v v = 0,   1<=|c|<=H, 1<=|m_v|<=H, S subset E\{w},
        |S| <= maxsupp.
    (A genuine Freiman relation has small support; large-support relations are
    not the obstruction we test. maxsupp=3 catches all single/pair/triple ties,
    which is the relevant 'low-height relation' notion for primitive k=8 sets.)
    Returns the list of strangers (empty list = every element bounded-related).
    Fast: enumerate small-support subsets and bounded coeffs, no reachable set.
    """
    E = sorted(set(E))
    others_all = [v for v in E if v != 0]
    strangers = []
    coeffs = [m for m in range(-H, H + 1) if m != 0]
    for w in E:
        if w == 0: continue
        rest = [v for v in E if v != w and v != 0]
        found = False
        # support size 1,2,3 among rest
        for r in range(1, maxsupp + 1):
            if found: break
            for S in itertools.combinations(rest, r):
                if found: break
                for combo in itertools.product(coeffs, repeat=r):
                    base = sum(m * v for m, v in zip(combo, S))
                    # need c*w = -base for some c in [1,H] (or [-H,-1])
                    for c in range(1, H + 1):
                        if c * w == -base or -c * w == -base:
                            found = True; break
                    if found: break
        if not found:
            strangers.append(w)
    return strangers

def main():
    print("=" * 78)
    print("SHARP FREIMAN DICHOTOMY — LRC(14) (opus-2026-06-19)")
    print("=" * 78)

    random.seed(7771)

    # ---- generate a broad k=8 corpus including HIGH-doubling sets ----
    corpus = set()
    # random moderate spread
    for _ in range(5000):
        sp = random.randint(8, 40)
        E = sorted(random.sample(range(sp + 1), 8))
        E = tuple(e - E[0] for e in E)
        corpus.add(E)
    # Sidon-like / large-doubling: geometric-ish, powers, random sparse large spread
    for _ in range(3000):
        sp = random.randint(40, 90)
        E = sorted(random.sample(range(sp + 1), 8))
        E = tuple(e - E[0] for e in E)
        corpus.add(E)
    # explicit Sidon set (B_2): 0,1,3,7,12,20,30,44 (Mian-Chowla)
    corpus.add(tuple([0,1,3,7,12,20,30,44]))
    # powers of 2
    corpus.add(tuple([0,1,2,4,8,16,32,64]))
    # 2D GAPs (small doubling, dim 2)
    for q in range(3, 12):
        for l1 in range(2, 5):
            for l2 in range(2, 5):
                if l1*l2 != 8 and l1*l2 < 8: continue
                P = sorted({a*q + b for a in range(l1) for b in range(l2)})
                if len(P) >= 8:
                    corpus.add(tuple(P[:8]))

    rows = []
    for E in corpus:
        E = list(E)
        if not primitive(E) or is_full_AP(E): continue
        sig = Fraction(doubling_size(E), 8)
        gap = min_gap2_size(E)
        gsize = gap[0] if gap else None
        strg = has_bounded_relation_for_each(E, 2)
        Lv, _ = L_y(E, 8)
        rows.append((E, float(sig), gsize, len(strg), float(Lv)))

    print(f"\n[corpus] {len(rows)} primitive non-AP k=8 sets analyzed")

    # Q1: sigma vs min GAP-size ratio
    print("\n[Q1] sigma vs smallest proper 2-dim GAP size (rho = gapsize/k):")
    print("     sigma-bin | #sets | #with-2dGAP | median gapsize | max gapsize | #strangers(H2)")
    bins = {}
    for E, sig, gsize, ns, Lv in rows:
        b = round(sig * 2) / 2  # 0.5 bins
        bins.setdefault(b, []).append((gsize, ns, Lv))
    for b in sorted(bins):
        vals = bins[b]
        withgap = [g for g, ns, L in vals if g is not None]
        nstr = sum(1 for g, ns, L in vals if ns > 0)
        med = sorted(withgap)[len(withgap)//2] if withgap else None
        mx = max(withgap) if withgap else None
        print(f"     sig~{b:4.1f} | {len(vals):5d} | {len(withgap):11d} | {str(med):14s} | {str(mx):11s} | {nstr}")

    # Q1b: BOUNDED-GAP dichotomy. Fix C; a set is 'GAP-bounded' if min 2-dim GAP
    # size <= C*k. Otherwise it should have a stranger. Sweep C and check.
    print("\n[Q1b] BOUNDED-GAP vs STRANGER dichotomy (the real cut). k=8.")
    print("     For threshold C: GAP-bounded = min2dGAPsize<=C*8. Else expect stranger.")
    for C in [2, 3, 4, 5, 6, 8, 10]:
        Ck = C * 8
        gapbounded = [r for r in rows if r[2] is not None and r[2] <= Ck]
        notbounded = [r for r in rows if not (r[2] is not None and r[2] <= Ck)]
        nb_no_strg = [r for r in notbounded if r[3] == 0]
        maxL_nb = max((r[4] for r in notbounded), default=0.0)
        print(f"     C={C:2d} (size<= {Ck:3d}): GAP-bounded={len(gapbounded):5d}  not-bounded={len(notbounded):5d}  "
              f"of-not-bounded-with-NO-stranger={len(nb_no_strg):5d}  maxL_y(not-bounded)={maxL_nb:.5f}")

    # Q2: third pocket = small sigma, NO 2dGAP, AND no stranger
    print("\n[Q2] THIRD POCKET search (small sigma & no 2-dim GAP & no H2 stranger):")
    SIG_SMALL = 4.0  # 'small doubling' threshold to probe
    third = [(E, sig, gsize, ns, Lv) for E, sig, gsize, ns, Lv in rows
             if sig <= SIG_SMALL and gsize is None and ns == 0]
    print(f"     sets with sigma<={SIG_SMALL}, no 2dGAP, no stranger: {len(third)}")
    for E, sig, gsize, ns, Lv in third[:30]:
        print(f"       {E} sigma={sig:.3f} L_y={Lv:.5f}")
    # also: small sigma & no 2dGAP (regardless of stranger)
    nogap_smallsig = [(E, sig, ns, Lv) for E, sig, gsize, ns, Lv in rows
                      if sig <= SIG_SMALL and gsize is None]
    print(f"     sets with sigma<={SIG_SMALL} and no detected proper 2-dim GAP: {len(nogap_smallsig)}")
    for E, sig, ns, Lv in nogap_smallsig[:20]:
        print(f"       {E} sigma={sig:.3f} strangers={ns} L_y={Lv:.5f}")

    # Q2b: of those with no 2dGAP, do they have a stranger? (dichotomy completeness)
    no2d = [(E, sig, ns, Lv) for E, sig, gsize, ns, Lv in rows if gsize is None]
    no2d_nostr = [t for t in no2d if t[2] == 0]
    print(f"\n[Q2b] DICHOTOMY COMPLETENESS over ALL sets:")
    print(f"      sets with no proper 2-dim GAP found      : {len(no2d)}")
    print(f"      of those, also NO H2 stranger (residue)  : {len(no2d_nostr)}")
    if no2d_nostr:
        print(f"      residue examples (re-test with H=3):")
        still = []
        for E, sig, ns, Lv in no2d_nostr[:50]:
            s3 = has_bounded_relation_for_each(E, 3)
            tag = "STILL RESIDUE" if len(s3) == 0 else f"stranger@H3:{s3}"
            if len(s3) == 0: still.append(E)
            if no2d_nostr.index((E,sig,ns,Lv)) < 15:
                print(f"        {E} sigma={sig:.3f} -> {tag}")
        print(f"      residue surviving H=3: {len(still)}")

    # Q3: L_y vs sigma (dimension penalty visible)
    print("\n[Q3] L_y by sigma bin (dimension penalty: higher sigma -> lower L_y):")
    print("     sigma-bin | #sets | max L_y | mean L_y | cap_8=0.38153")
    for b in sorted(bins):
        Ls = [L for g, ns, L in bins[b]]
        print(f"     sig~{b:4.1f} | {len(Ls):5d} | {max(Ls):.5f} | {sum(Ls)/len(Ls):.5f}")
    globalmaxL = max(r[4] for r in rows)
    argm = max(rows, key=lambda r: r[4])
    print(f"     GLOBAL max L_y over non-AP = {globalmaxL:.5f} at {argm[0]} (sigma={argm[1]:.3f})")
    print(f"     margin to cap_8 = {CAP[8]-globalmaxL:.5f}")

if __name__ == "__main__":
    main()
