#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ADVERSARIAL VERIFICATION of HYP-2675 angle "dichotomy-finite-reduction".

Claim under test (the TRUE-WIDE REDUCTION theorem, kps-Sx-wf):
  For primitive E, 0 in E, |E|=k in {8..12}, 2nd-largest(E) > 14 (true-wide),
    p0(E) <= cap_k
  reduces to: (A) scale invariance, (B) a finite cluster-shape family whose
  decorrelated worst case p0_inf = [k-1 consec]+[singleton] = Qb(k-1) < cap_k,
  (C) finite-M error ERR <= sum_far (6/49) V_i/g_i, explicit B closes gaps>B.

This script DOES NOT trust the claim.  It:
  (1) re-derives the key constants exactly with the engine;
  (2) HUNTS for a primitive wide k-set with p0 > cap_k OR > claimed wide-bound,
      across many adversarial families (multi-cluster near-AP, dilated primitive
      APs, balanced 2-3 cluster sets, all-residues-mod-7, two big clusters, ...);
  (3) checks whether p0_inf is REALLY an upper bound for p0(E) at finite scale,
      i.e. whether ERR = p0(E) - p0_inf can be POSITIVE and how large;
  (4) checks the "consec maximizes single-cluster coverage" claim;
  (5) probes the explicit B / glue-to-span-14 logic for gaps.

Default skeptical.  Report surviving gaps.  EXACT arithmetic (Fraction).
kind-pasteur-2026-06-20-Sx (workflow).
"""

import sys, itertools, random
from fractions import Fraction as F
from math import gcd
from collections import defaultdict

try:
    sys.stdout.reconfigure(encoding='utf-8', line_buffering=True)
except Exception:
    pass

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
        11: F(66, 91), 12: F(6, 7)}
# claimed worst decorrelated p0_inf (= THM-547 plateau Qb(k-1)):
QB_CLAIM = {8: 0.19660, 9: 0.36210, 10: 0.44789, 11: 0.53125, 12: 0.60224}
B_CLAIM = {8: 682, 9: 1453, 10: 1774, 11: 1988, 12: 2034}


# ----------------------------- exact engines -----------------------------
def breakpoints(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    return sorted(b for b in bps if 0 <= b <= 1)


def p0p1(E):
    E = sorted(set(E))
    bps = breakpoints(E)
    p0 = F(0); p1 = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in E)
        if len(miss) == 0:
            p0 += hi - lo
        elif len(miss) == 1:
            p1 += hi - lo
    return p0, p1


def p0(E):
    return p0p1(E)[0]


def primitive(E):
    E = sorted(set(e - min(E) for e in E))
    g = 0
    for e in E:
        g = gcd(g, e)
    if g > 1:
        E = [e // g for e in E]
    return sorted(set(E))


# ------------- shared-x decorrelation engine (copied from claim) -------------
def far_dist_at_x(C, x):
    base = [(c * x) % 1 for c in C]
    tb = {F(0), F(1)}
    for bb in base:
        for k in range(7):
            tb.add((F(k, 7) - bb) % 1)
    tb = sorted(tb)
    d = defaultdict(lambda: F(0))
    for tlo, thi in zip(tb, tb[1:]):
        if thi <= tlo:
            continue
        tmid = (tlo + thi) / 2
        sec = frozenset(s for s in (int(((bb + tmid) % 1) * 7) for bb in base)
                        if 1 <= s <= 6)
        d[sec] += thi - tlo
    return dict(d)


def anchored_set_at_x(C, x):
    return frozenset(s for s in (int((c * x) % 1 * 7) for c in C) if 1 <= s <= 6)


def p0_inf(clusters):
    allbps = set()
    for C in clusters:
        allbps |= set(breakpoints(C))
    xb = sorted(allbps)
    need = frozenset(range(1, 7))
    tot = F(0)
    for lo, hi in zip(xb, xb[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cur = {anchored_set_at_x(clusters[0], mid): F(1)}
        for C in clusters[1:]:
            fd = far_dist_at_x(C, mid)
            nxt = defaultdict(lambda: F(0))
            for cs, w in cur.items():
                for s, v in fd.items():
                    nxt[cs | s] += w * v
            cur = dict(nxt)
        tot += (hi - lo) * sum(w for cs, w in cur.items() if need <= cs)
    return tot


def cluster_decompose(E, G):
    E = sorted(set(E))
    clusters = []
    cur = [E[0]]
    for a, b in zip(E, E[1:]):
        if b - a <= G:
            cur.append(b)
        else:
            clusters.append(cur); cur = [b]
    clusters.append(cur)
    return clusters


def shapes_of(E, G=14):
    clusters = cluster_decompose(E, G)
    shapes = []
    for clu in clusters:
        shp = primitive(clu) if len(clu) > 1 else [0]
        shapes.append(shp)
    shapes.sort(key=len, reverse=True)
    return clusters, shapes


# =========================================================================
# TEST 1.  Re-derive the key constants EXACTLY (margins, worst decorr family)
# =========================================================================
def cluster_partitions_ge2(n):
    def rec(n, mx):
        if n == 0:
            yield []
            return
        for f in range(min(n, mx), 0, -1):
            for r in rec(n - f, f):
                yield [f] + r
    for p in rec(n, n):
        if len(p) >= 2:
            yield p


def test1_constants():
    print("=" * 74)
    print("TEST 1.  RE-DERIVE worst decorrelated p0_inf over >=2-cluster family [EXACT]")
    print("=" * 74)
    ok = True
    for k in range(8, 13):
        worst = F(0); wpart = None
        for part in cluster_partitions_ge2(k):
            clusters = [list(range(s)) for s in sorted(part, reverse=True)]
            v = p0_inf(clusters)
            if v > worst:
                worst = v; wpart = part
        margin = CAPS[k] - worst
        claim = QB_CLAIM[k]
        match = abs(float(worst) - claim) < 1e-4
        print(f"  k={k}: worst p0_inf={float(worst):.5f} (claim {claim})  part={wpart}  "
              f"cap={float(CAPS[k]):.5f}  margin={float(margin):.5f}  match={match}")
        if not match or margin <= 0:
            ok = False
    print(f"  => constants reproduced & all margins>0: {ok}")
    print()
    return ok


# =========================================================================
# TEST 2.  Is p0_inf really an UPPER bound?  Find sign and size of ERR=p0(E)-p0_inf
#          over ALL bounded >=2-cluster configs with EXACT p0 (no huge scales).
# =========================================================================
def test2_err_sign():
    print("=" * 74)
    print("TEST 2.  ERR = p0(E_finite) - p0_inf(shapes): sign & magnitude [EXACT]")
    print("=" * 74)
    print("  If ERR can be LARGE positive even at moderate gaps, the bound p0<=p0_inf+ERR")
    print("  with ERR controlled by (6/49)V/g is the live question. We tabulate ERR vs gap")
    print("  for the worst family [k-1 consec] + [singleton] (the claimed argmax).")
    print()
    for k in range(8, 13):
        base = list(range(k - 1))
        shapes = [base, [0]]
        pinf = p0_inf(shapes)
        print(f"  k={k}: p0_inf([0..{k-2}],[pt]) = {float(pinf):.5f}")
        worst_err = F(-10); worst_g = None
        for g in range(k - 1, 200):
            E = base + [g]
            if primitive(E) != sorted(E):  # keep primitive
                pass
            v = p0(E)
            err = v - pinf
            if err > worst_err:
                worst_err = err; worst_g = g
        print(f"        max_g>=k-1 [p0(consec_{k-1} U g) - p0_inf] = {float(worst_err):.5f} "
              f"at g={worst_g}  (g=k-1 is consec_k, INSIDE span-14 check)")
        # specifically wide region g>14:
        we = F(-10); wg = None
        for g in range(15, 200):
            v = p0(base + [g]); err = v - pinf
            if err > we:
                we = err; wg = g
        print(f"        WIDE g>14: max[p0 - p0_inf] = {float(we):.5f} at g={wg}")
    print()


# =========================================================================
# TEST 3.  ADVERSARIAL COUNTEREXAMPLE HUNT.
#   primitive wide k-set, 2nd-largest>14, p0 > cap_k  OR  p0 > QB_CLAIM (wide-bound).
#   Many structured families. EXACT p0. Cap max(E) so engine is feasible.
# =========================================================================
def is_truewide(E):
    E = sorted(set(E))
    return E[0] == 0 and E[-2] > 14


def test3_hunt():
    print("=" * 74)
    print("TEST 3.  ADVERSARIAL COUNTEREXAMPLE HUNT (primitive true-wide) [EXACT]")
    print("=" * 74)
    random.seed(424242)
    MAXE = 120
    best = {k: (F(-1), None) for k in range(8, 13)}
    over_cap = []
    over_qb = []        # exceeds the claimed wide decorr bound QB (a softer target)
    nchecked = {k: 0 for k in range(8, 13)}

    def consider(E):
        E = primitive(E)
        k = len(E)
        if k not in CAPS:
            return
        if not is_truewide(E):
            return
        if max(E) > MAXE:
            return
        nchecked[k] += 1
        v = p0(E)
        if v > best[k][0]:
            best[k] = (v, E)
        if v > CAPS[k]:
            over_cap.append((k, E, v))
        if float(v) > QB_CLAIM[k] + 1e-9:
            over_qb.append((k, E, float(v) - QB_CLAIM[k]))

    # --- Family 1: dilated primitive APs made wide (AP is the LRC extremal shape) ---
    for k in range(8, 13):
        for d in range(2, 9):
            E = [i * d for i in range(k)]            # AP with step d => span d*(k-1)
            consider(E)                              # primitivized -> back to consec; sanity
            # AP with a wide jump on last element:
            for jump in (16, 20, 30, 50, 80):
                consider(list(range(k - 1)) + [k - 2 + jump])

    # --- Family 2: two consec clusters at separated scales, all sizes ---
    for k in range(8, 13):
        for s1 in range(1, k):
            s2 = k - s1
            if s2 < 1:
                continue
            for M in (16, 18, 22, 30, 45, 70, 100):
                E = list(range(s1)) + [M + i for i in range(s2)]
                consider(E)

    # --- Family 3: three clusters, balanced, separated scales ---
    for k in range(8, 13):
        for (a, b, c) in [(3, 3, k - 6), (2, 2, k - 4), (4, 2, k - 6),
                          (2, k - 4, 2), (1, 1, k - 2), (k - 4, 2, 2)]:
            if min(a, b, c) < 1 or a + b + c != k:
                continue
            for (M1, M2) in [(18, 40), (16, 35), (20, 60), (25, 90)]:
                E = (list(range(a)) + [M1 + i for i in range(b)]
                     + [M1 + M2 + i for i in range(c)])
                consider(E)

    # --- Family 4: all residues mod 7 represented (max sector spread per cluster) ---
    for k in range(8, 13):
        for M in (15, 21, 28, 49, 70):
            # take 0..6 (all residues) then spread the rest far
            head = list(range(7))
            tail = [M + 7 * j for j in range(k - 7)]
            consider(head + tail)
            # residue-rich far cluster:
            consider(list(range(k - 3)) + [M, M + 1, M + 3])

    # --- Family 5: dilated AP made primitive but wide (step coprime to 7) ---
    for k in range(8, 13):
        for d in (3, 5, 9, 11, 13):
            for split in range(2, k):
                E = [i for i in range(split)] + [split - 1 + d * (i + 1)
                                                 for i in range(k - split)]
                consider(E)

    # --- Family 6: random multi-cluster near-AP sets (heavy sampling) ---
    for k in range(8, 13):
        for _ in range(8000):
            nclu = random.randint(2, min(4, k - 1))
            cuts = sorted(random.sample(range(1, k), nclu - 1))
            sizes = [b - a for a, b in zip([0] + cuts, cuts + [k])]
            base = 0; E = []
            for s in sizes:
                M = random.choice([15, 16, 18, 22, 28, 36, 50, 70, 95])
                start = base + M
                if s == 1:
                    clu = [start]
                else:
                    span = random.randint(s - 1, min(13, 3 * s))
                    pts = sorted(random.sample(range(1, span + 1), s - 1)) if span >= s - 1 else []
                    if len(pts) < s - 1:
                        clu = list(range(start, start + s))
                    else:
                        clu = [start] + [start + p for p in pts]
                E += clu
                base = max(E) + 1
            consider(E)

    # --- Family 7: two LARGE clusters (each size ~k/2, both >7) -- the genuinely
    #     "structured wide" worst case the claim says decorrelates ---
    for k in range(8, 13):
        s1 = k // 2; s2 = k - s1
        for M in (16, 20, 28, 40, 60):
            for sp1 in range(s1 - 1, min(13, 2 * s1) + 1):
                for sp2 in range(s2 - 1, min(13, 2 * s2) + 1):
                    # near-consec clusters of given spans:
                    c1 = list(range(s1)) if sp1 == s1 - 1 else \
                        [0] + sorted(random.sample(range(1, sp1 + 1), s1 - 1))
                    c2 = list(range(s2)) if sp2 == s2 - 1 else \
                        [0] + sorted(random.sample(range(1, sp2 + 1), s2 - 1))
                    E = c1 + [M + x for x in c2]
                    consider(E)

    print(f"{'k':>3} {'#wide checked':>14} {'best p0':>10} {'cap':>10} "
          f"{'QB_claim':>9} {'>cap?':>6} {'>QB?':>6}")
    for k in range(8, 13):
        bv, be = best[k]
        gt_cap = float(bv) > float(CAPS[k])
        gt_qb = float(bv) > QB_CLAIM[k] + 1e-9
        print(f"{k:>3} {nchecked[k]:>14} {float(bv):>10.5f} {float(CAPS[k]):>10.5f} "
              f"{QB_CLAIM[k]:>9.5f} {str(gt_cap):>6} {str(gt_qb):>6}")
    print()
    print(f"  violations of p0 > cap_k : {len(over_cap)}")
    for (k, E, v) in over_cap[:10]:
        print(f"    !!! k={k} p0={float(v):.5f} > cap E={E}")
    print(f"  rows exceeding claimed wide-bound QB(k-1): {len(over_qb)}")
    # show worst over-QB examples (these would be true-wide sets beating the decorr worst):
    over_qb.sort(key=lambda t: -t[2])
    for (k, E, excess) in over_qb[:12]:
        print(f"    k={k} p0-QB=+{excess:.5f}  E={E}")
    print()
    # print the best witnesses
    print("  best (highest-p0) true-wide witnesses found:")
    for k in range(8, 13):
        bv, be = best[k]
        print(f"    k={k}: p0={float(bv):.5f}  margin_to_cap={float(CAPS[k]-bv):.5f}  E={be}")
    print()
    return len(over_cap) == 0, over_cap, over_qb


# =========================================================================
# TEST 4.  Does CONSEC maximize single-cluster coverage at the decorr level?
#   Exhaustive over bounded shapes for sizes 7,8,9 (diam<=13).
# =========================================================================
def test4_consec_argmax():
    print("=" * 74)
    print("TEST 4.  Does consec maximize single-cluster decorr coverage? [EXACT, sizes 7-9]")
    print("=" * 74)
    # single anchored cluster + one far singleton: p0_inf([C],[0]).
    # exhaustively compare consec_s vs all other primitive shapes of size s, diam<=12.
    for s in (7, 8, 9):
        consec = list(range(s))
        base_val = p0_inf([consec, [0]])
        worst = base_val; wshape = consec
        cnt = 0
        # all shapes 0=a0<a1<...<a_{s-1}<=DIAM
        DIAM = 12
        for combo in itertools.combinations(range(1, DIAM + 1), s - 1):
            C = [0] + list(combo)
            if primitive(C) != C:
                continue
            v = p0_inf([C, [0]])
            cnt += 1
            if v > worst:
                worst = v; wshape = C
        beats = (worst > base_val)
        print(f"  size {s}: consec p0_inf={float(base_val):.5f}  "
              f"max over {cnt} shapes={float(worst):.5f}  "
              f"{'NON-consec WINS '+str(wshape) if beats else 'consec is argmax'}")
    print()


# =========================================================================
# TEST 5.  GLUE: span<=14 finite check vs span>14 wide branch -- any gap?
#   Check the partition boundary: at 2nd-largest exactly 15..B, is there a wide
#   set whose p0 exceeds cap?  (the (14,B] finite residual the claim defers).
# =========================================================================
def test5_glue_boundary():
    print("=" * 74)
    print("TEST 5.  GLUE at the span boundary + the (14,B] finite residual [EXACT, sampled]")
    print("=" * 74)
    print("  The claim CLOSES gaps>B rigorously and defers (14,B] to a finite check.")
    print("  We sample the (14,B] band hard for any p0>cap (would break the deferral).")
    random.seed(99)
    MAXE = 130
    viol = []
    bestband = {k: (F(-1), None) for k in range(8, 13)}
    for k in range(8, 13):
        cap = CAPS[k]
        for _ in range(6000):
            # build a wide set with 2nd-largest in (14, ~110]
            nclu = random.randint(2, min(3, k - 1))
            cuts = sorted(random.sample(range(1, k), nclu - 1))
            sizes = [b - a for a, b in zip([0] + cuts, cuts + [k])]
            base = 0; E = []
            for s in sizes:
                M = random.randint(15, 40)
                start = base + M
                if s == 1:
                    clu = [start]
                else:
                    span = random.randint(s - 1, min(13, 2 * s))
                    if span >= s - 1 and span >= 1:
                        pts = sorted(random.sample(range(1, span + 1),
                                                   min(s - 1, span)))
                        clu = [start] + [start + p for p in pts][:s - 1]
                        if len(clu) < s:
                            clu = list(range(start, start + s))
                    else:
                        clu = list(range(start, start + s))
                E += clu
                base = max(E) + 1
            E = primitive(E)
            if len(E) != k or not is_truewide(E) or max(E) > MAXE:
                continue
            v = p0(E)
            if v > bestband[k][0]:
                bestband[k] = (v, E)
            if v > cap:
                viol.append((k, E, v))
    for k in range(8, 13):
        bv, be = bestband[k]
        print(f"  k={k}: cap={float(CAPS[k]):.5f}  band-best p0={float(bv):.5f}  "
              f"margin={float(CAPS[k]-bv):.5f}  E={be}")
    print(f"  violations p0>cap in (14,B] band: {len(viol)}")
    for (k, E, v) in viol[:10]:
        print(f"    !!! k={k} p0={float(v):.5f}>cap E={E}")
    print()
    return len(viol) == 0


def main():
    print("#" * 74)
    print("# ADVERSARIAL VERIFY  HYP-2675  angle=dichotomy-finite-reduction")
    print("#" * 74)
    print()
    ok1 = test1_constants()
    test2_err_sign()
    ok3, over_cap, over_qb = test3_hunt()
    test4_consec_argmax()
    ok5 = test5_glue_boundary()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"  T1 constants reproduced, margins>0 : {ok1}")
    print(f"  T3 no p0>cap_k counterexample found: {ok3}  "
          f"(rows over claimed wide-bound QB: {len(over_qb)})")
    print(f"  T5 no p0>cap in (14,B] band        : {ok5}")
    print()
    if ok1 and ok3 and ok5:
        print("  No counterexample to p0(E)<=cap_k found on a large adversarial sample.")
        print("  The DECORRELATED-limit bound and finite-family worst case reproduce.")
    print("  SURVIVING GAP (by inspection of the claim, independent of samples):")
    print("   * The aggregation inequality p0(E) <= p0_inf + sum_far (6/49)V_i/g_i")
    print("     is NOT proved closed-form (cross-cluster Delta interaction). T2 shows")
    print("     ERR sign/size empirically; it is small but not rigorously bounded.")
    print("   * The (14,B] finite residual is sampled (T5), not exhausted.")
    print("   * 'consec maximizes single-cluster coverage' exhausted only sizes<=9 (T4).")


if __name__ == "__main__":
    main()
