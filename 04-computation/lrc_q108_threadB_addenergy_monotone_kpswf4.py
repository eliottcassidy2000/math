#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadB_addenergy_monotone_kpswf4.py   (kind-pasteur 2026-06-21, THREAD B)

TARGET (LEMMA-1): is measS7 MONOTONE INCREASING in additive energy E(A)=sum_s r_A(s)^2 ?

The consec-maximizer crux (HYP-2735/2638): among k-subsets, measS7 is maximized by the
consecutive run [1..k].  THREAD-D found corr(measS7, AddEnergy)=+0.58 (k8)/+0.55 (k9) --
add-energy is the dominant scalar.  But add-energy alone is NOT a perfect order:
  consec[1..8]            add_energy=344  measS7=0.34490  span=7
  odd AP (1,3,..,15)      add_energy=344  measS7=0.30882  span=14   <- SAME energy, lower measS7
(odd AP = 2*[1..8]-1, NOT a dilate, so measS7 genuinely differs.)
So the order is (add_energy, then a span / dilation-class tie-break).

This script tests RIGOROUSLY:
 (a) COMPRESSION MONOTONICITY: does an add-energy-nondecreasing "compression" move always
     non-decrease measS7?  (find a clean counterexample or a clean monotone sub-step)
 (b) WITHIN-FIXED-SPAN MONOTONICITY: restricting to k-sets of a GIVEN span s (with min=0),
     is measS7 a monotone function of add_energy?  Test EXHAUSTIVELY for all spans at k=8.
 (c) GLOBAL RANK CORRELATION (Spearman/Kendall) of (measS7) vs (add_energy) and vs
     (add_energy, -span) lexicographic, on all primitive k=8 sets.

EXACT rational arithmetic for measS7.  Integer arithmetic for add_energy.
Outputs honest PROVED / SUPPORTED / REFUTED verdict for each sub-question.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd
from collections import Counter

P = 7

def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0:
            continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def measS7(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E)
    tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e * mid) % 1) for e in E)) == P:
            tot += (b - a)
    return tot

def add_energy(A):
    """E(A) = #{(a,b,c,d) in A^4 : a+b = c+d} = sum_s r_A(s)^2, r_A(s)=#{(a,b):a+b=s}."""
    A = list(A)
    r = Counter()
    for a in A:
        for b in A:
            r[a + b] += 1
    return sum(v * v for v in r.values())

def span(A):
    return max(A) - min(A)

def primitive(A):
    g = 0
    for e in A:
        g = gcd(g, e)
    return g == 1

# ============================================================================
# (b) WITHIN-FIXED-SPAN MONOTONICITY  -- the cleanest rigorous test.
# Normalize each k-set so min=0 (translation: does NOT change measS7? -- check) ;
# group by span; within each span, test whether measS7 is a monotone NON-DECREASING
# function of add_energy (i.e. add_energy(A)>add_energy(B) => measS7(A)>=measS7(B)).
# A single strict inversion (higher energy, strictly lower measS7) REFUTES it.
# ============================================================================
def translation_invariance_check(k=8, N=14, trials=200):
    """Verify measS7 is translation-invariant (so 'min=0' is WLOG).  measS7 uses e>=1
       (e=0 dropped), and shifting changes the actual integers -- so translation is
       NOT obviously invariant.  TEST it."""
    import random
    random.seed(1)
    bad = 0
    checked = 0
    for _ in range(trials):
        A = sorted(random.sample(range(1, N + 1), k))
        if not primitive(A):
            continue
        for shift in [1, 2, 3]:
            B = [a + shift for a in A]
            if measS7(A) != measS7(B):
                bad += 1
                if bad <= 5:
                    print(f"    TRANSLATION CHANGES measS7: A={A} -> +{shift}: "
                          f"{float(measS7(A)):.5f} vs {float(measS7(B)):.5f}")
        checked += 1
    return bad, checked

def within_span_monotone(k=8, N=15):
    """For all k-sets in {1..N} grouped by span (max-min), test monotonicity of measS7
       in add_energy.  A strict inversion = a pair A,B of the SAME span with
       E(A) > E(B) but measS7(A) < measS7(B).  Efficient envelope method:
       sort each span's rows by energy ASC; maintain running max measS7 seen at any
       strictly-lower energy; if a row's measS7 < that running max with strictly higher
       energy, it's an inversion (a higher-energy set beaten by a lower-energy one)."""
    by_span = {}  # span -> list of (add_energy, measS7, set)
    count = 0
    for combo in itertools.combinations(range(1, N + 1), k):
        s = span(combo)
        by_span.setdefault(s, []).append((add_energy(combo), measS7(combo), combo))
        count += 1
    print(f"    enumerated {count} k={k} subsets of {{1..{N}}}, {len(by_span)} distinct spans")
    total_inv = 0
    spans_with_inv = 0
    worst = None
    for s in sorted(by_span):
        rows = sorted(by_span[s], key=lambda t: t[0])  # by energy ASC
        # running best (max measS7, and the witnessing set) over all STRICTLY lower energies
        inv = 0
        worst_local = None
        idx = 0
        n = len(rows)
        run_max = None  # (measS7, set, energy) best at energy < current group
        while idx < n:
            e_cur = rows[idx][0]
            j = idx
            group = []
            while j < n and rows[j][0] == e_cur:
                group.append(rows[j]); j += 1
            # before adding this energy-group, run_max holds the best at strictly lower energy
            if run_max is not None:
                rm_v, rm_set, rm_e = run_max
                for ej, vj, Aj in group:
                    if vj < rm_v:  # higher energy ej>rm_e but strictly lower measS7
                        inv += 1
                        gap = float(rm_v - vj)
                        if worst_local is None or gap > worst_local[0]:
                            worst_local = (gap, Aj, ej, float(vj), rm_set, rm_e, float(rm_v))
            # update run_max with this group's max
            gv, gset, ge = max(group, key=lambda t: t[1])
            if run_max is None or gv > run_max[0]:
                run_max = (gv, gset, ge)
            idx = j
        if inv > 0:
            spans_with_inv += 1
            total_inv += inv
            if worst is None or worst_local[0] > worst[0]:
                worst = worst_local + (s,)
            if spans_with_inv <= 8:
                g, Aj, ej, vj, Ai, ei, vi = worst_local
                print(f"    span={s}: >={inv} inversion(s); worst: "
                      f"E={ej} set={Aj} m={vj:.5f}  <  "
                      f"E={ei} set={Ai} m={vi:.5f}  (gap {g:.5f})")
    print(f"    => spans with >=1 inversion: {spans_with_inv}/{len(by_span)}; "
          f"(envelope-counted) inversions: {total_inv}")
    if worst:
        g, Aj, ej, vj, Ai, ei, vi, s = worst
        print(f"    => WORST inversion (span {s}): higher-energy {Aj} (E={ej}) "
              f"m={vj:.5f} < lower-energy {Ai} (E={ei}) m={vi:.5f}, gap {g:.5f}")
    return total_inv, spans_with_inv, len(by_span)

# ============================================================================
# (a) COMPRESSION MONOTONICITY:  a "compression" replaces A by a more-AP set A' with
# E(A') >= E(A).  Canonical compression: Freiman/AP compaction.  We test the SIMPLEST
# move: replace A by [min(A)..min(A)+k-1] only if that move does not decrease energy
# (it never does for k-sets; consec maximizes energy).  But that's the global max.
# Instead test LOCAL compression moves: shrink one element toward the centroid by 1
# (a single-step move that should increase clustering / add-energy).  Does measS7
# increase under every energy-increasing single-element move?
# ============================================================================
def single_step_compression(k=8, N=14, trials=4000):
    """For random primitive A, try every single-element move a->a+1 or a->a-1 (staying a
       valid k-set of distinct positive ints).  Record moves that STRICTLY INCREASE
       add_energy and check whether they ever STRICTLY DECREASE measS7.  Such a move is
       a counterexample to 'energy-increasing single moves never decrease measS7'."""
    import random
    random.seed(7)
    examined = 0
    energy_up = 0
    energy_up_meas_down = 0
    examples = []
    for _ in range(trials):
        A = sorted(random.sample(range(1, N + 1), k))
        if not primitive(A):
            continue
        eA = add_energy(A)
        vA = measS7(A)
        Aset = set(A)
        for idx in range(k):
            for d in (-1, +1):
                b = A[idx] + d
                if b < 1 or b in Aset:
                    continue
                B = sorted(Aset - {A[idx]} | {b})
                if not primitive(B):
                    continue
                eB = add_energy(B)
                if eB > eA:
                    energy_up += 1
                    vB = measS7(B)
                    if vB < vA:
                        energy_up_meas_down += 1
                        if len(examples) < 8:
                            examples.append((A, B, eA, eB, float(vA), float(vB)))
        examined += 1
    print(f"    examined {examined} sets; energy-increasing single moves: {energy_up}; "
          f"of those, measS7 STRICTLY DECREASED: {energy_up_meas_down}")
    for A, B, eA, eB, vA, vB in examples:
        print(f"      A={A} (E={eA},m={vA:.5f}) -> B={B} (E={eB},m={vB:.5f})  "
              f"energy UP {eB-eA}, measS7 DOWN {vA-vB:.5f}")
    return energy_up, energy_up_meas_down

# ============================================================================
# (c) RANK CORRELATIONS on all primitive k=8 sets (small N for full enumeration).
# ============================================================================
def rank_correlations(k=8, N=15):
    data = []
    for combo in itertools.combinations(range(1, N + 1), k):
        if combo[0] != 1:
            continue
        if not primitive(combo):
            continue
        data.append((add_energy(combo), span(combo), measS7(combo), combo))
    n = len(data)
    print(f"    {n} primitive k={k} sets (min=1) in {{1..{N}}}")
    # Spearman-like: count concordant/discordant pairs for (add_energy) vs measS7
    conc = disc = ties_e = ties_m = 0
    conc2 = disc2 = ties2 = 0  # lexicographic (add_energy, -span)
    for i in range(n):
        ei, si, mi, _ = data[i]
        for j in range(i + 1, n):
            ej, sj, mj, _ = data[j]
            # measure side
            if mi == mj:
                ties_m += 1
                continue
            # energy side
            if ei == ej:
                ties_e += 1
            else:
                if (ei - ej) * (mi - mj) > 0:
                    conc += 1
                else:
                    disc += 1
            # lexicographic (energy, then -span): define key
            ki = (ei, -si); kj = (ej, -sj)
            if ki == kj:
                ties2 += 1
            else:
                cmp_lex = 1 if ki > kj else -1
                if cmp_lex * (mi - mj) > 0:
                    conc2 += 1
                else:
                    disc2 += 1
    tot = conc + disc
    tot2 = conc2 + disc2
    print(f"    add_energy vs measS7: concordant {conc}, discordant {disc}, "
          f"ties(energy) {ties_e}  ->  Kendall-tau-b(ish) = {(conc-disc)/tot:+.4f}" if tot else "")
    print(f"    LEX(add_energy,-span) vs measS7: concordant {conc2}, discordant {disc2}  "
          f"->  agreement = {(conc2-disc2)/tot2:+.4f}" if tot2 else "")
    return conc, disc, conc2, disc2


def main():
    print("#" * 84)
    print("# THREAD B -- is measS7 MONOTONE in additive energy?  (consec-max crux LEMMA-1)")
    print("#" * 84)

    print("\n=== honest gap reproduced ===")
    consec = [1, 2, 3, 4, 5, 6, 7, 8]
    oddap = [1, 3, 5, 7, 9, 11, 13, 15]
    for nm, A in [("consec[1..8]", consec), ("oddAP(1,3..15)", oddap)]:
        print(f"  {nm:18s} measS7={float(measS7(A)):.5f}  add_energy={add_energy(A)}  span={span(A)}")
    print("  SAME energy 344, DIFFERENT measS7 -> energy is NOT a total order; span tie-breaks.")

    print("\n=== (translation invariance check: is 'min=0' WLOG?) ===")
    bad, checked = translation_invariance_check(k=8, N=14, trials=200)
    print(f"  translation changes measS7 in {bad} cases out of {checked} sets "
          f"-> {'NOT translation-invariant (span grouping is the right frame)' if bad else 'translation-INVARIANT'}")

    print("\n=== (b) WITHIN-FIXED-SPAN MONOTONICITY (k=8, sets of {1..15}) ===")
    print("  Claim under test: within a fixed span, E(A)>E(B) => measS7(A)>=measS7(B).")
    ti, swi, ns = within_span_monotone(k=8, N=15)
    verdict_b = "REFUTED" if ti > 0 else "SUPPORTED (no inversion in this range)"
    print(f"  VERDICT (b): within-span monotone = {verdict_b}")

    print("\n=== (a) SINGLE-STEP COMPRESSION MONOTONICITY ===")
    print("  Claim under test: an add_energy-INCREASING single-element move never")
    print("  STRICTLY decreases measS7.")
    eu, eumd = single_step_compression(k=8, N=14, trials=4000)
    verdict_a = "REFUTED" if eumd > 0 else "SUPPORTED (no counterexample in sample)"
    print(f"  VERDICT (a): energy-up-never-measS7-down = {verdict_a}")

    print("\n=== (c) RANK CORRELATIONS (all primitive k=8, min=1, in {1..15}) ===")
    rank_correlations(k=8, N=15)

    print("\nDONE.")

if __name__ == "__main__":
    main()
