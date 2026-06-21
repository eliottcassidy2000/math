#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_consec_max_adversarial_audit_macmini.py   (mac-mini 2026-06-21, THREAD D audit)

ADVERSARIAL AUDIT of the consec-max claim and the additive-energy lead (HYP-2735).

The ledger claims:
  (C1) consec [1..k] is the UNIQUE GLOBAL measS7-maximizer over ALL primitive k-sets
       of positive integers (exhaustive only proved for small N: k=8 N<=16, k=9 N<=15).
  (C2) additive energy AE = sum_s r_E(s)^2 is the dominant scalar (corr +0.58/+0.55),
       but is NECESSARY-NOT-SUFFICIENT (consec and AP-d2 have EQUAL AE=344 yet
       different measS7).

ADVERSARIAL QUESTIONS:
  Q1. Does WIDENING the span N produce a beater of consec at k=8? (push N as far as
      feasible; the global claim needs SOME a-priori span bound, currently absent.)
  Q2. Is there a CLEAN a-priori span bound proving wide sets cannot beat consec? Test:
      measS7(E) as a function of span = max(E)-min(E). Is consec (minimal span) really
      the top, AND does p0 DECAY with span at fixed k? (monotone-in-span would give the
      a-priori bound that closes C1 to a finite check.)
  Q3. The AE tie-break: among MAX-additive-energy primitive k-sets, is consec the unique
      measS7-max? And is the SECONDARY order really "minimal span"? Test directly: rank
      max-AE sets by measS7 and by span; check the correlation is monotone.

Reuses the EXACT measS7 (Fraction Lebesgue measure) from lrc_q108_consec_maximizer.
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
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def measS7(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e * mid) % 1) for e in E)) == P:
            tot += (b - a)
    return tot

def add_energy(E):
    # AE = #{(a,b,c,d) in E^4 : a+b = c+d} = sum_s r(s)^2, r(s)=#{(a,b): a+b=s}
    r = Counter()
    for a in E:
        for b in E:
            r[a + b] += 1
    return sum(v * v for v in r.values())

def primitive(combo):
    g = 0
    for e in combo: g = gcd(g, e)
    return g == 1

def main():
    print("#" * 80)
    print("# ADVERSARIAL AUDIT: consec-max global claim + additive-energy lead")
    print("#" * 80)

    consec8 = list(range(1, 9))
    pc8 = measS7(consec8)
    print(f"\nconsec8 = {consec8}  measS7 = {pc8} = {float(pc8):.6f}  AE={add_energy(consec8)}")

    # --- Q1: widen the span. Search all primitive 8-subsets with min=1 up to N. ---
    print("\n=== Q1: widen span N at k=8, hunt a beater of consec ===")
    for N in [16, 18, 20, 22, 24]:
        best = pc8; bestE = tuple(consec8); beat = 0; total = 0
        for combo in itertools.combinations(range(2, N + 1), 7):
            E = (1,) + combo
            if not primitive(E):
                continue
            total += 1
            v = measS7(E)
            if v > pc8 + Fr(1, 10**12):
                beat += 1
                if v > best:
                    best = v; bestE = E
        flag = "NO BEATER" if beat == 0 else f"!!! {beat} BEATERS, top={bestE} v={float(best):.6f}"
        print(f"  N={N:2d}: searched {total:7d} primitive sets (min=1)  ->  {flag}")

    # --- Q2: is p0 monotone-decreasing in span at fixed k? (the a-priori bound) ---
    print("\n=== Q2: measS7 vs span at k=8 (does p0 decay as the set spreads?) ===")
    # For each span s, find the MAX measS7 over primitive 8-sets with max-min=s, min=1.
    span_best = {}
    for N in [20]:
        for combo in itertools.combinations(range(2, N + 1), 7):
            E = (1,) + combo
            if not primitive(E):
                continue
            s = max(E) - min(E)
            v = measS7(E)
            if s not in span_best or v > span_best[s][0]:
                span_best[s] = (v, E)
    print("  span : max measS7 over primitive 8-sets with that span (min=1, max<=20)")
    for s in sorted(span_best):
        v, E = span_best[s]
        print(f"   span={s:2d}: max_p0={float(v):.6f}  witness={E}")
    spans = sorted(span_best)
    vals = [float(span_best[s][0]) for s in spans]
    is_mono = all(vals[i] >= vals[i + 1] - 1e-12 for i in range(len(vals) - 1))
    print(f"  -> monotone-decreasing in span? {is_mono}  (min span={spans[0]} is consec; "
          f"max-p0 at span={spans[vals.index(max(vals))]})")

    # --- Q3: AE tie-break. Among the TOP-AE primitive 8-sets, rank by measS7 & span. ---
    print("\n=== Q3: additive-energy tie-break audit at k=8 (min=1, max<=18) ===")
    rows = []
    for combo in itertools.combinations(range(2, 19), 7):
        E = (1,) + combo
        if not primitive(E):
            continue
        rows.append((add_energy(E), measS7(E), max(E) - min(E), E))
    maxAE = max(r[0] for r in rows)
    top = sorted([r for r in rows if r[0] >= maxAE - 0], key=lambda r: -r[1])
    print(f"  max AE among these sets = {maxAE}; #sets achieving it = {len(top)}")
    print("  top-AE sets ranked by measS7 (AE, measS7, span, E):")
    for ae, v, sp, E in top[:8]:
        print(f"    AE={ae} p0={float(v):.6f} span={sp} E={E}")
    # within max-AE class: is min-span == max-p0 ?
    if top:
        by_span = sorted(top, key=lambda r: r[2])
        best_p0 = max(top, key=lambda r: r[1])
        min_span = by_span[0]
        print(f"  within max-AE: min-span set = {min_span[3]} (span {min_span[2]}, p0 {float(min_span[1]):.6f})")
        print(f"                 max-p0   set = {best_p0[3]} (span {best_p0[2]}, p0 {float(best_p0[1]):.6f})")
        print(f"  tie-break 'min span => max p0' holds within max-AE class? "
              f"{min_span[3] == best_p0[3]}")

    # --- Q3b: GLOBAL spearman-like check: does (AE, -span) order measS7 perfectly? ---
    print("\n=== Q3b: does the order (AE desc, span asc) reproduce the measS7 order? ===")
    rows_sorted_true = sorted(rows, key=lambda r: -r[1])  # by measS7 desc
    rows_sorted_pred = sorted(rows, key=lambda r: (-r[0], r[2]))  # AE desc, span asc
    # count inversions: pairs where pred says A>B but measS7 says A<B
    # cheap: check top-10 agreement
    print("  rank | TRUE (by measS7)                    | PRED (by AE desc, span asc)")
    for i in range(min(10, len(rows))):
        t = rows_sorted_true[i]; p = rows_sorted_pred[i]
        print(f"   {i:2d}  | p0={float(t[1]):.5f} AE={t[0]} sp={t[2]} {t[3]}  | "
              f"p0={float(p[1]):.5f} AE={p[0]} sp={p[2]} {p[3]}")
    # how often does the predicted #1 equal the true #1?
    print(f"  PRED #1 == TRUE #1 (consec)? {rows_sorted_pred[0][3] == rows_sorted_true[0][3]}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
