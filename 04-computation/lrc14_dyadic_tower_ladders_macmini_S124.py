#!/usr/bin/env python3
"""HYP-7970 referee — the dyadic tower read of the two bottom ladders
(mac-mini-2026-07-19-S124).

Content:
  (D)  DYADIC SET IDENTITY: F_m(N) := {1..N-2, N, m(N-1)} satisfies, for odd N,
         F_m(N) = odds(N)  ∪  2·K_m((N-1)/2),      K_c(h) := {1..h-1, c·h},
       and recursively K_m(6) = odds(5) ∪ 2·K_m(3): the F-ladder is the AP
       with its single defect pushed down the 2-adic tower and stretched at
       the bottom.  (Trivial arithmetic; the value of the identity is that
       opus-S395's ladder M(F_m(13)) = m/(12m+5) is kind-pasteur's half-level
       K-ladder M(K_c(h)) = c/(ch+1) LIFTED, tying the GW-acceleration /
       binder-gate / determinant-slack vocabularies to the 2-adic tower.)

  (A)  CONTROLS (instrument-validation gate, kps rule: rediscover knowns
       from inside the instrument before any negative counts):
       F_2(13)=1/14, F_3(13)=3/41, F_4(13)=4/53, K_2(13)=2/27, K_3(13)=3/40,
       K_4(13)=4/53.

  (B)  THE FRACTIONAL CROSSING MEMBER F_{5/2}(13) = {1..11,13,30} (the set
       exists though m=5/2 is not integral; slack 2m-5 = 0 there).
       PREDICTION registered before the run: M = 3/41 with straddle (11,30).

  (C)  THE 4/55 SINGLE-FAR CANDIDATES: within {1..11,13,x} the only pairs
       summing to 55 need x ∈ {42,44,46,48,50,52,54}.  Exact M for each —
       independent verification of the binder-gate closure ("4/55 needs a
       non-single-far realizer") + the first flip-spectrum values there.

  (E)  ODD-LAYER DEFECT SWEEP (the ρ-antisymmetric flip direction from
       HYP-7960): for base ∈ {AP-14, GW-14}, replace one odd o by a new odd
       o' ∈ {15..29}: exact M of all replacements.  Instrument check: none
       may be 1/14 (THM-1142 tight locus = exactly two families) and any
       value in (1/14, 3/41) would be a sensation (expect none).

  (F)  THE MOD-6 MECHANISM PROBE: F_2(13) (ties, N=13 ≡ 1 mod 6) vs
       F_2(9) = {1..7,9,16} (N=9 ≢ 1 mod 6, kps-c86: does NOT tie 1/10):
       exact M + straddle of both — the failure mechanism, dyadically read.

Exactness: M(W) = max_t min_w ||w t|| is attained at a breakpoint of the
lower envelope of the clearance tents, i.e. at some t = k/q with
q ∈ {u+v} ∪ {|u-v|} ∪ {2w}; we scan all of them with integer arithmetic.
"""

from fractions import Fraction as F
from math import gcd
from itertools import combinations

ONE14 = F(1, 14)
GAP_HI = F(3, 41)


def exact_M(W):
    """Exact M(W) plus a maximizer t* and the straddle data at t*."""
    W = sorted(W)
    mods = set()
    for u, v in combinations(W, 2):
        mods.add(u + v)
        if u != v:
            mods.add(v - u)
    for w in W:
        mods.add(2 * w)
    mods.discard(0)
    best_n, best_d = 0, 1          # best clearance as fraction best_n/best_d
    best_t = (0, 1)
    for q in sorted(mods):
        for k in range(1, q):
            m = q                   # min_w min(wk mod q, q - wk mod q)
            for w in W:
                r = (w * k) % q
                if q - r < r:
                    r = q - r
                if r < m:
                    m = r
                    if m == 0:
                        break
            if m * best_d > best_n * q:
                best_n, best_d, best_t = m, q, (k, q)
    M = F(best_n, best_d)
    k, q = best_t
    t = F(k, q)
    # straddle: active runners at t*, and active pairs with (u+v)t* ∈ ℤ
    act = [w for w in W if abs_dist(w * t) == M]
    pairs = [(u, v) for u, v in combinations(act, 2) if ((u + v) * k) % q == 0]
    return M, t, act, pairs


def abs_dist(x: F) -> F:
    fr = x - (x.numerator // x.denominator)
    return min(fr, 1 - fr)


def dyadic_word(ws):
    if not ws:
        return ""
    s = "".join("O" if w % 2 else "E" for w in sorted(ws))
    ev = [w // 2 for w in sorted(ws) if w % 2 == 0]
    return s + ("|" + dyadic_word(ev) if ev else "")


def report(name, W, expect=None):
    M, t, act, pairs = exact_M(W)
    pstr = []
    for (u, v) in pairs:
        s = u + v
        D = M * s
        slack = 14 * D - s if D.denominator == 1 else None
        par = "both-odd" if (u % 2 and v % 2) else ("both-even" if (u % 2 == 0 and v % 2 == 0) else "mixed")
        pstr.append(f"({u},{v}) s={s} D={D} slack={slack} {par}")
    ok = ""
    if expect is not None:
        ok = "  MATCH" if M == expect else f"  *** MISMATCH (expected {expect})"
    flag = ""
    if ONE14 < M < GAP_HI:
        flag = "  <<< IN-GAP VALUE (1/14, 3/41) !!!"
    if M == ONE14:
        flag = "  [floor exactly]"
    print(f"  {name:<28} M = {str(M):>8}  t* = {t}  active={act}")
    print(f"      straddles: {'; '.join(pstr) if pstr else 'none at t*'}{ok}{flag}")
    return M


def main():
    # (D) set identities
    print("== (D) dyadic set identities ==")
    for m in range(2, 11):
        FmN = set(range(1, 12)) | {13} | {12 * m}
        odds = {1, 3, 5, 7, 9, 11, 13}
        Km6 = set(range(1, 6)) | {6 * m}
        assert FmN == odds | {2 * x for x in Km6}, m
        Km3 = {1, 2} | {3 * m}
        assert Km6 == {1, 3, 5} | {2 * x for x in Km3}, m
    print("  F_m(13) = odds ∪ 2·K_m(6) = odds ∪ 2·(odds5 ∪ 2·K_m(3))  for m=2..10  OK")

    print("\n== (A) controls (instrument-validation gate) ==")
    report("F_2(13) = GW {1..11,13,24}", list(range(1, 12)) + [13, 24], F(1, 14))
    report("F_3(13) {1..11,13,36}", list(range(1, 12)) + [13, 36], F(3, 41))
    report("F_4(13) {1..11,13,48}", list(range(1, 12)) + [13, 48], F(4, 53))
    report("K_2(13) {1..12,26}", list(range(1, 13)) + [26], F(2, 27))
    report("K_3(13) {1..12,39}", list(range(1, 13)) + [39], F(3, 40))
    report("K_4(13) {1..12,52}", list(range(1, 13)) + [52], F(4, 53))

    print("\n== (B) the fractional crossing F_{5/2}(13) = {1..11,13,30} ==")
    print("  PREDICTION (registered pre-run): M = 3/41, straddle (11,30)")
    report("F_{5/2}(13) {1..11,13,30}", list(range(1, 12)) + [13, 30])

    print("\n== (C) the 4/55 single-far candidates {1..11,13,x}, x pairs to 55 ==")
    for x in (42, 44, 46, 48, 50, 52, 54):
        M = report(f"x={x}", list(range(1, 12)) + [13, x])
        if M == F(4, 55):
            print("      *** 4/55 REALIZED — SENSATION, recheck everything ***")
    print("  (binder-gate closure predicts none hits 4/55)")

    print("\n== (E) odd-layer defect sweep (ρ-antisymmetric flips) ==")
    for base_name, base in (("AP-14", list(range(1, 14))),
                            ("GW-14", list(range(1, 12)) + [13, 24])):
        odds_in = [o for o in base if o % 2 == 1]
        vals = {}
        for o in odds_in:
            for op in range(15, 30, 2):
                if op in base:
                    continue
                Wd = sorted(set(base) - {o} | {op})
                M, t, act, pairs = exact_M(Wd)
                vals[(o, op)] = M
                assert M != ONE14, f"new tight family?! {Wd}"
                if ONE14 < M < GAP_HI:
                    print(f"      IN-GAP: {base_name} -{o}+{op}: M={M}  <<<")
        lo = min(vals.values())
        hi = max(vals.values())
        lo_k = [k for k, v in vals.items() if v == lo]
        below = sum(1 for v in vals.values() if v < ONE14)
        print(f"  {base_name}: {len(vals)} defect families; M range [{lo}, {hi}]; "
          f"min at {lo_k}; below-floor count = {below}")

    print("\n== (E2) even-layer defect sweep + the Hamming-1 flip spectrum ==")
    from collections import Counter
    for base_name, base in (("AP-14", list(range(1, 14))),
                            ("GW-14", list(range(1, 12)) + [13, 24])):
        evens_in = [e for e in base if e % 2 == 0]
        spec = Counter()
        ingap = []
        for e in evens_in:
            for ep in range(14, 62, 2):
                if ep in base:
                    continue
                Wd = sorted(set(base) - {e} | {ep})
                M, t, act, pairs = exact_M(Wd)
                spec[M] += 1
                assert M >= ONE14, f"below floor?! {Wd} M={M}"
                if M == ONE14 and set(Wd) not in ({*range(1, 14)},
                                                  {*range(1, 12), 13, 24}):
                    print(f"      NEW TIGHT?! {base_name} -{e}+{ep}")
                if ONE14 < M < GAP_HI:
                    ingap.append((e, ep, M))
        lo = min(spec)
        print(f"  {base_name}: {sum(spec.values())} even-defect families "
              f"(window e' ≤ 60); min M = {lo}; in-gap: {ingap if ingap else 'NONE'}")
        top = sorted(spec.items())[:8]
        print(f"      spectrum bottom: {[(str(v), c) for v, c in top]}")
    print("  TAIL (exact, all e' > 60): for the (12→e') / (24→e') defect the")
    print("  12-duty forces e' ∈ 12Z (else t=1/12 gives M = 1/12): the full")
    print("  infinite tail is the F-ladder m/(12m+5) ≥ 3/41; other-even defects")
    print("  keep their duties inside the window — windowed data honest as such.")
    report("control {1..11,13} (12 speeds)", list(range(1, 12)) + [13])

    print("\n== (E3) cross-parity defect sweeps (completing windowed Hamming-1) ==")
    for base_name, base in (("AP-14", list(range(1, 14))),
                            ("GW-14", list(range(1, 12)) + [13, 24])):
        spec = Counter()
        ingap = []
        tights = []
        # even -> odd'
        for e in [x for x in base if x % 2 == 0]:
            for op in range(15, 30, 2):
                if op in base:
                    continue
                Wd = sorted(set(base) - {e} | {op})
                M, *_ = exact_M(Wd)
                spec[M] += 1
                assert M >= ONE14, (base_name, e, op, M)
                if M == ONE14:
                    tights.append((e, op))
                if ONE14 < M < GAP_HI:
                    ingap.append((e, op, M))
        # odd -> even'
        for o in [x for x in base if x % 2 == 1]:
            for ep in range(14, 62, 2):
                if ep in base:
                    continue
                Wd = sorted(set(base) - {o} | {ep})
                M, *_ = exact_M(Wd)
                spec[M] += 1
                assert M >= ONE14, (base_name, o, ep, M)
                if M == ONE14:
                    tights.append((o, ep))
                if ONE14 < M < GAP_HI:
                    ingap.append((o, ep, M))
        lo = min(spec)
        print(f"  {base_name}: {sum(spec.values())} cross-parity defect families; "
              f"min M = {lo}; floor-hits: {tights if tights else 'none'}; "
              f"in-gap: {ingap if ingap else 'NONE'}")
        print(f"      spectrum bottom: {[(str(v), c) for v, c in sorted(spec.items())[:8]]}")

    print("\n== (F) the mod-6 mechanism probe ==")
    report("F_2(13) (N≡1 mod 6, ties)", list(range(1, 12)) + [13, 24])
    report("F_2(9) = {1..7,9,16} (N≢1)", list(range(1, 8)) + [9, 16])
    print("  (n=10 floor is 1/10 = 0.1; F_2(9) straddle shows the non-tie mechanism)")

    print("\nDone.")


if __name__ == "__main__":
    main()
