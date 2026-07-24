#!/usr/bin/env python3
"""
LRC(14): SGC'(13) and the tight classification on the k=3 (three-perturbation) family.
mac-mini-2026-07-24-S170 (Opus).

Family: S = ({1..13}\{i,j,l}) u {w1,w2,w3}, primitive, 13 distinct speeds.
GOAL: (a) no such S has gap in (1/14, 3/41);  (b) the only tight ones (gap=1/14) are AP and GW.

BOUNDED SEARCH (from the stranger lemmas, theta = 3/41):
  multi-stranger lemma: if f_C >= theta on an interval of length delta and k strangers all satisfy
  w_i >= 1/delta, then each bad set has measure <= 4*theta*delta, so 4*k*theta < 1 => some tau survives
  => gap >= theta.   theta=3/41: k=3 gives 36/41<1, k=2 gives 24/41<1, k=1 gives 12/41<1.
  Hence, with w1<w2<w3:
     w1 < B1 = 1/delta(10-core)          [else all three decouple]
     w2 < B2 = 1/delta(10-core u {w1})   [else the remaining two decouple]
     w3 < B3 = 1/delta(11-core u {w2})   [single-stranger lemma]
  All deltas are exact rationals; the region is therefore finite and DERIVED, not assumed.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
import itertools, time, importlib.util

spec = importlib.util.spec_from_file_location(
    "m", "04-computation/lrc14_sgc_prime_single_perturbation_theorem_macmini_S169.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

THETA = F(3, 41)
ONE14 = F(1, 14)
QMAX = 60
QK = [(q, k) for q in range(2, QMAX + 1) for k in range(1, q // 2 + 1) if gcd(k, q) == 1]
Qs = np.array([q for q, _ in QK])
Ks = np.array([k for _, k in QK])
TH_f = float(THETA)


def core_min_arr(C):
    m = Qs.copy()
    for v in C:
        r = (v * Ks) % Qs
        m = np.minimum(m, np.minimum(r, Qs - r))
    return m


def main():
    t0 = time.time()
    AP = tuple(range(1, 14))
    GW = tuple(list(range(1, 12)) + [13, 24])
    tights, band = set(), []
    exact_calls = 0
    processed = 0
    for (i, j, l) in itertools.combinations(range(1, 14), 3):
        C10 = [x for x in range(1, 14) if x not in (i, j, l)]
        d10 = M.longest_safe_interval(C10, THETA)
        if d10 <= 0:
            continue
        B1 = ceil(1 / d10)
        for w1 in range(1, B1 + 1):
            if w1 in C10:
                continue
            C11 = sorted(C10 + [w1])
            if len(C11) != 11:
                continue
            d11 = M.longest_safe_interval(C11, THETA)
            if d11 <= 0:
                continue
            B2 = ceil(1 / d11)
            for w2 in range(w1 + 1, B2 + 1):
                if w2 in C11:
                    continue
                C12 = sorted(C11 + [w2])
                if len(C12) != 12:
                    continue
                d12 = M.longest_safe_interval(C12, THETA)
                if d12 <= 0:
                    continue
                B3 = ceil(1 / d12)
                processed += 1
                w3s = np.array([w for w in range(w2 + 1, B3 + 1) if w not in C12], dtype=np.int64)
                if w3s.size == 0:
                    continue
                cm = core_min_arr(C12)
                r = (np.outer(w3s, Ks)) % Qs
                d = np.minimum(r, Qs - r)
                lb = (np.minimum(cm[None, :], d) / Qs[None, :]).max(axis=1)
                for w3 in w3s[lb < TH_f - 1e-12]:
                    V = tuple(sorted(set(C12 + [int(w3)])))
                    if len(V) != 13 or reduce(gcd, V) != 1:
                        continue
                    exact_calls += 1
                    g, t = M.exact_gap(list(V))
                    if g == ONE14:
                        tights.add(V)
                    elif ONE14 < g < THETA:
                        band.append((g, V, t))
        if (i, j, l)[0] != getattr(main, "_last", None):
            main._last = (i, j, l)[0]
            print(f"  ... i={i} done, elapsed {time.time()-t0:.0f}s, "
                  f"(triple,w1,w2) processed={processed:,}, exact={exact_calls}", flush=True)
    print(f"\nDONE in {time.time()-t0:.0f}s")
    print(f"  (triple,w1,w2) nodes processed: {processed:,};  exact gap evaluations: {exact_calls:,}")
    print(f"  sets with gap in (1/14, 3/41): {len(band)}")
    for g, V, t in sorted(band)[:10]:
        print(f"     gap={g} V={V} tau={t}")
    print(f"  TIGHT sets (gap = 1/14): {len(tights)}")
    for V in sorted(tights):
        lab = "AP" if V == AP else ("GW" if V == GW else "*** NEW TIGHT SET ***")
        print(f"     {lab:>22}  {V}")
    if not band:
        print("\n  ==> k=3 SGC' THEOREM: no 3-perturbation set has gap in (1/14, 3/41).")
    if tights <= {AP, GW}:
        print("  ==> k=3 TIGHT CLASSIFICATION: only AP and GW (no new tight sets).")


if __name__ == "__main__":
    main()
