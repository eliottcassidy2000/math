#!/usr/bin/env python3
"""
LRC@14 via the PROVEN n=7 case: the even-fold lever ("the 7 impossibility").
opus-2026-06-01-S554 (remote-control).

n=14 = 2*7 and LRC(7) is THEOREM (Barajas-Serra).  Exact bridge: an even speed
v=2u satisfies ||v t|| = ||u*(2t)||.  So with fold(S)={v/2 : v in S even} and
s=2t, the even runners' collar at time t equals the n=7-collar of fold at time s:
    g_even(t) = min_{u in fold} ||u s|| = g_fold(2t).
Hence
    M14(S) = max_t min(g_fold(2t), g_odd(t))  <=  max_s g_fold(s) = M(fold).      (*)

LEMMA (rigorous). If |fold(S)| <= 6 then M(fold) >= 1/7 by LRC(7); the even
runners are simultaneously >= 1/7 (>= 1/14) safe on E_good = {t : g_fold(2t) >=
1/7}, a set of POSITIVE measure.  So for such S the whole of LRC(14) is the ODD
residual:  is there t in E_good with g_odd(t) >= 1/14 ?

Antipodal structure of the residual: the two doubling-preimages of s, namely
s/2 and (s+1)/2, give every ODD runner ANTIPODAL positions (differ by 1/2), so
each odd runner is unsafe (< 1/14) at AT MOST ONE preimage.  Thus a single
even-good s yields a full n=14 witness UNLESS the odds "split" (some unsafe at
s/2, others at (s+1)/2).  LRC(14)[e<=6]  <=  "some even-good s has no odd-split".

This script:
 (1) verifies (*) exactly over many configs;
 (2) tests the structural law: do tight n=14 configs have e=6 and is fold n=7-TIGHT?
 (3) tests the fold reduction: over e<=6 configs, is there always an even-good s
     (g_fold(s) >= 1/7) whose better preimage is a full witness?  Characterise
     the odd-split failures.
"""

from fractions import Fraction
from math import gcd


def collar_candidates(V):
    """Exact candidate times for max-collar of speed set V (peaks + crossings)."""
    c = set()
    for v in V:
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    for k in range(abs(d) + 1):
                        c.add(Fraction(k, d) % 1)
    return c


def M_collar(V):
    """Exact max-collar M(V) = max_t min_i ||v_i t||."""
    best = Fraction(0)
    for t in collar_candidates(V):
        mn = min(min((v * t) % 1, 1 - (v * t) % 1) for v in V)
        if mn > best:
            best = mn
    return best


def collar_at(V, t):
    return min(min((v * t) % 1, 1 - (v * t) % 1) for v in V)


def fold(S):
    return tuple(sorted(v // 2 for v in S if v % 2 == 0))


def odds(S):
    return tuple(v for v in S if v % 2 == 1)


def primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in S)) if g else tuple(sorted(S))


# --------------------------------------------------------------------------
def verify_inequality(samples):
    print("== (1) Inequality  M14(S) <= M(fold(S))  (exact) ==")
    bad = 0
    for S in samples:
        f = fold(S)
        if not f:
            continue
        m14 = M_collar(S)
        mf = M_collar(f)
        if m14 > mf:
            bad += 1
            print(f"   VIOLATION: S={S} M14={m14} > M(fold)={mf}")
    print(f"   checked {len(samples)} sets, violations of M14<=M(fold): {bad}")
    print()


def tight_fold_structure():
    print("== (2) Tight n=14 configs: even-count and n=7 fold structure ==")
    AP = tuple(range(1, 14))
    Vstar = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
    for name, S in [("AP", AP), ("V* (AP,12->24)", Vstar)]:
        f = fold(S)
        e = len(f)
        m14 = M_collar(S)
        mf = M_collar(f)
        print(f"   {name}: S={S}")
        print(f"      even count e={e}, fold (halved evens)={f}")
        print(f"      M14={m14} (=1/14? {m14==Fraction(1,14)}), "
              f"M(fold) n=7-collar={mf} (=1/7? {mf==Fraction(1,7)})")
        print(f"      odds={odds(S)}")
    print()


def fold_reduction(samples):
    print("== (3) Fold reduction (e<=6): even-good s with full-witness preimage? ==")
    seven = Fraction(1, 7)
    fourteen = Fraction(1, 14)
    n_e6 = 0
    recovered = 0
    failures = []
    for S in samples:
        f = fold(S)
        if len(f) == 0 or len(f) > 6:
            continue
        n_e6 += 1
        m14 = M_collar(S)               # ground truth (LRC: should be >=1/14)
        # candidate even-good s: collar candidates of fold with g_fold(s) >= 1/7
        good_s = [s for s in collar_candidates(f)
                  if collar_at(f, s) >= seven]
        # plus crossings that define intervals -- sample interval midpoints too
        cs = sorted(collar_candidates(f))
        for i in range(len(cs)):
            a = cs[i]; b = cs[(i + 1) % len(cs)]
            length = (b - a) if b > a else (b - a + 1)
            mid = (a + length / 2) % 1
            if collar_at(f, mid) >= seven:
                good_s.append(mid)
        ok = False
        for s in good_s:
            for t in (s / 2 % 1, (s + 1) / 2 % 1):
                # even safe automatically (g_fold(2t)=g_fold(s)>=1/7); check odds
                od = odds(S)
                if not od or min(min((v * t) % 1, 1 - (v * t) % 1)
                                 for v in od) >= fourteen:
                    if collar_at(S, t) >= fourteen:   # full check (paranoia)
                        ok = True
                        break
            if ok:
                break
        if ok:
            recovered += 1
        else:
            failures.append((S, m14))
    print(f"   e<=6 configs: {n_e6}")
    print(f"   recovered a full witness from an even-good preimage: {recovered}")
    print(f"   fold-reduction FAILURES (witness exists but not via this route): "
          f"{len(failures)}")
    for S, m14 in failures[:15]:
        print(f"        fail: S={S}  M14={m14} (>=1/14? {m14>=fourteen})  "
              f"e={len(fold(S))}")
    print()
    return failures


def gen_samples():
    import random
    rng = random.Random(20260601)
    out = [tuple(range(1, 14)),
           (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)]
    # structured
    out.append(tuple([1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]))  # Fib
    out.append(tuple(2 * k + 1 for k in range(13)))                         # odd AP
    # random primitive, capped speeds to keep exact M fast
    made = 0
    seen = set()
    while made < 220:
        S = primitive(tuple(rng.sample(range(1, 40), 13)))
        if S in seen:
            continue
        seen.add(S)
        out.append(S)
        made += 1
    return out


if __name__ == "__main__":
    samples = gen_samples()
    verify_inequality(samples)
    tight_fold_structure()
    fold_reduction(samples)
