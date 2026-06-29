#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 floor via 2-ADIC PARITY DESCENT (mac-mini-2026-06-29-S4).

Connects the two creative angles: PARITY (Z2, t->1-t) and 2-ADIC DESCENT (t->2t).

EXACT RECURSION.  For a speed set S, split S = O u E (odd / even), let S' = E/2.
Since ||2s' t|| = ||s'(2t)|| and u=2t is measure-preserving 2-to-1,
   lonely(S) = lonely(O) cap lonely(E),     meas(lonely E) = meas(lonely S').
Define the PARITY CORRELATION
   rho(S) = meas(lonely S) / [ meas(lonely O) . meas(lonely S') ].
Recurse on S' until no evens remain.  Unrolling:
   meas(lonely S) = PROD_j rho_j  .  PROD_j meas(lonely O_j),
O_j the odd part at level j.  Each meas(lonely O_j) > 0 (odd set lonely at t=1/2,
>= cap_{|O_j|}).  So the floor reduces to:  rho_j >= c > 0 uniformly.

KEY SIMPLIFICATION.  The odd/even resonance lattice is 2Z (not 14Z), so the
Cauchy-Schwarz bound on rho uses the 2-SHEET count
   N2_O(t) = #{a in {0,1}: t + a/2 is O-safe}  (only 2 sheets!), giving
   rho(S) >= 1 - CV(N2_O) . sqrt((1-m')/m'),  m' = meas(lonely S').
This is far simpler than the 14-sheet count.  This script verifies the exact
recursion, sweeps rho_j over covering S and all levels, and tests the 2-sheet
CV bound.
"""
from __future__ import annotations
import functools, math, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F

W = F(1, 14)


def danger_intervals(p):
    p = abs(int(p))
    if p == 0:
        return [(F(0), F(1))]
    ivs = []
    for k in range(p + 1):
        lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
        if hi > lo:
            ivs.append((lo, hi))
    return ivs


def union_intervals(lists):
    ivs = sorted(iv for lst in lists for iv in lst)
    if not ivs:
        return []
    out = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo > out[-1][1]:
            out.append([lo, hi])
        else:
            out[-1][1] = max(out[-1][1], hi)
    return [(a, b) for a, b in out]


def complement(intervals):
    out = []; cur = F(0)
    for lo, hi in intervals:
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        out.append((cur, F(1)))
    return out


def lonely(S):
    if not S:
        return [(F(0), F(1))]
    return complement(union_intervals([danger_intervals(s) for s in S]))


def measure(intervals):
    return sum((hi - lo for lo, hi in intervals), F(0))


def m_lonely(S):
    return measure(lonely(S))


def intersect(a, b):
    out = []
    for (lo1, hi1) in a:
        for (lo2, hi2) in b:
            lo = max(lo1, lo2); hi = min(hi1, hi2)
            if hi > lo:
                out.append((lo, hi))
    return out


def shift_half(intervals):
    """shift by 1/2 mod 1."""
    res = []
    for lo, hi in intervals:
        a = lo + F(1, 2); b = hi + F(1, 2)
        if b <= 1:
            res.append((a if a < 1 else a - 1, b if b <= 1 else b - 1))
        else:
            res.append((a - 1 if a >= 1 else a, F(1) if a < 1 else b - 1))
            if a < 1:
                res.append((F(0), b - 1))
    return union_intervals([res])


def descend_levels(S):
    """Yield (level, O_j, S_j, S_{j+1}) until evens vanish."""
    cur = tuple(sorted(set(abs(int(x)) for x in S if x)))
    levels = []
    while cur:
        O = tuple(x for x in cur if x % 2 == 1)
        E = tuple(x for x in cur if x % 2 == 0)
        nxt = tuple(sorted(set(x // 2 for x in E)))
        levels.append((O, cur, nxt))
        if not E:
            break
        cur = nxt
    return levels


def is_covering(S):
    Sset = set(abs(int(s)) for s in S if s)
    return all(any(s % q == 0 for s in Sset) for q in range(2, 15))


def random_covering(rng, k=13, maxspeed=48):
    needed = [4, 8, 9, 5, 7, 11, 13]
    S = set()
    for q in needed:
        S.add(q * rng.randint(1, max(1, maxspeed // q)))
    while len(S) < k:
        S.add(rng.randint(1, maxspeed))
    S = tuple(sorted(set(list(S)[:k])))
    return S if (len(S) == k and is_covering(S)) else None


def cv2_two_sheet(O):
    """CV^2 of the 2-sheet count N2_O(t)=#{a in {0,1}: t+a/2 is O-safe}."""
    LO = lonely(O)
    mO = measure(LO)
    if mO == 0:
        return None, mO
    C = measure(intersect(LO, shift_half(LO)))   # both sheets safe
    EN = 2 * mO
    EN2 = 2 * mO + 2 * C        # E[(1_0+1_half)^2] = E1_0+E1_half + 2 E[1_0 1_half]
    Var = EN2 - EN * EN
    return Var / (EN * EN), mO


def main():
    print("=" * 86)
    print("LRC14 floor via 2-adic PARITY DESCENT (mac-mini-S4)")
    print("=" * 86)
    rng = random.Random(2718)

    samples = []
    tries = 0
    while len(samples) < 30 and tries < 8000:
        tries += 1
        S = random_covering(rng)
        if S and S not in samples:
            samples.append(S)

    # ---- verify exact recursion + product identity, collect rho_j ----
    print("\n[1] Exact recursion  meas(lonely S) = PROD rho_j . PROD meas(lonely O_j):")
    all_rho = []
    all_mO = []
    bad = 0
    depths = []
    for S in samples:
        levels = descend_levels(S)
        depths.append(len(levels))
        prod_rho = F(1); prod_mO = F(1)
        rho_list = []
        for (O, Sj, Snext) in levels:
            mSj = m_lonely(Sj); mO = m_lonely(O); mNext = m_lonely(Snext)
            if mO * mNext == 0:
                rho = None
            else:
                rho = mSj / (mO * mNext)
                prod_rho *= rho
                rho_list.append(float(rho))
                all_rho.append(float(rho))
            prod_mO *= mO
            all_mO.append(float(mO))
        check = (prod_rho * prod_mO == m_lonely(S))
        if not check:
            bad += 1
    print(f"   product identity exact for {len(samples)-bad}/{len(samples)} covering samples")
    print(f"   depth (levels) range: {min(depths)}..{max(depths)}")
    print(f"   rho_j over ALL levels/samples: min={min(all_rho):.4f}  max={max(all_rho):.4f}  "
          f"mean={sum(all_rho)/len(all_rho):.4f}")
    print(f"   meas(lonely O_j) (odd-set base measures): min={min(all_mO):.4f}  "
          f"(>0 always: {min(all_mO)>0})")

    # ---- the 2-sheet CV bound on rho ----
    print("\n[2] 2-sheet Cauchy-Schwarz bound  rho >= 1 - CV(N2_O).sqrt((1-m')/m'):")
    print(f"   {'level O':<26}{'m_O':>7}{'CV2_2sheet':>11}{'m_next':>8}{'rho_actual':>11}{'bound':>9}")
    worst_bound = 1.0
    cnt = 0
    for S in samples[:14]:
        for (O, Sj, Snext) in descend_levels(S):
            if not O:
                continue
            cv2, mO = cv2_two_sheet(O)
            mNext = m_lonely(Snext)
            mSj = m_lonely(Sj)
            if cv2 is None or mNext == 0:
                continue
            rho = mSj / (mO * mNext)
            bound = 1 - math.sqrt(float(cv2)) * math.sqrt(float((1 - mNext) / mNext)) if mNext > 0 else None
            worst_bound = min(worst_bound, bound)
            cnt += 1
            if cnt <= 14:
                print(f"   {str(O):<26}{float(mO):>7.3f}{float(cv2):>11.4f}{float(mNext):>8.3f}"
                      f"{float(rho):>11.4f}{bound:>9.4f}")
    print(f"   worst 2-sheet bound over sampled levels: {worst_bound:.4f}")

    # ---- compare resonance: parity split vs nothing ----
    print("\n[3] min rho_j (parity floor) vs the original R-14Q floor magnitude:")
    print(f"   parity-descent rho floor (min over all levels/samples) = {min(all_rho):.4f}")
    print(f"   => meas(lonely S) >= (min rho)^depth . PROD caps  -- positive if rho floor>0")
    print(f"   depths<= {max(depths)}, so worst product factor (min rho)^{max(depths)} = "
          f"{min(all_rho)**max(depths):.4e}")

    print("\n" + "=" * 86)
    print("If min rho_j > 0 uniformly and meas(lonely O_j)>0 (odd caps), the 2-adic")
    print("parity descent PROVES the covering floor recursively.  The 2-sheet CV bound")
    print("is the clean per-level certificate (2 sheets, not 14).")
    print("=" * 86)


if __name__ == "__main__":
    main()
