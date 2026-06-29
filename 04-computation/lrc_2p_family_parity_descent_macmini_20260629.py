#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC for the DIHEDRAL family n = 2p: validate the 2-adic parity descent
recursively, and test the p mod 4 (Brouwer vs Borsuk-Ulam) distinction.
(mac-mini-2026-06-29-S5)

n = 2p = |D_p| (dihedral order).  Threshold 1/n, danger comb ||s t|| < 1/n
(width 2/n = 1/p), p sectors = the p-gon vertices zeta_p^k.  The complement
x->-x is the p-gon reflection: an AUTOMORPHISM iff p=1 mod4, an ANTI-automorphism
(self-converse, free Z_2) iff p=3 mod4 (kps S31av).

THM-580 parity descent (general n): split S = O u E, S' = E/2,
   lonely(S) = lonely(O) cap lonely(E),  meas(lonely E) = meas(lonely S'),
   meas(lonely S) = PROD_j rho_j . PROD_j meas(lonely O_j).

TEST across the family:
 - n=6  (p=3, =3mod4): LRC(6) is KNOWN/proven -> the descent route MUST validate it.
 - n=10 (p=5, =1mod4): Brouwer/SOS side.
 - n=14 (p=7, =3mod4): the target (Borsuk-Ulam side).
 - n=22 (p=11,=3mod4): next p=3mod4.
Report: does meas(lonely S)>0 for all covering S? is min rho_j > 0 uniformly?
does the p mod 4 split show in rho_j or the floor structure?
"""
from __future__ import annotations
import functools, math, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F


def make_ops(n):
    """Return danger/lonely machinery for threshold 1/n."""
    W = F(1, n)

    def danger(p):
        p = abs(int(p))
        if p == 0:
            return [(F(0), F(1))]
        ivs = []
        for k in range(p + 1):
            lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
            if hi > lo:
                ivs.append((lo, hi))
        return ivs

    def uni(lists):
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

    def comp(intervals):
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
        return comp(uni([danger(s) for s in S]))

    def meas(ivs):
        return sum((hi - lo for lo, hi in ivs), F(0))

    return lonely, meas


def descend(S):
    cur = tuple(sorted(set(abs(int(x)) for x in S if x)))
    out = []
    while cur:
        O = tuple(x for x in cur if x % 2 == 1)
        E = tuple(x for x in cur if x % 2 == 0)
        nxt = tuple(sorted(set(x // 2 for x in E)))
        out.append((O, cur, nxt))
        if not E:
            break
        cur = nxt
    return out


def is_covering(S, n):
    Sset = set(abs(int(s)) for s in S if s)
    return all(any(s % q == 0 for s in Sset) for q in range(2, n + 1))


def rand_cov(rng, n, maxs):
    k = n - 1
    needed = []
    # prime powers <= n force divisibility
    for q in range(2, n + 1):
        # keep prime powers
        qq = q; pf = {}
        m = qq; d = 2
        while m > 1:
            while m % d == 0:
                pf[d] = pf.get(d, 0) + 1; m //= d
            d += 1
        if len(pf) == 1:
            needed.append(q)
    S = set()
    for q in needed:
        S.add(q * rng.randint(1, max(1, maxs // q)))
    while len(S) < k:
        S.add(rng.randint(1, maxs))
    S = tuple(sorted(set(list(S)[:k])))
    return S if (len(S) == k and is_covering(S, n)) else None


def main():
    print("=" * 80)
    print("LRC dihedral family n=2p: 2-adic parity descent validation (mac-mini-S5)")
    print("=" * 80)
    rng = random.Random(1729)

    for n in (6, 10, 14, 22):
        p = n // 2
        pmod4 = p % 4
        kind = "Brouwer/SOS (p=1mod4)" if pmod4 == 1 else "Borsuk-Ulam (p=3mod4)"
        known = "KNOWN/proven" if n <= 7 else "OPEN"
        lonely, meas = make_ops(n)
        maxs = {6: 18, 10: 28, 14: 40, 22: 56}[n]
        # sample covering (n-1)-sets
        samples = []
        tries = 0
        target = 60 if n <= 14 else 25
        while len(samples) < target and tries < 12000:
            tries += 1
            S = rand_cov(rng, n, maxs)
            if S and S not in samples:
                samples.append(S)
        if not samples:
            print(f"\nn={n} (p={p}, {kind}, {known}): no covering samples found")
            continue
        # descent stats
        all_rho = []; all_mO = []; min_floor = F(2); zero = 0; depths = []
        for S in samples:
            lv = descend(S)
            depths.append(len(lv))
            mS = meas(lonely(S))
            if mS == 0:
                zero += 1
            min_floor = min(min_floor, mS)
            for (O, Sj, Snext) in lv:
                mO = meas(lonely(O)); mN = meas(lonely(Snext)); mSj = meas(lonely(Sj))
                if O:
                    all_mO.append(float(mO))
                if mO * mN != 0:
                    all_rho.append(float(mSj / (mO * mN)))
        print(f"\nn={n} (p={p}, {kind}, LRC {known}):  {len(samples)} covering sets, "
              f"depth {min(depths)}..{max(depths)}")
        print(f"   meas(lonely S): min={float(min_floor):.5f}  (any zero? {'YES!' if zero else 'no'})")
        print(f"   rho_j: min={min(all_rho):.4f}  mean={sum(all_rho)/len(all_rho):.4f}  "
              f"frac<0.7={sum(1 for r in all_rho if r<0.7)/len(all_rho):.3f}")
        print(f"   odd-base meas(lonely O_j): min={min(all_mO):.4f}  (>0: {min(all_mO)>0})")

    print("\n" + "=" * 80)
    print("Validation: if min meas(lonely S)>0 for the KNOWN n=6 (and the descent")
    print("rho_j/odd-caps are bounded), the descent route is sound; the SAME structure")
    print("at n=14 then differs only in scale.  Watch whether rho_j min depends on p mod4.")
    print("=" * 80)


if __name__ == "__main__":
    main()
