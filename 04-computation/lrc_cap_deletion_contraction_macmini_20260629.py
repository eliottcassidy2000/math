#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The LRC cap is a measure-valued CLAIM A: a deletion-contraction with a
factor-2 (sigma-symmetry) and a per-tooth ("cycle") sum.  (mac-mini-2026-06-29-S7)

CLAIM A (tournament core, THM-070):  H(T) - H(T-v) = 2 * sum_{C odd cycle through v} mu(C).
LRC ANALOG (this script): for cap(S) = meas(lonely(S)) = meas{t: ||s t||>=1/n for all s in S},
   cap(S\{s}) - cap(S) = meas( lonely(S\{s}) cap D_s )      [the conditional danger of s]
                       = 2 * sum_{k=0}^{s-1} mu_LRC(s,k)     [factor 2 from sigma:t->-t; teeth = cycles]
where D_s = s "teeth" (arcs ||s t||<1/n at t ~ k/s), tooth_k the k-th, and
   mu_LRC(s,k) = meas( lonely(S\{s}) cap tooth_k cap [0,1/2) ).
The teeth k=0..s-1 are the s cyclic positions of s (the "odd cycles through s"); the factor 2
is the complement Z_2 (lonely cap D_s is sigma-symmetric).  Iterating from cap(empty)=1 gives the
cap as 1 - sum of conditional dangers = a chain rule = the measure-valued OCF / Redei analog.

We VERIFY the deletion-contraction, the factor 2, the tooth decomposition, and the chain rule.
"""
from __future__ import annotations
import functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F


def make(n):
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

    def tooth(p, k):
        """the k-th tooth of D_p: arc around t=k/p (k=0..p), clipped to [0,1)."""
        lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
        return [(lo, hi)] if hi > lo else []

    def uni(lsts):
        ivs = sorted(iv for L in lsts for iv in L)
        if not ivs:
            return []
        out = [list(ivs[0])]
        for lo, hi in ivs[1:]:
            if lo > out[-1][1]:
                out.append([lo, hi])
            else:
                out[-1][1] = max(out[-1][1], hi)
        return [(a, b) for a, b in out]

    def comp(ivs):
        out = []; cur = F(0)
        for lo, hi in ivs:
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

    def inter(a, b):
        out = []
        for (l1, h1) in a:
            for (l2, h2) in b:
                lo = max(l1, l2); hi = min(h1, h2)
                if hi > lo:
                    out.append((lo, hi))
        return out

    def clip_half(ivs):
        return inter(ivs, [(F(0), F(1, 2))])

    return danger, tooth, lonely, meas, inter, clip_half


def main():
    print("=" * 84)
    print("LRC cap deletion-contraction = measure-valued Claim A (mac-mini-S7)")
    print("=" * 84)

    for n in (6, 8, 10, 14):
        danger, tooth, lonely, meas, inter, clip_half = make(n)
        # a covering-ish S for this n (just need a few speeds to test the identity)
        S = tuple(range(1, n))  # {1,...,n-1}
        print(f"\n--- n={n}, S={S} ---")
        capS = meas(lonely(S))
        # pick a speed s to delete
        for s in (n - 1, n // 2 if (n // 2) in S else 3):
            if s not in S:
                continue
            Sm = tuple(x for x in S if x != s)
            capSm = meas(lonely(Sm))
            Lm = lonely(Sm)
            # deletion-contraction RHS
            marg = meas(inter(Lm, danger(s)))
            dc_ok = (capSm - capS == marg)
            # factor 2 (sigma symmetry): marg == 2 * meas(lonely(Sm) cap D_s cap [0,1/2))
            half = meas(clip_half(inter(Lm, danger(s))))
            f2_ok = (marg == 2 * half) or abs(float(marg - 2 * half)) < 1e-12
            # tooth decomposition: marg == sum_k meas(lonely(Sm) cap tooth_k)
            tooth_sum = sum((meas(inter(Lm, tooth(s, k))) for k in range(s + 1)), F(0))
            tooth_ok = (marg == tooth_sum)
            print(f"  delete s={s}: cap(S\\s)-cap(S) = {capSm-capS} ; marginal danger = {marg}  "
                  f"[DC {'OK' if dc_ok else 'FAIL'}]")
            print(f"      factor-2 (sigma): marg = 2*half? half={half}  [{'OK' if f2_ok else 'FAIL'}]  "
                  f"; tooth-sum decomposition [{'OK' if tooth_ok else 'FAIL'}], #teeth={s+1}")
            # per-tooth mu_LRC (on the half-circle)
            mus = [meas(clip_half(inter(Lm, tooth(s, k)))) for k in range(s + 1)]
            nz = [(k, m) for k, m in enumerate(mus) if m != 0]
            print(f"      mu_LRC(s,k) nonzero teeth (half-circle): {[(k, float(m)) for k,m in nz]}")

        # chain rule: cap(S) = 1 - sum of conditional dangers along a peel order
        order = list(S)
        chain = F(1); prefix = []
        for s in order:
            Lp = lonely(tuple(prefix))
            chain -= meas(inter(Lp, danger(s)))
            prefix.append(s)
        print(f"  chain rule: 1 - sum(conditional dangers) = {chain} ; cap(S) = {capS}  "
              f"[{'OK' if chain == capS else 'FAIL'}]")

    print("\n" + "=" * 84)
    print("Reframe: the LRC cap satisfies a deletion-contraction cap(S\\s)-cap(S) =")
    print("2 * sum_teeth mu_LRC -- the measure-valued analog of Claim A H(T)-H(T-v)=2 sum_C mu(C).")
    print("The s teeth are the 'odd cycles through s'; the 2 is the complement Z_2 (sigma).")
    print("The chain rule cap = 1 - sum(conditional dangers) is the measure-valued OCF/Redei.")
    print("=" * 84)


if __name__ == "__main__":
    main()
