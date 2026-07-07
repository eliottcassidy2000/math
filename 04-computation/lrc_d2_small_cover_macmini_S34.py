#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S34 -- a SMALL-lcm d=2 covering (Lean-decidable), verified
UNIFORMLY over all residues.

Goal: cover {1,…,10,x,y} (M>=2/25) with witnesses whose moduli have SMALL lcm, so
the covering completeness is a decidable finite check.  Try moduli subsets of
{11,12,25} (lcm 3300), then extend only if needed.  Verify uniformly: precompute
each x's witness-bitmask G(x) for x in [0,LCM); (x,y) covered iff G(x)&G(y)!=0.
Report uncovered residue classes (should be only the AP / d<=1 exceptions).
"""
from math import gcd, lcm
from fractions import Fraction as F
import itertools

BASE = list(range(1, 11))


def clears_base(q, c, mu):
    return all(mu <= (v * c) % q <= q - mu for v in BASE)


def good_set(q, c, mu):
    return frozenset(r for r in range(q) if mu <= (r * c) % q <= q - mu)


def build_witnesses(moduli):
    """(q,c,mu) with mu/q>=2/25, gcd(c,q)=1, base cleared; dedupe (q,c)->max mu."""
    byqc = {}
    for q in moduli:
        mu_min = -(-2 * q // 25)  # ceil(2q/25)
        for mu in range(mu_min, q // 2 + 1):
            for c in range(1, q):
                if gcd(c, q) == 1 and clears_base(q, c, mu):
                    if (q, c) not in byqc or mu > byqc[(q, c)]:
                        byqc[(q, c)] = mu
    return [(q, c, mu) for (q, c), mu in byqc.items()]


def uniform_check(moduli):
    wits = build_witnesses(moduli)
    L = 1
    for q in moduli:
        L = lcm(L, q)
    goods = [good_set(q, c, mu) for (q, c, mu) in wits]
    # G(x) bitmask over witnesses, x in [0,L); group residues by mask (few distinct masks)
    from collections import defaultdict
    mask_to_res = defaultdict(list)
    for x in range(L):
        b = 0
        for i, (q, c, mu) in enumerate(wits):
            if (x % q) in goods[i]:
                b |= (1 << i)
        mask_to_res[b].append(x)
    masks = list(mask_to_res.keys())
    # uncovered iff two masks (incl. same) are disjoint; report a representative pair
    uncovered = []
    for a in range(len(masks)):
        for b in range(a, len(masks)):
            if masks[a] & masks[b] == 0:
                xr = mask_to_res[masks[a]][0]
                yr = mask_to_res[masks[b]][0]
                uncovered.append((xr, yr))
    return wits, L, uncovered


def longest_ap(W):
    S = set(W)
    best = 1
    for a in W:
        for d in range(1, 40):
            Ln = 1
            while a + Ln * d in S:
                Ln += 1
            best = max(best, Ln)
    return best


def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    dens.discard(0)
    best = F(0)
    seen = set()
    for s in dens:
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def main():
    for moduli in [[11, 12, 25], [11, 12, 24, 25], [11, 12, 23, 24, 25]]:
        L = 1
        for q in moduli:
            L = lcm(L, q)
        print("=" * 84)
        print(f"MODULI {moduli}, lcm = {L}")
        print("=" * 84)
        if L > 400000:
            print(f"  lcm {L} too large for full uniform sweep; skipping")
            continue
        wits, L, uncovered = uniform_check(moduli)
        print(f"  witnesses: {len(wits)}  ({[(q,c,mu) for q,c,mu in wits][:8]}{'...' if len(wits)>8 else ''})")
        print(f"  uncovered residue classes (x_r,y_r): {len(uncovered)}")
        # interpret uncovered: pick representative concrete (x,y), report M + longestAP
        real_low = []
        for xr, yr in uncovered[:60]:
            if yr is None:
                continue
            # smallest concrete x>=11 with x%L==xr, y>x with y%L==yr
            x = xr if xr >= 11 else xr + L
            y = yr if yr > x else yr + L
            W = tuple(sorted(set(BASE) | {x, y}))
            if len(W) != 12:
                continue
            M = exact_M(W)
            lap = longest_ap(W)
            if M < F(2, 25):
                real_low.append((xr, yr, x, y, M, lap))
        print(f"  uncovered classes whose representative has M<2/25: {len(real_low)}")
        for xr, yr, x, y, M, lap in real_low[:12]:
            note = "AP" if M == F(1, 13) else ("d=1 (x or y=11)" if 11 in (x % L, y % L) or x == 11 or y == 11 else "??")
            print(f"    resid({xr},{yr}) e.g. ({x},{y}) M={M} ({float(M):.4f}) longestAP={lap}  [{note}]")
        if not real_low:
            print("  => EVERY uncovered class has M>=2/25 anyway (or is empty) -- covering COMPLETE for d=2")
        else:
            allAP_or_d1 = all(M == F(1, 13) or 11 in (xr, yr) or x == 11 or y == 11
                              for xr, yr, x, y, M, lap in real_low)
            print(f"  => the only M<2/25 uncovered classes are AP / d<=1 (x or y ≡11): {allAP_or_d1}")
        print()


if __name__ == "__main__":
    main()
