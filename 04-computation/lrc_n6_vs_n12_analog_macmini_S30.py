#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S30 -- WHY is N=6 nonempty but N=12 empty?  The exact
kps-S39 order-3 analog (dilated AP + 2 borders), scaled N=6 -> N=12, with the
OVERSHOOT mechanism (which resonance beats the intended gap value).

N=6 member = {1,5,6,11,16,17} = core {1,6,11,16} (d=5, N-2=4 terms)
             + inner border 5 (=core[1]-1) + outer border 17 (=core[-1]+1),
             M = 5/33 = s/(Ns+k) with s=5, k=3 (order 3), spacing d = s = 5.

The EXACT analog at N=12: core = (N-2)=10-term dilated AP spacing d, + inner
border core[1]-1 + outer border core[-1]+1.  For each spacing d, compute M and,
if it OVERSHOOTS the gap, the competing denominator q' and the resonance
(pairwise sum/diff/double) that carries it -- the opus-S119/HYP-4592 overshoot.
"""
import itertools
from fractions import Fraction as F
from math import gcd


def _dens(W):
    d = set()
    for v, w in itertools.combinations(W, 2):
        d.add(v + w)
        if v != w:
            d.add(abs(v - w))
    for v in W:
        d.add(2 * v)
    d.discard(0)
    return d


def exact_M_arg(W):
    best = F(0)
    arg = None
    seen = set()
    for s in _dens(W):
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
                arg = t
    return best, arg


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(set(v // g for v in W))) if g > 1 else tuple(sorted(set(W)))


def analog_family(N, d):
    """core = (N-2)-term dilated AP {1, 1+d, ..., 1+(N-3)d}, + inner border
    (core[1]-1) + outer border (core[-1]+1)."""
    core = [1 + i * d for i in range(N - 2)]
    inner = core[1] - 1
    outer = core[-1] + 1
    W = primitive(tuple(sorted(set(core) | {inner, outer})))
    return W if len(W) == N else None


def resonance_label(W, q):
    """which pairwise sum / difference / double equals q (the overshoot carrier)."""
    labels = []
    for v, w in itertools.combinations(W, 2):
        if v + w == q:
            labels.append(f"{v}+{w}")
        if abs(v - w) == q:
            labels.append(f"|{v}-{w}|")
    for v in W:
        if 2 * v == q:
            labels.append(f"2*{v}")
    return labels[:3]


def main():
    LO6, HI6 = F(1, 7), F(2, 13)
    LO12, HI12 = F(1, 13), F(2, 25)
    print("=" * 88)
    print("N=6 reference: the kps-S39 order-3 analog for each spacing d")
    print("=" * 88)
    print(f"  gap (1/7, 2/13) = ({float(LO6):.4f}, {float(HI6):.4f})")
    for d in range(2, 8):
        W = analog_family(6, d)
        if W is None:
            continue
        M, t = exact_M_arg(W)
        ingap = LO6 < M < HI6
        tag = "  <-- IN GAP (order-3 member)" if ingap else ""
        print(f"  d={d}: W={list(W)} M={M} ({float(M):.4f}) t={t}{tag}")

    print()
    print("=" * 88)
    print("N=12: the SAME analog (10-term dilated AP + 2 borders) -- does any land in gap?")
    print("=" * 88)
    print(f"  gap (1/13, 2/25) = ({float(LO12):.5f}, {float(HI12):.5f});  order-3 targets:")
    print(f"    5/63={float(F(5,63)):.5f} (d=5,k=3), 4/51={float(F(4,51)):.5f} (d=4,k=3)")
    for d in range(2, 10):
        W = analog_family(12, d)
        if W is None:
            print(f"  d={d}: (degenerate / not 12 distinct primitive speeds)")
            continue
        M, t = exact_M_arg(W)
        q = M.denominator
        ingap = LO12 < M < HI12
        if ingap:
            print(f"  d={d}: W={list(W)} M={M} ({float(M):.5f}) t={t}  <== IN GAP!!!")
        else:
            pos = "ABOVE" if M >= HI12 else "below"
            carr = resonance_label(W, q) if M >= HI12 else []
            carrs = f"  overshoot q={q} via {carr}" if carr else (f"  q={q}" if M >= HI12 else "")
            print(f"  d={d}: M={M} ({float(M):.5f}) [{pos} gap] t={t}{carrs}")
            print(f"        W={list(W)}")

    print()
    print("=" * 88)
    print("THE MECHANISM: why the intended order-3 value is beaten at N=12")
    print("=" * 88)
    print("  For each N=12 analog that overshoots, the intended gap value is s/(12s+k);")
    print("  report the intended value, the ACTUAL M, and the overshoot margin.")
    for d in range(3, 8):
        W = analog_family(12, d)
        if W is None:
            continue
        M, t = exact_M_arg(W)
        # intended order-3 value: s=d, k=3 => d/(12d+3)
        intended = F(d, 12 * d + 3)
        beats = M > intended
        print(f"  d={d}: intended (order-3) {intended} ({float(intended):.5f}) in-gap="
              f"{LO12 < intended < HI12}; ACTUAL M={M} ({float(M):.5f}); "
              f"M {'BEATS' if beats else '<='} intended by {float(M-intended):+.5f}")


if __name__ == "__main__":
    main()
