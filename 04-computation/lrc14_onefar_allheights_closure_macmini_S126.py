#!/usr/bin/env python3
"""One-far all-heights gap closure — finite half (mac-mini-S126).

THE CLOSURE: [plateau-dodging lemma, k=1: any 12-subset S ⊂ {1..13} has
M(S) ≥ 1/13 (settled LRC(13)); plateau at level 3/41 has |I| ≥ 4/6929;
one far x ≥ 297 dodges: M(S∪{x}) ≥ 3/41]  +  [THIS SCAN: every one-far
family with 14 ≤ x ≤ 296, NO filters, exhaustive].

If the scan shows no family with M ∈ (1/14, 3/41), then NO one-far family
at ANY height has M in the open gap — the first shape class closed at all
heights for the gap question (modulo the lemma's independent read).
Expected known full-scan hit: GW = {1..11,13}∪{24} at exactly 1/14 (tight).
"""

from fractions import Fraction as F
from itertools import combinations
import sys

sys.path.insert(0, "04-computation")
from lrc14_dyadic_tower_ladders_macmini_S124 import exact_M

ONE14 = F(1, 14)
GAP_HI = F(3, 41)
XMAX = 296


def core_tables(S, qmax):
    core = [None, None]
    for q in range(2, qmax + 1):
        row = [q] * q
        for u in S:
            um = u % q
            for k in range(q):
                r = (um * k) % q
                d = q - r if q - r < r else r
                if d < row[k]:
                    row[k] = d
        core.append(row)
    return core


def main():
    qmax = 13 + 2 * XMAX
    full_hits = []
    n = 0
    for a in range(1, 14):
        S = [u for u in range(1, 14) if u != a]
        core = core_tables(S, qmax)
        for x in range(14, XMAX + 1):
            n += 1
            exited = False
            for q in range(2, qmax + 1):
                row = core[q]
                xm = x % q
                for k in range(1, q):
                    d = row[k]
                    r = (xm * k) % q
                    dx = q - r if q - r < r else r
                    if dx < d:
                        d = dx
                    if 41 * d >= 3 * q:
                        exited = True
                        break
                if exited:
                    break
            if not exited:
                W = sorted(set(S) | {x})
                M, t, act, pairs = exact_M(W)
                tag = ("TIGHT (known iff GW)" if M == ONE14 else
                       "IN-GAP !!!" if ONE14 < M < GAP_HI else
                       "BELOW FLOOR?!" if M < ONE14 else "above")
                full_hits.append((W, M, tag))
                print(f"  full-scan: removed {a}, x={x}: M = {M}  => {tag}")
    print(f"\n{n} one-far families scanned exhaustively (x ≤ {XMAX}).")
    ingap = [h for h in full_hits if ONE14 < h[1] < GAP_HI]
    below = [h for h in full_hits if h[1] < ONE14]
    tights = [h for h in full_hits if h[1] == ONE14]
    gw = sorted(set(range(1, 12)) | {13, 24})
    assert not below, below
    assert all(W == gw for W, _, _ in tights), tights
    print(f"tights found: {len(tights)} (exactly GW: "
          f"{all(W == gw for W, _, _ in tights)})")
    if not ingap:
        print("IN-GAP: NONE.  With the k=1 plateau lemma (x ≥ 297 ⟹ M ≥ 3/41),")
        print("the ONE-FAR STRATUM has no family with M ∈ (1/14, 3/41) at ANY")
        print("height — the first all-heights shape-class closure of the gap.")
    else:
        for W, M, _ in ingap:
            print(f"  IN-GAP: {W}: M = {M}")


if __name__ == "__main__":
    main()
