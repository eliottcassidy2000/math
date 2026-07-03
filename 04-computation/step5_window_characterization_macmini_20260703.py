#!/usr/bin/env python3
"""
STEP 5 characterization (mac-mini-2026-07-03-S25): the far-peel (far_peel_lonely, kps) closes CoveringFarLonely
for w > threshold = #pieces(base) / (3 * length(base)) [from hbig: #pieces*4h < (1-2h)*length*w, h=1/14].
Step 5 = the residual 22 < w <= threshold. For a covering family, peel the LARGEST far runner w; classify:
 - step 3 (peel closes): w > threshold(base of the other 12).
 - step 5 (residual): 22 < w <= threshold.
Question: are step-5 families BOUNDED (finite census extension) or COMPRESSED/large-magnitude (need THM-608
renormalization)? Compute #pieces + length of the base good region and the threshold for real families.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def good_pieces_and_length(base, h=1/14, Ng=400000):
    """approximate #pieces (connected safe intervals) and measure of {t: min_i ||b_i t||>=h}."""
    inside = False; pieces = 0; cnt = 0
    for k in range(Ng):
        t = (k + 0.5)/Ng
        safe = min(nd(b*t) for b in base) >= h
        if safe:
            cnt += 1
            if not inside: pieces += 1; inside = True
        else:
            inside = False
    return pieces, cnt/Ng

if __name__ == "__main__":
    rng = random.Random(5)
    print("STEP 5 characterization: peel largest far w; threshold = #pieces/(3*length(base)). step3 vs step5.")
    print("=" * 96)
    print(f"{'family (sorted)':>30} {'w=max':>7} {'#pieces':>8} {'len(base)':>10} {'threshold':>10} {'route':>8}")
    n3 = n5 = 0; step5_maxspeeds = []
    for _ in range(40):
        sp = sorted(set(rng.sample(range(1, 400), 13)))
        if len(sp) != 13 or not is_covering(sp) or gcd_all(sp) != 1 or not any(v > 22 for v in sp):
            continue
        w = sp[-1]; base = sp[:-1]
        pieces, length = good_pieces_and_length(base)
        thr = pieces / (3 * length) if length > 0 else float('inf')
        route = "step3" if w > thr else "step5"
        if route == "step3": n3 += 1
        else: n5 += 1; step5_maxspeeds.append(w)
        if n3 + n5 <= 16:
            print(f"{str(sp[:6])+'..':>30} {w:>7} {pieces:>8} {length:>10.4f} {thr:>10.1f} {route:>8}")
    print(f"\nover the sample: step3 (peel closes) = {n3}, step5 (residual) = {n5}")
    if step5_maxspeeds:
        print(f"step-5 max-speeds: min={min(step5_maxspeeds)}, max={max(step5_maxspeeds)} -- these are the")
        print(f"   COMPRESSED/comparable families (largest runner not far above the base threshold).")
    print("\n=> step 5 = the compressed case (w not >> base). Bounded per family, but the base magnitude is")
    print("   unbounded, so step 5 is NOT a fixed finite window -- it is exactly the compressed/near-equal")
    print("   residual handled by THM-608 (renormalization, near-equal cluster) + the bounded-magnitude census.")
    print("   Peeling the LARGEST helps only for DOMINANT far runners; comparable-magnitude families need THM-608.")
