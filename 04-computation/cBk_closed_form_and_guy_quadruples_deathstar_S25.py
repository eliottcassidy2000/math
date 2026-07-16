#!/usr/bin/env python3
"""death-star-2026-07-16-S25 (HYP-7063/THM-906): (1) THE c_B(k) CLOSED FORM referee;
(2) near-AP quadruples vs Guy's conjecture / A000241.

PART 1 — THEOREM (c_B(k) corner form): for runners v with unique primitive relation
k ⊥ v (Σkᵢvᵢ = 0, no sub-pair relation), section boxes Bᵢ ⊆ Z₇ (intervals E_i = union of
[c/7,(c+1)/7)), the box-hit remainder plateau equals
    c_B(k) = Σ_{j≠0} ∏ᵢ Ê_i(j kᵢ)
           = −(1/(24·∏ᵢkᵢ)) · Σ_{sections cᵢ∈Bᵢ} Σ_{ε∈{0,1}⁴} (−1)^{Σε} B₄({Σᵢ kᵢ(cᵢ+εᵢ)/7})
with B₄(x) = x⁴ − 2x³ + x² − 1/30.  [Derivation: 1̂_{[a,b)}(m) = (e(−ma)−e(−mb))/(2πim);
∏ over i; (2πij)⁴ = (2π)⁴j⁴; Σ_{j≠0} e(jx)/j⁴ = −(2π)⁴B₄({x})/24.]
Referee: exact hit-measure strata (boxeph-S33 machinery) along fixed-relation ladders —
R must converge to c_B(k); plus the height-scaling law and codex's maximizer sectors.

PART 2 — additive quadruples vs Guy Z(n) = (1/4)⌊n/2⌋⌊(n−1)/2⌋⌊(n−2)/2⌋⌊(n−3)/2⌋
(A000241 = conjectured cr(K_n)): linear count Q_lin(n) = #{4-subsets a<b<c<d: a+d=b+c},
circular count on Z_q (parallel-chord pairs), near-AP perturbation response; the exact
relationship, honestly.
"""
from fractions import Fraction as Fr
from math import gcd, comb
from itertools import combinations
import sys, time

def B4(x):
    return x**4 - 2 * x**3 + x**2 - Fr(1, 30)

def frac(x):
    return x - (x.numerator // x.denominator) if isinstance(x, Fr) else x - int(x)

def cB_closed(k, boxes):
    """the corner form; k = 4 ints (nonzero), boxes = lists of sections in Z7."""
    tot = Fr(0)
    prodk = 1
    for ki in k:
        prodk *= ki
    from itertools import product as iproduct
    for secs in iproduct(*boxes):
        for eps in iproduct([0, 1], repeat=4):
            arg = Fr(0)
            for ki, ci, ei in zip(k, secs, eps):
                arg += Fr(ki * (ci + ei), 7)
            arg = arg - (arg.numerator // arg.denominator)  # {arg}
            sgn = (-1) ** sum(eps)
            tot += sgn * B4(arg)
    return -tot / (24 * prodk)

def hit_measure(vs, boxes):
    bps = sorted(set(Fr(kk, 7 * v) for v in vs for kk in range(7 * v + 1)))
    tot = Fr(0)
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        ok = True
        for v, B in zip(vs, boxes):
            if int((v * mid % 1) * 7) not in B:
                ok = False; break
        if ok:
            tot += bps[i + 1] - bps[i]
    return tot

def strata_R(vs, boxes):
    full = hit_measure(vs, boxes)
    prod = Fr(1)
    for B in boxes:
        prod *= Fr(len(B), 7)
    p2 = Fr(0)
    n = len(vs)
    for i in range(n):
        for j in range(i + 1, n):
            pij = hit_measure([vs[i], vs[j]], [boxes[i], boxes[j]])
            dev = pij - Fr(len(boxes[i]) * len(boxes[j]), 49)
            rest = Fr(1)
            for kk in range(n):
                if kk != i and kk != j:
                    rest *= Fr(len(boxes[kk]), 7)
            p2 += dev * rest
    return full - prod - p2

def part1():
    print("PART 1: c_B(k) closed form vs exact strata remainder (fixed-relation ladders)")
    # relation k = (1,1,-1,-1): v = (a, b, c, a+b-c) families, no sub-pair relations
    k = (1, 1, -1, -1)
    boxes = [[0], [2], [4], [6]]
    pred = cB_closed(k, boxes)
    print(f"  k={k} boxes={boxes}: c_B closed = {pred} ≈ {float(pred):.6e}")
    for (a, b, c) in [(8, 13, 10), (16, 27, 20), (31, 54, 41), (61, 108, 82), (121, 215, 163)]:
        v = (a, b, c, a + b - c)
        if len(set(v)) < 4: continue
        R = strata_R(list(v), boxes)
        print(f"    v={v}: R = {float(R):.6e}  (R/c_B = {float(R/pred) if pred else float('nan'):.3f})")
        sys.stdout.flush()
    # height-2 relation k = (1,1,-2,... ) e.g. v1+v2 = 2v3 (AP-type!) with 4th runner free-ish:
    # use k = (1,1,-2,0)? relation must involve all four for the quadruple stratum; take
    # k = (1,1,-1,-1) vs k = (2,1,-2,-1): v = (a,b,c,(2a+b-2c)):
    k2 = (2, 1, -2, -1)
    pred2 = cB_closed(k2, boxes)
    print(f"  k={k2} boxes={boxes}: c_B closed = {pred2} ≈ {float(pred2):.6e}")
    for (a, b, c) in [(9, 14, 11), (17, 30, 21), (33, 58, 41), (65, 114, 81)]:
        d = 2 * a + b - 2 * c
        v = (a, b, c, d)
        if d <= 0 or len(set(v)) < 4: continue
        R = strata_R(list(v), boxes)
        print(f"    v={v}: R = {float(R):.6e}  (R/c_B = {float(R/pred2) if pred2 else float('nan'):.3f})")
        sys.stdout.flush()
    # wide boxes (codex sectors): B = {0,1,2} type
    boxesW = [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
    predW = cB_closed(k, boxesW)
    print(f"  k={k} boxes=[0-2]^4: c_B closed = {predW} ≈ {float(predW):.6e}")
    for (a, b, c) in [(16, 27, 20), (61, 108, 82)]:
        v = (a, b, c, a + b - c)
        R = strata_R(list(v), boxesW)
        print(f"    v={v}: R = {float(R):.6e}  (R/c_B = {float(R/predW) if predW else float('nan'):.3f})")

def Zguy(n):
    return Fr((n // 2) * ((n - 1) // 2) * ((n - 2) // 2) * ((n - 3) // 2), 4)

def part2():
    print("\nPART 2: additive quadruples vs Guy Z(n) / A000241")
    print(f"  {'n':>3} {'Qlin':>6} {'Z(n)':>7} {'Qlin/Z':>7} {'Qcirc(Zn)':>9} {'C(n,4)':>7}")
    for n in range(5, 16):
        S = list(range(1, n + 1))
        Qlin = sum(1 for q in combinations(S, 4) if q[0] + q[3] == q[1] + q[2])
        # circular on Z_n: unordered pairs-of-disjoint-pairs with equal sum mod n
        Qc = 0
        for q in combinations(range(n), 4):
            a, b, c, d = q
            # pairings: {a,b},{c,d} / {a,c},{b,d} / {a,d},{b,c}
            for (p1, p2) in [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]:
                if (p1[0] + p1[1]) % n == (p2[0] + p2[1]) % n:
                    Qc += 1
        Z = Zguy(n)
        print(f"  {n:>3} {Qlin:>6} {str(Z):>7} {float(Qlin/Z):>7.3f} {Qc:>9} {comb(n,4):>7}")
    # near-AP response: kernel 13-sets
    print("  near-AP 13-sets: quadruple counts (v1+v2 = v3+v4 as 4-subsets)")
    import random
    rnd = random.Random(7)
    AP = list(range(1, 14))
    def qcount(S):
        return sum(1 for q in combinations(S, 4) if q[0] + q[3] == q[1] + q[2])
    print(f"    AP {{1..13}}: {qcount(AP)}   (C(13,4) = {comb(13,4)})")
    for tag, S in [("deep well", list(range(1, 13)) + [182]),
                   ("near-AP (GW-ish)", [1,2,3,4,5,7,8,9,10,11,12,13,24]),
                   ("random", sorted(rnd.sample(range(1, 200), 13)))]:
        print(f"    {tag}: {qcount(S)}")

if __name__ == "__main__":
    t0 = time.time()
    part1()
    part2()
    print(f"[total {time.time()-t0:.1f}s]")
