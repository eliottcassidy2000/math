#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S13 -- HYP-4392: CLOSING THE TIGHT-SIDE STRICT LIFT-RIGIDITY.

residue_pinning_13 (GREEN in the corpus): a tight (M=1/13) 12-family with no
multiple of 13 has residues EXACTLY {1..12} mod 13, injectively.  The remaining
tight-side piece: from that + M=1/13 conclude S is a DILATED AP c*{1..12}; and
since c*{1..12} is primitive only for c=1, PRIMITIVE tight => S = {1..12}.

This investigates:
 (1) is the tight locus EXACTLY the dilated APs?  (search primitive tight
     families != {1..12}; expect NONE)
 (2) THE ESCAPE WITNESS: for a lift of the AP (residues {1..12} mod 13, but not
     the AP), WHERE does M exceed 1/13?  Characterize the witness t = a/q --
     is it at a bounded, structured modulus?
 (3) HEIGHT-BOUNDEDNESS: does 'residues {1..12} mod 13 + M=1/13' force bounded
     height?  If yes, 'tight => AP' reduces to a FINITE check (closeable).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def exact_M_argmax(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); bt = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best, bt = mv, t
    return best, bt

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

FLOOR = F(1, 13)

def part1_tight_locus():
    print("=" * 78)
    print("PART 1: is the tight locus EXACTLY the dilated APs? (primitive tight != AP?)")
    print("=" * 78)
    random.seed(13)
    AP = tuple(range(1, 13))
    # dilated APs c*{1..12} -- confirm tight, confirm primitive only c=1
    print("  dilated APs c*{1..12}:")
    for c in (1, 2, 3, 5, 7):
        W = tuple(c * i for i in range(1, 13))
        M, t = exact_M_argmax(W)
        print(f"    c={c}: M={M} {'= 1/13 TIGHT' if M==FLOOR else ''}  primitive={primitive(W)==W or gcd(*W)==1}  gcd={gcd(*W)}")
    # search for PRIMITIVE tight families != {1..12}
    found = []
    tested = 0
    fl = float(FLOOR)
    # (a) families with residues {1..12} mod 13 (lifts), primitive
    for _ in range(60000):
        perm = list(range(1, 13)); random.shuffle(perm)
        W = [perm[i] + 13 * random.randint(0, 3) for i in range(12)]
        W = primitive(tuple(sorted(set(W))))
        if len(W) != 12 or W == AP:
            continue
        tested += 1
        if float_M(W) < fl + 1e-7:      # candidate tight
            M, t = exact_M_argmax(W)
            if M == FLOOR:
                found.append(W)
    print(f"  primitive residue-{{1..12}} lifts tested: {tested}")
    print(f"  primitive TIGHT families != {{1..12}}: {len(found)}")
    for W in found[:6]:
        print(f"    *** {list(W)}")
    if not found:
        print("    NONE -- primitive tight locus = {1..12} uniquely (on this search).")
        print("    (dilated APs are the only tight families; primitive => c=1 => the AP)")

def part2_escape_witness():
    print()
    print("=" * 78)
    print("PART 2: the ESCAPE WITNESS -- where does a lift of the AP exceed 1/13?")
    print("=" * 78)
    AP = list(range(1, 13))
    random.seed(131)
    witness_q = {}
    margins = []
    for _ in range(1200):
        # single-element lift: replace one element by +13k
        W = list(AP)
        idx = random.randrange(12)
        W[idx] += 13 * random.randint(1, 3)
        W = primitive(tuple(sorted(set(W))))
        if len(W) != 12:
            continue
        M, t = exact_M_argmax(W)
        if M > FLOOR:
            witness_q[t.denominator] = witness_q.get(t.denominator, 0) + 1
            margins.append(M)
    print(f"  single-element AP-lifts with M > 1/13: {len(margins)}")
    if margins:
        mn = min(margins)
        print(f"  MIN margin over lifts: {mn} = {float(mn):.5f} (all strictly > 1/13={float(FLOOR):.5f}? "
              f"{'YES' if mn > FLOOR else 'NO'})")
        print(f"  escape-witness denominator q distribution (top):")
        for q in sorted(witness_q, key=lambda x: -witness_q[x])[:12]:
            print(f"    q={q}: {witness_q[q]}")

def part3_height_bound():
    print()
    print("=" * 78)
    print("PART 3: HEIGHT-BOUNDEDNESS -- does M stay > 1/13 + delta as lift height grows?")
    print("(if M - 1/13 is bounded BELOW away from 0 for non-AP lifts, tight=>AP is finite)")
    print("=" * 78)
    AP = list(range(1, 13))
    random.seed(132)
    by_height = {}
    for _ in range(9000):
        W = list(AP)
        nlift = random.randint(1, 4)
        maxk = 0
        for _i in range(nlift):
            idx = random.randrange(12); k = random.randint(1, 5)
            W[idx] += 13 * k; maxk = max(maxk, k)
        W = primitive(tuple(sorted(set(W))))
        if len(W) != 12:
            continue
        fm = float_M(W)
        band = "1-2" if maxk <= 2 else ("3-4" if maxk <= 4 else "5-8")
        by_height.setdefault(band, []).append(fm)
    print("  min (M - 1/13) by lift-height band (nonzero AP-lifts):")
    fl = float(FLOOR)
    for band in ("1-2", "3-4", "5-8"):
        if band in by_height:
            gaps = [m - fl for m in by_height[band]]
            print(f"    height {band}: count={len(gaps)}, MIN(M-1/13)={min(gaps):.5f}, "
                  f"median={sorted(gaps)[len(gaps)//2]:.5f}")
    print("  => if MIN(M-1/13) stays bounded away from 0 across heights, the tight")
    print("     locus is height-isolated: no lift approaches tightness => tight=>AP finite.")

if __name__ == "__main__":
    part1_tight_locus()
    part2_escape_witness()
    part3_height_bound()
