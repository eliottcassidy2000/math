#!/usr/bin/env python3
"""
lrc14_angleC_compaction_opus_0620.py   (opus-2026-06-20)

ANGLE C, sharper:  HOLE-COMPACTION (CA defect annihilation).

OBSERVATION from the criticality run: at EVERY residual span N, the worst (max measS7)
shape with |E|=k is the one whose holes are pushed to the EDGES, leaving a long
consecutive block -- and the global worst over all spans is consec_k itself.  This is
the opposite of "spread the holes": measS7 is MAXIMISED by CLUSTERING the kept clocks
(equivalently clustering the holes at one end).

CONJECTURED COMPACTION LEMMA (CA "defect annihilation"):
   Among all E with |E|=k and 0 in E, measS7 is maximised by E = consec_k = {0..k-1}.
   More locally: any "compaction move" that slides a clock to fill an interior hole
   (reducing the spread / merging two blocks) does NOT decrease measS7.

We TEST a rigorous LOCAL move that would chain to the global statement:

  MOVE T (left-compaction by one): pick the leftmost hole position p (smallest p>=1 with
  p not in E but some e>p in E). Let q = smallest clock > p. Replace q by p, i.e. slide
  that clock left into the hole.  E' = E \ {q} U {p}.  This strictly reduces sum(E) by q-p>0.
  Iterating MOVE T from any E reaches consec_k.  If measS7(E') >= measS7(E) for EVERY
  such move, then consec_k is the global max (chain of >=).

  We ALSO test the stronger BLOCK-SLIDE: slide an entire trailing block left to close a gap.

If MOVE T is monotone (measS7 non-decreasing) we have a PROOF SKELETON: a finite, local,
checkable monotone-move lattice with consec_k at the top.  We probe whether the move is
ALWAYS non-decreasing or has exceptions (recall C2: span-monotonicity is FALSE, so naive
'smaller spread is bigger' is false -- compaction is a DIFFERENT, finer move; test it).

EXACT Fractions; exhaustive over bounded boxes; stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        res = set(int(7*e*xm) % 7 for e in E)
        if len(res) == 7: total += x1 - x0
    return total

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def left_compaction(E):
    """Yield E' for the single-step left-compaction move: slide the first clock that sits
    immediately after a hole into that hole. Generate ALL such single-slides (every clock
    that has an empty position immediately to its left, slide it left by 1)."""
    E = sorted(set(E)); S = set(E)
    moves = []
    for e in E:
        if e-1 >= 0 and (e-1) not in S:
            Ep = sorted((S - {e}) | {e-1})
            moves.append(tuple(Ep))
    return moves

def main():
    print("="*96)
    print("ANGLE C HOLE-COMPACTION (defect annihilation)  (opus-2026-06-20)")
    print("="*96)

    # ---- TEST 1: is the single-step left-slide (e -> e-1 into a hole) monotone up? ----
    print("(1) single-step left-slide monotonicity:  E' = E with one clock slid left by 1")
    print("    into an adjacent hole.  Check measS7(E') >= measS7(E) (compaction increases).")
    for k in [8, 9, 10]:
        ck = CAP[k]
        viol = 0; tested = 0; worst_drop = F(0); worst_pair = None
        # exhaustive over all k-subsets of {0..N} containing 0, bounded box
        Nmax = {8:13, 9:13, 10:13}[k]
        for body in itertools.combinations(range(1, Nmax+1), k-1):
            E = (0,)+body
            mE = measS7(E)
            for Ep in left_compaction(E):
                if 0 not in Ep:    # keep 0 pinned (slide of clock 1 -> 0 keeps 0; clock at 0 won't slide)
                    pass
                mEp = measS7(Ep)
                tested += 1
                if mEp < mE:        # compaction DECREASED measS7 -> violation
                    viol += 1
                    drop = mE - mEp
                    if drop > worst_drop:
                        worst_drop = drop; worst_pair = (E, Ep, mE, mEp)
        print(f"  k={k}: tested {tested} slides; compaction-DECREASES (violations) = {viol}")
        if worst_pair:
            E,Ep,mE,mEp = worst_pair
            print(f"        worst violation: {E}({float(mE):.5f}) -> {Ep}({float(mEp):.5f}) "
                  f"DROP {float(worst_drop):.5f}")
    print()

    # ---- TEST 2: does consec_k dominate ALL k-subsets (the crux), exhaustively? ----
    print("(2) DIRECT crux re-verify: consec_k = argmax measS7 over ALL k-subsets, box max(E)<=13")
    for k in [8, 9, 10]:
        consec_val = measS7(tuple(range(k)))
        best = F(-1); bestE = None; ties = 0; beats = 0
        for body in itertools.combinations(range(1, 14), k-1):
            E = (0,)+body
            v = measS7(E)
            if v > best: best, bestE = v, E
            if v > consec_val: beats += 1
            if v == consec_val and E != tuple(range(k)): ties += 1
        print(f"  k={k}: consec measS7={float(consec_val):.6f}  global max in box={float(best):.6f} "
              f"at {bestE}  #strictly-beating-consec={beats}  #ties={ties}")
    print()

    # ---- TEST 3: COMPACTION ORDER -- is measS7 monotone along the "merge gaps" partial order? ----
    # partial order: E <= E' if E' is obtained from E by left-slides (E' more compact).
    # We test a STRONGER claim that would give the proof: measS7 is monotone w.r.t. the
    # number of gaps AND their total length jointly.  Specifically test the "fill any single
    # interior hole by importing the next clock" reduces toward consec and never decreases.
    print("(3) BLOCK-SLIDE: slide entire trailing block left to close the FIRST gap.")
    print("    E = [0..a-1] gap [g..] ... ; move the whole tail left by gap size.")
    for k in [8,9,10]:
        viol = 0; tested = 0; worst = F(0); wp=None
        Nmax = {8:13,9:13,10:13}[k]
        for body in itertools.combinations(range(1, Nmax+1), k-1):
            E = list((0,)+body); S=sorted(E)
            # find first gap
            gap_at = None
            for i in range(len(S)-1):
                if S[i+1] > S[i]+1:
                    gap_at = i; break
            if gap_at is None:
                continue  # already consec
            g = S[gap_at+1]-S[gap_at]-1  # gap size
            # slide tail (S[gap_at+1:]) left by g
            Ep = tuple(S[:gap_at+1] + [x-g for x in S[gap_at+1:]])
            mE = measS7(tuple(S)); mEp = measS7(Ep)
            tested += 1
            if mEp < mE:
                viol += 1
                if mE-mEp > worst: worst = mE-mEp; wp=(tuple(S),Ep,mE,mEp)
        print(f"  k={k}: block-slide tested {tested}; DECREASES = {viol}")
        if wp:
            S,Ep,mE,mEp = wp
            print(f"        worst: {S}({float(mE):.5f}) -> {Ep}({float(mEp):.5f}) DROP {float(worst):.5f}")
    print("\nDONE.")

if __name__ == "__main__":
    main()
