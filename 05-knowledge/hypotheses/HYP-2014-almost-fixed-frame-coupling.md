---
id: HYP-2014
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S533
related:
  - HYP-2012
  - HYP-2011
---

# HYP-2014: 2^floor(n/2) | A000568(n) (pair-flip factor); clean frame-factorization only n<=4; the breakdown = the LRC coupling

**VERIFIED (`almost_fixed_frame_s533.py`):**
1. `2^floor(n/2) | A000568(n)` for n=3..12; quotients 1,1,3,7,57,430,11971,...
   v_2(A000568(n)) = floor(n/2) for odd n and even n<8, but EXCEEDS it (extra 2) for
   even n=8,10,12 -> the pair-flip factor is a lower bound on the 2-part.
2. n=3,4: quotient=1 -> a SINGLE fixed frame + floor(n/2) pair-flips realizes ALL
   iso classes (A000568(3)=2^1, A000568(4)=2^2). The user's picture, exactly true.
3. n=5: 12=3x4 but NO partition into 3 disjoint pair-flip blocks exists (136 of 256
   frames give a full 4-block, all overlapping); min cover=5 frames, Hamming-spread
   4/8. So neither a fixed nor an almost-fixed frame realizes the factorization.

**INTERPRETATION:** the frame cannot stay (almost) fixed for n>=5 because flipping
one independent pair changes how the others sit on the shared frame -> the pairs
are COUPLED. This is the same fact as: the n=4 parity law fails for n>=6 because the
floor(n/2) channels stay independent and couple (S531/S532). 'Frame won't stay
fixed' = 'pairs coupled' = 'no single-congruence parity law' = LRC hard.

**CLAIM:** the 2^floor(n/2) divisibility reflects a free (Z/2)^floor(n/2) shadow on
the iso-class set even when no fixed frame displays it; the multi-channel program is
to MEASURE the coupling = the deviation of the true loneliness/iso structure from
the free 2^floor(n/2) product, channel by channel.

**OPEN:** (A) prove 2^floor(n/2) | A000568(n) in general (and the even-n>=8 extra 2);
is it a known 2-adic result? (B) define the coupling tensor (pairwise frame
interaction) and relate its vanishing to the inside-debt vanishing (S531).

**Files:** `04-computation/almost_fixed_frame_s533.py` (+.out). Reflection:
`07-reflections/almost-fixed-frame-and-the-coupling-s533.md`.
