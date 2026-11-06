# Scalar unit-start synthesis lead

Status: FINITE-EXACT, independently audited aggregate result. No eventual threshold proof. The additional recipe packaging also passed independent review.

## Inheritance

The physical O/Z/J laws and complete unit-start grammar are inherited from Session203 rooted-cover closure and its Session198/200 source-classification dependency. Retain only (h,R,e,d), initialized at (1,2,0,1) for the rooted directed edge:

- O: (h,2R,0,d).
- Z: (h+R-e,2R-d,d,d).
- J: (h+R-e,h+R,h-e,h).

For positive-start actual sources, R>d>=e, so h increases strictly under Z/J. Its increment is even by tournament parity. O preserves h. Removing terminal O does not change the scalar target.

## Complete target pruning

At target bound M, a Z/J child is kept exactly when h+R-e<=M. An O child is needed only when h+2R<=M: if this fails, no subsequent first non-O letter can remain below M. Until that letter, further O operations only double R. Thus no represented count at most M is lost. Equal four-states can be merged because all continuations depend only on those coordinates.

The production exploration uses `while cur`, never a guessed order cap. The run empties at order21 for M=50000. Independently, deleting trailing O from any nonempty recipe leaves a last Z/J child, which is strong: the positive-start root reaches every old vertex, every old nonroot points to the newborn, and the newborn points to that root. Busch's exponential strong minimum also implies finiteness at every fixed M. Neither finite termination nor coverage of a long finite interval is an all-order occupancy proof.

## Signal

At M=50000 the support has24971 odd values. The29 missing values are

7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

Thus every odd2649..49999 appears in this physical grammar. This does not say the full tournament spectrum excludes those29 values.

## Reproduction and actual witnesses

    python3 verify.py
    g++ -O2 -std=c++17 audit/independent.cpp -o independent
    ./independent
    python3 audit/literal.py
    g++ -O2 -std=c++17 witnesses.cpp -o witnesses
    ./witnesses
    python3 check_witnesses.py

WITNESSES.tsv supplies one literal O/Z/J recipe for every supported count. The recipe checker constructs each actual tournament, validates every edge orientation, and certifies its count by the inherited all-input transfer. It does not claim literal subset-DP recounting of every large witness. The independent literal control separately recounts all364 recipes through order7, including their root-start count. Running `python3 check_witnesses.py 2649` additionally prints that target's literal outgoing rows and root.

The fixed witness generator uses a64-bit recipe encoding and stops with a nonzero error before that encoding could overflow. It contains no fixed search-depth cutoff and has ample capacity for the actual maximum order21. The Python scalar enumerator and independent C++ enumerator exhaust their frontiers, without a guessed depth cap.

The frozen dependency provides the four-coordinate formulas and complete unit-start source grammar. The independent audit records its provenance. Scalar completeness for the unit-start subclass uses this inherited source theorem; scalar existence of each listed recipe needs only the actual physical operations and exact transfer law.

## Remaining all-order obligation

Every odd h>=2649 having a unit-start rooted tournament is an OPEN sufficient lemma, supported but not proved by this finite experiment. No induction, interval closure, or infinite invariant family is supplied. The29 missing values are exclusions from the unit-start subclass only, not exclusions from the unrestricted tournament spectrum. Nothing above49999 is asserted by this certificate.
