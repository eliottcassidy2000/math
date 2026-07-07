---
id: THM-648
title: BLUE SELF-LOOP LINES EXIST ONLY AT EVEN n — at odd n no grid-symmetric tiling has its complement tiling in the same isomorphism class (proof via THM-644's ρ-score relation + THM-646's score-complement law); at even n the necessary condition is s(1) = (n−2)/2 (necessary, NOT sufficient: 4 of 24 candidates at n=6)
status: PROVED (odd-n impossibility, half page below; verified n=3,5,7 exhaustively: 0 loops). Even-n side: the s(1) = (n−2)/2 condition PROVED necessary, shown NOT sufficient (n=6: 24 candidates, 4 loop-endpoints); existence at n=4,6 by census (1 and 2 loops). Settles THM-643-C2 = klein-S161-C1.
source: mac-mini-2026-07-07-S48 (HYP-5047; owner: work the shaped targets)
depends_on:
  - THM-644   # opus-S139: gridsym ⟺ ρ(i)=n+1−i is an anti-automorphism
  - THM-646   # mac-mini-S47: the score-complement law s + s' = c
related:
  - THM-643   # C2 was posed there (and as klein-S161 C1)
  - THM-647   # monad-S8: Anti-Rédei (the fiber-side parity companion)
---

# THM-648 — Blue self-loops only at even n

**Statement.** A *blue self-loop* is a line {t, flip(t)} with t grid-symmetric and both
endpoints in the same isomorphism class. For ODD n there are none. For EVEN n, any blue
self-loop's tilings satisfy s(1) = (n−2)/2 (base-path labeling); this is necessary but
not sufficient.

**Proof (odd-n impossibility).** Let t be grid-symmetric with labeled scores s, and let
a = s(1).
1. By THM-644, ρ(i) = n+1−i is an anti-automorphism of t's tournament, so
   s(ρ(v)) = (n−1) − s(v) for every v; in particular s(n) = n−1−a, and the multiset
   identity **{ (n−1) − s(v) : v } = { s(v) : v }** holds.
2. By THM-646, the flip partner's scores are s'(v) = c(v) − s(v) with
   c = (n−2, n−1, …, n−1, n). Hence, as multisets,
   {s'} = {(n−1) − s(v)} with the v=1 entry lowered by 1 and the v=n entry raised by 1
   = {s} − {n−1−a, a} + {n−2−a, a+1} (using step 1 twice; both removed values are
   genuinely present: a = s(1), n−1−a = s(n)).
3. If flip(t) lies in the same class, {s'} = {s}, forcing the exchanged pairs to agree:
   {n−1−a, a} = {n−2−a, a+1}. The match n−1−a = n−2−a is impossible, so
   n−1−a = a+1 and a = n−2−a, both giving **a = (n−2)/2 — an integer only for even n**. ∎

**Verification (exhaustive, n=3..7).** The ρ-score relation holds on every grid-symmetric
tiling (2/4/16/64/512 of them); the step-2 bookkeeping identity holds on all; blue
self-loop counts are 0/1/0/2/0 with every even-n loop satisfying s(1) = (n−2)/2; at n=6,
24 grid-symmetric tilings satisfy the condition but only 4 are loop endpoints (the two
loops sit at classes with score multisets {2,2,2,3,3,3} — the near-regular SC classes).

**Remarks.** (i) The proof is pure score arithmetic: the two "defective" entries of the
complement vector c (at the base path's endpoints) are exactly what obstructs odd n —
the same endpoint asymmetry that makes c(1) = n−2 and c(n) = n. (ii) Via THM-646(iii),
the even-graph shadow offers a second mechanism view: at odd n the flip translates by
K* = K_n (the full even graph), at even n by K_n minus the consecutive matching — the
deficiency is what leaves room for class-preservation. (iii) Combined with THM-643/644/
647 and klein-S163's rigidity, the strict blue/black structure now has: all parities,
all masses, the self-loop law, Anti-Rédei, and line-rigidity — the mod-2 program closed.

**Files.** `04-computation/gn_blue_selfloop_parity_macmini_S48.py` (+ `.out`).
