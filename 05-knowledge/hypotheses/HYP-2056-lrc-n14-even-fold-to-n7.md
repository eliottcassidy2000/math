---
id: HYP-2056
status: PARTIAL — rigorous reduction for e<=6 + computational (127/127); open residual
source: opus-2026-06-01-S554
related:
  - HYP-2055
  - HYP-2052
  - THM-369
---

# HYP-2056: LRC@14 reduces to the proven n=7 case via the even-fold (the "7 lever")

## Setup

`n=14=2*7`, LRC(7) proven. Even speed `v=2u` ⇒ `||v t||=||u·(2t)||`. With
`fold(S)={v/2: v∈S even}`, `s=2t`: `g_even(t)=g_fold(2t)`, so
`M14(S) ≤ M(fold(S))` (verified exactly, 0/224 violations).

## Claims

- **(rigorous) `min(#even,#odd) ≤ 6`** (13 is odd) ⇒ LRC(7) always engages on the
  minority parity.
- **(rigorous) e ≤ 6:** `M(fold) ≥ 1/7` (LRC(7)), so even runners are `≥1/7`-safe on
  a positive-measure `E_good={t: g_fold(2t)≥1/7}`. The two doubling-preimages
  `s/2,(s+1)/2` put each ODD runner at antipodal points (differ by 1/2) ⇒ each odd
  is unsafe (`<1/14`) at ≤1 preimage. Hence **LRC(14)[e≤6] ⟸ some even-good `s` has
  no odd-split.**
- **(computational) the reduction succeeds 127/127** on e≤6 configs: an even-good
  preimage that is a full n=14 witness always exists.
- **(structure)** both known tight configs (AP, V*) have e=6 and the same 7 odd
  runners `{1,3,5,7,9,11,13}`. AP's fold `{1..6}` IS the n=7 AP (n=7-tight); V*'s
  fold `{1,2,3,4,5,12}` is n=7-loose (2/13). So tightness at 14 does NOT need the
  fold n=7-tight — the odd coupling carries it.

## Open / handoff

1. Prove "no odd-split over `E_good`" for e≤6 ⇒ would make LRC(14)[e≤6] a theorem
   conditional only on the proven LRC(7).
2. The e≥7 (o≤6) branch: dodging ≥7 even runners needs LRC(8+) — find an odd-side
   fold or a covering argument.

## Status

Not a proof of LRC@14. A rigorous reduction of ~half the configs to LRC(7) plus a
crisp odd-residual, verified computationally, isolating the difficulty into the
odd antipodal split.

**See:** `07-reflections/lrc-n14-even-fold-to-the-proven-n7-case-s554.md`,
`04-computation/lrc_n14_even_fold_to_n7_s554.py` (+.out); proven LRC(7); S520
(n=7 menu collapse); HYP-2055 (V*), HYP-2052.
