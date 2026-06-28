---
id: HYP-3161
title: The score→iso compression boundary (bijective n≤4, fails n=5) IS the k=8 LRC hard row — the cap dip turns on exactly at |P|=13-k crossing 4→5; and the dip's gK8 content -9S₃+6S₄ splits by PARITY into the EVEN biquadratic face (a+b/a·b, S70) and the ODD Worpitzky face (a^b/b^a, codex HYP-3147, dominant)
status: VERIFIED (score→iso n≤4 bijective/n=5 fails; dip onset at |P|=4→5; parity split, odd dominates) + SYNTHESIS (the k=8 dip = even biquadratic + odd Worpitzky; the compression-trick). Not a proof.
source: mac-mini-2026-06-27-S71
merges:
  - HYP-3147   # codex n=3 edge-flip Worpitzky kernel = the ODD/antisymmetric face
  - HYP-3132   # mac-mini k=8 biquadratic resolvent = the EVEN/symmetric face
related:
  - HYP-3150   # function-compression guardrail
  - HYP-3151   # executable Worpitzky/function-compression scout
  - HYP-3152   # Lee-Yang/Galois correction
  - HYP-3153   # Lee-Yang/Worpitzky/quartic packet follow-on
  - HYP-3199   # n=4 Einheit/minimality chart
  - HYP-3122   # φ⁴ cumulants (the even/odd cumulants = the parity faces)
  - THM-577    # cap = the symmetric face, exact for |P|<=3 (compression works)
  - THM-062    # deformed Eulerian (Worpitzky)
  - OPEN-Q-108
external: Worpitzky's identity; Eulerian numbers A(3,k)=[1,4,1]
---

# HYP-3161 — The score-compression boundary is the k=8 hard row, split by parity

## Invariants as arc-cube functions; the function quartet by parity
An invariant is a function on `(Z/2)^{C(n,2)}`. The owner's four operations split an edge's endpoint pair:
`a+b, a·b` SYMMETRIC (the score/H, even, order-blind — "cannot see an arc flip", codex); `a^b, b^a` ASYMMETRIC
(the orientation `i→j`, odd, order-aware). This is the even⊕odd / cut⊕cycle split.

## The compression (VERIFIED, `lrc_arc_cube_compression_parity_macmini_S71.py`)
The SCORE (commutative `a+b` face) **determines the iso class for n≤4** (`n=3,4` bijective) and **FAILS at n=5**
(12 iso vs 9 scores). So at n≤4 the iso class has a low-degree coordinate chart: the owner's **scheme 2 is a
clean Klein-four slice** (the 2-arc slice closes into a group), while **scheme 1 (3 arcs) is the same data
over-coordinatized** (a magma). S72/HYP-3152 corrects the broader claim: the full flip action is a transformation
monoid with an absorbing apex arc, not a V4 group. The trick: **compress to the gauge where the symmetric face is a
complete invariant ⟹ finite-degree computability.**

HYP-3199 sharpens this at n=4: the fixed-path `a,b,c` table is a high-multiplicity cover with a deletable
`c`/`S`-bulk coordinate, while the partial-score `x,y` chart is the exact Einheit section.  The score boundary is
therefore not just about class counts; it is about when the class chart has a minimal section plus named deletion
sidecar.

## The payoff: k=8 IS the compression boundary
The cap dip turns on exactly at the n=4→5 boundary (`|P|=13−k`):
```
 k:    13  12  11  10    9        8
 |P|:   0   1   2   3    4        5
 dip:   0   0   0   0  1/4004  1081/76440
```
`k=8 ⟺ |P|=5` — the quintic level where `score→iso` fails. **The cap dip is the failure of the commutative
compression**: `|P|≤3` (k≥10) the symmetric pair-Pascal face is complete (`dip=0`, THM-577); `|P|=4` (k=9) it
begins to leak; `|P|=5` (k=8) the antisymmetric (orientation) content is irreducible — the binding row.

## The parity merge (VERIFIED for consec_8)
The gK8 dual's higher-order term `−9S₃ + 6S₄` splits by parity:
- **even `+6S₄`** = the symmetric face = the **biquadratic resolvent** `u⁴−5u²+4` (S70/HYP-3132), solvable;
- **odd `−9S₃`** = the antisymmetric face = the **Worpitzky/ordered** content (codex HYP-3147), **dominant**
  (`|odd|/|even| ≈ 3.15`), a sum of n=3 edge-flip oscillations (eigenvalue **−1/3**) with Eulerian descent
  weights `A(3,k)=[1,4,1]` (`x³=C(x+2,3)+4C(x+1,3)+C(x,3)`).

**Improved argument:** bound the k=8 dip = [even biquadratic coefficient, solvable degree-2, S70] + [odd
Worpitzky-weighted sum of n=3 `−1/3` oscillations, codex]. The owner's `a+b,a·b` vs `a^b,b^a` is exactly this
parity decomposition; the dominant (odd) face is the content the score-compression cannot see.

## Next
1. Make the parity decomposition of the k=8 dip rigorous: dip = even(biquadratic) + odd(Worpitzky).
2. Bound the odd part as a Worpitzky-Eulerian-weighted sum of n=3 `−1/3` edge-kernel oscillations (codex's
   HYP-3147 Next Test #3) + the even part by the biquadratic (S70). Together ⟹ close the k=8 dip.
