# LRC14 A000568 Edge Sandwich

Anchors: HYP-3133, HYP-3129, HYP-3128, HYP-3125, HYP-3124, T1200,
LTI-261, LTT-159, OPEN-Q-108.

The useful part of the user's observation is the shift.  The numbers `12` and
`56` are A000568 values one vertex above the edge-witness census layer:

```text
edge m=4: sector 10 < A000568(5)=12 < paired child 16
edge m=5: sector 20 < A000568(6)=56 < paired child 80
edge m=6: sector 35 < A000568(7)=456 < paired child 632
```

This is exactly the kind of hidden middle quotient we have been missing.  The
sector word is only equinumerosity: it counts how many outside vertices land in
the four edge sectors.  The paired child deck is equidecomposability: delete
the tail, delete the tip, keep both cards.  A000568 is neither one.  It is the
unrooted one-extra-vertex extension shadow, the equidistribution layer between
local sector sizes and two-ended recursive decomposition.

For LRC14 this does not replace the HYP-3129 SPEC certificate.  HYP-3129 is
stronger and more proof-facing: exact low Fourier modes plus a Parseval tail
already give a certified positive floor on tested multi-far rows.  The new
middle quotient should instead stratify the finite constant chase.  If a row is
bad at the sector level, ask whether its A000568 extension shadow already
selects the correct resonance class.  If not, escalate to the paired tail/tip
child deck.  If that still fails, name the resonance-lattice debt or route to
H7/Asano/Lee-Yang/phi4/Cech sidecars.

The guardrail is HYP-3128.  Naive Asano across the crowded `R` tail fails, so
this A000568 layer is not a zero-free shortcut.  It is a controlled-forgetting
diagnostic: sector word, then free extension shadow, then paired endpoint
children.  The proof-facing field to add is:

```text
a000568_extension_shadow
```

placed between `edge_tail_tip_sector_word` and the two deletion-child Rprime
ratios in HYP-3125's edge-floor packet.

Post-rebase, HYP-3131 independently made the finite-window warning explicit:
the A000568 sandwich/tameness pattern holds in the apex-7 window and breaks at
`n=8`.  That is the right scope.  The field is meant for finite LRC14
few-apex row stratification, not as an infinite tournament-enumeration law.
