---
id: HYP-2626
title: LRC(14) prime-mask/coimage transfer seam
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S19
depends_on:
  - HYP-2624
  - HYP-2617
  - HYP-2619
  - THM-539
related:
  - HYP-2627
  - HYP-2625
  - HYP-2608
  - THM-538
  - THM-540
  - HYP-2623
  - OPEN-Q-108
---

# HYP-2626 - LRC(14) Prime-Mask/Coimage Transfer Seam

## Claim

The hidden recurrence suggested by the LRC14 support-six residual is not a
one-dimensional recurrence in the number of runners.  It is a finite transfer
through two quotients:

```text
prime-mask transfer of low-height wall supports
-> LRC14 unit seam (Z/14Z)^*
-> projective mod-7 coimage class
-> signed repeated-residue tail.
```

The exact seam identity is elementary:

```text
(Z/14Z)^* = {1,3,5,9,11,13} -> F_7^* = {1,2,3,4,5,6}.
```

Thus the HYP-2617 quotient by scalar multiplication in `F_7^*` is exactly the
coimage of the 14-runner unit action.  The mod-7 coimage atlas is not an
arbitrary residue trick; it is the quotient forced by the `14=2*7` clock.

## Computation

Script:

- `04-computation/lrc14_prime_mask_coimage_transfer_codex_s19.py`
- output: `05-knowledge/results/lrc14_prime_mask_coimage_transfer_codex_s19.out`

The script reruns the height `<=2` one-large wall enumeration from HYP-2624
and records the prime mask of the large wall speed `M` for primes `{2,3,5,7}`.
It then asks how much nonzero coimage signed mass is addressable when only
apex masks from a chosen finite mask set are allowed.

Support-level census:

```text
k   supports  classes  M divisible by 7  support touches 7
8     272007      147              38219            168865
9     177648      147              24507            115641
10     68124      108              10293             32479
```

The prime-mask transfer table is:

```text
k=8:  mask {} already hits 46/46 nonzero classes, 100% signed mass.
k=9:  mask {} already hits 79/79 nonzero classes, 100% signed mass.
k=10: mask {} hits 73 classes, 72.120496% signed mass.
k=10: mask {2,3,5} still hits only 73 classes, 72.120496% signed mass.
k=10: mask {7} hits 85 classes, 84.229179% signed mass.
```

So the old mod-30 slice is not the live transfer coordinate in the fixed LRC14
support-six residual.  Adding primes `2,3,5` to the apex mask does nothing at
`k=10`; adding the LRC prime `7` accounts for the extra `12` wall-addressed
classes and raises the addressed mass from `72.120496%` to HYP-2624's
`84.229179%`.

Minimal apex-mask antichains:

```text
k=8:  ({})  -> 46 classes, 15.765012 signed mass
k=9:  ({})  -> 79 classes, 12.821615 signed mass
k=10: ({})  -> 73 classes,  7.7789776 signed mass
k=10: ({7}) -> 12 classes,  1.3060527 signed mass
```

## Repeated-Tail Character Split

After prime-mask transfer and height-2 wall addressing, the same `31` k=10
tail-only classes remain.  The dominant repeated packet has a multiplicative
character split over `F_7^*`.

For `4+2` classes `(1,1,1,1,a,a)`:

```text
a=2,4   chi_7(a)=+1   |S_9|=0.23891209
a=3,5,6 chi_7(a)=-1   |S_9|=0.17201670
```

For `4+1+1` classes `(1,1,1,1,a,b)`, the signed mass collapses to a small list
of signatures in

```text
(chi(a), chi(b), chi(ab), chi((a-1)(b-1))).
```

The largest entries have `|S_9|=0.076451868`; two smaller entries have
`|S_9|=0.0095564835`.  Thus the residual is not merely "repeated residues" in
a combinatorial sense.  It is a repeated-root character-sum packet.

## Interpretation

This explains why mod `30` looked real and why it was not enough.  In the
max-min spectrum, the corrected THM-539/HYP-2623 picture still shows the mask
mechanism clearly: the family `F(k,a)` gets its confirmed `a=3,4` dips by
making the large speed `a(k-1)` kill clocks whose denominators divide `k-1`,
while the attempted `a=5` primorial continuation collapses to the floor rather
than producing an unbounded dip.  In the fixed LRC14 support-six tail, the live
seam is instead the unit action of `(Z/14Z)^*` and the prime `7` coordinate
that becomes invisible in the mod-7 coimage relation.

The slogan is:

```text
prime masks route the wall ledger;
unit seams quotient the residue address;
coimage characters carry the signed tail.
```

This is compatible with HYP-2624.  The prime-mask/coimage transfer does not
delete the final tail; it says exactly where the final theorem must live:
a signed cotangent/Dedekind estimate for repeated-root packets
`(1,1,1,1,a,b)` with explicit quadratic-character cases.

It is also complementary to the concurrent HYP-2625/THM-540 modular-recurrence
thread.  HYP-2625 supplies the small-part address rows
`{2,3} -> {2,3,5} -> {2,3,5,7}` and explains mod `210` as the squarefree
divisor-profile interface to the support-six mod-7 tail.  HYP-2626 asks which
of those addresses actually survive the height-2 wall/coimage quotient.  The
answer is asymmetric: the mod-30 address is inert for the fixed k=10 signed
coimage wall coverage, while the `{7}` address is live and contributes the
extra `12` wall-addressed classes.  So HYP-2625 gives the recurrence address
space; HYP-2626 gives the signed coimage transfer occupancy inside it.

HYP-2627 gives the same warning from the complete-graph side.  The raw
Harary-Hill `K_14` product has profile `1260 -> rad 210`, while the divided
crossing value `315` loses prime `2`.  The match is exact at the known
THM-523 measure champion:
`15/36-2/5-1/70-1/504=1/(2*1260)`, doubled by symmetry to `1/1260`.
This matches the coimage-transfer discipline here: quotient only after the
proof-relevant squarefree address has been retained.

## Tournament Analysis

The session explicitly challenges runner vertices again.  Candidate vertices
included runners, gaps, raw supports, large speeds, prime masks, mod-7 residue
tuples, unit orbits, coimage classes, and proof obligations.  The quotient that
preserves the LRC14 support-six predicate is:

```text
unit_seam_coimage
> prime_mask_transfer
> height2_wall_classes
> repeated_tail_packet
> signed_dedekind_tail
> raw_supports
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1]
```

The quotient destroys exact witness times and row identities.  It preserves
the support-six signed-tail address, which is the relevant predicate after
THM-538 and the finite wall ledgers.

## HYP-2630 Update

HYP-2630 completes the Euler-copy re-indexing test suggested by HYP-2629.  The
copy ledger gives exact-period packet capacity, but it is uniform over
`F_7^*`.  For raw `q=1260`, exact top-period packets give `48` copies per unit
residue and the full `{2,3,5,7}` mask gives `96` copies per unit residue.
Consequently copy capacity is identical inside the nonzero `4+2` row
`(1,1,1,1,a,a)`.

The observed split is instead the quadratic-character phase:

```text
chi_7(a)=+1: |S_9|=0.23891209
chi_7(a)=-1: |S_9|=0.17201670
```

So the remaining repeated-root theorem should retain the `chi_7` channel
explicitly.  Raw prime masks and raw Euler-copy capacities are necessary
addresses, but not the final separator.

## HYP-2632 Update

HYP-2632 refines the repeated-root target.  The `4+2` packet is genuinely
Legendre:

```text
2*S_9(1,1,1,1,a,a)/U = -43 - 7*chi_7(a),  a=2..6.
```

However, the `4+1+1` signature listed above is incomplete if used alone.  The
finite kernel has a hidden affine zero lane:

```text
(0,2), (3,6), (4,5)  <=>  a+b=2 mod 7  <=>  S_9=0.
```

Off that lane, the remaining high/low split is controlled by

```text
Q(a,b)=ab*(1+3(a+b))-1,
S/U=8 iff chi_7(Q)=+1,
S/U=1 iff chi_7(Q)=-1.
```

So the repeated-root cotangent/Dedekind theorem should use a `chi_7` plus
affine-line packet table, not just the four multiplicative characters.

## Status

Partially confirmed by exact enumeration of height `<=2` one-large wall
supports and exact recomputation of the HYP-2617 coimage masses.  LRC(14)
remains open.

Next steps:

```text
1. Convert the 4+2 character split into an explicit cotangent/Dedekind bound
   for (1,1,1,1,a,a), separated by chi_7(a).
2. Do the same for 4+1+1 packets using the character signatures above.
3. Check whether multi-large low-height walls add new prime-mask antichains or
   only recycle the same {7}/repeated-root packet.
4. Fold the resulting signed tail theorem into HYP-2608's wide-spread bound.
```
