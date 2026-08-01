---
id: THM-2995
title: "Projected k2 band 1580--1599: prefix-tail and translated-capacity closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The exact
  projected k=2 scalar atlas on 1580<=z_1<=1599 has 26 rows. Twenty-four
  cheap rows have 11,978 finite packets, all closed by located order-two
  torsion. A large finite row has 29,372 packets and contracts losslessly
  to 5,507 prefix--denominator translated-capacity tests, with minimum
  slack six. The exceptional row's nonpositive lane has 955 translated-
  capacity tests with minimum slack one; its positive lane has 250,647
  zero-high packets and 997 one-high projective classes, all closed. The
  band is empty and, with THM-2980, the projected cap is z_1<=1579. This is
  not LRC(14) and removes no profile from the 165-row ledger.
source: codex-lrc14-k2-band-1580-1599-closure-2026-07-31
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2980-projected-k2-1600-1679-zero-high-and-high-order-phase-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - MISTAKE-334
  - MISTAKE-338
script:
  - 04-computation/lrc14_j7_k2_scalar_band_1580_1599_thm2995.py
  - 04-computation/lrc14_j7_k2_1580_1599_cheap_finite_torsion_thm2995.py
  - 04-computation/lrc14_j7_k2_1580_1599_translated_capacity_thm2995.py
  - 04-computation/lrc14_j7_k2_z1586_exceptional_positive_prefix_projective_thm2995.py
  - 04-computation/lrc14_j7_k2_1580_1599_integrated_promotion_audit_thm2995.py
output:
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1580_1599_thm2995.out
  - 05-knowledge/results/lrc14_j7_k2_1580_1599_cheap_finite_torsion_thm2995.out
  - 05-knowledge/results/lrc14_j7_k2_1580_1599_translated_capacity_thm2995.out
  - 05-knowledge/results/lrc14_j7_k2_z1586_exceptional_positive_prefix_projective_thm2995.out
  - 05-knowledge/results/lrc14_j7_k2_1580_1599_integrated_promotion_audit_thm2995.out
script_sha256:
  - 520957b782560b1c86aa7baadf4928034a21cdd49bed0ff33b956d534263c7ab
  - d9de1cd6f5c8d1f79adb9ee9e14f5d2d5bc0f96f762f81e1411b6397c6280967
  - c5ac80964eec09fe5b75f9d2cf62954d19258136fed1ff233b128030f0305427
  - 16bea2aa2acdda1f2bd0dde6417e300a5a5a0febec6381304bce7fd6a34db50c
  - 0c9804f0d4fad9525ecfbd9c0ccac6956fc57e38711cf3653c02432523cc0b19
output_sha256:
  - 01577b566a51d37c96c2fed783f133e27403ad3566c3df002aa5f129883f2ba4
  - caec7c86b4a0556dcee7a683d383cc0addecc574791b9056f2d91bf579afb6e6
  - 12e96d131dfc5bd9cc1f602aebc90dd56be3cd0dbc55c88adad82410f9128c8d
  - d037df44cec3d18302dfd9b2be044d51307c3eba7d408f52082356dd1e940f15
  - 6ffc5196ee33b283a5e190486932a2973b803ba8a6befb61b0a8cbc04bbe6c1a
hash_basis: LF-normalized bytes
---

# THM-2995 -- projected k2 band 1580--1599: prefix-tail and translated-capacity closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement and exact scalar inheritance

The THM-2941 direct scalar engine leaves exactly `26` projected `k=2`
survivor rows on

```text
1580<=z_1<=1599.                                          (1)
```

Their height distribution is

```text
1581:3, 1586:11, 1588:1, 1590:3, 1594:2, 1595:2, 1599:4. (2)
```

Every admissible packet on every row is incompatible with the complete-body
safe carrier. Consequently `(1)` is empty. Combined with the proved THM-2980
closure of `1600<=z_1<=1679` and the higher band above it, the projected cap is

```text
z_1<=1579.                                                (3)
```

This is a projected scalar-atlas descent only. It is not LRC(14), changes no
entry of the 165-profile ledger, and asserts no uniform theorem below 1580.

The scalar companion checks all `60,060` candidate body-height rows. Its
profile and survivor digests are

```text
profile   0f02d797f9a3bc39871d4a4ccf72c42525ca3ca4109773469b8fe95b780780d3
survivor  aab25b18f0bb9e2cb33eaa378176e720e1375c19891be9c97c67371f72ae2636. (4)
```

All scalar comparisons are exact rational comparisons.

## 2. Twenty-four cheap finite rows

On a positive amplitude ray the contribution `A/z` decreases strictly with
height. THM-2980's global top-three cutoff therefore gives a finite exhaustive
candidate list whenever its residual `eta` is positive. On 24 rows this list
has at most 378 entries. Exact four-suffix enumeration gives

```text
admissible packets  11,978
failures                  0
order-two witnesses 11,978.                              (5)
```

For each packet the first puncture already works: puncture `z_1`, fold the
complete carrier safe for the other four labels modulo the punctured
denominator, and select the recorded order-two pair. Its two phases have
circular separation `1/2>1/7`, so no common translation puts both inside one
strict-open danger arc.

## 3. The large finite row and the prefix--tail contraction

The only positive-cutoff row too large for the cheap bank is

```text
z_1=1586, E=(1,8,10,12,13,14), L=152880,
eta=2179445/164752242291>0.                               (6)
```

It has 12,588 scalar candidates and 29,372 admissible ordered suffixes
`i<j<k<ell`. Puncture the tail `ell`. The fixed-safe carrier depends only on
the leading triple `(i,j,k)`, while the translated danger predicate for the
tail depends only on its denominator `d`; height and primitive unit may be
forgotten only after that predicate is universal. The exact quotient has

```text
leading triples                              249
leading-triple x tail-denominator tests    5,507.        (7)
```

For the carrier projection `S` modulo `d`, THM-2984 gives the sharp arbitrary-
translation danger capacity

```text
kappa(d)=ceil(d/7).                                      (8)
```

Every test satisfies `|S|>kappa(d)`. The minimum slack is six, attained at

```text
leading labels (1612,2004,2340), d=7,
|S|=7, kappa(d)=1.                                      (9)
```

Thus every fibre of the prefix--denominator quotient contracts. This is a
genuine holotopy contraction of the packet incidence set: the preserved
predicate is universal translated escape; the discarded height and unit
coordinates carry no remaining obstruction. At slack zero the contraction
would be unlawful and THM-2984's affine interval-orbit sidecar would be needed.

## 4. Exceptional row: nonpositive rays

The unique nonpositive-cutoff row is

```text
z_1=1586, E=(1,10,11,12,13,14), L=840840.               (10)
```

The sign audit of MISTAKE-338 is load-bearing. Exact best-two, mixed, and
two-high gaps show that a packet containing a zero or negative suffix has
exactly three positive companions, all below the high floor. There are 35,697
positive lows but only five admissible low triples. Folding their fixed
carriers over the 191 nonpositive denominators gives

```text
5*191=955 tests, failures 0, minimum slack 1.            (11)
```

The minimum occurs at the low triple `(1612,1800,2340)` and `d=2`, where
`|S|=2` and `kappa(d)=1`. The unprojected cell-count shortcut proves only
889/955 cases; residue projection is essential. For residue zero and later
height representatives, THM-2980's full-turn argument gives projected safe
mass at least `3/7>25/91`. This closes the entire nonpositive tail without
reversing its monotonicity.

## 5. Exceptional row: positive projective completion

The all-positive part of `(10)` splits into zero-high and exactly-one-high
families. Two exact inequalities rule out two high labels, so the split is
exhaustive.

For a zero-high packet, order the four lows and puncture the fourth. The same
prefix--denominator quotient gives

```text
zero-high packets                 250,647
leading triples                       710
prefix--denominator tests          19,587
coarse cell-count passes           19,576
exact residue folds                    11
failures                                0
minimum exact-fold slack               10.              (12)
```

The minimum exact witness is prefix `(1612,1650,1800)`, denominator 12,
cardinality 12 and `kappa(12)=2`.

For exactly one high label, retain its denominator and primitive unit. The
projective census is

```text
one-high classes          997
small-torsion classes     982
full-projection classes    15
relevant unit instances   355.                           (13)
```

The 15 residual carriers have full cyclic projection; their denominators are
`11,13,143`, five classes each. For every relevant unit, choose the
unit-dependent two-cell chord of THM-2980. The minimum full-projection margin
is 59 and the minimum selected-pair margin is one. Strict-open danger makes
the equality endpoint safe. Hence every positive projective ray closes.

## 6. Exact verification and independent hostile audit

Run all four companions both normally and with `-O`, directing each output to
a scratch path, and compare with the corresponding stored transcript:

```text
python 04-computation/lrc14_j7_k2_scalar_band_1580_1599_thm2995.py
python 04-computation/lrc14_j7_k2_1580_1599_cheap_finite_torsion_thm2995.py
python 04-computation/lrc14_j7_k2_1580_1599_translated_capacity_thm2995.py
python 04-computation/lrc14_j7_k2_z1586_exceptional_positive_prefix_projective_thm2995.py
```

The independent audit checked the exhaustion directions rather than accepting
the totals: the top-three cutoff for finite rows, the two-high and mixed-shape
inequalities on `(10)`, the reversed monotonicity of negative rays, and the
height-free primitive-unit typing of projective classes. It rederived the
arbitrary-translation bound `(8)` and checked that equality would not justify
the prefix contraction. The exact `d=28` hostile from THM-2980/MISTAKE-334
still defeats the smaller centered threshold, while the order-seven/order-
eight strict-open controls separate the valid and invalid torsion boundaries.

Normal, optimized, and stored transcripts agree byte-for-byte. The four
script and four output LF-normalized hashes are exactly those in the
frontmatter. The scalar output has 7,777 LF bytes; the remaining outputs have
2,637, 1,046, and 1,339 LF bytes. No floating comparison enters a proof gate.
As an additional publication audit, reversing only the four canonical path
literals makes every script byte-identical to its independently replayed
normal/optimized source, and reversing the one printed scalar-output path makes
each stored transcript byte-identical to that replay. Thus the namespace move
changes no executable or mathematical branch.

### Independent monolithic set-join cross-audit

A fifth verifier reruns the scalar atlas monolithically over all
`3,003*20=60,060` body-height rows and recovers exactly the same `26`
survivors and height distribution `(2)`.  It then joins only certified
terminal records from the five component referees.  The comparison is

```text
monolithic atlas rows       26
pinned atlas rows           26
component terminal rows     26
duplicates                   0
missing                      0
extra                        0.
```

The terminal partition is nine `EMPTY`, thirteen `SCALAR-EMPTY`, and four
`SECTION-CLOSURE-BELOW`; those four section rows match the four translated
complete-section witnesses one-for-one.  This independently reconstructs the
atlas and exact atlas-to-component set join, while deliberately not claiming
an independent implementation of every component argument.  Ordinary and
optimized runs are byte-identical to the stored transcript; the fifth
source/output hashes are recorded in the frontmatter.

**QED.**
