---
id: THM-3087
title: "z234 two-high common-torsion and finite-height reduction"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-lrc-z234-boundary-scout-2026-08-01
audit: >
  An independent hostile audit first reconstructed the common-Q pointwise
  capacity certificate and its exact 1,045/654 mask split.  A second immutable
  pass rederived the complete-cell discrepancy inequality, both rational
  height cutoffs, the three exceptional loads and strict margins, and the
  108,966,498-packet finite census; it independently replayed normal and
  optimized modes against stored output and verified the script, output,
  bank, dependency, and semantic hashes.  The audit confirms that this is a
  finite reduction only and makes no projected closure, cap, or ledger claim.
depends_on:
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1094-exact-two-comb-component-theorem
related:
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-3033-projected-k3-z246-to-z244-descent-and-z243-high-floor-addendum
script: 04-computation/lrc14_j7_k3_z234_two_high_common_torsion_finite_height_thm3087.py
output: 05-knowledge/results/lrc14_j7_k3_z234_two_high_common_torsion_finite_height_thm3087.out
script_sha256: 440077a9fda3f9f5cb8898659d0e5006b6e2529b6bbb8c12ced37d7a67e5c69a
output_sha256: cd3059f6a83b2b5cd63176ea359ec76bd3f46fe75e144cef95752f6b462f5bff
semantic_sha256: e7f27276343813afa2b34471e8d3b98a28f32491b2be49204e6ca5acc3bc144d
hash_basis: LF-normalized bytes
---

# THM-3087 -- z234 two-high common-torsion and finite-height reduction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3078 leaves four projected `k=3`, `z_1=234` rows at the two-high
boundary.  The present theorem isolates the genuinely scalar-unbounded part
of that boundary and makes it finite.  It also closes a majority of its
denominator masks by an exact common-torsion capacity test.

This is a **finite reduction**, not a row closure.  It makes no ledger or
projected-cap decrement.

## 1. Exact scalar-unbounded mask bank

Let `L` be the body period, `H=floor(13L/132)+1` the high floor, and let
`delta(z)` denote the exact suffix surplus.  In a two-high packet, one of the
three suffix slots is occupied by a low literal label and two are occupied by
high labels.  The scalar threshold for the two high labels is nonpositive
exactly in the following two THM-3078 bodies:

```text
E=(1,5,9,11,12,14):
  L=194040, H=19111, low=243,
  low excess = 34729/119189070,
  886 denominator masks;

E=(1,9,10,11,12,14):
  L=194040, H=19111, low=260,
  low excess = 3307/4414410,
  813 denominator masks.                                  (1)
```

Thus the exact scalar-unbounded sector has

```text
886+813=1699 masks.                                       (2)
```

The tracked mask bank is

```text
05-knowledge/results/
  lrc14_j7_k3_z234_two_high_mask_bank_thm3087.tsv
```

with LF hash

```text
b3d424991e816bd5ddc0540554483cd3bc40efbaf44a3f3064197c337ae6ca8d. (3)
```

It is the sorted text projection of the pinned THM-3078 residual screen.
The companion can rebuild it from the THM-3078 checkpoint fingerprint using
`--extract-bank`; no untracked scalar search horizon enters (1).

## 2. Common-torsion two-band certificate

Fix one mask with high denominators `d_1,d_2|L`, put

```text
Q=lcm(d_1,d_2),   kappa_i=ceil(d_i/7),                  (4)
```

and let `J` be the set of complete `1/L` body cells which are safe for both
the fixed label `234` and the mask's lone low label.  Write

```text
S={j mod Q:j in J}.                                      (5)
```

At a fixed normalized within-cell coordinate `x`, a high label of denominator
`d_i` occupies an open unit-oriented cyclic band of at most `kappa_i`
residues modulo `d_i`.  Its lift to `Z/QZ` therefore occupies at most

```text
kappa_i Q/d_i                                           (6)
```

residues.  We deliberately grant the two bands independent translations and
orientations.  Hence the strict inequality

```text
|S| > kappa_1 Q/d_1 + kappa_2 Q/d_2                    (7)
```

ensures that for every `x` some complete body cell is safe from **both** high
labels.  Thus the exact four-drift projection is the full circle,

```text
P_(E,Z)=T.                                               (8)
```

The three aligned tails have open-union mass at most `36/91` by THM-1166,
so (8) contradicts the cover criterion THM-2941 `(25g)`.

The exact census is

```text
body                         masks   (7) closes   remains
(1,5,9,11,12,14)              886       663        223
(1,9,10,11,12,14)             813       382        431
---------------------------------------------------------
total                         1699      1045        654.   (9)
```

The positive slacks are `1575..17434` and `700..13482` respectively.
The worst failures are `-80678` and `-80124`.  Each body has exactly one
equality mask, the high-denominator pair `(2,2)`.  All `1045` successes use
the exact support (5); the coarse cardinality lower bound closes none.

## 3. A uniform physical height contract

The `654` masks left by (9) still contain infinitely many literal high-label
pairs at the scalar level.  We now make that sector finite without confusing
raw safe mass with the projected cover criterion.

Choose any complete cell `I` from `J`, so `mu(I)=1/L`.  On the normalized
cell, the three aligned tails occupy at most `36/91`.  THM-1094's exact
one-comb interval discrepancy gives, for every integer high label `z`,

```text
L mu(I intersect D_z) <= 1/7 + 6L/(49z).                (10)
```

If both high labels are at least `N`, a cover of `I` would require

```text
1 <= 36/91 + 2/7 + 12L/(49N).                          (11)
```

The right side is strictly below one once

```text
N > 1092L/1421 = 4324320/29.                            (12)
```

Consequently the smaller high label is at most

```text
149114.                                                  (13)
```

For each of the two bodies, the companion then checks the exact weak endpoint
criterion for every integer

```text
19111 <= z <= 149114.                                   (14)
```

Except at

```text
z=97020, 129360, 145530,                                (15)
```

there is a complete body/`234`/low cell wholly safe from `D_z`.  On that cell
only the other high comb and the three aligned combs remain.  By (10), the
other high label must satisfy

```text
w <= 13L/49 = 51480.                                    (16)
```

The three exceptional heights are not loopholes.  Their least normalized
loads on an exact clean cell, their far-label cutoffs, and the strict margin
obtained already at `w=z` are

```text
z       least load  cell   far cutoff       margin at w=z
97020      1/7      14850  2162160/29       47/637
129360     3/28     14850  2882880/43       435/2548
145530     2/21     14851  324324/5         388/1911.    (17)
```

Every cutoff in (17) is strictly smaller than its exceptional `z`.  Since
`z` was chosen as the smaller high label, (17) is impossible.  Combining
(13)--(17) yields the exact finite contract

```text
19111 <= z,w <= 51480.                                  (18)
```

This argument is physically typed on an actual complete interval.  It does
not replace the two high combs by scalar masses and it retains the exact
high-ray law

```text
(z+L) delta(z+L)=z delta(z).                            (19)
```

## 4. Exact finite remainder

The companion imposes (18), the mask denominators, distinctness when the two
denominators agree, and the exact scalar inequality using (19).  After the
`1045` common-torsion closures, the remainder is

```text
body                         masks   ordered literal pairs
(1,5,9,11,12,14)              223       36,200,669
(1,9,10,11,12,14)             431       72,765,829
---------------------------------------------------------
total                          654      108,966,498.      (20)
```

Fourteen masks have no labels at all in the finite box (nine in the first
body, five in the second).  In every nonempty mask all raw denominator-
compatible ordered pairs in (18) satisfy the scalar inequality; the number
in (20) is therefore both the raw and exact scalar count.

## 5. Positive and hostile controls

For each of the four literal scalar maximizers displayed by THM-3078, the
common-torsion test is already strict.  It proves `P_(E,Z)=T`, with slacks

```text
502, 9504, 13482, 8566,                                  (21)
```

and independently reproduces the four positive raw safe masses recorded in
THM-3078.  This closes those **four literal assignments only**.  It does not
close their rows, because another assignment in the same residual body can
have a different projected section.

The first failed implication is therefore

```text
finite height or positive raw drift-safe mass
  ==> mu(P_(E,Z)) > 36/91 for every remaining packet.   (22)
```

What is missing is the exact placement of the three-aligned tail union
`U_A` inside the four-drift projection `P_(E,Z)`.  The `108,966,498` packets
in (20) are a finite input to that problem, not a closure certificate.

## 6. Exact evidence and scope

The companion performs explicit checks rather than optimization-sensitive
Python assertions.  It verifies the mask bank, all `1699` common-modulus
supports, all `260008` whole-cell height tests, the exception loads and
margins, the full high-ray recurrence used by (20), and the four literal
controls.  Normal and `python -O` output agree byte-for-byte with the stored
transcript.

**Scope.**  This theorem closes `1045` scalar-unbounded denominator-mask
families and reduces the other `654` to (20).  It does **not** prove
`mu(P_(E,Z))>36/91` for those packets, close any of the four THM-3078 rows,
change the projected cap, or decrement the LRC ledger.

QED.
