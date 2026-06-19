---
id: HYP-2633
title: LRC(14) two-large reciprocal coupling - the character kernel must be paired with residue-lift equidistribution
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S24
depends_on:
  - HYP-2632
  - HYP-2630
  - HYP-2626
  - HYP-2614
  - THM-538
related:
  - HYP-2619
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2633 - LRC(14) Two-Large Reciprocal Coupling

## Claim

HYP-2632 gives the finite `d=9` repeated-residue character kernel, but the
next proof obligation is not finished by that finite table alone.  The signed
kernel must be coupled to the actual reciprocal hyperplane sum

```text
sum_{n_i != 0, 7 not| n_i, sum e_i n_i=0}
    C_9(n mod 7) / (n_1...n_6).
```

Scratch exact shell sums on representative two-large supports show the key
hazard: finite coimage-kernel lanes and actual reciprocal-lift lanes can have
different signs at bounded height, because denominator signs and relation
resonances weight residue classes non-uniformly.

The next theorem should therefore be:

```text
finite chi_7/affine kernel cancellation
+ residue-lift equidistribution / summation by parts
=> two-large reciprocal tail bound.
```

This is a sharper version of the HYP-2614 cotangent/Dedekind target.  The
finite packet table supplies the signs; the missing analytic step is to prove
that the relation lattice samples the lanes evenly enough after the finite
low-height ledger is removed.

## Computation

Script:

- `04-computation/lrc14_two_large_reciprocal_coupling_codex_s24.py`
- output: `05-knowledge/results/lrc14_two_large_reciprocal_coupling_codex_s24.out`

The script compares HYP-2632 packet weights with exact cumulative reciprocal
relation-lattice sums through height `H=16`.  The finite packet is computed in
the HYP-2632 unit

```text
U = 147/(16*pi^6).
```

The reciprocal lift is enumerated over exact integer supports:

```text
42_QR_a2      (1,2,8,9,15,22)   finite packet (1,1,1,1,2,2)  -> -25U
42_QR_a4      (1,4,8,11,15,22)  finite packet (1,1,1,1,4,4)  -> -25U
42_NQR_a3     (1,3,8,10,15,22)  finite packet (1,1,1,1,3,3)  -> -18U
411_high_23   (1,2,3,8,15,22)   finite packet (1,1,1,1,2,3)  ->  +8U
411_low_26    (1,2,6,8,15,22)   finite packet (1,1,1,1,2,6)  ->  +1U
411_zero_02   (1,2,7,8,15,22)   finite packet (1,1,1,1,0,2)  ->   0U
```

At `H=16` the actual cumulative reciprocal lifts are:

```text
case             finite S/U   lift signed       abs/signed   guardrail
42_QR_a2             -25      +0.002676143676       20.128   finite/lift sign flip
42_QR_a4             -25      -0.000130513735      432.469   same sign, huge envelope
42_NQR_a3            -18      +0.002424138026       24.703   finite/lift sign flip
411_high_23           +8      +0.000178063082      357.576   same sign, huge envelope
411_low_26            +1      -0.001230813632       50.455   finite/lift sign flip
411_zero_02            0      +0.000411593461      155.020   finite-zero lift leak
```

The strongest warning is that `42_QR_a2` and `42_QR_a4` have the same finite
HYP-2632 weight `-25U`, but their bounded reciprocal lifts have opposite signs
at `H=16`.  The affine-zero finite packet also has a nonzero bounded lift.
Thus the finite Dedekind table is the correct coefficient layer, but it is not
yet the analytic tail bound.

## Refined Proof Target

The theorem-shaped target is now:

```text
finite chi_7/affine/Q Dedekind packet
+ finite low-height wall deletion
+ residue-lift equidistribution on relation lattices
+ signed summation by parts inside additive frequency shells
=> two-large reciprocal tail bound using -108U+54U, not 162U.
```

A useful analytic form is to write each additive-frequency shell as a Stieltjes
sum over relation-lattice height, prove bounded cumulative residue-lift
discrepancy after wall deletion, and then apply Abel summation to the
`1/prod(n_i)` denominator.

## HYP-2634 Follow-Up

HYP-2634 explains the first opposite-sign pair from this file.  In the seed
family `S_a=(1,a,8,a+7,15,22)`, the finite QR packets `a=2` and `a=4` both
have weight `-25U`, but only `a=4` has extra negative low-height relation
motifs with defects

```text
2a - 8 = 0,
7a - 28 = 0.
```

The HYP-2633 lift lemma should therefore split into two pieces:

```text
finite low-height relation-defect zero sieve
+ residual residue-lift equidistribution / Abel summation.
```

Equidistribution should not be asked to absorb exact low-height resonances
that can be isolated in a finite wall ledger.

## Tournament Analysis

Candidate vertices included runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, additive Fourier modes,
exact-period packets, matroid circuits, and proof obligations.

Chosen Hamiltonian path:

```text
residue_lift_equidistribution
> additive_frequency_summation
> finite_chi_affine_kernel
> low_height_wall_deletion
> exact_period_squarefree_packets
> boundary_face_cancellation
> raw_support_values
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1,1]
hamiltonian_paths = 1
```

The quotient preserves the signed two-large reciprocal-tail predicate after
low-height wall deletion.  It destroys exact witness times and individual
low-height relation identities.

## Status

Partially confirmed as a guardrail and next-lemma locator.  LRC(14) remains
open.  HYP-2633 narrows the remaining proof obligation: prove the lift
equidistribution / summation-by-parts lemma that lets the HYP-2632 finite
`chi_7`/affine/Q cancellation survive the actual reciprocal hyperplane sum.
