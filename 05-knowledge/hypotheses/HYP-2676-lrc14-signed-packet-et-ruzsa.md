---
id: HYP-2676
title: LRC(14) signed Erdos-Turan packet / Ruzsa-model scout
status: OPEN; exact computational evidence, proof obligation sharpened
source: codex-2026-06-20-S48
depends_on:
  - HYP-2675
  - HYP-2674
  - THM-546
  - HYP-2638
  - HYP-2639
related:
  - HYP-2632
  - HYP-2662
  - HYP-2667
  - HYP-2671
  - HYP-2672
  - HYP-2653d
  - OPEN-Q-108
---

# HYP-2676 - Signed Erdos-Turan Packet / Ruzsa Model

## Claim Being Tested

The same-sign one-missed-sector packet obstruction from HYP-2674 is not a
generic Erdos-Turan discrepancy failure.  It appears to live in compressed
additive models: AP, near-AP, shell-full dyadic blocks, doubled-odd packets,
and the boundary collar rows already isolated by HYP-2675.

The proposed split is:

```text
low additive growth / small model
  -> finite Ruzsa/Freiman model classification of ++++++ packet alignments

large additive growth / high model dimension
  -> signed packet cancellation or small packet mass, measured before
     absolute values and paired with THM-546's BV one-far estimate.
```

This is deliberately narrower than a scalar `w*Delta_w` constant.  THM-546
already proves the gapped BV estimate

```text
|Delta_w(E',w)| <= kappa * V(E') / (pi^2 w).
```

HYP-2676 asks what remains in the ungapped/bounded-near-plateau regime:
classify the finite same-sign packet models, then prove signed cancellation
outside them.

## Computation

Script:

- `04-computation/lrc14_signed_packet_et_ruzsa_codex_s48.py`
- output: `05-knowledge/results/lrc14_signed_packet_et_ruzsa_codex_s48.out`

The script keeps exact Fraction arithmetic for:

- `Delta_w = p0(E' union {w}) - Phi(E')`;
- the packet telescope
  `w*Delta_w = sum_runs [G0(w*b - s/7) - G0(w*a - s/7)]`;
- packet sign word by missed sector `s=1..6`;
- sector-absolute and run-absolute cancellation ratios;
- additive profiles: `|E'+E'|`, excess, `K2`, `K3`, additive energy,
  longest AP/run, and squarefree profile.

The only floating display is the THM-546 comparison constant
`kappa=1.856901...`.

## Exact Findings

Named rows confirm the packet story:

```text
B13 same-sign packet:
  E'=(0,1,2,4,6,7,8,10), w=12
  Delta_w=997/5880
  sign=++++++
  sumset excess=5, K2=5/2
  |sum|/sector_abs=1, |sum|/run_abs=997/1237

HYP-2671 dyadic block:
  E'=(0,1,2,4,8,12,16,20), w=24
  Delta_w=457/3920
  sign=++++++
  excess=9, K2=3
  |sum|/run_abs=457/607

HYP-2672 doubled-odd tail:
  E'=(0,1,2,4,8,14,26,34), w=38
  Delta_w=483281/5761028
  sign=++++++
  excess=15, K2=15/4
  still below the k=9 comparison margin by 7647117/158428270
```

The counterweight is the wide/third-pocket row:

```text
KPS third pocket:
  E'=(0,3,5,16,28,30,33), w=35
  Delta_w=1171/452760
  sign=++-+--
  excess=13, K2=26/7
  |sum|/run_abs=1171/15473
```

This is the useful contrast: it has comparable additive width to the dangerous
rows, but the signed run-level packet already cancels by an order of magnitude.

## B14 Near-Speed Bank

The B14 near-speed bank scans primitive
`E'={0}+7-subsets of [1,14]`, `w=max(E')+1..max(E')+4`.

Bucketed by sumset excess:

```text
AP:       count=4,     ++++++=1,    max Delta=445/8232
near_AP:  count=376,   ++++++=58,   max Delta=5347/30870
small:    count=12220, ++++++=1465, max Delta=155831/905520
medium:   count=1124,  ++++++=107,  max Delta=3071/26754
```

The top positive rows are all `++++++` and live in low-to-small excess
patterns:

```text
(0,2,4,6,7,8,9,10), w=12:
  Delta=5347/30870, excess=3, K2=9/4

(0,4,6,8,10,11,12,14), w=16:
  Delta=155831/905520, excess=6, K2=21/8

(0,1,4,6,8,10,12,14), w=16:
  Delta=2357/13720, excess=7, K2=11/4
```

This supports a Ruzsa-model interpretation: low-growth rows should be
normalized into finite cyclic/AP/GAP packets before any scalar discrepancy
bound is attempted.

## Dyadic/Bertrand Scale Ledger

For the dyadic family

```text
E_s={0,1,2,4,8,3s,4s,5s}
```

the `w=6s` spike remains finite at `s=4`, with
`Delta=457/3920`.  Later dyadic blocks have much smaller maxima:

```text
s=12..23:  best Delta=2539/64680
s=24..47:  best Delta=3221/94080
s=48..95:  best Delta=127/3675
s=96..160: best Delta=3607/113680
```

The analogy to Bertrand/Chebyshev is only methodological: local bias can
persist on a finite scale, but a useful proof needs a scale reset ledger.  Here
the reset is not "there is a prime in every interval"; it is that every dyadic
block after the finite packet has either smaller same-sign efficiency or enough
signed cancellation to stay far below the margin.

## Integration With Prior Work

HYP-2638 supplies the finite Freiman `3k-4` pocket: if
`|E+E| <= 3k-4`, translate/dilate to a finite normal form and verify exact
sector values there.

HYP-2639 supplies the warning: additive energy and sumset size alone are too
coarse.  A valid quotient must preserve summand shell, visibility, and sign.
HYP-2676 therefore records packet signs, not just `K2` and energy.

HYP-2632 supplies the finite signed-kernel analogy: Legendre selectors and
affine zero lanes can reduce absolute packet mass before analytic tail bounds
are applied.  The present packet signs are the one-far analogue of that lesson.

HYP-2662 and `lrc14_far_delta_galois_phase_codex_s38.py` supply a reusable
phase algebra for future refinements:

```text
G0(y - s/7)
  = trace average
  + chi_7(s) * QR/NQR quadratic channel
  + intra-quadratic residual.
```

HYP-2667 and `lrc14_p1_tax_packet_frontier_codex_s41.py` supply the existing
phase-packet contribution machinery: positive/negative packet mass, QR/NQR
counts, and the B13 `Delta^+/p1=997/2562 < 2/5` frontier.  HYP-2676 is the
global `Delta_w` counterpart of that local `p1`-tax packet ledger.

The per-sector Koksma script is also a guardrail: only exact singleton-missed
runs telescope correctly.  Any future signed ET proof should use

```text
w*Delta_w = sum_s sum_{exact-{s} runs} [G0(wb-s/7)-G0(wa-s/7)],
```

not a full "sector-s-missed" telescope.

THM-546 supplies the rigorous gapped estimate.  HYP-2676 is the complementary
ungapped finite-model target.

## Next Proof Obligation

1. Define an exact finite "same-sign packet model" class:
   AP/near-AP, shell-full dyadic block, doubled-odd packet, and HYP-2675
   boundary collar variants.
2. Use Ruzsa modeling/Freiman normalization to make the low-growth class finite.
3. Prove a signed Erdos-Turan packet estimate for the complement, with a
   high-growth or high-dimension condition forcing either a non-`++++++` sign
   word or a run-absolute cancellation ratio bounded away from 1.
4. Feed this back into HYP-2675: direct `p0` margins close wide rows, while
   packet signs handle bounded near-plateau peels.

No LRC(14) proof is claimed here.

## Tournament Analysis

Vertices are packet proof obligations / additive model profiles, not runners.
The pairwise observable is the share of the THM-546 BV budget consumed by the
exact signed `Delta_w`.  The stored tournament is transitive on the top twelve
vertices (`0` directed 3-cycles), with Hamiltonian path beginning:

```text
nonshell_warning
> B13_same_sign_packet
> B14 small-excess ++++++ rows
> HYP2671_dyadic_block
> HYP2675_true_wide_base
> HYP2672_doubled_odd_tail
```

The quotient preserves the far-peel discrepancy and additive model class; it
destroys the full endpoint/state-word address, so any theorem must reattach the
HYP-2648 state-word and HYP-2675 direct-`p0` ledgers before final assembly.
