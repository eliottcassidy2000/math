---
id: HYP-2636
title: LRC(14) block-frequency transfer - the two-large reciprocal tail should be bounded as a core/pair channel norm
status: OPEN
source: codex-2026-06-19-S24
depends_on:
  - HYP-2632
  - HYP-2630
  - HYP-2624
  - HYP-2617
  - HYP-2614
related:
  - HYP-2608
  - THM-538
  - OPEN-Q-108
---

# HYP-2636 - LRC(14) Block-Frequency Transfer

## Claim

HYP-2632's finite `chi_7` / affine / `Q` packet should lift to the exact
support-six reciprocal tail only after the hyperplane sum is grouped by exact
additive frequency.  For a model two-large face

```text
c1*n1+c2*n2+c3*n3+c4*n4 + A*x + B*y = 0,
```

do not estimate the six reciprocal variables independently.  First write the
tail as a finite transfer over exact channels `s`:

```text
T_{A,B}(H) = sum_s < Core_s(u,v), Pair_s^{A,B}(u,v) >_{u,v in F_7^*}.
```

Here `Core_s` is the four-variable signed reciprocal table with the two pair
residues `(u,v)` left open, and `Pair_s` is the exact two-large table over
`A*x+B*y=-s`.  This keeps both pieces of information that the earlier absolute
Minkowski envelope destroyed:

1. the exact additive relation channel `s`;
2. the full `6 x 6` mod-7 residue matrix carrying HYP-2632's signed packet.

The theorem-shaped target is therefore a channelwise Cauchy/Schur estimate:

```text
sum_s ||Core_s||_2 ||Pair_s^{A,B}||_2 <= acceptable support-six tail margin
```

after finite low-height wall deletion, rather than a raw six-dimensional
absolute harmonic bound.

## Evidence

Script:

- `04-computation/lrc14_block_frequency_transfer_codex_s24.py`
- output: `05-knowledge/results/lrc14_block_frequency_transfer_codex_s24.out`

At truncation `H=24`, ambient dimension `d=9`, the script enumerates
two four-core transfer tables, each with `3,111,696` terms, and compares three
envelopes:

- `block L2`: `sum_s ||Core_s||_2 ||Pair_s||_2`;
- `block L1`: the entrywise product of already-summed blocks;
- `raw atom abs`: the old envelope that takes absolute values before summing
  inside each channel.

Selected rows:

```text
case              packet  active s   signed        L2/sig  L1/sig  raw/sig
4+2 QR             -25U       13   +8.4796e-03       2.6     1.11     21.4
4+2 NQR            -18U       11   +8.8629e-03      2.72     1.05     21.7
4+1+1 high           8U      191   -1.4990e-03      38.7     3.34     91.2
4+1+1 low            1U       49   -5.3594e-04      54.7     3.45      205
4+1+1 zero           0       187   +8.5335e-05       459     18.9     1420
spread 4+2 QR      -25U       97   +2.0594e-04       127     16.7      302
spread zero          0      1273   -3.9505e-04       145     16.1      185
```

The zero-lane row is the clearest signal.  The fully raw absolute mass is about
`1420` times the signed mass, but the block `L1` envelope drops this to `18.9`.
The same-residue spread core `(1,8,15,22)` is less symmetric and has many more
active `s` channels, but it still keeps block `L1/signed` near `14-17` while
the raw envelope is `185-302`.  This does not prove the tail bound, but it
shows that the large absolute mass seen in HYP-2614 is at least partly a
bookkeeping artifact of taking absolute values before the exact transfer
channel and residue matrix have been summed.

## Pair-Line Form

For fixed `A,B`, the pair block is one-dimensional on each exact channel.  Let
`g=gcd(A,B)`, `a=A/g`, `b=B/g`, and `c=-s/g`.  If `g` does not divide `s`, the
channel is empty.  Otherwise solutions of

```text
a*x + b*y = c
```

are

```text
x = x0 + b*t,
y = y0 - a*t.
```

For `c != 0`,

```text
1/(x*y) = (b/c)/(x0+b*t) + (a/c)/(y0-a*t).
```

So each nonzero pair channel is a pair of arithmetic harmonic sums, with
residue restrictions modulo `7`.  The `s=0` channel is the separate line
through the origin and has quadratic decay.  This is exactly the place where a
cotangent/Dedekind estimate can keep HYP-2632's character phase visible.

## Connection To KPS S11

The incoming KPS S11 empty-sector scripts keep the whole distribution
`p_t=P(N=t)` rather than only `p_0`; the later `p_6` maximality scout sharpens
this to an exact all-missed/AP correlation endpoint.  The lesson matches this
block transfer: the useful object is not the scalar shadow but the retained side
channel.  On the moment side, the side channel is the empty-sector distribution.
On the support-six tail side, it is the exact additive channel `s` together
with the `6 x 6` residue matrix.

This suggests a common proof grammar for LRC(14):

```text
retain packet/distribution -> apply signed low-degree or channel estimate
-> scalarize only after cancellation is visible.
```

## Tournament Analysis

Candidate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, speed residues, cover arcs, additive
frequency channels, pair-line lattices, Fourier modes, residue matrices,
cotangent/Dedekind shells, and proof obligations.

Chosen Hamiltonian path:

```text
core_pair_block_transfer
> additive_frequency_packet
> affine_Q_selector
> signed_channel_cauchy_bound
> successive_minima_tail_count
> blind_pair_residue_matrix
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1]
```

The quotient preserves the exact relation channel, the pair residues, and the
finite signed packet phase.  It destroys runner labels, the ordering of the
six coordinates, and the full untruncated reciprocal lattice.  The challenged
assumption is that the residual Minkowski lemma should count six free
coordinates first; the two-large face has a natural exact core/pair transfer
before any lattice-volume envelope is introduced.

## Status

Open.  This is evidence and a proof route, not a proof of LRC(14).  The scout
currently tests normalized and same-residue spread four-core/two-large models at
finite height `H=24`.  The next step is to turn the pair-line formula into a
uniform Dedekind/cotangent bound and splice it into HYP-2608 after the finite
height-2 wall deletion.
