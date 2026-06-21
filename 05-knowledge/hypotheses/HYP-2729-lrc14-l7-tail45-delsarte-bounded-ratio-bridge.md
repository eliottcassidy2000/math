# HYP-2729: LRC14 L7 should be a tail45-to-Delsarte bounded-ratio bridge

Status: active proof target and namespace claim, codex-2026-06-21-S72.

## Claim

HYP-2728's generated `tail45 = q_5 + 5q_6` strip is not itself a row-level
L7 inequality.  It is a generated-word compatibility gate for normalized atom
moves.  The current proof route should hand this gate to the relation-code /
Delsarte layer (HYP-2726), and only then to the KPS L7 bounded-ratio window
`f_2/f_1 in (1, 2.15)`.

The proposed bridge is:

1. **Generated-word gate.**  Factorial boundary packets give
   `q_0 = sum_j (-1)^j W_j`, while HYP-2728 excludes the cheap abstract
   q0-hiding directions by the finite generated `tail45` strip.
2. **Packet classification.**  Surviving packets are classified by the
   relation-code / Delsarte observables: low-support relation spectrum,
   Krawtchouk-positive dual packets, and miss-count law `p_t`.
3. **L7 analytic residue.**  The only non-finite input should be the KPS
   bounded-ratio window: a finite resonance atlas for small rational
   `f_2/f_1 in (1,2.15)` plus a non-resonant two-dimensional
   Erdos-Turan-Koksma constant.

## Challenged Assumption

Tournament vertices should not be runners or arcs in this bridge.  The useful
vertices are proof obligations or ratio/resonance channels.  The quotient
preserves cap margin, maximum row `p_0`, and Delsarte packet type, but it
destroys runner ownership and exact sector ancestry; those must be restored
inside the finite atlas.

Planned tournament observable: for two ratio channels, channel `A` beats
channel `B` when its worst sampled cap margin is thinner, with ties broken by
larger raw `p_0`, then smaller denominator, then lexicographic label.  This
turns the L7 scout into a risk tournament on proof obligations, not a
tournament on runners.

## Exact Evidence

The scout now exists:

`04-computation/lrc14_l7_tail45_bounded_ratio_codex_s72.py`

Output:

`05-knowledge/results/lrc14_l7_tail45_bounded_ratio_codex_s72.out`

It scanned `3111` exact primitive bounded-ratio rows after skipping `237`
nonprimitive rows, with `0` cap violations.  The thinnest margins by `k=8..12`
were:

- `k=8`: `1471/5880`, row `(0,1,2,3,4,5,20,28)`.
- `k=9`: `62911/300300`, row `(0,2,4,6,8,10,12,25,28)`.
- `k=10`: `16/91`, row `(0,1,2,3,4,5,6,7,21,35)`.
- `k=11`: `424/1911`, row `(0,1,2,3,4,5,6,7,8,21,35)`.
- `k=12`: `4891/17640`, row `(0,1,2,3,4,5,6,7,8,9,35,42)`.

The KPS named ratio `28/25` with even base has

`p0=299/1050`, `cap_9-p0=62911/300300`, `tail45=607/14700`,
`U4=4793/14700`.

The ratio-channel tournament is transitive under the risk key:

`5/3 > 7/4 > 7/6 > 5/4 > 6/5 > 2/1 > 7/5 > 4/3 > 9/5 > 8/7 > 3/2 > 15/7 > 28/25 > 8/5`.

Raw row `tail45` disagrees with margin risk on `27/91` channel edges.  This is
evidence for the sequential handoff, not for a direct raw `tail45` L7 theorem.

## Formalization Target

`THM-562` records the proof skeleton:

`01-canon/theorems/THM-562-lrc14-l7-bounded-ratio-handoff-skeleton.md`

The finite formal layer is the exact cell law: breakpoints `a/(7e)`, constant
sector words between adjacent breakpoints, and rational interval sums for
`p_t(E)`.  The two missing mathematical premises are a finite rational
resonance atlas and a non-resonant two-dimensional Erdos-Turan-Koksma bound
after the HYP-2726 Delsarte readout.

## Remaining Missing Pieces

Any positive evidence remains only evidence until the finite resonance atlas
and non-resonant discrepancy lemma are written.  The next concrete target is a
single-channel periodic proof, likely ratio `5/3`, because it is the sampled
risk leader.
