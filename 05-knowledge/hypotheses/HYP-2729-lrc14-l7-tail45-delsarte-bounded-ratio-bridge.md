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

## Missing

This stub reserves HYP-2729 for the next exact scout:

`04-computation/lrc14_l7_tail45_bounded_ratio_codex_s72.py`

The scout should compute exact missed-count laws for bounded-ratio rows,
include the KPS ratio `28/25` witness, report cap margins and row-level
`U4/tail45` diagnostics, and output a tournament fingerprint on ratio
channels.  Any positive evidence must remain labeled as evidence until the
finite resonance atlas and non-resonant discrepancy lemma are written.
