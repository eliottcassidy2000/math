---
id: HYP-3144
title: "LRC14 Worpitzky Pair-Function Three-Edge Quotient"
status: RESERVED
date: 2026-06-27
owner: codex-2026-06-27-S274
tangent: T1209
technique: LTI-270
tournament_technique: LTT-168
scripts:
  - 04-computation/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.py
results:
  - 05-knowledge/results/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.out
anchors:
  - HYP-3142
  - HYP-3143
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3137
  - HYP-3134
  - HYP-3124
  - THM-084
  - T1209
  - LTI-270
  - LTT-168
---

# HYP-3144: LRC14 Worpitzky Pair-Function Three-Edge Quotient

Status: RESERVED / executable scout pending; not a proof.

This lane claims the user's three-edge prompt as a finite guardrail for the
LRC14 generating-function packet.  The model object is the quotient of
labelled three-vertex tournaments by score class:

- `T = (0,1,2)`, the transitive class;
- `C = (1,1,1)`, the cyclic class.

Flipping one directed edge gives a two-node, three-edge quotient kernel.  From
`C`, every edge flip exits to `T`.  From `T`, two edge flips stay in the
transitive class and one edge flip enters `C`.  The expected quotient matrix is

```text
        to T  to C
from T    2     1
from C    3     0
```

The same kernel appears for three coin flips after quotienting by "all same"
versus "two-to-one mix": the two-to-one class has two self flips and one exit
to all-same; all-same has three exits.  The planned scout will verify this
kernel exactly, compute its stationary distribution/eigenmodes, and compare it
with the `6:2` labelled tournament split.

Post-rebase placement: HYP-3143/S276 now owns the n=4 packet-subbasis
exact-order audit.  This K3 lane is deliberately smaller: it studies the
two-state edge-flip kernel before the n=4 lower-order leakage issue appears,
then hands its order-sensitive sidecar fields upward to the HYP-3143 packet
basis test.

## Pair-Function Guardrail

The prompt's functions

```text
a+b, a*b, a^b, b^a
```

split into two quotient types.  The sum and product survive the unordered pair
quotient `{a,b}`.  The two exponentials require an ordered pair sidecar:
`a^b` and `b^a` generally differ and are swapped by reversing the edge.

LRC14 reading: scalar PGF coefficients, score class, and Bravais-flat residue
summaries behave like `a+b` or `a*b` only when the target predicate is constant
on quotient fibers.  Conditional fiber-PGF data from HYP-3140, edge tip/tail
data from HYP-3141, and boundary/sector leakage data from HYP-3139 behave like
`a^b` and `b^a`: their order is part of the proof object.  This hypothesis
will test the small K3 kernel as a named `ordered_pair_exponent_sidecar`
warning before scalarizing LRC generating functions.

## Worpitzky Transfer

The old forward-edge polynomial thread (THM-084 and the `F(T,x)` variable)
already says Worpitzky coefficients are a graded refinement of OCF: ascents
along Hamiltonian paths retain more information than the single value `H`.
For this lane, Worpitzky is not decorative vocabulary.  The finite question is:

```text
When does an ascent/edge-flip PGF remember only a score-class scalar, and when
does it need the ordered edge/function sidecar?
```

The K3 quotient should become the smallest exact diagnostic.  Its two modes
should separate a symmetric count mode from an antisymmetric order-sensitive
mode, echoing the Worpitzky/Eulerian ascent basis and warning against reducing
the HYP-3140 fiber-PGF curve to a single value too early.

## Planned Signals

- `three_edge_flip_kernel`
- `coin_flip_kernel_isomorphism`
- `pair_function_order_word`
- `ordered_pair_exponent_sidecar`
- `score_class_scalar_survival`
- `worpitzky_ascent_payload`
- `quotient_transition_multiplicity`
- `edge_flip_role`
- `fiber_pgf_order_loss_alarm`
- `tip_tail_commutator_shadow`

## Challenged Assumption

Do not assume the tournament vertices are runners or raw arcs.  This scout will
use quotient classes, pair-function carriers, edge-flip roles, Worpitzky ascent
payloads, and LRC proof obligations as possible vertices.  The quotient
preserves the class-size/stationary-count predicate, but destroys which edge is
the unique `T -> C` exit and whether a pair function is evaluated as `a^b` or
`b^a`.

## Next Pull

Implement the scout, store its output, then upgrade this file from RESERVED to
EVIDENCE with:

1. exact K3 labelled/quotient flip counts;
2. exact coin-flip quotient isomorphism;
3. pair-function survival table;
4. Worpitzky `n=3` identity check;
5. Tournament Analysis over proof carriers rather than runners.
