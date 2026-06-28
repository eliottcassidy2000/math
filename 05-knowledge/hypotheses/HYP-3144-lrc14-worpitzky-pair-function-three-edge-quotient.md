---
id: HYP-3144
title: "LRC14 Worpitzky Pair-Function Three-Edge Quotient"
status: EVIDENCE
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

Status: EVIDENCE / exact K3 quotient scout; not a proof.

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

## Scout Readout

Script:
`04-computation/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.py`.
Stored output:
`05-knowledge/results/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.out`.

Exact labelled K3 enumeration gives class sizes `T=6`, `C=2`.  Single-edge
flips give raw labelled counts:

```text
        to T  to C
from T   12     6
from C    6     0
```

Dividing by source-class size gives the quotient kernel:

```text
        to T  to C
from T    2     1
from C    3     0
```

The normalized Markov kernel has stationary distribution `{T:3/4,C:1/4}`,
matching the labelled `6:2` class split, and eigenvalues `1` and `-1/3`
(`3` and `-1` before normalization).  Edge roles split exactly as:

```text
adjacent_order_edge_returns_to_T: 12
cycle_edge_breaks_to_T: 6
long_source_sink_edge_exits_to_C: 6
```

Thus the transitive score class has one order-sensitive long-edge exit and two
self-class adjacent-edge flips.  The score-class quotient preserves the
transition multiplicities but forgets which edge is the unique exit.

The three-coin quotient is isomorphic after identifying `T <-> mix` and
`C <-> same`:

```text
        to mix  to same
from mix    2        1
from same   3        0
```

## Worpitzky / PGF Warning

The forward-edge PGF audit is the sharpest evidence:

```text
class T aggregate_F=(1,4,1)
  distinct F: (1,0,0)^1, (0,1,0)^4, (0,0,1)^1

class C aggregate_F=(1,4,1)
  distinct F: (1,2,0)^1, (0,2,1)^1
```

So both score classes have the same aggregate Worpitzky/Eulerian-looking
payload `(1,4,1)`, while their state-level PGF curves differ.  This is the
smallest executable witness for the user's "single value is never the whole
PGF curve" rule.  The `n=3` Worpitzky identity

```text
x^3 = C(x,3) + 4 C(x+1,3) + C(x+2,3)
```

checks for `x=0..7`, but the scout warns not to identify that ascent basis
with raw score-class counts.  It is the coordinate system for order-sensitive
function payloads.

## Pair-Function Readout

Over samples `(2,3),(2,4),(4,2),(3,5),(5,3)`, `a+b` and `a*b` are invariant
under swap in `5/5` cases.  `a^b` and `b^a` are swap-invariant only in `2/5`
sample cases, exactly the accidental equality cases around `(2,4)` and
`(4,2)`.  The proof packet therefore needs:

```text
unordered_pair_scalar_allowed: a+b, a*b
ordered_pair_sidecar_required: a^b, b^a
```

LRC14 transfer: HYP-3140's `F_R(y)` / `F_R,Q(y)`, HYP-3141's tip/tail edge
witness, HYP-3139's reflection-block leakage, and HYP-3143's packet order
audit should all carry `pair_function_order_word` before scalarization.

## Tournament Analysis

Tournament Analysis used proof carriers rather than runners:

```text
pairwise_observable =
  weighted retained proof payload:
  exactness, order, curve, LRC transfer, quotient guard,
  HYP-3143 compatibility, minimality
binary_gauge =
  A->B iff weighted payload score(A)>score(B);
  ties use lexical Hamiltonian path
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
selected_path =
  ordered_pair_exponent_sidecar
  -> three_edge_flip_kernel
  -> worpitzky_ascent_payload
  -> packet_subbasis_order_audit
  -> fiber_pgf_conditional_curve
  -> tip_tail_commutator_shadow
  -> symmetric_sum_product_quotient
  -> score_class_scalar
  -> raw_single_value
```

## New Hypotheses

1. The normalized `-1/3` antisymmetric K3 mode is a reusable order-loss
   contraction signature for one-coordinate quotienting in HYP-3140 fiber
   PGFs.
2. HYP-3143's n=4 lower-order leakage is the K3 unique long-edge exit unfolded
   by one more free bit.
3. LRC14 scalar extremality should be phrased as a legal function-evaluation
   theorem: symmetric functions may quotient; ordered exponent-like functions
   require a sidecar or named debt.
4. Worpitzky ascent payloads are the right coordinate for detecting when a PGF
   equality at `x=1` hides an order-sensitive curve difference.

## Challenged Assumption

Do not assume the tournament vertices are runners or raw arcs.  This scout will
use quotient classes, pair-function carriers, edge-flip roles, Worpitzky ascent
payloads, and LRC proof obligations as possible vertices.  The quotient
preserves the class-size/stationary-count predicate, but destroys which edge is
the unique `T -> C` exit and whether a pair function is evaluated as `a^b` or
`b^a`.

## Next Pull

Attach `pair_function_order_word`, `ordered_pair_exponent_sidecar`,
`three_edge_flip_kernel`, `worpitzky_ascent_payload`, `edge_flip_role`,
`fiber_pgf_order_loss_alarm`, and `tip_tail_commutator_shadow` to the active
HYP-3140/HYP-3141/HYP-3139/HYP-3143 packet rows.  Then test whether the
`-1/3` antisymmetric mode or the class-level `aggregate_F=(1,4,1)` collision
predicts which HYP-3140 fiber-PGF coefficient inequalities need a sidecar.
