---
id: HYP-3144
title: "LRC14 Worpitzky Pair-Function Three-Edge Quotient"
status: EVIDENCE / exact K3 quotient scout plus HYP-3147/S71/HYP-3150/3151/3152/HYP-3200 integration; not a proof
date: 2026-06-27
owner: codex-2026-06-27-S274
tangent: T1209
technique: LTI-270
tournament_technique: LTT-168
scripts:
  - 04-computation/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.py
results:
  - 05-knowledge/results/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.out
reflections:
  - 07-reflections/worpitzky-pair-function-three-edge-quotient-codex-20260627.md
anchors:
  - HYP-3142
  - HYP-3200
  - HYP-3161
  - HYP-3199
  - HYP-3160
  - HYP-3153
  - HYP-3152
  - HYP-3151
  - HYP-3150
  - HYP-3147
  - HYP-3145
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

Status: EVIDENCE / exact K3 quotient scout plus HYP-3147/S71/HYP-3150/3151/3152 integration; not a proof.

This lane turns the user's three-edge prompt into a finite guardrail for the
LRC14 generating-function packet.  The model object is the quotient of
labelled three-vertex tournaments by score class:

- `T = (0,1,2)`, the transitive class;
- `C = (1,1,1)`, the cyclic class.

Flipping one directed edge gives a two-node, three-edge quotient kernel.  From
`C`, every edge flip exits to `T`.  From `T`, two edge flips stay in the
transitive class and one edge flip enters `C`.  The class-level multiplicity
matrix is

```text
        to T  to C
from T    2     1
from C    3     0
```

The same kernel appears for three coin flips after quotienting by "all same"
versus "two-to-one mix": the two-to-one class has two self flips and one exit
to all-same; all-same has three exits.  The scout verifies this kernel exactly,
computes its stationary distribution/eigenmodes, and compares it with the
`6:2` labelled tournament split.

Post-rebase placement: HYP-3143/S276 now owns the n=4 packet-subbasis
exact-order audit.  This K3 lane is deliberately smaller: it studies the
two-state edge-flip kernel before the n=4 lower-order leakage issue appears,
then hands its order-sensitive sidecar fields upward to the HYP-3143 packet
basis test.

Later mainline work sharpened the lane rather than replacing it.  HYP-3147
recasts the same two-state kernel with a cyclic coin reference and names the
minority-edge gate inside the transitive fiber.  S71 then links the same
order-sensitive face to the k=8 correction: score-to-iso compression is
bijective through n=4 and fails at n=5, exactly where the bounded-core cap dip
turns on; the odd Worpitzky `-9S3` part dominates the even biquadratic `+6S4`
part.  HYP-3150 packages that as a completed factor-through audit, HYP-3151
executes the same audit as a function-compression scout and proves the n=4
`x=a OR c, y=b OR c` compression has no affine `GF(2)^3 -> GF(2)^2`
replacement, and the KPS Lee-Yang lambda scout adds the root-curve warning
that scalar PGF values are not enough.  HYP-3152 corrects the group-theoretic
language around the two-bit picture: the useful slice can look Klein-four, but
the full flip action is a transformation monoid with an absorbing apex arc.
The solvability claim should therefore be stated as bounded degree `<=4` /
Galois `<=S4`, not as a V4 law.  HYP-3153 reserves the executable packet that
should fuse these lanes, and HYP-3160 sharpens the k=8 node: the even
biquadratic face is a degree-two variance/total-covariance extremality, while
the odd `-9S3` Worpitzky face is the non-associative sidecar.  The follow-up
S31ai correction keeps the covariance target but drops two tempting shortcuts:
the associativity-defect ratio is not a universal `1/7`, and the odd part is
not explained by a removable anchor residue.
HYP-3200 makes that correction exact on the bounded bank
`E={0} union A`, `A subset {1,...,14}`, `|A|=7`: consec has
`Sigma kappa_3/S3=407891843/2855269200`, differs from `1/7` by
`-3757/2855269200`, and no row has exact ratio `1/7`.  The surviving
proof-facing result is that consec maximizes `Sigma kappa_2` among primitive
rows (`0/3431`).
HYP-3199 supplies the prompt-exact n=4 companion: the fixed-path `a,b,c` cube
is an abundant cover whose `S` fiber has five representatives, while the
partial-score `(0,1,1,2)` chart with live `x=a,y=b` is the exact Einheit
section.  The class-preserving map `x=a OR c, y=b OR c` is small but lossy, so
representation abundance must carry minimality/deletion sidecars before it can
be used as proof evidence.

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
records the small K3 kernel as a named `ordered_pair_exponent_sidecar`
warning before scalarizing LRC generating functions or HYP-3145/HYP-3149
filler-core quotients.  HYP-3199 adds the adjacent n=4 rule: if a quotient uses
the fixed-path cover, it must also record whether `c` is filler/deletable and
whether the exact `x,y` Einheit section is available.

## Worpitzky Transfer

The old forward-edge polynomial thread (THM-084 and the `F(T,x)` variable)
already says Worpitzky coefficients are a graded refinement of OCF: ascents
along Hamiltonian paths retain more information than the single value `H`.
For this lane, Worpitzky is not decorative vocabulary.  The finite question is:

```text
When does an ascent/edge-flip PGF remember only a score-class scalar, and when
does it need the ordered edge/function sidecar?
```

The K3 quotient is the smallest exact diagnostic.  Its two modes separate a
symmetric count mode from an antisymmetric order-sensitive mode, echoing the
Worpitzky/Eulerian ascent basis and warning against reducing the HYP-3140
fiber-PGF curve to a single value too early.

This also sharpens the bounded-core meta-point from
HYP-3132/HYP-3142/HYP-3150/HYP-3152/HYP-3160/HYP-3200.  The hard LRC14 core stays below the
Abel-Ruffini wall only if the proof packet keeps enough order to choose the
right branch.  The k=8 even face is the solvable biquadratic/covariance side;
the odd Worpitzky/minority-edge face is the orientation sidecar that prevents a
false scalar compression.  The algebraic ceiling is the degree-four/S4 ceiling,
and the proof target is now more concrete: prove the even total-covariance /
excess co-emptiness extremality, then bound the odd Worpitzky/non-associative
residual with its sidecars intact and without assuming a now-refuted universal `1/7` or
anchor law.

## Signals

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

Do not assume the tournament vertices are runners or raw arcs.  This scout uses
quotient classes, pair-function carriers, edge-flip roles, Worpitzky ascent
payloads, and LRC proof obligations as possible vertices.  The quotient
preserves the class-size/stationary-count predicate, but destroys which edge is
the unique `T -> C` exit and whether a pair function is evaluated as `a^b` or
`b^a`.

## Next Pull

Attach `pair_function_order_word`, `ordered_pair_exponent_sidecar`,
`three_edge_flip_kernel`, `worpitzky_ascent_payload`, `edge_flip_role`,
`fiber_pgf_order_loss_alarm`, and `tip_tail_commutator_shadow` to the active
HYP-3140/HYP-3141/HYP-3139/HYP-3143/HYP-3145/HYP-3149 packet rows.  Then test
whether the `-1/3` antisymmetric mode, the class-level `aggregate_F=(1,4,1)`
collision, the S71 odd `-9S3` dominance, or the KPS off-circle root variance
predicts which HYP-3140 fiber-PGF coefficient inequalities need a sidecar.
