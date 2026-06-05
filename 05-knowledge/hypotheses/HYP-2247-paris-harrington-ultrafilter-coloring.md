---
id: HYP-2247
status: OPEN recursion/filter-rank hypothesis with exact pair-coloring miniature
source: codex-2026-06-05-S672
related:
  - HYP-2246
  - HYP-2245
  - HYP-2244
  - HYP-2243
  - HYP-2241
  - HYP-2232
  - HYP-2189
---

# HYP-2247: Paris-Harrington Ultrafilter Coloring Recursion

## Claim

Paris-Harrington turns the HYP-2245 ultrafilter/metagraph picture into a
recursion-rank problem.

The useful finite translation is:

```text
coloring of tuple atoms
  = side-choice on a Boolean atom cube
bad coloring
  = side choice with no relatively-large homogeneous filter trace
outer extension
  = add a new ambient point and extend the side choice coherently
least PH witness stage
  = first height where every coherent bad branch dies
```

So the metagraph is still not literally the ultrafilter.  It is the base over
which side choices must descend.  Paris-Harrington adds the missing derivative
coordinate:

```text
not only upper/lower membership, but extension rank:
how many coherent bad children remain?
```

This is the repo-facing form of the user's prompt.  The black/blue upper/lower
lines resemble a principal ultrafilter on the divisor lattice of `210`; the
Paris-Harrington largeness condition asks whether such side choices can be
continued forever after the initial segment is named.  The theorem says no in
the standard finite world, while the witness function is too fast for PA to
control uniformly.

External anchors used for the theorem shape: Paris and Harrington's 1977 paper
[`A Mathematical Incompleteness in Peano Arithmetic`](https://philpapers.org/rec/PARAMI);
MathWorld's summary of the domination statement for the PH witness function
[`Paris-Harrington Theorem`](https://mathworld.wolfram.com/Paris-HarringtonTheorem.html);
and the proof-length/witness-function discussion in
[`Proof lengths for instances of the Paris-Harrington principle`](https://www.sciencedirect.com/science/article/pii/S0168007217300052).

## Exact Pair-Coloring Miniature

S672 computes the small diagonal instance:

```text
x = 2
color pairs of [N] into 2 colors
seek homogeneous H of size 3 with |H| > min(H)
```

The exact bad-coloring counts are:

| N | pair atoms | relatively-large triples | bad colorings |
|---:|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 1 | 0 | 2 |
| 3 | 3 | 1 | 6 |
| 4 | 6 | 4 | 18 |
| 5 | 10 | 10 | 12 |
| 6 | 15 | 19 | 0 |

Thus the tiny least forced stage is `N=6`.

The outer-extension tree is more informative than the count:

```text
N=1: {2: 1} bad children
N=2: {3: 2}
N=3: {3: 6}
N=4: {0: 6, 1: 12}
N=5: {0: 12}
```

At `N=4`, only the middle edge-count shell extends:

```text
edge count 2: 3 bad nodes, all dead
edge count 3: 12 bad nodes, all extend once
edge count 4: 3 bad nodes, all dead
```

At `N=5`, the surviving middle shell `edge_count=5` has no child.

This is the exact toy pattern behind the ultrafilter metaphor:

```text
upper/lower side choice gives a blue/black split,
but the PH pressure lives in the derivative child-count coordinate.
```

For this tiny case, even the scalar `edge_count` separates extendability at
`N=4`; that is a small-size accident, not the general method.  The general
method is to retain extension rank/profile as a side channel.

## Ternary Pressure Scout

S672 also samples:

```text
x = 3
color triples into 3 colors
seek homogeneous H of size 4 with |H| > min(H)
```

Random colorings avoid all targets frequently through `N=9`; in `2000` random
trials the avoid counts are:

```text
N=4..9: 1921, 1679, 1131, 499, 130, 25
N=10..12: 0, 0, 0
```

But heuristic local repair still finds explicit bad colorings for every
`N=10..20` tried.  Therefore the random pressure is not evidence of a forced
stage.  It is evidence of the correct qualitative obstruction:

```text
the branch can keep pushing coherence into the tail long after random
colorings stop doing so.
```

That tail-escape feature is the finite shadow of the fast-growing
Paris-Harrington witness function.

## Relation To HYP-2245 And HYP-2246

HYP-2245:

```text
literal ultrafilters live on Boolean tiling/address cubes;
tournament isomorphism quotients leak side membership.
```

HYP-2246:

```text
endpoint tournament enumeration repairs quotient leaks with
Phi(T)=multiset_v(iso(T-v), L/M/U side(outdeg(v))).
```

HYP-2247:

```text
coloring recursion repairs another leak by adding extension rank:
Phi_bad(C)=profile of coherent bad children under outer extension.
```

In short:

```text
upper/lower filter side
  -> L/M/U owner address
  -> extension-rank derivative profile
```

The third layer is where recursion theory enters.

## Tournament Analysis

Vertices are transfer routes:

1. `PH_initial_segment_filter`
2. `bad_coloring_outer_extension_tree`
3. `endpoint_half_filter_trace`
4. `metagraph_ultrafilter_sheaf`
5. `LRC_owner_carry_filter_rank`
6. `unit_distance_spine_bulk_filter`
7. `ordinary_Ramsey_shadow`
8. `raw_coloring_count`
9. `generic_ultrafilter_CH_analogy`

Pairwise observable:

```text
exact_small_evidence,
retains_address,
outer_extension_rank,
quotient_descent_fit,
LRC_transfer,
overclaim_risk_control
```

Switch: majority.  Tie Hamiltonian path: the listed priority order.

Fingerprints:

- `score_histogram={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

Top route:

```text
PH_initial_segment_filter
  -> bad_coloring_outer_extension_tree
  -> endpoint_half_filter_trace
  -> LRC_owner_carry_filter_rank
```

## LRC14 Transfer

For LRC14, the proposed PH-style bad tree is:

```text
node:
  owner/carry fiber over the C=27 shell with no lonely witness
outer extension:
  coherent +27 carry/owner lift
relative-largeness gate:
  early named residue/owner section that cannot be deferred to the tail
rank:
  number/profile of coherent bad child fibers after the HYP-2241
  owner-private deletion filter is attached
```

The concrete next proof target is:

```text
After adding the owner-private deletion bit from HYP-2241,
the bad-child rank strictly decreases on every local carry extension
that preserves the visible Res_27 floor shadow.
```

This would turn "owner/carry labels repair the quotient" into a termination
certificate rather than just a no-leak classification.

## Unit-Distance Transfer

For unit-distance construction search, use the same shape:

```text
initial segment:
  unit Hamiltonian spine / named unit-distance owner section
tail:
  bulk edges and late slab/cap coordinates
bad node:
  construction shadow that postpones every required unit-spine witness
rank:
  spine-owner deletion profile plus remaining extension count
```

This refines the earlier spine/bulk split: edge count alone is the raw scalar;
the proof-relevant coordinate is whether bad constructions still have coherent
children after the spine-owner section is named.

## Challenged Assumption

Do not assume tournament-analysis vertices must be runners, arcs, or coloring
classes.  In this session the useful vertices are:

```text
tuple atoms,
bad-coloring nodes,
outer-extension children,
initial-segment gates,
owner/carry filters,
unit-spine owners,
deletion side states,
proof obligations.
```

Preserved predicate:

```text
existence/nonexistence of a relatively-large homogeneous trace,
or its LRC/unit-distance analogue.
```

Destroyed information:

```text
raw labels, exact colors, and late-tail coordinates once the quotient is taken.
```

The repair is an extension-rank side channel.

## Next Tests

1. Define the LRC14 bad-child rank on the HYP-2241 AP/Vstar/2AP carry fibers
   and test whether owner-private deletion makes rank strictly decrease.
2. Extract a canonical `Phi_bad` trace for tournament endpoint enumeration:
   `half-filter trace + child-count profile`.
3. Build a unit-distance spine-owner rank over the `n=21` Moser carrier:
   compare unit-spine-preserving extensions against bulk-only extensions.
4. Replace the ternary heuristic scout with a SAT/CP exact lower-bound search
   for small `x=3` PH bad colorings.
