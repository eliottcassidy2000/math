---
id: HYP-2894
status: NEGATIVE GUARDRAIL / realizability obstruction for a naive unital map
source: codex-2026-06-22-S107
tags: [lrc14, unital, pair-slots, realizability, exact-tiling, goddyn-wong, tournament-analysis]
related:
  - HYP-2893
  - HYP-2892
  - HYP-2891
  - HYP-2890
  - HYP-2889
  - THM-560
  - OPEN-Q-108
results:
  - 04-computation/lrc14_unital_pair_slot_realizability_codex_s107.py
  - 05-knowledge/results/lrc14_unital_pair_slot_realizability_codex_s107.out
external_prompts:
  - https://mathworld.wolfram.com/Unital.html
  - https://mathworld.wolfram.com/ClebschGraph.html
  - https://mathworld.wolfram.com/TruncatedOctahedralGraph.html
---

# HYP-2894: the q=3 unital is not a canonical AP8 pair-slot design

HYP-2892 noticed the numerology

```text
q=3 unital points = q^3+1 = 28 = C(8,2),
block size = q+1 = 4,
lambda = 1.
```

That makes the Hermitian unital a real tight frame candidate for the k=8
pair-slot layer.  S107 tests the most canonical possible strengthening:

> Are the unital blocks realizable as a natural `S8`-invariant block system on
> the 28 edge slots of `K8`?

The answer is no.

## Exact check

Script:

- `04-computation/lrc14_unital_pair_slot_realizability_codex_s107.py`
- output: `05-knowledge/results/lrc14_unital_pair_slot_realizability_codex_s107.out`

The script enumerates the `S8` orbits of four-edge subsets of `K8`.  There are
11 orbit types.  For each orbit it computes:

```text
orbit size
component signature
point replication if the orbit is included
lambda on adjacent edge-pairs
lambda on disjoint edge-pairs
```

It then checks every union of orbit types.  No union has:

```text
63 blocks,
edge-point replication 9,
adjacent-edge pair lambda 1,
disjoint-edge pair lambda 1.
```

Equivalently, no `S8`-invariant union realizes a `2-(28,4,1)` design on the
edges of `K8`.

The script also records the AP8 pair packets:

```text
sum class size histogram  = {1:4, 2:4, 3:4, 4:1}
diff class size histogram = {1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1}
equal-sum pair count      = 22
equal-diff pair count     = 56
```

So the natural AP pair-slot packets are highly nonuniform, unlike the unital
lambda-1 frame.  A literal symmetric identification would erase exactly the
observer-relative structure that THM-560 and S41 say the proof must keep.

## Meaning after THM-560

Incoming THM-560 proves the structured difference-closed exact tilers are
exactly the AP-dilates `d*{1,...,13}`.  It also isolates Goddyn-Wong
`{1,...,11,13,24}` as a sporadic tight row, explained by a balanced
gap-plus-collision residue pattern rather than by a complete residue system.

Therefore HYP-2894 redirects HYP-2892:

```text
wrong target:
  identify the q=3 unital with canonical AP8/K8 pair slots

live target:
  use the q=3 unital as a weighted tight frame after choosing a labelled
  category-1 map tied to AP-dilates plus the Goddyn-Wong sporadic locus
```

This is a useful negative result.  It prevents the proof search from promoting
an elegant observer-blind design into a coverage theorem.  The unital can still
average residual pair packets, but the labelling/weighting must be supplied by
the exact-tiling geometry, not by `S8` symmetry.

## Tournament Analysis

Vertices are proof-carrier candidates:

```text
category1_exact_tiling_AP_plus_GW,
AP8_sum_difference_packets,
unital_q3_weighted_pair_frame,
scalar_additive_energy,
S8_invariant_unital_block_system.
```

Pairwise observable: priority tuple

```text
(LRC faithfulness, tight-locus specificity, incidence strength, canonicality).
```

The resulting Hamiltonian path is:

```text
category1_exact_tiling_AP_plus_GW
  > AP8_sum_difference_packets
  > unital_q3_weighted_pair_frame
  > scalar_additive_energy
  > S8_invariant_unital_block_system.
```

Assumption challenged: the most symmetric finite design on the right number of
points should be the missing LRC object.  The check says no.  The useful object
must preserve the observer-relative AP/Goddyn-Wong tight-locus data, while the
unital supplies only a secondary pair-frame once that labelling exists.

## Next proof obligation

Define an explicit labelled map from the category-1 tight-locus atlas to pair
slots:

```text
AP-dilate residue bijection packets
  plus Goddyn-Wong one-gap/one-collision packets
  -> weighted pair-frame coordinates
  -> HYP-2890 residual leak blocks.
```

Then test whether each labelled block has nonpositive residual leak after
subtracting the positive same-frequency `A*` packet.  If that works, the
remaining analytic part is still HYP-2636/HYP-2887 tail and curl cancellation.
