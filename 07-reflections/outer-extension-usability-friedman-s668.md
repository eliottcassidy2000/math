# Outer Extension Usability and Friedman Finite Trees (S668)

Prompt: outer extension usability theorem, outer extension embedding theorem,
tree of 3, recursion theory, infinite sequence of finite threes, n colors,
homeomorphic embeddings, isomorphisms, micro incompleteness, Harvey Friedman.

The move that clicked is this:

```text
embedded maximality tells us what a maximum is;
outer-extension usability tells us whether the quotient survives the next
ambient growth operation.
```

I did not find the exact phrase "outer extension" in the Friedman PDFs opened
this session, so I am treating it as our repo phrase.  It means: grow the
ambient while preserving the internal bonds.

## Friedman Anchors

The official Friedman pages point to exactly the right neighborhood:

- reverse mathematics and recursion theory;
- reverse mathematics of homeomorphic embeddings;
- internal finite tree embeddings;
- finite trees requiring large cardinals;
- invariant maximality and concrete mathematical incompleteness;
- the 2025 invariant maximality chapters with usability and reversals.

The finite-tree paper is especially resonant.  It studies insertion domains and
rules for adding new higher vertices into finite trees.  That is a literal
outer-extension operation.  The invariant maximality chapters add the second
ingredient: a maximal object is not just maximal internally; it must be stable
under the invariance/usability requirement.

## The Toy

The S668 script builds colored rooted trees and one-step outer extensions by
adding a colored leaf.  The target predicate is whether the extended tree
contains the rooted homeomorphic color chain:

```text
0 -> 1 -> 2
```

This is intentionally small: an infinite sequence of finite threes is too big
for a session artifact, but the first finite three already shows the failure
mode.

With `3` colors and trees through `5` nodes:

- `1788` trees;
- `3204` one-step extension rows;
- `295` extensions contain the `0->1->2` pattern.

Coarse quotients leak:

- `size_height`: `6` mixed fibers;
- `color_hist`: `10` mixed fibers;
- `frontier`: `46` mixed fibers;
- `outer_address`: `55` mixed fibers.

Embedding/downset profiles repair:

- `small_embedding_profile`: `0` mixed fibers;
- `full_embedding_profile`: `0` mixed fibers;
- `address_plus_small_embedding`: `0` mixed fibers.

The surprising bit is that the insertion address alone is not enough.  You need
an embedding profile: which small probes are already represented in the
extended object.

## LRC Translation

This maps cleanly onto the current LRC14 route.

| Tree toy | LRC14 |
|---|---|
| internal tree | visible `Res_27` shell |
| outer leaf insertion | `+27` carry or owner-route lift |
| tree-of-3 predicate | floor-vs-strict obstruction |
| embedding downset | proof-obligation/owner/carry profile |
| usable quotient | no mixed floor/strict fibers |

S666 already found the first small version.  Visible `Res_27` data leaks.  Cheap
pairs leak.  Owner cover counts almost work but still leak in two named cases.
The private-owner deletion bit gives zero mixed fibers.  In this new language,
that bit is a tiny embedding profile of the proof-obligation tree.

So the next theorem should not merely say "owner bit works locally."  It should
say:

```text
Every allowed outer carry extension either changes the bounded
proof-obligation embedding profile, pays positive tax, or is a globally
coherent scalar floor lift.
```

## Micro Incompleteness

The micro-incompleteness warning is not that our toy is independent of ZFC.  It
is much smaller and more operational:

```text
finite windows can look settled while the uniform extension bound is the real
theorem.
```

That is the proof-design lesson from Friedman-style finite combinatorics.  A
finite bad sequence is cheap.  A uniform theorem preventing all bad sequences is
where the hidden proof strength lives.

For LRC, the corresponding danger is clear.  We can separate all local carry
rows through a finite radius, but the real proof requires a uniform statement
over coherent carry spaces.  The right bounded experiments are still valuable,
but only if each one discovers a reusable embedding address.

## Next Move

Build an LRC14 outer-extension embedding profile:

- probes: small D/U/N owner fragments, endpoint protector fragments, pair-pinch
  fragments, carry-cocycle fragments;
- extension operation: coherent `+27` carry and HYP-2165 owner-route lifts;
- target predicate: floor versus strict;
- success criterion: zero mixed fibers after adding the smallest embedding
  profile, except scalar coherent floor lifts.

That would be the real successor to S666 and S667.
