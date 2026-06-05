---
id: HYP-2236
status: OPEN synthesis / exact small finite evidence through n=6
source: codex-2026-06-05-S660
related:
  - HYP-2237
  - HYP-2235
  - HYP-2234
  - HYP-2233
  - HYP-2228
  - HYP-2187
  - HYP-2176
  - HYP-2171
---

# HYP-2236: Tournament Deck Collisions Are Repaired by Derivative Side Channels

## Claim

Tournament deck collisions should be treated as marked-projection failures.
The ordinary vertex-deleted deck forgets the boundary/owner label of the
missing vertex.  The useful repair is not just "add more scalar data"; it is a
paired deletion derivative:

```text
card T-v  +  how v attached to that card.
```

The cleanest derivative found in S660 is:

```text
(card isomorphism type, deleted vertex outdegree).
```

This paired deleted-score deck resolves every checked tournament deck
collision through `n=6`.

## Computation

S660 imports the S658 tournament utilities, recursively generates tournament
isomorphism classes through `n=6`, and compares repair channels:

```text
full_deck
Kelly subtournament counts
scalar = (H, score sequence, c3, SCC sizes)
scalar_deck
full_deck + H
full_deck + score sequence
full_deck + full scalar
unpaired deleted-score multiset
paired (card, deleted score)
paired (card, c3 loss)
paired (card, H loss)
paired (card, deleted score, c3 loss)
paired all derivatives
S658 loss deck
card-scalar all derivatives.
```

The exact class counts are:

```text
n=3: 2
n=4: 4
n=5: 12
n=6: 56.
```

## Main Findings

The ordinary full deck still collides through `n=6`:

```text
n=3: full_deck buckets=1, max bucket=2
n=4: full_deck buckets=3, max bucket=2
n=5: full_deck buckets=11, max bucket=2
n=6: full_deck buckets=52, max bucket=2, 4 colliding buckets.
```

Kelly-style subtournament counts have the same collision strength as the full
deck in this range, so they are deck-visible but do not repair the projection.

The scalar channel is much coarser at `n=6`:

```text
scalar buckets=40, max bucket=4, 14 colliding buckets, 30 colliding classes.
```

The important negative result is that global scalar repairs do not resolve the
`n=6` deck collisions:

```text
full_deck + H:        still 4 colliding buckets
full_deck + score:    still 4 colliding buckets
full_deck + H+score:  still 4 colliding buckets
full_deck + scalar:   still 4 colliding buckets.
```

Nor does the unpaired deleted-score multiset:

```text
unpaired deleted score: still 4 colliding buckets at n=6.
```

But the paired deleted-score deck resolves every checked collision:

```text
paired (card, deleted score):
  n=3: max bucket=1
  n=4: max bucket=1
  n=5: max bucket=1
  n=6: max bucket=1.
```

The stronger derivative repairs also resolve all checked cases:

```text
paired (card, score, c3 loss): max bucket=1 through n=6
S658 loss deck:                max bucket=1 through n=6
card-scalar all derivatives:   max bucket=1 through n=6.
```

## Collision Anatomy

At `n=6`, all four full-deck collision pairs are converse pairs, have the same
`H`, score sequence, `c3`, and SCC shape, and are separated by the paired
score/`c3`-loss derivative:

```text
0x18 vs 0x1c: H=23, c3=6, scores=(1,2,2,3,3,4), converse, flip distance 1
0x92 vs 0x95: H=31, c3=6, scores=(1,2,2,3,3,4), converse, flip distance 2
0xb1 vs 0xb5: H=29, c3=6, scores=(1,2,2,3,3,4), converse, flip distance 1
0x158 vs 0x159: H=43, c3=8, scores=(2,2,2,3,3,3), converse, flip distance 1.
```

This explains why raw scalar repairs fail: the collision lives in the
card-to-boundary assignment, not in the global totals.

## Deck-Visible Checks

The expected Kelly-style `c3` formula holds:

```text
sum_v c3(T-v) = (n-3) c3(T)
```

for all checked classes through `n=6`.  The sum of card Hamiltonian counts is
deck-visible but is not a replacement for the paired boundary derivative.

## Tournament Analysis

Vertices are repair channels.  The pairwise observable is

```text
(separates_checked_collisions, deck_visibility,
 side_channel_transfer, simplicity).
```

The repair-channel tournament is transitive:

```text
paired_deleted_score
> paired_score_c3_loss
> full_deck
> kelly_subcounts
> full_deck_plus_H
> full_deck_plus_score
> paired_H_loss
> scalar.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
hamiltonian_paths=1.
```

## Cross-Repo Synthesis

This is the deck version of the repo's larger side-channel principle:

```text
ordinary deck       -> paired deleted score
LRC Res_27 quotient -> owner/carry/Cprime labels
unit-distance graph -> point-deletion degree/gain/direction support
union-closed freq   -> set pressure
finite-field Kakeya -> pinned/concurrency/owner labels
LRC wild strategy  -> no-leak deletion-derivative atlas
CH cardinal shadow  -> forcing/model side channel.
```

The same pattern repeats: scalar equality is often a cache, not a proof
object.  The useful object is a derivative or boundary label attached to the
projection fiber.

## Work Queue

- **S660a:** Extend to `n=7` using canonical augmentation or a nauty-style
  automorphism-pruned generator, not raw permutation canonicalization.
- **S660b:** Prove or refute that paired deleted-score data can be recovered
  naturally from some enriched ordinary deck operation.
- **S660c:** Build the analogous point-deletion derivative deck for finite
  unit-distance graphs: card graph plus lost degree, direction support, and
  unit-spine damage.
- **S660d:** Compare paired deleted-score with THM-410 long-edge derivatives
  and S652 square-node substitution blocks.

## Honest Status

No general tournament reconstruction theorem is claimed.  The exact evidence
is through `n=6`.  A naive raw canonicalization attempt for `n=7` was too slow
for this session, and the script records that the next extension should use
canonical augmentation or automorphism pruning.

**See:** `04-computation/tournament_deck_derivative_s660.py`;
`05-knowledge/results/tournament_deck_derivative_s660.out`;
`07-reflections/tournament-deck-derivative-carrier-s660.md`; HYP-2237,
HYP-2235, HYP-2234, HYP-2228, HYP-2187, HYP-2176, HYP-2171.
