---
id: HYP-2234
status: OPEN synthesis / small finite carrier evidence
source: codex-2026-06-05-S658
related:
  - HYP-2233
  - HYP-2235
  - HYP-2228
  - HYP-2227
  - HYP-2197
  - HYP-2187
  - HYP-2176
---

# HYP-2234: Frontier Problems Need Small Predicate-Preserving Carrier Atlases

## Claim

The HYP-2233 frontier lanes become useful only when they are converted into
small predicate-preserving carrier atlases.  The atlas should name:

```text
finite carrier,
predicate preserved by the projection,
side channel destroyed by the projection,
small repair/lift that restores useful information.
```

S658 tests this on three lanes: tournament reconstruction decks,
union-closed frequency carriers, and finite-field Kakeya/Falconer incidence
toys.

## Tournament Reconstruction Deck Lane

The script recursively generates tournament isomorphism classes through `n=6`
by vertex insertion, then compares four projections:

```text
scalar = (H, sorted score sequence, c3, SCC sizes)
full_deck = multiset of vertex-deleted card isomorphism types
scalar_deck = multiset of scalar(card)
loss_deck = scalar_deck plus sorted (H-H_card, c3-c3_card, deleted score)
```

Results:

```text
n=3 classes=2
  scalar max bucket 1, full_deck max bucket 2, scalar_deck max bucket 2,
  loss_deck max bucket 1.
n=4 classes=4
  scalar max bucket 1, full_deck max bucket 2, scalar_deck max bucket 2,
  loss_deck max bucket 1.
n=5 classes=12
  scalar max bucket 2, full_deck max bucket 2, scalar_deck max bucket 2,
  loss_deck max bucket 1.
n=6 classes=56
  scalar max bucket 4, full_deck max bucket 2, scalar_deck max bucket 2,
  loss_deck max bucket 1.
```

So raw vertex-deck and scalar projections collide in the small tournament
atlas, but the `H/c3/deleted-score` loss profile resolves every checked
collision.  This does not yet prove the loss profile is a pure deck invariant;
it identifies the first side-channel lift to audit.  Reconstruction progress
should therefore ask which loss components are reconstructible from cards and
which require extra owner/endpoint labels.

## Union-Closed Frequency Lane

The script enumerates nontrivial union-closed families up to coordinate
permutation for ground size `m<=4`, excluding the family containing only the
empty set.  It compares sorted coordinate-frequency vectors with two retained
set-side channels:

```text
size distribution = counts of sets by cardinality
set pressure(A) = sum_{B in F} |A union B|
```

Results:

```text
m=1 canonical families=2, failures=0, tight=1
m=2 canonical families=8, failures=0, tight=3
m=3 canonical families=36, failures=0, tight=6
m=4 canonical families=366, failures=0, tight=11
```

For `m<=3`, the sorted frequency vector separates all canonical families.  At
`m=4`, it becomes lossy:

```text
frequency buckets:        297 buckets, max bucket 4, 58 colliding buckets
frequency + size buckets: 336 buckets, max bucket 3, 22 colliding buckets
frequency + pressure:     366 buckets, max bucket 1, 0 colliding buckets
```

The half-frequency witness is easy in this range, but the closure structure
first starts hiding at `m=4`.  The cheap repair is not another element scalar;
it is a set-side pressure channel measuring how sets union with the whole
family.

## Finite-Field Kakeya/Falconer Lane

For each odd prime `p in {3,5,7}`, the script enumerates all choices of one
affine line in every direction of `F_p^2`.  There are `p^(p+1)` choices.  It
records the union size distribution, minimum-size line-choice sets, point-line
multiplicity profiles, and sampled distance/pinned-distance fingerprints.

Results:

```text
p=3: 81 choices, min size 7, min count 72
p=5: 15625 choices, min size 17, min count 3000
p=7: 5764801 choices, min size 31, min count 16464
```

The minimum sizes follow the visible odd-prime pattern

```text
min_size = (p^2 + 2p - 1)/2
```

in the checked cases.  Minimum multiplicity profiles are stable:

```text
p=3: ((1,3),(2,3),(3,1))
p=5: ((1,6),(2,9),(3,2))
p=7: ((1,9),(2,19),(3,3))
```

The sampled minimum sets also carry distance side channels.  For `p=5`, sampled
minimum sets realize all `5` quadratic-distance values, with pinned counts
either all `5` or one point pinned at `4`.  For `p=7`, sampled minimum sets
realize all `6` nonzero distance values, with pinned counts mostly all `6`.

This is the right finite toy because the Kakeya predicate is direction
coverage, while Falconer-style side information is distance-fiber coverage.
The line-choice carrier preserves the first and lets the second be measured.

## Cross-Lane Tournament Analysis

Vertices are the three worked lanes:

```text
tournament_decks
union_closed_frequency
finite_field_incidence
```

The pairwise observable is

```text
(finite_exactness, ambiguity_reduction,
 side_channel_visibility, near_term_actionability).
```

The route tournament is transitive:

```text
union_closed_frequency
> tournament_decks
> finite_field_incidence
```

with `score_hist={0:1,1:1,2:1}`, `directed_3cycles=0`, and one Hamiltonian
path.  The ordering is not a statement of mathematical importance; it is a
work-priority ranking for immediate finite carrier progress.  Union-closed
families win because `m=4` exposes an exact scalar failure and an exact cheap
repair.  Tournament decks are next because the loss-profile repair is strong
but still needs a deck-invariance audit.  Finite-field incidence is third
because it gives clean exact data but the next theorem likely needs known
finite-geometry context.

## Work Queue

- **S658a:** Prove or refute that the tournament `loss_deck` components are
  reconstructible from the ordinary vertex-deck, or isolate the smallest
  owner-label needed.
- **S658b:** Classify the `m=4` union-closed frequency collisions and ask
  whether set pressure is a canonical closure-side sufficient statistic beyond
  `m=4`.
- **S658c:** Extend the finite-field incidence atlas to affine equivalence
  classes of minimum Kakeya line-choice sets and compare their distance fibers.
- **S658d:** Use the same "scalar projection -> side-channel repair" template
  on the HYP-2233 Hadwiger-Nelson finite unit-distance lane.

## Honest Status

No external famous problem is solved here.  The progress is that three
frontier lanes now have small exact carriers and concrete next lemmas.  The
common pattern is the repo's recurring scalar-shadow lesson: raw counts,
decks, frequency vectors, and direction coverage become useful only after the
destroyed side channel is named and reattached.

S659 extends the finite-field Kakeya/Falconer lane as HYP-2235, adding a
dedicated distance/pinned/unit ledger, exact `p<=7` minima, and an LRC `n=14`
no-leak transfer statement.  Tournament reconstruction decks and union-closed
frequency carriers remain open lanes under this atlas.

**See:** `04-computation/frontier_small_carriers_s658.py`;
`05-knowledge/results/frontier_small_carriers_s658.out`;
`07-reflections/frontier-small-carrier-atlases-s658.md`;
`04-computation/finite_field_kakeya_falconer_s659.py`;
`05-knowledge/results/finite_field_kakeya_falconer_s659.out`;
`07-reflections/finite-field-kakeya-falconer-lrc-carrier-s659.md`; HYP-2235,
HYP-2233, HYP-2228, HYP-2227, HYP-2197, HYP-2187, HYP-2176.
