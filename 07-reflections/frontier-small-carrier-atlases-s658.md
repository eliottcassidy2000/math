# S658 Reflection: Frontier Small-Carrier Atlases

Prompt: spend a long session working on some HYP-2233 missed-frontier lanes.

I worked three lanes because they are small enough to compute honestly:
tournament reconstruction decks, union-closed frequency carriers, and
finite-field Kakeya/Falconer incidence toys.

The shared lesson is sharper than expected:

```text
important problem -> small carrier -> scalar collision -> side-channel repair.
```

## Tournament Decks

The tournament deck lane produced the most surprising local result.  Through
`n=6`, raw vertex-decks by card isomorphism type still collide:

```text
n=3 full_deck max bucket 2
n=4 full_deck max bucket 2
n=5 full_deck max bucket 2
n=6 full_deck max bucket 2
```

Scalar invariants also collide, with `n=6` having a scalar bucket of size `4`
for `(H, score sequence, c3, SCC sizes)`.  But adding a sorted loss profile

```text
(H(T)-H(T-v), c3(T)-c3(T-v), deleted vertex score)
```

to the scalar card-deck resolves every checked bucket through `n=6`.

This is not yet a reconstruction theorem, because the loss profile may include
side data not directly visible from the ordinary deck.  But it names the next
proof obligation: determine which parts of the loss profile are deck-visible
and which require an owner/endpoint lift.  That is exactly the repo's
side-channel style of progress.

## Union-Closed Families

For union-closed families up to coordinate permutation:

```text
m=1: 2 canonical families
m=2: 8
m=3: 36
m=4: 366
```

Frankl's half-frequency witness has zero failures in this range.  The
interesting event is not the witness; it is the first scalar collision.  Sorted
coordinate-frequency vectors separate all canonical families through `m=3`,
then fail at `m=4`:

```text
freq only:        297 buckets, max bucket 4, 58 colliding buckets
freq + size dist: 336 buckets, max bucket 3, 22 colliding buckets
freq + pressure:  366 buckets, max bucket 1, 0 colliding buckets
```

The repair is set-side, not element-side:

```text
pressure(A) = sum_{B in F} |A union B|.
```

That is a good sign.  It says the frequency witness lives on elements, but the
closure predicate lives on sets.  A proof route should not try to squeeze
closure out of frequency alone; it should retain a closure-pressure channel.

## Finite-Field Incidence

For `F_p^2`, I enumerated all choices of one affine line in every direction:

```text
p=3: 81 choices, min Kakeya size 7
p=5: 15625 choices, min Kakeya size 17
p=7: 5764801 choices, min Kakeya size 31
```

The checked minimum sizes fit

```text
(p^2 + 2p - 1)/2.
```

Minimum line-choice carriers have stable point-line multiplicity profiles:

```text
p=3: ((1,3),(2,3),(3,1))
p=5: ((1,6),(2,9),(3,2))
p=7: ((1,9),(2,19),(3,3))
```

Distance fibers are the Falconer side channel.  Sampled minimum sets for `p=5`
realize all `5` quadratic-distance values.  For `p=7`, they realize all `6`
nonzero distance values in the samples.  So the line-choice carrier preserves
Kakeya direction coverage while making distance coverage measurable.

## Priority

Cross-lane Tournament Analysis ranks immediate next work as:

```text
union_closed_frequency
> tournament_decks
> finite_field_incidence.
```

That ranking is about near-term repo traction, not external importance.
Union-closed wins because the scalar failure and repair are both exact at
`m=4`.  Tournament decks are close behind because the loss-profile repair is
very strong but needs an invariance audit.  Finite-field incidence is clean but
will probably need more finite-geometry context before it becomes theorem
material.

The next useful move is not to chase a famous statement directly.  It is to
turn each side-channel repair into a lemma:

- loss-profile deck visibility;
- closure-pressure sufficiency for union-closed collisions;
- affine classification of minimum Kakeya line-choice carriers;
- then repeat the same method on finite unit-distance coloring.

Artifacts: HYP-2234, T727, `04-computation/frontier_small_carriers_s658.py`,
and `05-knowledge/results/frontier_small_carriers_s658.out`.
