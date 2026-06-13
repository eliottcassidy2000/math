# S660 Reflection: Tournament Deck Derivative Carrier

Prompt: think creatively and explore deeply among many tangents the tournament
deck angle; see how it applies and connects.

The deck angle clicked when I stopped treating reconstruction as a yes/no
question and started treating the deck as a projection.  A vertex-deleted deck
forgets the missing vertex's boundary label.  So the natural repairs are not
more global totals.  They are derivatives:

```text
what changed when this vertex was removed?
```

## Main Computation

S660 compares repair channels for tournament isomorphism classes through
`n=6`.  The ordinary full deck still collides:

```text
n=3: max full-deck bucket 2
n=4: max full-deck bucket 2
n=5: max full-deck bucket 2
n=6: 4 full-deck collision buckets, max bucket 2
```

Global scalar repair is not enough at `n=6`.  Adding `H`, score sequence, or
the full scalar `(H, score, c3, SCC)` to the full deck still leaves the same
four collision buckets.

The clean repair is:

```text
(card isomorphism type, deleted vertex outdegree).
```

That paired deleted-score deck separates every checked class through `n=6`.
Unpaired deleted scores do not.  The pairing matters.  The information is not
"what degrees existed" but "which deleted boundary degree belongs to which
card."

This is pleasantly close to the repo's owner/carry language.  The card is the
quotient; the deleted score is the owner label on the missing boundary.

## Collision Shape

At `n=6`, every full-deck collision pair in the atlas is a converse pair.  Each
pair has the same `H`, the same score sequence, the same `c3`, and the same SCC
shape.  So the collision is not living in familiar scalar territory.  It lives
in the assignment of boundary data to cards.

The pairs are:

```text
0x18  vs 0x1c   H=23, c3=6, scores=(1,2,2,3,3,4)
0x92  vs 0x95   H=31, c3=6, scores=(1,2,2,3,3,4)
0xb1  vs 0xb5   H=29, c3=6, scores=(1,2,2,3,3,4)
0x158 vs 0x159  H=43, c3=8, scores=(2,2,2,3,3,3)
```

They are tiny flip-distance changes, but the full deck cannot see which card
gets which missing-vertex score.

## Tangent Web

This gives a compact dictionary:

```text
tournament deck: card + deleted score
LRC n=14: Res_27 atom + owner/carry/Cprime
unit distance: point-deleted graph + lost degree/gain/direction support
union-closed: frequency vector + set pressure
finite-field Kakeya: direction cover + pinned/concurrency labels
CH: cardinal shadow + model/generic side channel
```

The same shape keeps recurring.  The scalar quotient is useful only as a
cache.  The proof object is the projection plus a fiber derivative.

In LRC terms, the deck is like a row after reducing modulo `C=27`.  The
deleted score is like the carry coordinate: a small label that says how the
object reattaches to the quotient.  A raw residue row or raw deck can look
identical while the boundary assignment differs.

In unit-distance terms, the analogue is not just the graph deck.  It is:

```text
delete point p,
record the remaining unit graph,
record lost degree,
record lost direction support,
record whether the unit spine or dense core was damaged.
```

That is exactly the S614/S617/S622/S623 pressure story in a deck dialect.

In union-closed terms, `set pressure(A)=sum_B |A union B|` is a deletion
derivative for the closure operation.  Frequency alone names how often
elements appear; pressure names how a set interacts with the whole closure
system.  The tournament deck's paired deleted score is the same sort of
interaction label.

Incoming S659's finite-field Kakeya/Falconer work makes the same point:
direction coverage and distance support saturate early, so the proof-relevant
labels are pinned distances, concurrency multiplicity, and owner labels.

Incoming S661's LRC14 wild-strategy claim is the same pattern one level closer
to the main target: the no-leak theorem should not trust scalar wall/support
data until owner, carry, and deletion-derivative labels are attached.

## Why This Matters

The deck angle reframes reconstruction as a general method:

1. Choose a projection that is natural but lossy.
2. Find the first collision bucket.
3. Ask what local derivative splits it.
4. Check whether that derivative is intrinsic, deck-visible, or needs marking.
5. Export the derivative to other problems as a side-channel candidate.

For the tournament lane, the next theorem target is not "all tournaments are
reconstructible" in one leap.  It is:

```text
paired deleted score resolves full-deck collisions through n=6;
does this remain true through n=7, and is there a natural deck-visible way to
recover the pairing?
```

The naive `n=7` canonicalization was too slow, so the next computational step
should use canonical augmentation or a nauty-style automorphism-pruned
generator.

Artifacts: HYP-2236, T729, `04-computation/tournament_deck_derivative_s660.py`,
and `05-knowledge/results/tournament_deck_derivative_s660.out`.
