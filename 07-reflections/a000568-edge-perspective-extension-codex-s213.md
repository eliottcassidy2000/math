# A000568 edge-perspective extension - S213

S211 answered the count question: the shifted equality stops at `n=6` because
`P(5)=48` while `U(6)=56`.  The extra S213 observation is that the next rung is
not vague.  A rooted 5-perspective plus the new observer's incident word is
exactly an ordered-pair perspective on a 6-tournament.

The exact counts are:

```text
U(5)=12, P(5)=48, U(6)=56
rooted 5-perspective + incident word states = 1408
ordered-pair perspectives on U(6) = 1408
directed-edge perspectives on U(6) = 704
```

So the edge perspective is the old/new-role quotient of the full extension
state.  It is natural and dualistic, but it already forgets which endpoint was
the old root and which was the new observer.  That is harmless for abstract
tournament counting only when the sidecar is not being used to certify an
observer predicate.  In LRC it is not harmless unless endpoint-owner or
threshold sidecars restore the role.

The compact sector-deck result is the sharpest part.  Around an ordered pair
`(r,o)`, put every other vertex into one of four sectors
`(r beats x, o beats x)`.  At `n=6`, the multiset deck of sector sizes separates
`55/56` classes, and adding the internal tournament types in each sector still
separates only `55/56`.  The only collision is the converse pair with masks
`344` and `345`, both with score `(2,2,2,3,3,3)`, `c3=8`, `H=43`, and trivial
automorphism group.  Adding cross-sector orientation counts separates all
`56/56`.

That says the first missing coordinate is literally cross-sector orientation.
It is the tournament toy version of the controlled-forgetting rule in the LRC
packet work: a quotient can forget the root role, sector orientation, or owner
address only after a sidecar has retained, reconstructed, annihilated, or
routed that coordinate.

For the prompt's creative menu:

1. Edge-perspectives are the first honest lift after rooted nodes.
2. Cycle-perspectives should be used to explain the one converse/chirality
   collision.
3. Clique/insertion perspectives should treat observer insertion as a
   low-rank update, not as another node color.
4. Conflict perspectives should live on edge-sector/cycle incidence cells.
5. The LRC-facing object is still endpoint-owner packet sheaf: source
   threshold arc, incident word, cross-sector orientation, owner labels, and
   proof-obligation state.

The current best slogan:

```text
depth-2 node color recovers the rooted cache;
edge-sector cross orientation recovers the first extension sidecar.
```

This also completes the immediate HYP-3048 matrix-atlas pull: the four-sector
edge block matrix is not only a suggested encoding.  At the first failure,
sector size/internal block data gives `55/56`, and adding the cross-sector
orientation block gives `56/56`.  In observability-matrix language, the first
separating column is `cross_sector_orientation_word`.

The next exact test is whether the same sector-cross deck remains a useful
compressor at `n=7`, and whether the single converse-pair miss of the
size/internal deck is the first member of a stable chirality-defect family.
