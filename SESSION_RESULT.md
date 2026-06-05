# Session Result

## Task Chosen

I chose a tiny exhaustive sanity check connected to the local `TournamentH7`
thread: verify directly, for all labeled tournaments on `n <= 5` vertices, that
no tournament has exactly `H(T) = 7` Hamiltonian paths.

This is only a low-order verification of the stronger `H(T) != 7` theorem
recorded in `04-computation/lean/TournamentH7/README.md`, not a proof of the
general theorem.

## What I Did

I ran a transient brute-force Python enumeration. For each `n = 1, 2, 3, 4, 5`,
I enumerated all `2^(n choose 2)` labeled tournament orientations. For each
tournament I counted Hamiltonian paths independently by testing all `n!`
vertex permutations against the directed edges.

Tournament Analysis declaration:

- Pairwise observable: the orientation of each unordered vertex pair.
- Switch/gauge: the labeled bit mask over unordered pairs, with each bit
  choosing one of the two possible arc directions.
- Tie Hamiltonian path: not needed, since every pairwise comparison is already
  decisive in a tournament.
- Fingerprints recorded: full `H(T)` value histograms for each `n <= 5`.

Assumption challenge:

- I kept tournament vertices as the vertex set because the claim being checked
  is directly about Hamiltonian paths in `T`.
- Alternate vertices such as arcs, odd cycles, score gaps, residues, or proof
  obligations could preserve local orientation or obstruction data, but would
  destroy the direct one-to-one permutation test for `H(T)`.

## Concrete Result

No failures were found. The exhaustive distributions were:

```text
n=1, labeled_tournaments=1,    H=7 count=0, distribution={1: 1}
n=2, labeled_tournaments=2,    H=7 count=0, distribution={1: 2}
n=3, labeled_tournaments=8,    H=7 count=0, distribution={1: 6, 3: 2}
n=4, labeled_tournaments=64,   H=7 count=0, distribution={1: 24, 3: 16, 5: 24}
n=5, labeled_tournaments=1024, H=7 count=0, distribution={1: 120, 3: 120, 5: 240, 9: 240, 11: 120, 13: 120, 15: 64}
```

Thus the direct exhaustive check confirms `H(T) != 7` for every labeled
tournament with `n <= 5`.

## Confidence Note

Confidence is high for the bounded claim. The enumeration is exhaustive over
all labeled tournaments for `n <= 5`, and `H(T)` was counted directly from
permutations rather than inferred from the forbidden-value theorem.
