# Session Result

## Task Chosen

I chose a tiny exhaustive sanity check of the odd-cycle formula for labeled
tournaments on `n = 5` vertices.

At `n = 5`, no two odd directed cycles can be vertex-disjoint, so the conflict
graph independence polynomial has only the empty set and singleton odd-cycle
terms. Thus OCF specializes to

```text
H(T) = I(Omega(T), 2) = 1 + 2 * (c3(T) + c5(T)).
```

## What I Did

I ran a transient Python enumeration over all `2^10 = 1024` labeled
tournaments on five vertices. For each tournament I:

- counted Hamiltonian paths directly over all `5!` vertex permutations;
- counted directed 3-cycles and directed 5-cycles directly from cyclic orders;
- checked `H(T) = 1 + 2 * (c3(T) + c5(T))`;
- recorded small tournament fingerprints.

Tournament Analysis declaration:

- Pairwise observable: the orientation of each unordered vertex pair.
- Switch/gauge: the 10-bit labeled tournament mask, with each bit selecting one
  of the two possible arc orientations.
- Tie Hamiltonian path: not needed; the relation is already a tournament with no
  ties.
- Assumption challenge: I kept tournament vertices as the vertex set because the
  claim being checked is exactly about `H(T)` and odd directed cycles of `T`.
  Replacing vertices by arcs or cycle-objects would preserve some local conflict
  information but destroy the direct Hamiltonian-path count being tested.

## Concrete Result

No failures were found:

```text
n=5, labeled_tournaments=1024, failures=0
```

The Hamiltonian-path count distribution was:

```text
H=1:  120
H=3:  120
H=5:  240
H=9:  240
H=11: 120
H=13: 120
H=15: 64
```

Equivalently, the distribution of `c3+c5` was:

```text
0: 120
1: 120
2: 240
4: 240
5: 120
6: 120
7: 64
```

The sorted score-sequence fingerprint distribution was:

```text
(0, 1, 2, 3, 4): 120
(0, 1, 3, 3, 3): 40
(0, 2, 2, 2, 4): 40
(0, 2, 2, 3, 3): 120
(1, 1, 1, 3, 4): 40
(1, 1, 2, 2, 4): 120
(1, 1, 2, 3, 3): 240
(1, 2, 2, 2, 3): 280
(2, 2, 2, 2, 2): 24
```

## Confidence Note

Confidence is high for this bounded check. The enumeration is exhaustive over
all labeled `n=5` tournaments, and both sides were counted independently by
direct brute force rather than by using the claimed identity.
