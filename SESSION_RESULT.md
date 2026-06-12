# Session Result

## Task Chosen

I chose one tiny exhaustive sanity check of the odd-cycle collection formula
at `n=5`.  At this size, two vertex-disjoint odd cycles cannot coexist, so
the formula specializes to the directly checkable identity

```text
H(T) = 1 + 2 * (# directed 3-cycles + # directed 5-cycles).
```

This is small enough to verify over all `2^10 = 1024` labeled tournaments.

## What I Did

I ran a transient Python enumeration over all labeled tournaments on five
vertices.  For each tournament I computed:

- `H(T)` by direct permutation scan of all Hamiltonian paths,
- the number of cyclic triples,
- the number of directed 5-cycles, counted up to cyclic rotation, and
- the specialized OCF right-hand side above.

Tournament Analysis declaration:

- Pairwise observable: the orientation bit on each unordered vertex pair.
- Switch/gauge: the bit value gives the binary tournament edge.
- Tie Hamiltonian path: no tie handling is needed for tournaments; all
  distinct vertex pairs have exactly one directed edge.
- Fingerprints recorded: `H` histogram, cyclic-triple histogram,
  directed-5-cycle histogram, and failure count.

Assumption challenge:

- I kept tournament vertices as the original five vertices because this task
  tests the literal OCF statement.
- Alternate vertex sets such as arcs, odd cycles, proof obligations, or
  Fourier modes are useful for other quotient views, but would obscure this
  direct small-`n` identity.  The odd-cycle conflict graph is still present:
  at `n=5` its independent sets have size at most one.

## Concrete Result

The identity held for every labeled tournament on five vertices.

```text
n = 5
tournaments = 1024
failures = 0
H_hist = {1: 120, 3: 120, 5: 240, 9: 240, 11: 120, 13: 120, 15: 64}
c3_hist = {0: 120, 1: 120, 2: 240, 3: 240, 4: 280, 5: 24}
c5_hist = {0: 480, 1: 360, 2: 144, 3: 40}
max_c5 = 3
```

Thus the direct `n=5` specialization of OCF passes an exhaustive independent
check, including the contribution from directed 5-cycles.

## Confidence Note

Confidence is high for this bounded verification.  The enumeration is
exhaustive, the path count is by direct permutation scan, and the cycle counts
are independent of any repository helper scripts.
