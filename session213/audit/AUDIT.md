# Session213 independent audit: cycle-to-path transfer and endpoint equality

## Verdict and scope

The frozen files listed in `APPROVED_PACKETS.json` pass independent proof
review and the stated exact controls. The main new conclusions are a precise
restriction on the binary cycle compiler, and an all-tournament classification
of equality in the inherited endpoint half-bound. Neither conclusion proves
or refutes the unrestricted Hamiltonian-path spectrum conjecture.

No historical priority is claimed. The general attachment quadratic and
submodularity are inherited from Session203; one-apex recurrences, insertion
monotonicity, the cycle-cut bound, odd-cycle formula, endpoint half-bound, and
dominance-join multiplicativity are explicitly credited. The novel endpoint
classification upgrades the finite-only observation in the cited THM-4876
Audit D by supplying a missing all-order path witness.

## Proof review

### Exact comb path response

The recovered attachment identity partitions paths by whether the new vertex
is first, last, or internal. Removing an internal vertex produces an ordered
two-path cover, so the pair coefficient is exactly `-(K_ij+K_ji)`.

All evaluated comb coefficients were checked combinatorially, including the
two nonempty post-apex cases, the empty post-apex exceptions, endpoint
indicators, and the subtraction of empty cover bins. The bound
`H>=2^m+1` gives a genuine all-length cutoff. Independent naive quadratic
evaluation, rather than the author's Gray traversal, reproduces all 5,072
attained values through 65,535. The first missing odd target at least 23 is
35; 11 is also missing from this family. Neither is claimed globally absent.

The fixed-cycle/unbounded-path examples also prove the stated obstruction
to retaining every input as an induced core while bounding the output only
by its cycle count. Choosing one controlled input per cycle count is outside
that particular obstruction.

### Two independent sparsity proofs

The path-allocation proof counts two disjoint apex orders, each with three
choices for chain vertices below a last selected bit. This gives the stronger
bound `H >= (1+tau)2^m + 2V` and a geometric-series estimate on the number
of parameter presentations with `H<=X`.

The odd-cycle proof separately constructs a pair of vertex-disjoint odd
cycles. Its exact parity projection includes the `j=t=0` correction and
counts distinct collections, not rooted or ordered duplicates. It proves
the same exponent `alpha=log_3 2` with a larger constant.

Both estimates count all leading-zero padding and both apex orientations.
Counting presentations safely overcounts the numerical image. The resulting
`O(X^alpha)` statement is a bound, not an optimal exponent or an asymptotic
equivalence, and it concerns the defined comb family rather than every
two-apex tournament.

The finite-menu extension corollary requires a fixed finite number of
outputs per marked presentation and induced retention. The product-closure
addendum is also valid: dyadic summation makes the generator Dirichlet series
converge at every `beta>alpha`; the positive Euler product overcounts repeated
factorizations. Generator 1 is excluded from the product, finite extra seeds
are harmless, and density zero follows by choosing `alpha<beta<1`. The result
does not cover arbitrary connected interactions, an infinite independent
seed bank, or a new operation applied to a product.

### Arithmetic, deletion, and arc reversal

The sentinel core satisfies `H=2h+D(h-1)`. Its exact `2h+1` locus consists
of powers of two because every zero-before-one pair has positive weight.
The four-branch addendum correctly adds the forward-apex targets
`3*2^j+1` and no others for the specified deletion. This is distinct from
the image density of arbitrary one-apex counts.

The full reversed-comb readout uses the correct terminal count
`e_x1=(2-tau)epsilon_top`. The preliminary missing factor at `tau=0` was
repaired before freezing. The sharp upper deletion ray is stated with
`s>=1`, avoiding the false `g(1)=2` boundary. Arc reversal's positive,
exponential path increment and unique extremizers are family-specific;
no general arc-flip monotonicity follows.

### All-tournament endpoint equality

For positive endpoint count `e_a`, equality `H=2e_a+1` is equivalent to a
transitive prefix followed by the one-backedge tournament `R_m`, marked at
its penultimate vertex. The proof checks are:

1. The inherited injection gives `H>=2e_a+sum_{r in N+(a)}e_r`.
2. In a strong tournament, positive end counts at every vertex force a
   unique exit `a->r` and `e_r=1` under equality.
3. Deleting the forced final segment `a,r` identifies `e_r` with the path
   count of the remainder, so that remainder is transitive.
4. If `r` beats two vertices in that chain, the displayed spanning path
   ends in the remainder but is outside the injection image: deleting `a`
   creates a missing arc. Its validity includes the empty final-suffix case.
5. Strongness forces the sole exit of `r` to the first chain vertex. Strong
   component factorization then forces the entire earlier part to be
   transitive, since its multiplicative factor divides 1.

An independent equivalent recognizer is: `a` has sole outneighbor `r`,
`r` has outdegree one, and deleting `a,r` leaves a transitive tournament.
The unique exit of `r` identifies where the final strong block begins.

The positive-count assumption is essential. Transitive tournaments with
`e_a=0` have the degenerate equality `1=2*0+1`. Sinks give no positive
exception. The dual classification therefore makes a positive marked
endpoint count a power of two whenever `H=2c+1`. Applied to a sole-exit
cycle closure, this is an all-tournament obstruction to that endpoint-
preserving route, not an obstruction to an unrelated realizer of `2c+1`.

## Independent exact controls

- 511 one-apex cases, 127 sentinel cores, and 157 permutation controls.
- 504 complete comb variants through order ten, with 480 nonzero-code
  density checks, 270 literal disjoint-cycle certificates, and 56 complete
  permutation counts through order seven.
- 75 labeled attachment cores and all 1,098 cuts, checked against a
  permutation-and-cut two-cover oracle; 100 ordered coefficient entries.
- 262,136 independent full quadratic evaluations, matching the complete
  5,072-value support through 65,535 without Gray-code updates.
- 504 four-branch deletion cases and 126 exact `D=3` classifications.
- All 33,867 labeled tournaments through order six, with a separate C++
  endpoint DP and degree/remainder recognizer; the 456 and 6,880 inherited
  isomorphism-class representatives at orders seven/eight; and 1,981
  explicit paths outside the endpoint injection image.

The independent C++ equality counts agree with the author's doubled
orientation census after dividing the latter by two. Its supplementary
unlabeled universe is inherited and hash-pinned, not independently generated.
All normal/optimized author replays and independent replay receipts are
listed separately. Exact arithmetic is used throughout.

## Publication boundary

The audit made no edits to the stale shared mathematical checkout. Frozen
research bytes and source identities are recorded separately from publication
receipts. A proof approval does not by itself certify a remote push; the
session publisher verifies each checkpoint's remote bytes.
