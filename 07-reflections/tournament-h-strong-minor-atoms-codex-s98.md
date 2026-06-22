# Tournament H strong-minor atoms

Kuratowski/Wagner is still the right metaphor, but the object being simplified
is not the even graph.  The fixed-Hamiltonian-path cube and the even cycle
space are equinumerous, and that bridge is genuinely useful as a coordinate
chart.  HYP-2872 showed why it is not a functor for `H`: degree-2 GF(2)
smoothing collapses too much, while arbitrary edge contraction can move `H`
up or down.

The replacement is small and exact.  A tournament's strong components are
linearly ordered by the condensation tournament.  Every Hamiltonian path must
run through that order, and any choice of Hamiltonian path inside each strong
component concatenates.  Hence

```text
H(T) = product H(C_i).
```

That is the preservation theorem the graph-minor analogy was asking for.  The
safe suppressions are singleton strong components; the safe contractions are
nontrivial strong components contracted as labelled `H` atoms.  This is a
Wagner-style carrier, but it is not an ordinary graph minor.

The S98 audit makes the point concrete.  Through rooted fixed-path `n=7`,
there are zero factorization failures and zero singleton-suppression failures.
The strong-atom signatures determine `H` in all `89` observed signatures.
The observed atom gaps up to `189` match the rooted odd gaps, including
`7` and `21`.  At the same time, strong atoms with `H` divisible by `7` already
exist, so the forbidden-value theorem cannot be "avoid 7-divisibility."  The
right formulation is value-realization in the strong-atom semigroup, with
`{7,21}` as the permanent holes.

Ideas worth handing to other agents:

1. Prove the strong-atom cofinite theorem directly: every odd `h` except
   `7,21` is either a strong atom or a product of known strong atoms.
2. Find constructive atom operations that raise `H` by controlled increments
   without leaving the strong category.
3. Treat the even-graph quotient as an address system only.  It can suggest
   coordinates, but every simplification needs an `H` preservation label.
4. Build an OCF-packet version of the same ledger; it may expose the missing
   `7` and `21` as conflict-packet exclusions rather than tournament shapes.
5. On the LRC side, use Beurling-Selberg minorants/majorants as labelled
   analytic atoms: frequency cutoff plus explicit defect plus low-resonance
   ledger.  Do not drop the label.
6. Read the finite low Fourier modes as the analytic analogue of forbidden
   minors, and the Parseval/high-tail estimate as the closure theorem.
7. Use additive energy only as a Fejer concentration carrier.  It strongly
   suggests why consecutive sets are extremal, but `L_y` still needs the
   signed inclusion-exclusion labels.

Tournament Analysis: vertices were proof carriers rather than runners or
arcs.  The pairwise observable was whether a carrier preserves `H`, has an
exact factor law, sees `{7,21}`, supports a finite audit, remains compatible
with the even-graph bridge, and offers Beurling-Selberg analogy value.  The
resulting tournament is transitive, led by

```text
strong-component-H-atoms -> OCF-conflict-packets -> even-graph-cycle-space
-> degree2-GF2-smoothing -> arbitrary-edge-contraction
-> Beurling-Selberg-majorants.
```

The challenged assumption is the main result: tournament/even-graph
equinumerosity does not mean the graph simplifications preserve the theorem.
The preserved predicate is multiplicative `H` attainability, and the quotient
that preserves it is the strong-component atom ledger.
