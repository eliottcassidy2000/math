# Odd-layer conjugacy, the divisor monoid, and the harmonic atlas

**Research reflection / provenance, not a truth source.**  Exact claims are
routed to the audit-pending THM-3472 and its proved dependencies THM-3405 and
THM-3453.  LRC(14) remains open.

## The map that closes the typing gap

The inheritance board is unusually clean:

- closest proved mechanism: THM-3405's primitive-divisor/two-layer reduction;
- canonical near miss: THM-3469 used half-layer lower bounds before making the
  omitted fixed-zero layer explicit;
- hostile: at even `Q=8`, `ell -> 2ell+1` reaches only odd sheets;
- least-used sidecar: the active primitive divisor `Q`, not merely ambient
  modulus `q`.

For odd `Q`, the exact transport is

```text
primitive fixed-zero cover on Q
  -- ell -> 2ell+1; s -> 2s --> primitive half cover on Q
  -- divisor dilation ----------> literal half cover on odd q.
```

Every arrow preserves strict incidence, cover cardinality, and active gcd.
The final grading forgets which primitive divisor supplied the cover.  This
explains THM-3464 without overstating it: `Q=41` and `Q=123` are disjoint
primitive layers whose minima happen to agree at eight.

## A monoid statement, not a list statement

Divisor dilation composes multiplicatively.  Primitive mask packets are the
atoms; multiplication by an odd scale moves them through ambient moduli
without changing their grade.  The fixed-zero layer does not add a new atom
on odd moduli because the sheet conjugacy embeds it into the half layer.

This gives a precise version of the user's degree-graded monoid viewpoint:

```text
object       = primitive labelled mask packet,
operation    = divisor dilation,
grade        = cover cardinality,
unit data    = owner signs and odd sheet permutation,
forgotten    = primitive ancestry after taking the minimum.
```

The grade is not a complete invariant.  Equal grades may come from several
primitive layers, and a nonminimal representation may still be irredundant.
That last distinction is exactly what the THM-3469 private-sheet frontier now
tests.

## Why the carrier is stronger than a tournament

Here the source and target are labelled cover hypergraphs, and the map is an
actual sheet permutation.  It preserves every intersection order, not just
pairwise observables.  A tournament on owners would lose the map because
there is no intrinsic orientation between two masks.  XOR would also be too
small: multiplicities are irrelevant to union coverage but remain useful
sidecars for private-sheet and current questions.

The ternary word `uncovered/private/multiple` is a faithful coverage summary
only after the owner-by-sheet incidence matrix remains attached.  THM-3472
uses the full matrix isomorphism; it does not infer equality from a ternary or
pairwise shadow.

## Exact harmonic subsets

The odd cap-seven theorem turns four natural-number subsets into periodic
sets.  Their ambient harmonic coefficients are

```text
rank 4:   1/18,
rank 6:   4/45,
rank 7:   283376/6633315,
rank >7:  691712/2211105.
```

Thus “every subset of the naturals is a subset of the harmonic series” becomes
a theorem only when a distribution mechanism is supplied.  Here that
mechanism is a `576`-state CRT product with minimal annotated period
`729,664,650`.  Arbitrary subsets have no such coefficient for free.

## New frontiers

1. **Even equality.**  The odd proof fails at `Q=8`, but that does not prove
   `rho_ZMC(q)<rho_H(q)` for any even `q`.  Search for the first true rank
   separation, or construct a different even-layer transport.
2. **Irredundant versus minimum.**  Prove the exact private-sheet formulas for
   THM-3469's eight-owner family.  If every owner is private at all `k`, the
   family is irredundant even on rank-4/6/7 grades, cleanly separating local
   deletion minimality from global rank.
3. **Current-aware transport.**  The sheet permutation preserves Boolean
   incidence and zero cochain, but it does not manufacture a labelled
   endpoint current.  Add owner modes, boundary labels, and target activity
   before comparing relation-residue coefficients.
4. **Spectral closure.**  Equality of cover ranks is upstream of the LRC
   `7 tensor 13` bispectrum.  A bridge must preserve the current contraction,
   not merely the cover hypergraph.

The central lesson is that a small affine permutation can collapse an entire
omitted layer, but only because its divisor, parity, sign, and incidence
sidecars are all retained.
