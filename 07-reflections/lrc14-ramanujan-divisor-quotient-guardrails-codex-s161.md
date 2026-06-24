# LRC14 Ramanujan-Divisor Quotient Guardrails

The cleanest theorem shape from this pass is not a new numerical invariant.  It
is a rule about proof maps:

```text
A quotient is allowed to forget only labels that are irrelevant to the next
proof predicate, or labels that are explicitly recoverable by a certificate.
```

This is the common lesson behind the repo's recurring scalarization failures.
Irreducibility is not just degree.  A unital is not just point/block counts.
Faulhaber moments are not just positivity.  Pollock defects are not just a
polygonal number.  Unit-distance carriers are not just a norm-1 count unless we
know which norm layer is meant.  Tilings and solids are not just cell counts
unless the boundary-gluing unit is retained.  The Farey/perfect-number product
`ab=|E(K_{a,b})|` is useful, but it is not a graph-minor obstruction without
the incidence labels.

The divisor-function page made this sharper rather than merely metaphorical.
The scalar divisor functions are multiplicative, `sigma = Id * 1` under
Dirichlet convolution, and admit Ramanujan-sum expansions.  So even a scalar
`sigma_k(n)` secretly carries a primitive-root packet expansion.  In LRC14
language: divisor data becomes proof-relevant only after we decide whether the
Ramanujan/cyclotomic packet, endpoint owner, boundary equality, or K33 debt is
still visible.

The exact audit gives a finite warning.  On `2694` named and one-swap rows,
scalar divisor signatures have `138` mixed qdiv/safe-route fibers; unitary
divisor data reduces but does not remove the problem; Ramanujan speed and pair
traces still mix routes; exact-period packets nearly work but still identify AP
with positive `12->96`.  That last failure is the best clue: the Ramanujan
projector is the right language, but not yet the theorem object.  It needs the
open-vs-boundary label.

Thus HYP-2978 and HYP-2979 should merge as a two-layer object:

```text
exact-period Ramanujan projector
  + endpoint-owner / safe-measure / AP-GW boundary label
  + K33 or HYP-2908/THM-572 state-lift debt
```

The tournament analysis used quotient candidates as vertices rather than
runners.  The pairwise observable was not "which quotient is prettier," but
which quotient has fewer mixed fibers for the LRC route predicate.  The
transitive Hamiltonian path

```text
endpoint_measure > full_row > exact_period_packet > unitary_divisor >
scalar_divisor > ramanujan_pair > ramanujan_speed > qcover
```

should be read cautiously.  `endpoint_measure` wins only for the route predicate
we asked it to preserve.  It is not a proof of LRC14; it is a signpost saying
the proof object must recover exact open/boundary information after quotienting.
This lines up with HYP-2975 taut bridges and HYP-2970 endpoint-credit cycles.

The next concrete proof attempt should be a labelled packet theorem:

```text
Every non-AP/GW LRC14 residual, after exact-period Ramanujan projection and
endpoint-owner refinement, either has a positive strict-safe certificate or
routes to a K33/HYP-2908/THM-572 state-lift debt.
```

That formulation keeps the divisor/Ramanujan insight but avoids the trap of
treating multiplicative functions as final scalars.  Multiplicativity is the
irreducibility ledger.  Ramanujan sums are the phase ledger.  Endpoint labels
are the boundary ledger.  A proof must carry all three until the exact predicate
being proved no longer depends on them.

The concurrent S161 synthesis adds four useful inquiry leads.  Regular
A-function Ramanujan sums say that the selected divisor subset is itself part of
the quotient contract.  Busche-Ramanujan identities are the right template for
non-coprime product repairs when shared prime support has been forgotten.
Divisor-summatory hyperbola errors look like boundary-defect packets, so they
belong beside Pollock and tiling defects.  Euler pentagonal sigma recurrences
may offer sparse harmonic boundary moments for the Toeplitz/spectral stack.
