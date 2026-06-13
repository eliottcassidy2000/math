# HYP-2455 - Boundary-lift irreducibility transfer unifies polynomial, LRC, unit-distance, and code support gates

**Status:** OPEN synthesis; exact local anchors plus a reproducible carrier tournament.
**Source:** codex-2026-06-13.
**Companions:** HYP-2454, HYP-2453, HYP-2452, HYP-2451, HYP-2450,
HYP-2446, HYP-2445, HYP-2444, HYP-2301, HYP-2300, OPEN-Q-057.
**Computation:** `04-computation/boundary_lift_analogy_atlas_codex.py`;
stored output `05-knowledge/results/boundary_lift_analogy_atlas_codex.out`.

## Statement

Recent agents have produced several instances of the same proof grammar:

```text
visible scalar or boundary total  ->  hidden lift / support certificate
```

The working hypothesis is that this is not just a metaphor.  It is a portable
proof interface.  In each domain, the "reducible" object is the visible scalar
or product decomposition, while the proof-bearing object is the hidden lift:

```text
integer polynomial coefficients  -> integral convolution factor grid
LRC blocked denominator          -> runner/Pisano/divisor/owner support lift
unit-distance product count      -> non-product Moser-ring incidence fiber
triangular equalities            -> moment/fractional/Pell address lift
[72,36,16] weight enumerator     -> binary support-design incidence lift
```

Thus irreducibility should be read broadly as absence of a nontrivial hidden
factor/product/support lift, not merely as an algebraic property of a polynomial.
The transfer target is a common "lift feasibility" data model: boundary totals,
candidate hidden cells, local gates, surviving allocations, and proof owners.

## Anchors Already In The Repository

1. **Polynomial irreducibility.**  HYP-2452 makes reducibility literal:
   `a_k=sum_{i+j=k}b_i c_j`.  Stored degree-4/5 primitive scans have zero
   mismatches against symbolic reducibility.  This is the cleanest exact model.

2. **LRC14.**  HYP-2444 shows the plain scalar horizon is too coarse: Q27 covers
   all one-stranger rows, and the visible residual is the missing antipodal class
   `+-10 mod 27` plus the 13-clock.  HYP-2446 adds the operator warning: scalar
   local tests can lie both ways.

3. **Unit distance.**  OPEN-Q-057 and HYP-2300/HYP-2301 say the unit-distance
   threshold is irreducible in the product sense.  Products tie at `N=27` and
   first beat `3N` only at `N=32`; the known `N=28` crosser lives in the rank-4
   Moser ring rather than a product or rank-2 lattice shadow.

4. **Triangular towers.**  HYP-2453/HYP-2454 expose the scalar danger in a
   gentler toy model.  The p=1 and p=2 balances are exact, but p>=3 needs a
   fractional/moment address in the checked range.  The `78/90` row should be
   treated as a support address, not as stand-alone numerology.

5. **The `[72,36,16]` code frontier.**  Weight enumerator feasibility is only a
   boundary total.  Minimum-word supports and design/matroid incidence are the
   hidden lift.  The repeated `78` beacon is useful only if it becomes a
   constraint on that incidence lift.

## Main Analogy

The unit-distance "product family" is analogous to polynomial reducibility.
It is not bad; it is the structured factorized baseline.  The theorem pressure
comes exactly where the baseline cannot cross:

```text
polynomial:     reducible iff a nontrivial convolution lift exists
unit distance:  product-reducible iff a Cartesian/Minkowski factor explains the count
LRC:            scalar-reducible iff q-blocking forgets runner-owner support
code72:         enumerator-reducible iff weights forget support incidence
```

This suggests a concrete search discipline:

```text
1. Identify the scalar/product shadow.
2. Define the hidden lift variables.
3. Add local gates that prune lift feasibility.
4. Keep survivor ledgers instead of collapsing them to one number.
5. Prove the frontier by showing all lifts either realize a known easy structure
   or hit an obstruction.
```

For LRC14, this means multi-stranger rows should be represented like
factor-capture allocations: each blocked twist has tokens, runners, divisor
fibers, owner targets, and carries.  For unit distance, it means the `N=27/28`
frontier should be split into product-reducible fibers and Moser-irreducible
extension fibers.  For `[72,36,16]`, it means the search should move from
Gleason coefficients to binary incidence lift feasibility.

## Computation

The stored atlas compares seven carriers:

```text
polynomial_convolution_lift
lrc_q27_support_lift
unit_distance_moser_irreducible_fiber
triangular_moment_fractional_address
code72_support_design_lift
p_curvature_operator_ledger
product_quotient_diagonal_support_gate
```

Tournament Analysis uses carriers/proof obligations as vertices.  Pairwise
observable: majority comparison across

```text
exact_certificate, hidden_lift_visibility, irreducibility_signal,
lrc14_transfer, unit_distance_transfer, code72_transfer, computable_now,
scalar_warning.
```

The resulting fingerprint is intentionally nontrivial:

```text
score_hist={1: 2, 2: 1, 3: 2, 5: 1, 6: 1}
directed_3cycles=3
scc_sizes=[5, 1, 1]
hamiltonian_paths=9
leader=polynomial_convolution_lift
edge_flips_vs_scalar_warning=11/21
```

The edge flips matter.  If one ranks only by "scalar warning", the code and
operator ledgers look maximally urgent; the majority carrier tournament instead
keeps polynomial convolution and LRC Q27 support as the most useful engines
because they are exact and computable now.

## Immediate Research Program

1. **Common lift schema.**  Define a lightweight JSON/CSV schema for boundary
   totals, hidden cells, local gates, surviving allocations, and proof owners.

2. **Polynomial degree 6.**  Extend HYP-2452 with bounded ILP/SAT convolution
   lifts, using degree-4/5 as regression tests.

3. **LRC multi-stranger.**  Replace "q blocked" by an allocation ledger over
   twists, runners, Pisano classes, divisor fibers, carries, and Bprime owners.

4. **Unit-distance reducible/product split.**  Classify candidate `N=27/28`
   extensions by product-reducible versus Moser-irreducible fiber, mirroring
   reducible versus irreducible polynomial rows.

5. **Code72 support lift.**  Encode `[72,36,16]` minimum-word support incidence
   as a hidden lift problem over the `78/90` address.  The goal is not just to
   match `lambda_5=78`, but to decide whether the binary support lift can exist.

## Assumption Challenge

Do not assume the tournament vertices must be runners, arcs, coefficients, or
points.  This session explicitly considered polynomial rows, factor coefficients,
prime tokens, LRC runners, residues, Q27 divisors, unit-distance points, unit
directions, Moser ears, triangular endpoints, code supports, Fourier modes, and
proof obligations.  The chosen quotient uses carriers/proof obligations because
it preserves the predicate "does a hidden lift exist?" across the domains.

What it destroys: domain-specific geometry.  Therefore the analogy is only a
search and proof-interface generator until each domain reattaches its local side
channels.
