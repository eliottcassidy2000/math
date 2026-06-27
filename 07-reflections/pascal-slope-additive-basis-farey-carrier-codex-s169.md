# Pascal-Slope Additive Bases And Farey Operator Lanes

codex-2026-06-24-S169

## The Fibonacci row is a carrier, not a decoration

The row pattern from the prompt

```text
1,
1,
1+1,
1+2,
1+3+1,
1+4+3,
1+5+6+1,
...
```

is exactly the `d=2` Pascal-slope row system:

```text
row_2(n,k) = C(n-1-k,k),     sum_k row_2(n,k) = F_n.
```

The entries count packets with `k` long steps and no adjacent chosen positions.
Equivalently they count independent-set / tiling states before the row is
summed into a single Fibonacci number.  That means the old sequence picture
already had a quotient lesson inside it: the scalar Fibonacci term forgets
which row fiber, carry width, and packet support produced it.

Zeckendorf is the normal-form endpoint of that same row geometry.  The sparse
Fibonacci atoms become useful because the carry rule is confluent:

```text
F_i + F_{i+1} -> F_{i+2}
```

with the no-adjacent-digit condition as the canonical section.  In quotient
language, Zeckendorf is not "Fibonacci numbers are special"; it is "the
representation hypergraph has a terminating carry graph and a unique normal
form."

## Additive basis regimes

The S501 additive-basis pass becomes cleaner if arranged by representation
fiber geometry.

| Regime | Atom set | Fiber shape | Proof currency |
|---|---|---|---|
| Goldbach | primes | many two-atom fibers | local density, singular series, smoothing |
| Ternary Goldbach | primes | three-atom hyperedges | one extra smoothing dimension |
| Fermat polygonal | `s`-gonal numbers | bounded-arity residue absorber | finite local residue checks plus descent |
| Zeckendorf | Fibonacci numbers | unique normal-form fiber | carry confluence / no adjacent digits |
| Pascal-slope rows | tiling packets | finite row fibers before summation | carry width, row entropy, independent-set support |

These are not competing metaphors.  They are different points on a common
representation-hypergraph axis: abundance, smoothing dimension, bounded cover,
and normal form.

For LRC14, this suggests adding three additive-basis fields to labelled
packets:

```text
representation_entropy
local_residue_rank
carry_width
```

Goldbach wants high entropy and analytic lower bounds.  Fermat polygonal wants
small bounded rank and residue absorption.  Zeckendorf wants entropy collapse
through a carry normal form.  Pascal-slope rows are the finite audit surface
between these extremes.

## Farey operators as lanes

The Farey side should be read as a set of labelled lanes over the root packet
`p/q`, not as replacements for it.

```text
root:        p/q, q, exact threshold excess
sum lane:    p+q
product:     p*q
powers:      q^p, p^q
```

The sum lane is the additive-basis lane.  It lines up with Pascal-slope row
bookkeeping, the old `n+2` recursion shadow, and the minimum-loss transformation
on the LRC14 unit-excess chain:

```text
p/(14p-1):       q=14p-1,       p+q=15p-1.
```

The product lane is the incidence/coimage lane:

```text
p*q = |E(K_{p,q})|.
```

It is essential for the K33 story, but it is not the same kind of information
as the additive row.  For the unit-excess chain, `p*q=14p^2-p`, and the first
`p>=3` branch is exactly where K33-style state-lift debt becomes visible.
This is also where the earlier perfect-number/perfect-product theme belongs:
perfectness can mark a product fiber, but the proof coordinate is the labelled
factor/Farey/incidence packet, not the product value by itself.

The powers `q^p` and `p^q` are ordered scale stressors.  They help expose false
proofs that rank packets by magnitude, but they are not allowed to forget the
root Farey coordinate.

## The useful packet schema

The merged packet should carry:

```text
exact_M
q
Farey_excess
operator_lane
Pascal_slope_d
Pascal_row_vector
additive_basis_regime
representation_entropy
local_residue_rank
carry_width
Kpq_factor_fiber
power_shadow_flags
endpoint_owner_labels
Fejer/Ramanujan/Haar certificate labels
```

The quotient rule is the same as HYP-2990: a lane may forget only what is
fiber-constant, reconstructible, annihilated by a dual certificate, or routed
to a named residual sector.

## Tournament Analysis

Vertices:

```text
exact_pq_packet
pascal_slope_row
farey_sum_lane
zeckendorf_normal_form
fermat_polygonal_cover
ternary_goldbach_smoothing
goldbach_pair_fiber
farey_product_Kpq_lane
farey_power_stress_lane
```

Pairwise observable: carrier `A` beats carrier `B` when `A` preserves more of
the LRC packet predicate while exposing the additive-basis information needed
to classify the row family.

A conservative tie path is:

```text
exact_pq_packet
> pascal_slope_row
> farey_sum_lane
> zeckendorf_normal_form
> fermat_polygonal_cover
> ternary_goldbach_smoothing
> goldbach_pair_fiber
> farey_product_Kpq_lane
> farey_power_stress_lane
```

This path is not a theorem and not a scalar ranking.  It is a dependency order
for packet construction: exact scale first, additive row/carry information
second, product incidence only after the additive quotient has declared what it
forgets.

## LRC14 conjectural readout

The most promising concrete target is:

```text
AP/GW equality atoms = low-entropy carry-normal additive packets
unit-petal/C27 rows   = additive lane with marked second-gap residue transfer
K33 rows              = first product-incidence residuals not decided by p+q
covering rows         = bounded-cover residue absorbers with positive Haar front
F7                    = named residual where carry/product/endpoint cocycles do not close
```

This is still proof-interface work, not a proof.  Its value is that it gives
future scripts a schema for deciding whether a hard LRC14 row is being handled
by abundance, bounded cover, normal form, product incidence, or residual
cocycle debt.

## Next computations

1. Emit Pascal-slope row vectors and Zeckendorf carry graphs beside named
   LRC14 packet rows.
2. Add `additive_basis_regime`, `representation_entropy`, `local_residue_rank`,
   and `carry_width` to the HYP-2963 classifier output.
3. Compare Farey `p+q`, `p*q`, `q^p`, and `p^q` lanes only after grouping by
   exact `M`, endpoint owners, and packet route.
4. Test whether K33 rows are precisely the first product-lane fibers with
   nontrivial incidence debt after the sum lane has been fixed.
