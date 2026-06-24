# HYP-2999: Pascal-slope additive-basis Farey packet schema

**Status:** SYNTHESIS / proof-carrier guardrail, not a proof.

**Session:** codex-2026-06-24-S169

**Claim.** The user-specified Fibonacci decomposition

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

is the `d=2` Pascal-slope row system

```text
row_2(n,k) = C(n-1-k,k),       sum_k row_2(n,k) = F_n
```

and should be treated as the additive-fiber packet schema joining the repo's
Goldbach, ternary Goldbach, Fermat polygonal, Zeckendorf, Fibonacci, and
mutated-Farey threads.  This is a companion to HYP-2998: HYP-2998 gives the
computed Farey-Fibonacci additive-basis carrier, while HYP-2999 names the
extra packet fields and quotient guardrails the LRC classifier should retain.

## Carrier

Use a representation hypergraph:

```text
vertices = permitted atoms
hyperedges over N = atom multisets / packets summing to N
packet coordinates = entropy, local/residue rank, carry width, exact scale
```

This puts the old additive-basis regimes on one axis.

- **Goldbach:** abundant two-prime fibers; proof currency is smoothing,
  local singular series, and lower bounds for many representations.
- **Ternary Goldbach:** one more summand raises fiber dimension and gives
  enough smoothing room to defeat local obstruction.
- **Fermat polygonal:** bounded arity plus residue absorption; every integer
  is covered by a bounded number of `s`-gonal atoms.
- **Zeckendorf/Fibonacci:** the opposite extreme from Goldbach abundance:
  a sparse atom set plus a confluent carry rule gives a unique normal form
  with no adjacent Fibonacci digits.

The Pascal-slope rows are the finite bridge between bounded additive covers
and Zeckendorf normal forms.  The entry `C(n-1-k,k)` counts binary strings /
tilings / independent-set packets with `k` long steps and no adjacent selected
positions; summing over `k` gives the Fibonacci recurrence.

## Farey-operator relation

For an exact Farey/LRC packet `p/q`, do not replace the root coordinate by a
derived scalar.  Keep:

```text
root scale:        p/q, q, and Farey excess such as e=14p-q
additive lane:     p+q
product lane:      p*q
ordered power lane q^p and p^q
```

Interpretation:

- `p+q` is the least destructive additive ledger.  It is the lane closest to
  Pascal-slope, Goldbach/Fermat summand accounting, and the old `n+2` recursion
  shadow.
- `p*q` is the multiplicative/coimage ledger.  It remembers factor fibers and
  the complete bipartite incidence graph `K_{p,q}`, with `p*q=|E(K_{p,q})|`.
  Perfect-number and perfect-product analogies belong here as product-fiber
  guardrails; the product value is useful only with its factor and Farey labels.
- `q^p` and `p^q` are ordered magnitude stress tests.  They are useful for
  detecting false scalar orderings but are not proof coordinates unless their
  root packet and quotient kernel are retained.

The LRC14 unit-excess chain `p/(14p-1)` makes this especially sharp:
`q=14p-1`, `p+q=15p-1`, and `p*q=14p^2-p`.  The additive lane advances like a
small recursion ledger; the product lane crosses the `K33` incidence wall at
`p>=3`; the power lanes mostly amplify scale and order.

## Quotient guardrail

This carrier is allowed to forget a coordinate only if the LRC predicate is
constant on the forgotten fibers, the coordinate is reconstructible from the
retained packet, a dual certificate annihilates the loss, or the loss is routed
to a named residual sector.

In particular:

- A Fibonacci next-term shadow is not enough; keep the row vector and carry
  rule.
- A Goldbach count is not enough; keep local rank, residue obstruction, and
  smoothing regime.
- A polygonal cover is not enough; keep the side count and bounded arity.
- A Zeckendorf normal form is not enough; keep the carry graph that proves
  confluence.
- A Farey product or power is not enough; keep exact `p/q`, `q`, excess, and
  branch labels.

## Tournament Analysis

Vertices should be proof carriers, not runners:

```text
pascal_slope_row
goldbach_pair_fiber
ternary_goldbach_hyperedge
fermat_polygonal_cover
zeckendorf_carry_normal_form
farey_sum_lane
farey_product_Kpq_lane
farey_power_stress_lane
```

Pairwise observable: which carrier preserves the exact LRC packet predicate
while exposing additive-fiber structure.

Proof-gauge order, before any row-specific computation:

```text
exact p/q packet
> pascal_slope row vector / carry width
> farey sum lane
> zeckendorf normal form
> fermat polygonal bounded cover
> ternary goldbach smoothing
> goldbach abundance
> farey product Kpq lane
> power stress lane
```

This is not intended as a scalar theorem.  Its use is to decide what fields a
labelled packet schema must carry before comparing additive, multiplicative,
and power shadows.

## LRC14 use

Extend HYP-2963-style packet records with:

```text
additive_basis_regime
representation_entropy
local_residue_rank
carry_width
pascal_slope_row_id
farey_operator_lane in {root,sum,product,q^p,p^q}
product_factor_fiber / Kpq_state
```

Expected value: if a hard row survives Fejer/Ramanujan/Haar packet tests, this
schema can say whether it is a smoothing-rich additive fiber, a bounded-cover
residue absorber, a unique-carry normal form, or a product-incidence residual.

## Open hooks

1. Build a small script that emits the `d=2` Pascal-slope rows, Zeckendorf
   carry graph, and Farey operator lanes for named LRC14 rows.
2. Add `additive_basis_regime`, `carry_width`, and `local_residue_rank` to the
   packet classifier output.
3. Test whether AP/GW equality atoms are precisely low-entropy carry-normal
   additive packets after exact `M` and endpoint owners are retained.
4. Test whether K33 packets are exactly the first product-incidence residuals
   that cannot be read from the additive lane alone.

## Anchors

S501 Goldbach/polygonal/Zeckendorf additive-basis reflections; HYP-2523
Fibonacci/Zeckendorf bridge; HYP-2931 mutated Farey operator carriers;
HYP-2932/HYP-2934 Farey product/K33 incidence; HYP-2940 operator sequence
tournament; HYP-2963 labelled packet classifier; HYP-2998 Farey-Fibonacci
carrier; HYP-2990 quotient guardrail.
