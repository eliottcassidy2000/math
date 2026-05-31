# HYP-1820: Product-sum equations are operation-shadow critical pairs

**Status:** EXPLORATORY; exact for the finite shadow facts and small enumerations.

## Claim

The natural-number operation graph should be modeled as a two-colored gate
hypergraph, not as an ordinary directed graph.

- The additive gate `{x,y}->x+y` has simple one-input shadow `x->z` iff `x<z`,
  so its shadow is the transitive tournament.
- The multiplicative gate `{x,y}->xy` has simple one-input shadow `x->z` iff
  `x|z`, so its shadow is the divisibility DAG.
- Equations of the form `x_1+...+x_k=x_1...x_k` are precisely critical pairs
  where an additive fold and a multiplicative fold reach the same endpoint,
  after allowing unit padding.

For a seed `F=(a_1,...,a_m)` with all `a_i>1`, define

```text
P(F)=prod a_i,   S(F)=sum a_i,   D(F)=P(F)-S(F).
```

If `D(F)>=0`, then `D(F)` copies of `1` repair the additive fold:

```text
D(F)*1 + S(F) = P(F),
```

giving a product-sum solution of arity `k=m+D(F)`.

## Evidence

`natural_operation_digraphs_s365.py` verifies the finite one-shadow identities
and records the sparsity contrast. On `[100]`, the additive shadow has all
`4950` transitive comparisons, while the multiplicative shadow has only `382`
including unit edges and `283` excluding unit edges.

`NaturalOperationDigraphs.lean` now formalizes the shadow identities
axiom-free:

```text
AddShadow x z <-> x < z
MulUnitShadow x z <-> x | z        (positive z)
MulShadow x z <-> x | z and x < z  (positive x)
```

The binary gate collision is unique:

```text
x+y=xy  <=>  (x-1)(y-1)=1  <=>  (x,y)=(2,2).
```

The two-factor arity layer has the divisor formula

```text
r+a+b=ab,  k=r+2  <=>  (a-1)(b-1)=k-1.
```

The first arity where a multi-factor seed beats all two-factor seeds is `k=5`:
seed `(2,2,2)` gives product `8`, beating the two-factor minimum `9`.

The S366 extension checks all potentially minimal witnesses through
`k<=120`; strict multi-factor wins over the binary divisor layer occur in
`96/119` arities.

## Predictions

1. The product-sum number sequence is a minimal-endpoint problem in the
   additive transitive completion of the divisibility DAG.
2. Higher-factor improvements over the divisor layer should be explainable by a
   packing principle: several small compatible multiplicative atoms can beat one
   larger two-factor atom.
3. A two-colored operation complex with additive and multiplicative 2-cells
   should expose product-sum witnesses as low-dimensional critical cycles.

## See Also

- `04-computation/natural_operation_digraphs_s365.py`
- `04-computation/lean/TournamentH7/TournamentH7/NaturalOperationDigraphs.lean`
- `04-computation/natural_operation_modes_s366.py`
- `05-knowledge/results/natural_operation_digraphs_s365.out`
- `05-knowledge/results/lean_natural_operation_digraphs_s366.out`
- `05-knowledge/results/natural_operation_modes_s366.out`
- `07-reflections/natural-operation-digraphs-and-product-sum-s365.md`
- `07-reflections/formal-extension-investigation-operation-modes-s366.md`
- `04-computation/natural_numbers_s116i.py`
- `04-computation/understand_naturals_s90ay.py`
- `04-computation/sum_product_tournament.py`
