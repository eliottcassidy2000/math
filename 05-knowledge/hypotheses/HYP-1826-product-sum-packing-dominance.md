# HYP-1826: Product-sum minima are governed by factor-packing dominance

**Status:** EXPLORATORY; supported by exact enumeration through arity `120`.

## Claim

For product-sum equations

```text
x_1 + ... + x_k = x_1 * ... * x_k,
```

minimal endpoints are usually not controlled by the two-nonunit divisor layer.
They are controlled by a packing principle for small nonunit factors.

For a nonunit seed `F`, define

```text
P(F) = prod(F)
S(F) = sum(F)
D(F) = P(F)-S(F)
k(F) = |F| + D(F).
```

The seed gives a product-sum witness of arity `k(F)` by adjoining `D(F)` ones.
Thus the minimal product-sum endpoint at arity `k` is the minimum `P(F)` among
seeds satisfying

```text
P(F) - S(F) + |F| = k.
```

The conjectural dominance principle is: replacing one large binary factor by
several small factors often wins because the additive discount `S(F)-|F|`
increases faster than the endpoint cost.

## Evidence

Lean formalization in `NaturalOperationDigraphs.lean` proves the base layer:

```text
r + (a+1) + (b+1) = (a+1)(b+1)  <->  a*b = r+1.
```

It also proves the universal binary witness: `r` ones plus `2` and `r+2`
solve `sum=product`, giving `m(k)<=2k`.

The S366 computation enumerates all potentially minimal witnesses through
`k<=120`. Strict multi-factor wins over the best binary witness occur in
`96/119` arities.

Early examples:

```text
k=5:  (2,2,2),      P=8,  beats binary P=9
k=8:  (2,2,3),      P=12, beats binary P=16
k=12: (2,2,2,2),    P=16, beats binary P=24
k=21: (3,3,3),      P=27, beats binary P=30
k=27: (2,2,2,2,2),  P=32, beats binary P=42
```

Through `k<=120`, the minimal seed-length distribution is:

```text
length 2: 23
length 3: 52
length 4: 30
length 5: 11
length 6: 3
```

## Predictions

1. There is a finite dominance rule that prunes large factors in favor of
   bounded small-factor replacements.
2. The record jumps and plateaus of the minimal endpoint sequence can be
   explained by changes in the optimal multiplicative partition type.
3. A two-colored operation complex should interpret multi-factor wins as
   compatible packings of operation cells, analogous to OCF packings.

## See Also

- `04-computation/lean/TournamentH7/TournamentH7/NaturalOperationDigraphs.lean`
- `04-computation/natural_operation_modes_s366.py`
- `05-knowledge/results/lean_natural_operation_digraphs_s366.out`
- `05-knowledge/results/natural_operation_modes_s366.out`
- `07-reflections/formal-extension-investigation-operation-modes-s366.md`
