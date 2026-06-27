# LRC14 q=23 Drop/Add Haar Square

HYP-3038/S201 attacks the HYP-3032 squarefree `q=23` petal/covering residual
pair by putting it into the HYP-3031 Haar-tile repair dictionary.

The two residual rows are diagonal corners of the square:

```text
drop(10,13) add(20,26)    drop(10,13) add(16,24)
drop(8,12)  add(20,26)    drop(8,12)  add(16,24)
```

The diagonal corners have `M=2/23`; the off-diagonal corners open as `M=1/10`
and `M=1/8` q-witness rows.  Exact-M zeta is nonzero:

```text
zeta(M) = -47/920
```

But exact `M` still mixes the two diagonal routes.  Endpoint-owner strips split
them:

```text
petal external owners: 12:26x6,6:20x4
cover external owners: 2:16x6
```

Repair class:

```text
nested_refinement_to_q23_diagonal_then_owner_strip
```

Next pull: test the same double-pair drop/add square lemma across more
diagonal families and see whether every square either opens off diagonal,
descends by family q-diagonal, or exposes endpoint-owner strip data routing
petal, covering, K33, or F7/THM-572 debt.
