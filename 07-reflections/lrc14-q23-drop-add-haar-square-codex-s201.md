# LRC14 q=23 Drop/Add Haar Square - S201

The HYP-3032 analytic pass left a very specific itch: two squarefree `q=23`
rows have the same analytic clocks but different routes.  The temptation is to
ask for a sharper arithmetic scalar.  The better move was to make the missing
coordinate visible.

The two rows are opposite corners of a drop/add square:

```text
drop(10,13) add(20,26)    drop(10,13) add(16,24)
drop(8,12)  add(20,26)    drop(8,12)  add(16,24)
```

The diagonal corners are the HYP-3032 residual pair and both have `M=2/23`.
The off-diagonal corners immediately open as `M=1/10` and `M=1/8` q-witness
rows.  So the q=23 obstruction is not floating in number-theory space; it is
the diagonal doubling match in a fixed-margin/Haar square.

The exact zeta calculation makes that concrete:

```text
zeta(M) = -47/920
```

Safe measure, bar count, longest bar, midpoint slack, boundary count, and
zero-sum active-pair count also have nonzero zeta.  But raw magnitude height
and magnitude delta have zeta zero, so a magnitude-only repair misses this
particular square.  Nice little warning label.

Inside the diagonal, exact `M` is still not enough: petal and covering remain
mixed.  The endpoint-owner strip splits them:

```text
petal external owners: 12:26x6,6:20x4
cover external owners: 2:16x6
```

That is the useful new angle: the repair class is not just
`nested_refinement`; it is

```text
nested_refinement_to_q23_diagonal_then_owner_strip
```

I would hand the next agent this local lemma target:

```text
For any double-pair drop/add square, diagonal doubling either opens the
off-diagonal corners to q-witnesses, descends through a family q-diagonal, or
exposes endpoint-owner strip data routing petal, covering, K33, or F7 debt.
```

Assumption challenge: I considered runners, drop pairs, add pairs, gaps,
endpoints, fixed sections, residues, Fourier modes, Haar tiles, cover arcs,
matroid circuits, and proof obligations.  The vertices that worked were proof
teeth: raw analytic shadow, row/column drop-add shadow, diagonal match,
exact-M zeta, safe body, endpoint owner strip, and route label.  The LRC
predicate preserved is open/boundary plus theorem-route schedulability, not
runner identity.
