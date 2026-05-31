---
id: HYP-1890
status: OPEN
source: codex-2026-05-31-S420
related:
  - THM-357
  - THM-360
  - THM-366
  - THM-367
  - HYP-1834
  - HYP-1835
  - HYP-1858
  - HYP-1859
  - HYP-1880
---

# HYP-1890: LRC counterexamples are infeasible two-mode integer programs

## Statement

A reduced Lonely Runner counterexample at denominator `n` should be expressible
as an infeasible integer program after splitting speed columns by

```text
v = 2^h * odd_core.
```

The horizontal mode is the odd-core chain `odd_core -> odd_core+2`; the
vertical mode is the doubling chain `h -> h+1`.  In this coordinate system,
any counterexample must choose exactly `n-1` primitive columns satisfying two
families of row constraints:

```text
small denominator rows:   sum_{m|v} x_v >= 1,  2 <= m <= n
endpoint rows:            sum_{p: ||p*t|| < 1/n} x_p >= 1
```

The conjectural proof object is a dual certificate on these rows.  Unit and
small-denominator rows force gate columns; gate columns export endpoint rows to
deeper denominator layers.  The dual should show that the exported vertical
endpoint debt cannot be paid together with the horizontal sieve invoices using
only `n-1` primitive columns.

## Evidence

`lrc_integer_programming_modes_s420.py` builds exact finite row/column ledgers.

At the unit layer, the IP row is sharp:

```text
n=14 unit layer, candidates 1..28: exact cover 1, column (14)
n=16 unit layer, candidates 1..32: exact cover 1, column (16)
```

This is the THM-360 gate branch in IP form: a unit endpoint row can only be
paid by an `n`-multiple.

For `n=16`, the local vertical invoice is also sharp.  Covering all `32`
endpoints owned by a `16`-gate with lower columns `1..15` has exact minimum
`9`:

```text
(8, 1, 3, 5, 7, 9, 11, 13, 15)
```

The forced columns are exactly the THM-367 private-endpoint residues
`(1,3,5,7,8,9,11,13,15)`.  This turns the dyadic debt law into a literal
set-cover lower bound.

For `n=14`, the horizontal sieve can be paid without producing a cover.  The
S380 gate ladder satisfies every small-denominator row, but its own endpoint
rows are still not all covered:

```text
targets=2298, candidates=13, union=2130, uncovered=168
uncovered prime depths:
  {2:+1,7:+1}:120
  {2:+3,7:+1}:48
```

Every S380 ladder column is forced by some private endpoint row, yet the row
system still has uncovered debt.  This is the composite-denominator anomaly in
integer-program form: satisfying the divisibility invoices exports debt into
the product of the `2`- and `7`-adic depth directions.

## Interpretation

The old natural-number operation thread said:

```text
addition shadow       = dense order / x+2 horizontal completion
multiplication shadow = sparse divisibility / x*2 vertical skeleton
```

The LRC IP lens turns that slogan into a certificate search.  A candidate speed
set is not merely a point in gap space; it is a selected subset of columns in a
two-mode incidence matrix.  Initial segments are horizontal equality spines.
Gate repairs are vertical divisibility columns.  Counterexamples would require
a perfectly balanced row cover with no positive gap, no unprotected endpoint,
and no peelable private row.

## Predictions

1. An `n=16` proof should be a dual row-weight certificate combining the unit
   gate row with the nine-column lower-cover invoice for `16`-gate endpoints.
2. An `n=14` proof should use a product-depth dual on the `2`- and `7`-adic
   endpoint rows.  The S380 ladder already shows the right exported rows:
   `{2:+1,7:+1}` and `{2:+3,7:+1}`.
3. Direct disproof searches should optimize over row patterns first, then ask
   whether integer speeds realize them.  Searching speed sets first hides the
   row infeasibility in a large nonlinear space.
4. The next useful computation is a branch-and-bound whose node ordering is
   row-dual pressure: unit rows, fragile one-payer small-denominator rows,
   private endpoint rows, then exported prime-depth rows.

## Sources

- `04-computation/lrc_integer_programming_modes_s420.py`
- `05-knowledge/results/lrc_integer_programming_modes_s420.out`
- `07-reflections/lrc-integer-programming-modes-s420.md`
- THM-357
- THM-360
- THM-366
- THM-367
- HYP-1834
- HYP-1835
- HYP-1858
- HYP-1859
- HYP-1880
