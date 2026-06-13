---
id: HYP-1881
status: OPEN
source: codex-2026-05-31-S411
related:
  - THM-371
  - HYP-1833
  - HYP-1855
  - HYP-1859
  - HYP-1866
  - HYP-1867
  - HYP-1868
  - HYP-1880
  - HYP-1890
  - HYP-1891
  - HYP-1892
  - HYP-1893
  - HYP-1905
---

# HYP-1881: LRC denominator recursion splits into +2 columns and *2 rows

## Statement

Let an LRC threshold denominator be written as

```text
n = 2^r * b,  b odd.
```

The denominator recursion should be handled as a two-mode grid, matching the
old tournament-size recursion:

```text
top-row mode:  odd b -> b+2        unit/product-sum changes
row mode:      n -> 2n             quotient endpoint-debt transfer
```

The conjectural proof strategy is a two-dimensional induction:

```text
Prove +2 movement on the odd top row using unit skeleton and product-sum data.
Prove *2 movement down a fixed odd-core column using row-parent endpoint debt.
```

## Evidence

`lrc_column_row_modes_s411.py` verifies that the initial-segment LRC boundary
debt is exactly the unit skeleton:

```text
unprotected endpoints for {1,...,n-1} = phi(n),  3 <= n <= 40.
```

Along a fixed odd core `b`, the unit skeleton obeys:

```text
phi(2b) = phi(b)
phi(2^(r+1)b) = 2 phi(2^r b),  r >= 1.
```

So the first doubling `b -> 2b` is a seam: it creates nonunit quotient room but
does not create new unit witnesses.  Later doublings are honest doublings.  This
is the LRC analogue of the tournament column-family seam where an unmatched
vertex first gains a twin.

S453/THM-371 sharpens this into a residue-level statement: for odd `b`,
`U(2b) -> U(b)` is a bijection, with inverse given by the unique odd lift among
`a` and `a+b`.  Thus the first seam copies the unit skeleton into the odd
parity sheet rather than doubling it.  The newly available columns are
nonunit/quotient room, so using them should be charged as endpoint-debt export.

The same script computes largest-proper-divisor ladders.  For every even
denominator in the table,

```text
lpd(n) = n/2.
```

Thus the first visible quotient ladder is literally the row parent in the
2-adic grid.  Sample exact rows:

```text
n=14 = 2*7:   lpd=7,  gap/th=5/924,   debt=84,  product=5/11
n=16 = 16*1:  lpd=8,  gap/th=1/132,   debt=140, product=35/33
n=18 = 2*9:   lpd=9,  gap/th=1/176,   debt=176, product=1
n=20 = 4*5:   lpd=10, gap/th=1/325,   debt=180, product=36/65
n=24 = 8*3:   lpd=12, gap/th=11/4080, debt=264, product=121/170
```

Odd composites instead use an odd-divisor channel inside the top row:

```text
n=15: lpd=5, gap/th=7/660, debt=66, product=7/10
n=21: lpd=7, gap/th=3/455, debt=152, product=456/455
```

This separates the two LRC repair modes: even denominators descend to a row
parent, while odd composites pay a divisor debt within the top-row arithmetic.

## Interpretation

The old natural-number picture says every number is either:

```text
an odd top-row member reached by +2,
or an even row descendant reached by *2.
```

For tournaments, these are two recursion modes:

```text
n -> n+2         top-row staircase / Mode B motion
n -> 2n          blowup / row motion
```

For LRC, the same split says:

```text
odd +2 motion    changes the unit skeleton and product-sum factor packing
even *2 motion   keeps odd core and exports endpoint debt to n/2
```

This explains why `n=14`, `n=16`, and `n=18` are different laboratories:

```text
n=14: first prime-core seam, 7 -> 14
n=16: pure dyadic row lab, odd core 1
n=18: first square-odd-core seam, 9 -> 18
```

The S398 gap-debt product law and S399 Bruhat-Tits flow law are row-mode
phenomena, HYP-1868/HYP-1880 refine that row mode as Bruhat-Tits
frontier/descent, HYP-1890/HYP-1891 recast the same split as IP rows and adic
column transfer, and HYP-1892/HYP-1893 sharpen the seam between adelic
gap-debt and coarse/fine denominator regimes.  The S377/S378 product-sum target
work is top-row mode.  HYP-1881 asserts that a complete recursive proof needs
both.

## Predictions

1. Initial-segment boundary-only proofs should be formalized as a unit-skeleton
   theorem: the exposed endpoints are exactly the unit residues modulo `n`.
2. First-row seams `b -> 2b` should be disproportionately important: they open
   the divisor channel `b` without increasing `phi`.
3. Among one-divisor quotient ladders, `lpd(n)` should minimize gap pressure
   except where a genuine multifactor product-sum packing beats the binary
   divisor layer.
4. A future LRC feature extractor should include `(r,b)`, `phi(n)`, seam
   defect `2 phi(n/2)-phi(n)`, `lpd(n)`, row-parent debt, and product-sum core.
5. The `n=16` proof should be the pure row-mode proof; the `n=14` proof should
   be a first-seam proof; the `n=18` proof should test how odd square torsion
   rides on the first even row.

## Sources

- `04-computation/lrc_column_row_modes_s411.py`
- `05-knowledge/results/lrc_column_row_modes_s411.out`
- `07-reflections/lrc-column-row-modes-s411.md`
- `04-computation/lrc_tournament_first_doubling_seam_s453.py`
- `05-knowledge/results/lrc_tournament_first_doubling_seam_s453.out`
- `07-reflections/lrc-tournament-first-doubling-seam-s453.md`
