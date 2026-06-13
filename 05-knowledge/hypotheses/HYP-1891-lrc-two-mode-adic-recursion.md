---
id: HYP-1891
status: OPEN
source: codex-2026-05-31-S421
related:
  - THM-371
  - HYP-1831
  - HYP-1838
  - HYP-1854
  - HYP-1859
  - HYP-1866
  - HYP-1867
  - HYP-1868
  - HYP-1880
  - HYP-1890
  - HYP-1905
---

# HYP-1891: LRC splits into row-mode dyadic flow and column-mode odd-payload transfer

## Statement

Use the repo's natural-number grid

```text
n = 2^r * odd.
```

The odd roots form an `x+2` chain, and each odd root has an `x*2` vertical
chain.  Tournament history reads these as two modes:

```text
row mode:       n -> 2n        tournament blowup;
column mode:    odd -> odd+2   add-two / Mode-B top-row recursion.
```

The LRC proof should split the same way.

Row mode is p-adic: external denominator doubling and internal endpoint-depth
translation both move down a Bruhat-Tits cusp.

Column mode is odd-payload transfer: changing the odd root changes which
product-building trees carry endpoint debt.

## Evidence

For initial segments, boundary endpoint debt is exactly `phi(n)`.  Along a
vertical chain `n=2^r*m` with `m` odd,

```text
phi(m) -> phi(2m) -> phi(4m) -> ...
```

has the seam pattern

```text
phi(2m) = phi(m),
phi(2^(r+1)m) = 2 phi(2^r m) for r>=1.
```

Thus the first doubling keeps the boundary-debt count fixed, while all later
row steps double it.  This is the LRC analogue of the tournament row-seam in
the column-family notes: the first blowup is special.

S453 identifies the fixed-count mechanism: `U(2m) -> U(m)` is not a two-sheet
unit cover but a bijection, with inverse given by the unique odd lift among
`a` and `a+m`.  On the tournament side, the same seam is the first
twin-gaining event: the unmatched odd vertex gains a partner, so matching
pairs acquire the extra `+1` while unit residues refuse to double.

The first-even children of the odd `x+2` chain give the current LRC frontier:

```text
odd 7 -> n=14: gap/th=5/924, debt=84,  product=5/11
odd 9 -> n=18: gap/th=1/176, debt=176, product=1
```

Meanwhile `n=16` is not a column-neighbor of these rows.  It is the pure
vertical column over odd root `1`:

```text
1 -> 2 -> 4 -> 8 -> 16.
```

That explains why the `n=16` proof naturally wants dyadic Bruhat-Tits
frontier/harmonic-flow tools, while `n=14` and `n=18` want odd-payload
transfer tools.

## Proof Route

Separate the LRC proof technology into two lemmas.

```text
Row lemma:
  In a fixed odd-root column, repaired endpoint debt can move deeper in the
  p-adic cusp, but normalized Bruhat-Tits frontier mass remains positive.

Column lemma:
  Moving odd root m -> m+2 changes the odd payload tree.  Any apparent
  reduction in real gap must be paid by a transfer of endpoint debt into the
  new product-building coordinates.
```

Then the remaining proof architecture is a commutation problem.  A
counterexample would need to traverse both modes and make every real and
p-adic frontier vanish.  HYP-1866/HYP-1868 say row-mode vanishing is blocked by
frontier mass.  The missing piece is the column-mode transfer matrix on odd
payloads.

## Sources

- `04-computation/lrc_adic_column_modes_s421.py`
- `05-knowledge/results/lrc_adic_column_modes_s421.out`
- `07-reflections/lrc-two-mode-adic-recursion-s421.md`
- `04-computation/lrc_tournament_first_doubling_seam_s453.py`
- `05-knowledge/results/lrc_tournament_first_doubling_seam_s453.out`
- `07-reflections/lrc-tournament-first-doubling-seam-s453.md`
- `07-reflections/adic-column-families.md`
- HYP-1866.
- HYP-1867.
- HYP-1868.
