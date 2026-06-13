---
id: HYP-1905
status: OPEN
source: codex-2026-05-31-S453
related:
  - THM-371
  - THM-378
  - HYP-1881
  - HYP-1890
  - HYP-1891
  - HYP-1904
  - HYP-1910
  - HYP-1920
---

# HYP-1905: LRC first-even denominators are parity-seam export problems

## Statement

The critical component in the old odd `x+2` / even `x*2` grid is the first
doubling seam

```text
odd m -> 2m.
```

This step should not be treated as an ordinary doubling.  THM-371 says the
unit skeleton satisfies

```text
U(2m) ~= U(m),  phi(2m)=phi(m),
```

by choosing the unique odd lift of each unit, while later row steps double the
unit skeleton.  Thus the first-even LRC denominator opens new nonunit quotient
room without creating new unit witnesses.

The tournament analogue is the first twin-gaining event.  An odd tournament has
one unmatched vertex in a maximum matching; after doubling, that vertex gains a
twin and creates the extra pair.  THM-378 refines the SC blowup into a
double-round-robin voltage lift: every old edge becomes a `2 x 2` block where
the winner and loser take complementary perfect matchings.  In the zero-voltage
SC blowup, every base vertex becomes a strong/weak twin pair and all score
variation collapses to the universal near-regular split.

The proposed LRC proof move is:

```text
unit endpoints persist through odd lift;
repairs must use the new nonunit quotient room;
using that room exports endpoint debt to row-parent/product-depth layers.
```

Equivalently, read backward from the first-even child `2m` to its odd parent
`m`: the projection collapses two parity sheets, but the unit layer records
only the odd sheet.  Any proof that reduces `2m -> m` must carry a separate
debt ledger for the nonunit sheet it has just collapsed.

## Evidence

`lrc_tournament_first_doubling_seam_s453.py` prints the exact arithmetic seam:

```text
pairs(2m) - 2*pairs(m) = 1
2*phi(m) - phi(2m) = phi(m)
2*phi(2m) - phi(4m) = 0
```

For the relevant first-even rows:

```text
n=14 = 2*7:  lpd=7,  gap/th=5/924, debt=84,  product=5/11
n=18 = 2*9:  lpd=9,  gap/th=1/176, debt=176, product=1
n=22 = 2*11: lpd=11, gap/th=2/561, debt=260, product=520/561
```

These rows are not simply weaker versions of odd-root rows.  Their most visible
quotient ladder is exactly the row parent `m`.

The same script verifies the tournament-side collapse in small cases:

```text
SC blowup score for n=3 base: 2^3 3^3
SC blowup score for n=5 base: 4^5 5^5
```

independent of the base tournament.  Across all labelled `n=5` tournaments,
the SC blowup concentrates Hamiltonian path counts in `14937..15565`, while
lex blowup ranges across `1..15565`.

S501c sharpens this: signed double-round-robin doublings have sheet-flip gauge
classes classified by `binom(n-1,2)` triangle parity bits, matching the
fixed-base tiling cube dimension.

## Predictions

1. A proof for `n=14` should begin by treating `U(14)` as the odd-lift copy of
   `U(7)`, then charging every use of the new quotient room to exported
   endpoint debt.
2. `n=18` should be the next first-even test: it has the same seam mechanism as
   `n=14`, but the odd payload is square-rooted through `9`.
3. `n=16` should remain separate: it is not a first-even child but a pure
   dyadic row-lab over odd root `1`.
4. The right tournament analogy is SC blowup, not lex blowup.  Lex blowup
   magnifies hierarchy; SC blowup creates twin balance while hiding the old
   tournament in lane/cross-lane data.  LRC first-even repairs should likewise
   hide old unit data while exporting debt to a labelled second sheet.
5. Feature extractors should record a first-seam flag, odd-lift unit labels,
   row-parent ladder debt, and whether a candidate repair spends nonunit room
   or merely relabels the inherited unit skeleton.

## Sources

- `01-canon/theorems/THM-371-first-doubling-unit-pair-seam.md`
- `01-canon/theorems/THM-378-double-round-robin-vertex-doubling.md`
- `04-computation/lrc_tournament_first_doubling_seam_s453.py`
- `04-computation/double_round_robin_blowup_s501.py`
- `05-knowledge/results/lrc_tournament_first_doubling_seam_s453.out`
- `05-knowledge/results/double_round_robin_blowup_s501.out`
- `07-reflections/lrc-tournament-first-doubling-seam-s453.md`
- `07-reflections/double-round-robin-vertex-doubling-s501.md`
- `07-reflections/adic-column-families.md`
- `07-reflections/sc-blowup-and-twin-gaining.md`
- HYP-1881
- HYP-1891
