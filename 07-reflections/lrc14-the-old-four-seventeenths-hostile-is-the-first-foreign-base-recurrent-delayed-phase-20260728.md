# The old `4/17` hostile is the first foreign-base recurrent delayed phase

Status: **FINITE-EXACT SCOUT / STRUCTURAL REFRAME**, not a physical LRC
transition and not an LRC(14) exclusion.

## Trigger

[THM-2693](../01-canon/theorems/THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence.md)
separates a positive three-event lift from a fatal delayed base map

```text
y |-> {13y}.
```

The raw delayed word is

```text
W = D_(13^3) intersect D_14^c intersect D_27^c intersect D_40^c
    intersect D_53^c intersect D_66^c intersect D_(2*13^5)^c.
```

Its `13`-adic orbit has no four-step survivor.  The cheapest hostile probe is
therefore not another longer `13`-word, but a change of operation: for which
integer multiplier `b` does `y |-> {by}` first have an interior fixed point in
`W`?

## Exact verdict

The companion exact scan exhausts every fixed point `k/(b-1)` for
`2 <= b <= 34`.  No multiplier `2,...,17` works.  The first is

```text
b = 18,        y = 4/17 or 13/17.
```

Both points are strict interior points of reflected components of `W`.  At
`4/17` the delayed factors have exact distances

```text
target 13^3: 1/17,
units 14,27,40,53,66: 5/17,6/17,7/17,8/17,8/17,
high speed 2*13^5: 2/17.
```

Thus the target is dangerous and every guard is strictly safe.  The mechanism
is the small residue identity `13^3*4 = -1 (mod 17)`, together with
`18 = 1 (mod 17)`.  Every multiplier `b = 1 (mod 17)` fixes the same phase;
`18` is simply the first positive expanding one.

## The unexpected bridge to THM-789

This is not a new numerical phase.  It is exactly the hostile anchor in
[THM-789](../01-canon/theorems/THM-789-two-sheet-erosion-and-symmetric-return-packet.md):

```text
U_0 = {1,2,3,5,7,8,9,10,11,12},       t_0 = 4/17,
clearance numerators = (4,8,5,3,6,2,2,6,7,3).
```

THM-789 used it to show that a full symmetric return packet can remain trapped
inside one folded diamond.  The delayed-word scan finds the same phase from a
different direction: it is the **first recurrent phase after varying the base
dilation**.  The denominator `17` simultaneously:

1. makes the old ten-speed core uniformly deep;
2. makes `13^3` land one residue from zero;
3. keeps all six delayed guards at least two residues from zero; and
4. supports the smallest expanding congruence stabilizer `18`.

The shared object is therefore not merely the rational number `4/17`.  It is a
strict residue chamber stabilized by a dilation congruence.

## Connection audit

- **Source:** the THM-2693 raw guard word with its base-`13` odometer.
- **Target:** the THM-789 heavy-phase cell and return-packet obstruction.
- **Map:** replace the delayed evolution `{13y}` by `{by}`, then send an
  interior fixed point to its vector of runner clearances.
- **Preserved predicate:** target danger and strict safety of every listed
  guard; at `b = 18`, recurrence preserves the complete pointwise word.
- **Destroyed information:** the map forgets whether `b` comes from an actual
  LRC successor, endpoint transport, carry, or labelled row.
- **Needed sidecar:** a physical derivation of the induced multiplier and the
  accompanying affine phase for every factor, not merely a chosen circle map.
- **Cheapest decisive next test:** enumerate the exact slope/multiplier induced
  by every admissible two-event or heterogeneous-edge handoff already present
  in the THM-2614--2694 atlas, and ask whether any branch is congruent to `1`
  modulo `17` while transporting the full guard word.

## Concept-board comparison

- **Odometer:** `13` is nilpotent on this word at depth four; `18` has an
  honest recurrent component.  This isolates base dynamics, not word size, as
  the differing coordinate.
- **Toothpick self-similarity:** the successful phase is not centered in a
  `13`-adic tooth.  Its denominator is transverse to the tooth base, so a
  purely nested `13`-adic refinement cannot discover its recurrence.
- **H-drift:** an affine drift that only translates while retaining multiplier
  `13` cannot realize this fixed-point mechanism.  The linear coefficient of
  the return map is load-bearing.
- **Slope-seven face:** the current positive mixed face supplies heterogeneous
  geometry but no semantic transport.  Its next useful invariant is the
  induced multiplier modulo `17`, not only support or Bockstein class.
- **THM-789 global choice:** recurrence at the hostile anchor is not enough;
  THM-789 already proves that local persistence can be the obstruction.  Any
  use of `4/17` must transport alternatives between components, not merely
  stay at this component.

## Scope and stopping certificate

The scan is exhaustive only for fixed points of the abstract maps
`y |-> {by}` with `2 <= b <= 34` in the stated raw word.  It proves no physical
multiplier-`18` handoff, no endpoint chronology, no row exclusion, and no case
of LRC(14).  The exact positive signal is nevertheless sharp: among all
smaller expanding integer bases, none has even an interior fixed point.  Future
work should now search physical handoffs for the congruence class `1 mod 17`,
or prove that the admissible handoff semigroup avoids it.

## Reproduction

```bash
python3 04-computation/lrc14_delayed_word_foreign_slope_fixed_point_scout_20260728.py
python3 -O 04-computation/lrc14_delayed_word_foreign_slope_fixed_point_scout_20260728.py
```

The two transcripts are byte-identical.  SHA-256:

```text
script  a2b8e5fabba00e17fe5c50f40f8cfc0a5de99061bb2faf370a4c105a72640982
output  c9f1de6fc1016adcb614422b9069ed05318f837c60cdd3de42cb24fada75135c
```
