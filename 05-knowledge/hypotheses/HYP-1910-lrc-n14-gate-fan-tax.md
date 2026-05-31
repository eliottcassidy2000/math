---
id: HYP-1910
status: OPEN
source: codex-2026-05-31-S440
related:
  - THM-357
  - THM-360
  - THM-366
  - THM-365
  - HYP-1837
  - HYP-1839
  - HYP-1844
  - HYP-1880
  - HYP-1890
---

# HYP-1910: The n=14 proof is a gate fan-tax and product-depth debt certificate

## Statement

At LRC denominator `n=14`, any open-cover counterexample must contain a
`14`-gate by the unit endpoint obstruction.  But every `14`-gate carries a
local endpoint fan tax: if its `28` endpoints are covered only by lower columns
`1..13`, the exact minimum cover has size `8`:

```text
(7, 1, 2, 3, 5, 9, 11, 13)
```

The forced private-row columns are

```text
(1, 3, 5, 7, 9, 11, 13).
```

Thus a `14`-gate locally demands the six unit residue columns, the half-gate
`7`, and at least one even bridge before it can be globally useful.  The
conjectural `n=14` proof should charge every `14`-gate against this heptagonal
fan tax, then show the remaining `13`-speed budget cannot also pay the
product-depth endpoint debt exported by the gates.

## Evidence

`lrc_n14_creative_reframes_s440.py` audits seven canonical `n=14` rows through
modulus-cover, antipodal-fold, owner-debt, and product-depth lenses.

The local `14`-gate set-cover subproblem is exact:

```text
owner 14 endpoints, lower 1..13:
  exact lower-cover size = 8
  exact columns = (7, 1, 2, 3, 5, 9, 11, 13)
  forced columns = (1, 3, 5, 7, 9, 11, 13)
```

Each unit residue column has `2` private endpoints, while the half-gate `7`
has `12` private endpoints; all private endpoints sit at extra product depth
`{2:+1,7:+1}`.

The product-depth comparison shows why the gate route fails as a disproof
architecture.  The seven-ladder and S380 gate ladder have the same root-weighted
frontier mass but different denominator pressure:

```text
seven-ladder:
  unprotected = 84
  depths = {2:+0,7:+1}:60, {2:+2,7:+1}:24
  frontier_mass = 66/7
  denominator_pressure = 1092

S380 gate ladder:
  unprotected = 168
  depths = {2:+1,7:+1}:120, {2:+3,7:+1}:48
  frontier_mass = 66/7
  denominator_pressure = 4368
```

So the 14-multiple ladder does not erase debt.  It shifts the same frontier mass
one 2-adic step deeper and quadruples the denominator pressure.

Owner-wise endpoint debt identifies which speeds export the debt.  In S380 the
top exposed endpoint owners are:

```text
154 = 2*77: 48 labels
168 = 8*21: 48 labels
182 = 2*91: 48 labels
14,42,70,126: 24 labels each
```

This supports an owner-charge descent: endpoint debt should be charged to the
gate owner that exports it, not just to the speed set globally.

## Reframes

1. **Gate fan tax.**  A `14`-gate locally consumes a heptagonal fan:
   six unit residues, half-gate `7`, and an even bridge.
2. **Product-depth debt.**  Gate-heavy repairs push exposed endpoints into
   the product of the `2`- and `7`-adic depth directions.
3. **Antipodal CRT fold.**  The quotient `t ~ t+1/2` separates even columns
   from odd high/low columns; at `n=14`, this must be coupled to the mod-7
   unit cycle.
4. **Owner-charge descent.**  A proof should assign exported endpoint rows to
   owner speeds, then show every possible labelled endpoint cycle has a
   positive owner-charge divergence.

## Predictions

1. Any primitive `13`-speed branch containing a `14`-gate either lacks the
   local fan tax, has a small-denominator witness, or exports positive
   product-depth endpoint debt.
2. A branch-and-bound proof should branch first on the forced fan columns
   `(1,3,5,7,9,11,13)` around a chosen `14`-gate, not on raw speed magnitude.
3. The S380 ladder should be locally extremal for preserving frontier mass
   while pushing debt deeper; perturbations should either increase frontier
   mass, reopen a small modulus, or create positive gap sooner.
4. The final certificate should look like a Hall/Farkas dual over rows:
   unit rows force gates, gate fan rows force odd/half/even columns, and
   product-depth rows have positive owner-charge divergence.

## Sources

- `04-computation/lrc_n14_creative_reframes_s440.py`
- `05-knowledge/results/lrc_n14_creative_reframes_s440.out`
- `07-reflections/lrc-n14-creative-reframes-s440.md`
- THM-357
- THM-360
- THM-365
- THM-366
- HYP-1837
- HYP-1839
- HYP-1844
- HYP-1880
- HYP-1890
