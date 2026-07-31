# Independent hostile audit of the THM-2771 all-cell copy stratification

**Status:** `FINITE-EXACT`, scratch only.  This is an independent audit of
the proposed scope of reserved THM-2818.  It proves no row exclusion and no
LRC(14) conclusion.

## Verdict

The proposed all-cell decomposition survives.

The independent companion discovers, rather than supplies, the nonzero cell
set by scanning all `7*13=91` physical `(clock,target)` pairs.  Exactly
twenty-eight right-wing cells survive.  In every cell the coefficient is a
positive sum of equal-length, equal-weight, equal-content interval copies.
The only raw/effective mismatch occurs in the three target-twelve cells,
where the delayed terminal selector is an alternating `1,0,1,0,...` word on
half-step chains.  The effective copy counts are exactly

```text
2, 121, 265, 254.
```

The Bockstein is copywise additive in all twenty-eight cells.  The ancestry
claim is valid locally, cell by cell, but must not be promoted to a single
global prototype: the first global change is at `(clock,target)=(2,2)`.

## 1. Independence

The audit imports only four LF-hash-pinned canonical companions:

```text
lrc14_fully_marked_root_zero_target_profile_thm2749.py
  d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e

lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
  25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552

lrc14_root_zero_full_target_semantic_clutch_20260728.py
  208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8

lrc14_full_arm_orbit_path_sheet_audit_thm2791.py
  1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543
```

It does not import `.scratch/lrc_cofiber_copy_bridge/audit.py`,
`all_cells.py`, their stored output, or their expected tables.  In
particular, it does not encode:

- the twenty-eight cells;
- the raw/live/dead counts;
- the exceptional block lengths;
- either weight or ancestry digest;
- the common gcd or primitive copy factors; or
- any Bockstein residue.

Those quantities are reconstructed from the physical section, delayed
coefficient functional, literal ancestry walls, and canonical wing
coefficient functional.

## 2. The independently discovered cell table

Every copy has length

```text
26,444,880.
```

The compact count table is:

| clock | targets | cells | raw/live/dead per cell |
|---:|---|---:|---|
| 1 | `3,...,11` | 9 | `2/2/0` |
| 1 | `12` | 1 | `241/121/120` |
| 2 | `2,...,9` | 8 | `2/2/0` |
| 2 | `12` | 1 | `528/265/263` |
| 3 | `2,...,7,10,11` | 8 | `2/2/0` |
| 3 | `12` | 1 | `506/254/252` |

Thus there are twenty-five ordinary cells and three exceptional cells:

```text
25*(2 raw = 2 live + 0 dead),
241 raw = 121 live + 120 dead,
528 raw = 265 live + 263 dead,
506 raw = 254 live + 252 dead.
```

The stored output records all twenty-eight rows individually, including
weights, ancestry counts and digests, copy coefficients, primitive factors,
and Bockstein residues.

For every raw piece the audit evaluates source carry `12` and translated
target carry `6` independently.  The only outcomes are

```text
(c,c) with c=103,478,815,440,       or       (0,0).
```

There is no unequal source/target value, signed piece, or second nonzero
content.

## 3. Target twelve is an alternating half-step chain

The step is derived from the ambient period rather than inferred from the
claimed block table:

```text
h=T/(2*13^5)=401,080,680.
```

The script orders pieces by left endpoint and starts a new block precisely
when the next left endpoint is not `h` away.  All interblock jumps are

```text
20,054,034,000 = 50h.
```

It then labels a position by the independently evaluated delayed content.
Every discovered block is exactly

```text
1,0,1,0,...,
```

starting live.  The reconstructed profiles `(length,live,dead)` are

```text
clock 1: (145,73,72), ( 96,48,48),
clock 2: (143,72,71), (289,145,144), (96,48,48),
clock 3: (143,72,71), (289,145,144), (74,37,37).
```

Consequently the formerly opaque multipliers have the exact formulas

```text
121 = ceil(145/2)+ceil( 96/2),
265 = ceil(143/2)+ceil(289/2)+ceil(96/2),
254 = ceil(143/2)+ceil(289/2)+ceil(74/2).
```

This is a parity/chain mechanism, not cancellation.

## 4. Local ancestry and the global prototype change

For every cell the audit takes the hull from its first left endpoint to its
last right endpoint.  It constructs all `U` and `V` contributor walls from
the defining interval systems and verifies that no wall crosses either
cell hull.  Hence all pieces in that particular cell have identical literal
ancestry sets.  This proves local identity without enumerating a contributor
set at every one of the many target-twelve pieces.

There are three distinct wall-free chamber pairs but only two resulting
label prototypes:

```text
clock 1:
  |U|=966,606, |V|=28,534,
  digest=15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd;

clocks 2 and 3:
  |U|=966,574, |V|=28,534,
  digest=bc32242480e738421f7f9d7c1f300c09077d9c76ef4d1156083a87673b4d598c.
```

The first lexicographic departure from the clock-one prototype is `(2,2)`.
Direct set difference, not cardinality alone, gives

```text
U loses 32 labels and gains 0;
V loses 0 labels and gains 0.
```

The supplied path remains active on both sides.  Therefore:

> “ancestry-identical copies” is a within-cell statement.  It is false as a
> single global-prototype statement, even though clocks two and three happen
> to share their second prototype.

This scope distinction should remain explicit in THM-2818.

## 5. Independent coefficient and Bockstein audit

The two discovered copy values are

```text
clock 1:
  2,854,063,240,791,928,925,760;

clocks 2 and 3:
  2,853,968,755,527,296,447,040.
```

Their gcd and derived primitive factors are

```text
g0 = 5,905,329,039,529,920,
(copy values)/g0 = (483,303,483,287),
gcd(483,303,483,287)=1.
```

Both copy values and `g0` have exact `13`-adic valuation one.  Without
supplying a residue, the audit derives

```text
copy beta:
  clock 1 = 9,
  clocks 2 and 3 = 2.
```

For every cell it independently recomputes the complete seven-clock wing
vector and verifies

```text
total vector = live_count * individual-copy vector,
total beta   = live_count * copy beta mod 13.
```

The resulting totals are:

| cells | live multiplier | copy beta | total beta |
|---|---:|---:|---:|
| clock 1, ordinary | 2 | 9 | 5 |
| `(1,12)` | 121 | 9 | 10 |
| clocks 2/3, ordinary | 2 | 2 | 4 |
| `(2,12)` | 265 | 2 | 10 |
| `(3,12)` | 254 | 2 | 1 |

This verifies the Bockstein multipliers as literal positive-copy
multiplicities, including the three sporadic residues.

## 6. Incoming connection to the sharp eleven-label boundary

The newly installed THM-2809 candidate identifies

```text
L_*={2,3,...,12}
```

as the unique positive eleven-label residual after the exact label-zero and
label-one attachment obstructions.  Independently, the present cell census
gives the target supports

```text
S_1={3,4,5,6,7,8,9,10,11,12},
S_2={2,3,4,5,6,7,8,9,12},
S_3={2,3,4,5,6,7,10,11,12}.
```

Hence

```text
S_1 union S_2 union S_3 = L_*.
```

The right-cofiber support excludes exactly labels zero and one, the two
labels isolated by THM-2809's distinct anchor and delayed-digit
obstructions.  Moreover `t=12` is the unique target supported at all three
clocks whose raw/effective split is nontrivial; the other common targets
`3,4,5,6,7` are ordinary two-copy cells.

This is a precise set-theoretic coincidence, not yet a theorem-level map.
THM-2809 varies arbitrary THM-2640 label packets against the full marked
source, while this audit cuts the THM-2771 right wing after a fixed physical
section and coefficient functional.  The cheapest decisive bridge test is
to compare the two target-label section constructors before coefficient
evaluation and determine whether exclusion of labels zero and one is
literally the same factor failure.  Until that is done, the connection is
`VERIFIED` as support combinatorics and `OPEN` as a common mechanism.

## 7. Boundary

This audit supports the coefficient-copy and local-ancestry scope proposed
for THM-2818.  It does not construct a map from the THM-2806 common sheet to
the right cofiber.  In particular, it does not repair the already recorded
native-factor, carrier-mask, endpoint, target-convolution, or root-deck
failures.  Equal coefficient content and local ancestry remain too weak to
identify fully typed physical atoms.

Reproduction:

```bash
python3 .scratch/lrc_cofiber_all_cells_hostile_twist_20260728/audit.py
python3 -O .scratch/lrc_cofiber_all_cells_hostile_twist_20260728/audit.py
```
