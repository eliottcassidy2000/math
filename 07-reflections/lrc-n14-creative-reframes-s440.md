# LRC n=14 Creative Reframes S440

**Session:** codex-2026-05-31-S440
**Script:** `04-computation/lrc_n14_creative_reframes_s440.py`
**Stored output:** `05-knowledge/results/lrc_n14_creative_reframes_s440.out`

The best new object from this session is the **14-gate fan tax**.

The old split was:

```text
no 14-gate   -> unit endpoints a/14 survive
has 14-gate  -> unit layer is killed, but endpoint debt appears
```

S440 adds a local invoice between those two clauses.  A single `14`-gate has
`28` endpoints.  If a branch treats that gate as locally maximal, lower columns
`1..13` cover those endpoints only with exact minimum `8`:

```text
(7, 1, 2, 3, 5, 9, 11, 13)
```

The forced private-row columns are:

```text
(1, 3, 5, 7, 9, 11, 13)
```

So the local geometry is not vague.  The `14`-gate demands the six unit residue
columns, the half-gate `7`, and one even bridge.  This looks like the `n=14`
analogue of the `n=16` nine-cover, but with the CRT product `2 x 7` visible.

## Product-Depth Debt

The S380 ladder remains the cleanest negative lesson.  It pays every
small-denominator invoice and every one of its columns is forced by some
endpoint row, but it still leaves `168` endpoints exposed:

```text
{2:+1,7:+1}:120
{2:+3,7:+1}:48
```

Compared to the seven-ladder, it preserves the same frontier mass:

```text
frontier_mass = 66/7
```

but its denominator pressure jumps:

```text
seven-ladder pressure = 1092
S380 pressure         = 4368
```

That is a strong reframe: the 14-multiple ladder is not wasteful.  It is a
locally necessary way to push the same debt deeper into the `2`-adic direction
while carrying the `7`-adic payload.

## Antipodal CRT Fold

The quotient `t ~ t+1/2` did not create a new independent gap measure in the
sample: the pair-gap matched the ordinary gap ratio.  Still, it is useful
structurally.  Even speeds cover both sides of an antipodal pair the same way.
Odd speeds cover one side through a low band and the other through a high band.

For `n=14`, this means the even fold cannot be separated from the mod-7 unit
cycle.  A proof should couple:

```text
2-adic antipodal parity
mod-7 heptagonal residue fan
product-depth endpoint debt
```

## Owner-Charge Descent

Endpoint debt should be charged to owner speeds.  In the S380 row, the largest
exposed owners are:

```text
154 = 2*77
168 = 8*21
182 = 2*91
```

with `48` exposed endpoint labels each.  The smaller `14`, `42`, `70`, and
`126` owners carry `24` each.

This suggests the next proof object: a labelled endpoint cycle from THM-365
should carry an owner-charge potential.  Gates can move the charge to deeper
product-depth rows, but the primitive `13`-speed budget cannot make the charge
divergence vanish.

## Artifacts

- HYP-1910
- `04-computation/lrc_n14_creative_reframes_s440.py`
- `05-knowledge/results/lrc_n14_creative_reframes_s440.out`
