# LRC14 Safe-Component Stalk Descent (codex-2026-06-26-S193)

User asked for more creative LRC angles after the S188/S189 fiber-zipper and
carrier-fusion work.  The new angle was to reverse the direction of HYP-3026:
instead of adding more sidecars, ask whether the sidecars descend from a local
safe-component stalk.

## Move

The tested object is the stalk of the largest strict safe component:

```text
left endpoint owner residues
peak bottleneck owner residues
right endpoint owner residues
component length
local peak height
open/boundary status
```

This is deliberately not the full barcode and not the exact magnitude cocycle.
It forgets non-largest bars and asks whether one component carries enough
local topology/normal-fan information to split a hard automatic fiber.

## Result

`04-computation/lrc14_safe_component_stalk_codex_s193.py` filters the HYP-2963
bank to the AP/GW automatic word `MFCMMCCFFFCCC`, giving `639` packets.  The
stored result is
`05-knowledge/results/lrc14_safe_component_stalk_codex_s193.out`.

Readout:

```text
residue_terminal_fiber     mixed_route=27 max_mixed=30
owner_only_stalk           mixed_route=7  max_mixed=5
coarse_component_stalk     mixed_route=2  max_mixed=2
exact_component_stalk      mixed_route=0  max_mixed=0
magnitude_cocycle          mixed_route=0  max_mixed=0
```

The coarse residuals are only two size-2 open-route collisions:

```text
13->159 / 13->117
13->118 / 13->104
```

They mix `Q-WITNESS` with `COVERING-MOMENT`, not boundary with open.  The exact
local stalk resolves them.

## Proof-Carrying Interpretation

The exact stalk is more geometric than the magnitude cocycle:

```text
exact component length + peak height + owner residues
```

looks like a common descent source for:

```text
barcode shape
active normal fan
endpoint current
local Cech/topology data
Fejer interval anchor
```

If this survives beyond the target word, the HYP-3026 fusion signature can be
partly reorganized: some fields become derived local data rather than primitive
packet fields.

## Caution

This is not yet a full-bank theorem.  The exact stalk may be too sharp, and it
still uses exact rational interval data.  The valuable proof target is the
coarse-to-exact split: prove the two coarse residual patterns directly, then
stress the exact stalk on larger automatic fibers.

Assumption challenged: a proof carrier does not have to be a runner, an
automatic state, a global barcode, or exact magnitude.  A single local
safe-component germ may be the smaller object that explains why the larger
sidecar packet works.
