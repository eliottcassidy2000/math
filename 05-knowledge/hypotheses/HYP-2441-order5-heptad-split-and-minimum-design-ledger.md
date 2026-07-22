# HYP-2441 - The order-5 branch forces a heptad split and a minimum-design ledger

**Status:** OPEN synthesis / next attack. Codex 2026-06-11.

Assume a putative extremal Type II `[72,36,16]` code has an order-5
automorphism of type `5-(14,2)`. HYP-2439 forces the fixed projection to be
`e8+e8`, with one fixed coordinate in each `e8` block. Thus the fourteen
5-cycles split into two heptads:

```text
fixed point f0 + seven 5-cycles A1..A7,
fixed point f1 + seven 5-cycles B1..B7.
```

Inside each `e8` block, the seven tetrads through the fixed coordinate form a
Fano plane on the seven cycle-orbit coordinates. Lifting them gives:

- `7` fixed minimum words through `f0`;
- `7` fixed minimum words through `f1`;
- no fixed minimum word through both fixed points.

Since the extremal length-72 scalar gate has `A_16=249849`, the order-5 action
then gives the exact orbit ledger

```text
fixed minimum words = 14,
moving minimum words = 249849 - 14 = 249835 = 5 * 49967.
```

The design lambdas from the putative `5-(72,16,78)` minimum layer give immediate
residual counts:

```text
blocks through a fixed coordinate: 55522 - 7 = 55515 = 5 * 11103 moving blocks,
blocks through the pair {f0,f1}: 11730 - 0 = 11730 = 5 * 2346 moving blocks.
```

HYP-2430 gives an exact negative control on this ledger.  Even after all
fourteen fixed minimum blocks are frozen, the order-five-invariant incidence
map from moving `16`-subset orbits to `5`-subset orbits has kernel dimension at
least `823261001634556`.  Thus the displayed counts and every five-incidence
marginal remain radically non-identifying at the level of formal signed block
systems.  Binary linearity, self-orthogonality, and the `F_16` glue are
load-bearing; the trade does not itself construct a second simple design.

## Next Gate

Build the nonfixed-eigenspace completion problem over the primitive factor
`Phi_5`, i.e. over `F_16`, with this split-heptad fixed boundary:

1. enumerate Hermitian self-dual `[14,7]` `F_16` candidates up to the heptad
   stabilizer and the two Fano-plane structures;
2. enforce binary lift minimum distance `>=16`;
3. enforce the residual minimum-design orbit ledger above;
4. turn the surviving moves into a support-building Tournament Analysis.

The expected obstruction is no longer scalar. It should be a compatibility
failure between the `e8+e8` fixed Fano planes and the `F_16` nonfixed glue.

## Links

- Extends HYP-2439 and HYP-2440.
- Supports OPEN-Q-067.
