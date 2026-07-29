# Ranked-suffix bootstrap minimum-seed battery

**Status: FINITE-EXACT SCOPED RESEARCH NOTE.  No uniform theorem.**

This companion sharpens the four-root ranked-suffix scalar battery from a
fixed ordering into an exact target-set-selection problem.  Its closest
proved mechanism is THM-753's safe-peel recursion: both repeatedly delete
nonbinding objects and isolate an irreducible core.  THM-753 acts on a speed
family at complement tight points; the present object acts on finite
first-apex proof obligations.

The least-used prior sidecar is the tournament table in
`00-navigation/LRC-LENS-MAP.md`, which proposed insertion obligations,
residual length, and a heuristic decision order while requiring the literal
state and future-operation bank.  Here that heuristic becomes an exact
monotone activation system.  It is a directed threshold hypergraph, not a
tournament: activating one apex can require several earlier exclusions, so
no pair orientation contains the needed state.

## 1. Activation system

Fix one root and its finite hitting gate `A`.  For `a in A` and
`P subset A minus {a}`, let `T_a(P)` mean:

```text
the five largest individual coverages of the literal a-residual
over labels outside P union {a} have sum strictly below its mass.
```

Then `T_a(P)` is monotone in `P`: deleting more allowed labels cannot
increase the five largest values.  Starting from a seed `H subset A`,
repeatedly add every `a` with `T_a(P)` until reaching the least fixed point.
If it reaches all of `A`, the nonseed branches scalar-peel in any compatible
activation order.

The exact connection data are:

```text
source:       an apex branch a with excluded-prefix state P;
target:       scalar-terminal status T_a(P);
map:          P |-> the allowed residual top-five sum;
preserved:    the lossless partition of every putative six-cover;
lost:         deletion order if P and the apex rank are forgotten;
sidecar:      the marked active set P and literal residual.
```

This is compatible with THM-2894: the order is retained explicitly rather
than reconstructed from the unmarked residual semilattice.

## 2. Seed scope warning

**Seed membership is a workload target, not a branch proof.**  Declaring
`H` active only reserves those branches for nonscalar cap/flag arguments.
To finish a root, choose an order on `H` and certify seed `i` on the suffix
excluding only the earlier seeds.  It may not assume later seeds are absent.

Thus minimum seed size optimizes how many branches require deeper work; it
does not close those branches or the root.

## 3. Exact minimum seeds

The fixed decreasing-coverage order left `25` branches on the four-root
battery.  Exhaustive target-set search reduces the exact minimum seed sizes
to

```text
(4,3,5,3),
```

or `15` total nonscalar targets.  The **numbers of minimum seeds** are,
separately,

```text
(2,1,1,4).
```

The complete minimum families are:

```text
E=(2,8,9,10,11,13,14), K=19:
  (17,23,46,24), (17,23,46,29);

E=(1,3,9,10,11,12,14), K=20:
  (39,23,16);

E=(2,5,9,11,12,13,14), K=21:
  (16,23,19,40,46);

E=(2,3,4,5,6,7,8), K=13:
  (22,26,33), (22,18,36), (22,18,44), (22,18,33).
```

The first deterministic witnesses activate the remaining branches in round
sizes

```text
(2,10,1,2), (5,3,5,4), (1,3,3,3,1,1,1,1,2), (4,2,3,1).
```

In total `58` nonseed branches scalar-activate.  The verifier exhausts
`34,661` seed subsets and performs `574,940` exact activation-predicate
checks, including every subset smaller than the reported minima.

## 4. Interaction with the proved rank-one battery

The earlier complement-cap battery genuinely closes the rank-one apices

```text
19, 39, 16
```

on the first three roots.  The second and third unique minimum seeds contain
their certified apex.  Neither size-four minimum seed on the first root
contains `19`; the first successful size-five seed that does is

```text
(19,17,23,46,24).
```

A hostile control is decisive: using only the already certified rank-one
apex as seed activates **zero** further apices on every root.  The closures
remain `1/19`, `1/20`, and `1/21`, leaving `18`, `19`, and `20` hard
branches.  Therefore the prior rank-one certificates do not yet recompose
any whole root.  Future optimization must jointly consider seed size,
seed order, and cap/flag margins.

## 5. Exact controls and reproduction

For each of the `73` apex profiles, the complement-of-gate top five is
globally sealed by the strict discrepancy bound, then combined with every
other gate-label coverage.  The largest tail starts at `3030`.  Controls:

```text
16 root vector/scalar checks,
292 outside-gate vector/scalar checks,
1298 gate-label vector/scalar checks,
73 literal/direct residual reconstructions.
```

Ordinary and optimized runs are byte-identical.  The complete activation
system digest is

```text
43edbe3acf8ec1fbc66301d8c52500afe083c212db2242f76fc8279cf38d7124.
```

```text
04-computation/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.py
SHA-256 a7a77dc433b21d94a54524064ccd62e553ed67ae8a3d3364bf79c41c36849d04

05-knowledge/results/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.out
SHA-256 bdcf896d152d206b3ae77a3568609887190e1ba991909481769d7ad560a68835
```

The universe is only four of the `3432` seven-body roots.  This result does
not prove any seed branch, the seven-body/six-slot rung, or LRC(14).
