---
id: THM-1141
title: The clustered beat mechanism is refuted; its replacement actual-mean lemma is also refuted
status: CORRECTED NEGATIVE RESULT. Clustered phase sweep is too short to force alignment. The later `max gap >= 4/3 actual mean` target is false by THM-1148/MISTAKE-171. Sample telemetry remains descriptive; uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128c69, corrected by codex-2026-07-18-S74
depends_on: [THM-1140]
related: [THM-1137, THM-1147, THM-1148, MISTAKE-163, MISTAKE-164, MISTAKE-169, MISTAKE-171]
script: 04-computation/beat_structure_kps_S128c69.py (+ .out)
---

# THM-1141 -- two tempting four-comb mechanisms fail

## 1. The beat/alignment proposal fails

For two nearby killers, phase difference changes across a core component of
length `ell` by

```text
ell(k_j-k_i)
```

full turns.  On the sampled clustered quadruples this lies between `0.027`
and `0.54`, below the one full turn needed to force an alignment.  The
configuration is frozen across the component.  Alignment begins to sweep
only for multiplicatively spread killers, where THM-1137's exact interval
transfer is the appropriate proof mechanism.

This is a genuine negative conclusion about the proposed proof route.  It is
not a theorem that all clustered configurations are safe.

## 2. The measured gap surplus

The sampled frozen configurations nevertheless have long survivor gaps.  In
the old worst sample

```text
P={1,3,5,6,7,8,11,12},
killers=(371,374,377,379),
```

the four phases are genuinely spread, yet `7k4L=2.358`.  Across 500 sampled
clustered quadruples, the minimum observed `Lk4` was `0.3584`; none landed in
the proposed uniform-spread range.  These values are reconnaissance only.

The early interpretation was that distinct combs must therefore have

```text
longest gap >= (4/3)(actual mean survivor gap).           (1)
```

THM-1148 refutes (1) exactly on a maximal legal core interval:

```text
P={1,...,8},
J=[1/14,13/112],
killers=(108,109,110,111),
L/actual_mean=638/573<4/3.                               (2)
```

The desired metric conclusion still holds there with `7k4L=319/72>1`.
Tooth overlap lowers the number of components and raises their mean, so an
almost uniform gap word need not be a hard row.  The actual mean and the
crude union-bound baseline are different denominators.

## 3. Corrected structural lesson

The proof-bearing object is not scalar nonuniformity.  It is

```text
(tooth-overlap cluster count, labelled metric gap word,
 endpoint owners, wall slopes).
```

THM-1148 replaces (1) by three sound tools:

- a sharp four-residue multiplier Kakeya cone;
- the exact mass/component gate `Q4>0`;
- the corrected THM-1137 `Phi` transfer, with a `9/5` ratio corollary.

The old sample ratio `3.34` was measured against the uniform-interleaving
benchmark `m0=3/(7 sum k)`, not against the actual component mean.  To keep
the two effects separate, set

```text
D=L_max/mu_actual,
B=mu_actual/m0.
```

Then the sharp target is

```text
7 k4 L_max>1  iff  D B > (sum k)/(3 k4).               (3)
```

Endpoint dispersion controls `D`; tooth overlap, component mass, and
multiplier charts control `B`.  The exact row (2) succeeds through baseline
gain despite having `D<4/3`.

Their first infinite proof-method residual is `m(3,4,5,6)`, not the whole
clustered majority.  Uniform `r=5` remains open.

## Tournament audit

Ordering wall events produces a transitive tournament and loses the metric
surplus.  A useful orientation must place labelled residue gaps or interval
proof states at the vertices and retain overlap/endpoint data on the edges.
The naked runner tournament cannot distinguish the false statement (1) from
the valid metric cones.
