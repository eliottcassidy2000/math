---
id: THM-566
title: LRC14 covering sets have no uniform bounded-denominator rational witness
status: PROVED
date: 2026-06-22
session: codex-2026-06-22-S93
depends_on:
  - THM-523
  - THM-366
related:
  - HYP-2052
  - HYP-2865
  - THM-524
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_covering_bounded_denominator_obstruction_codex_s93.out
---

# THM-566: Covering Sets Have No Uniform Bounded-Denominator Witness

## Statement

For every integer `B >= 2` there is a primitive LRC14 covering 13-set `S_B`
such that no rational point `a/D` with `1 <= a < D <= B` is a level-`1/14`
lonely witness for `S_B`.

Equivalently, the THM-523 reduction to covering sets does not make a fixed
bounded-denominator witness theorem possible.

## Construction

Let

```text
L_B = lcm(1,2,...,B)
S_B = {1,2,...,11,13,84 L_B}.
```

Then `S_B` has 13 distinct positive speeds.  It is primitive because it contains
`1`.  It is a THM-523 covering set because:

- for each `q=2,...,11`, the speed `q` itself is present;
- `q=12` divides `84 L_B`;
- `q=13` is present;
- `q=14` divides `84 L_B`.

Now fix any `D <= B` and any numerator `a`.  Since `D | L_B`, also
`D | 84 L_B`.  Hence

```text
(84 L_B) * a / D
```

is an integer, so the last runner is exactly at the observer:

```text
||(84 L_B) a/D|| = 0 < 1/14.
```

Thus `a/D` cannot be a lonely witness for `S_B`.  This holds for every
`D <= B`, proving the claim.

## Consequences

This refutes the strong bounded-denominator formulation suggested by the
observed witness `17/41` for `{1,...,11,13,84}`.  The obstruction survives the
covering-set restriction: the examples above are already primitive covering
13-sets.

The result is compatible with the positive random evidence.  In bounded or
non-adversarial banks, small witnesses are common; the obstruction is divisor
loading.  For example the S93 scout found:

```text
{1,...,11,13,84}:        17/41
{1,...,11,13,84*6}:      22/53
{1,...,12,182}:          2/27
```

and on `S_m={1,...,11,13,84m}`, all `m<=5000` had a witness with denominator
at most `67`.  But for a proposed fixed bound `B`, the row with
`m=lcm(1,...,B)` kills every denominator up to `B` simultaneously.

## Proof Route Impact

The finite residue-atlas idea remains useful only in a scaled or conditioned
form:

1. bounded denominator for bounded-speed or non-divisor-loaded covering sets;
2. a denominator bound in terms of the least modulus not killed by any speed
   (the HYP-2052+ direction);
3. a hybrid with THM-524 binding-pair switches, HYP-2864 sheet-gcd quotients,
   or THM-565 finite-ruler sampling for divisor-loaded towers.

What is ruled out is the zero-analysis closure "prove one absolute `D <= B0`
for every covering set."

## Verification

The script

```text
04-computation/lrc14_covering_bounded_denominator_obstruction_codex_s93.py
```

prints exact divisibility certificates for `B=14,26,41,67,80`, verifies the
constructed rows are primitive covering sets, and then records local positive
scouts explaining why the false conjecture looked plausible.
