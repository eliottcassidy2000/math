---
id: THM-360
name: lrc-unit-endpoint-divisibility-filter
status: PROVED
date: 2026-05-30
session: codex-2026-05-30-S362
depends_on:
  - THM-357
  - THM-358
---

# THM-360: Unit Endpoint Divisibility Filter

## Statement

Let `V` be a `k`-speed Lonely Runner instance and set `n=k+1`.  For any unit
residue `a mod n`, the boundary point

```text
t = a/n
```

is strictly protected by a speed `v` only if

```text
v == 0 mod n.
```

Consequently, if no speed in `V` is divisible by `n`, then `t=1/n` is a lonely
witness.  Every full-open-cover counterexample must contain at least one speed
divisible by `n`.

Moreover, any speed divisible by `n` protects every unit endpoint `a/n`.

## Proof

A speed `v` strictly protects `a/n` exactly when

```text
||v a/n|| < 1/n.
```

Because `a` is a unit modulo `n`, multiplication by `a` permutes residue
classes modulo `n`.  If `v` is not divisible by `n`, then `v a mod n` is a
nonzero residue.  Therefore the circular distance from `v a/n` to an integer
is at least `1/n`, not strictly less.

Thus strict protection is possible only when `v a == 0 mod n`, equivalently
`v == 0 mod n`.

If `v` is divisible by `n`, then for every unit `a`,

```text
v a/n in Z,
```

so `||v a/n|| = 0 < 1/n`; such a speed protects every unit endpoint.

Taking `a=1` proves that if no speed is divisible by `n`, then `1/n` is not
strictly protected by any speed and in fact satisfies the closed Lonely Runner
inequalities for all speeds.  By THM-357, a full-open-cover counterexample
must protect every endpoint, so it must contain a speed divisible by `n`.

## Use

THM-358 says the initial-segment tight examples expose the nonzero unit
endpoints.  THM-360 says the first obstruction to a counterexample is
divisibility by `n`: before any higher endpoint quotient can matter, the speed
set must contain a multiple of `k+1`.

This is the first exact quotient-descent filter in HYP-1813.

## Related

- THM-357: endpoint-protection trichotomy.
- THM-358: initial-segment unit skeleton.
- HYP-1810: unit-boundary skeleton.
- HYP-1813: Bohr-boundary descent.
- `07-reflections/lonely-runner-bohr-descent-formal-session-s362.md`.
