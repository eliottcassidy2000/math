---
id: THM-366
name: lrc-small-denominator-divisibility-sieve
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S388
depends_on:
  - THM-357
  - THM-360
  - HYP-1850
results:
  - 05-knowledge/results/lrc_small_denominator_sieve_s388.out
---

# THM-366: LRC Small-Denominator Divisibility Sieve

## Statement

Let `V={v_1,...,v_{n-1}}` be a Lonely Runner speed set at threshold `1/n`.
Let `m` be an integer with `2 <= m <= n`, and let `a` be a unit modulo `m`.
At the point

```text
t = a/m,
```

a speed `v` can strictly forbid, or strictly protect, `t` only if

```text
m | v.
```

Equivalently, if no speed in `V` is divisible by `m`, then `t=a/m` is a
lonely witness in the closed Lonely Runner sense.  More precisely:

1. if `m < n`, then `t=a/m` is an open lonely witness, since every distance is
   strictly larger than `1/n`;
2. if `m = n`, then `t=a/n` is a boundary lonely witness.

Consequently, every full-open-cover counterexample at level `n` must contain,
for each `m in {2,...,n}`, at least one speed divisible by `m`.

## Proof

Fix `m <= n` and a unit residue `a mod m`.  A speed `v` strictly forbids
`t=a/m` exactly when

```text
|| v*a/m || < 1/n.
```

If `m` does not divide `v`, then `v*a` is nonzero modulo `m`, because
multiplication by a unit permutes residue classes modulo `m`.  Therefore the
circular distance from `v*a/m` to an integer is at least

```text
1/m.
```

Since `m <= n`, we have

```text
1/m >= 1/n.
```

Thus the strict inequality `||v*a/m|| < 1/n` is impossible unless `m | v`.

If no speed in `V` is divisible by `m`, every speed satisfies

```text
|| v_i*a/m || >= 1/m >= 1/n.
```

So `t=a/m` satisfies the closed Lonely Runner inequalities for all speeds.
When `m<n`, the inequality is strict because `1/m > 1/n`, so a small open
neighborhood of `t` is lonely.  When `m=n`, the point may lie on the boundary,
but it is still not in the open forbidden set.

The final claim follows by contrapositive.  If a speed set misses some
divisibility class `m in {2,...,n}`, the corresponding unit point `a/m` is not
strictly forbidden by any interval.  Hence the open forbidden union cannot be
the whole circle.

## Interpretation

THM-360 is the `m=n` unit-endpoint case of this theorem.  THM-366 promotes it
to a full small-denominator sieve:

```text
counterexample => every denominator m <= n is hit by a divisible speed.
```

This is a necessary condition, not a sufficient one.  The S388 computation
shows the distinction sharply: replacing a missing initial-segment speed `d`
by `lcm(n,d)` satisfies the sieve, but the audited `n=14,15,18` transfer rows
remain positive-gap, and the focused endpoint audits remain terminal-core
empty.

## Verification Record

`04-computation/lrc_small_denominator_sieve_s388.py` audits the theorem on
initial segments, known sporadic tight examples, `n=14` quotient ladders, and
`n=14,15,18` gate-transfer families.  The stored output records:

```text
initial n=14: missing=(14,), boundary t=1/14, unprotected=6, coreE=0
drop 13 add 14: missing=(13,), open t=1/13, gap/th=0.045455
drop 13 add 182: missing=-, sieve-complete, gap/th=0.065934, coreE=0
n14 seven-ladder: missing=-, sieve-complete, gap/th=0.005411, coreE=0
drop 17 add 306 at n=18: missing=-, gap/th=0.052288, coreE=0
```

The transfer tables separate two failure modes:

```text
miss a small denominator -> immediate rational witness;
cover it by an lcm gate   -> positive gap plus exported endpoint debt.
```

## Use

The theorem gives the first branch in any proof or disproof search:

```text
not sieve-complete -> done by a small-denominator witness;
sieve-complete     -> study endpoint debt, critical surplus, and labelled cores.
```

It also explains why the initial segment `{1,...,n-1}` is tight in exactly the
way it is: it covers every denominator below `n` and misses only the top
denominator `n`, leaving the unit boundary skeleton `a/n`.

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-360: unit endpoint divisibility filter.
- HYP-1850: modulus-sieve localization.
- HYP-1851: gap-floor/no-multiple-of-n at full measure.
- HYP-1844: quotient-debt export.
- `04-computation/lrc_small_denominator_sieve_s388.py`.
