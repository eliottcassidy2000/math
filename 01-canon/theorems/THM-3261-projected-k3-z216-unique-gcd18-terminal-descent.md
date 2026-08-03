---
id: THM-3261
title: "Projected-k3 z216 unique gcd18 terminal descent"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/z216-gcd18-terminal/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3251-projected-k3-z218-terminal-descent-and-cap216
related:
  - THM-3245-pointed-divisor-median-cubes-saturation-band-no-go-and-z219-supplier-support
  - THM-3259-charge-paired-modular-free-factor-lift-and-idempotent-localization-obstruction
script: 04-computation/lrc14_j7_k3_z216_unique_gcd18_terminal_descent_thm3261.py
output: 05-knowledge/results/lrc14_j7_k3_z216_unique_gcd18_terminal_descent_thm3261.out
script_sha256: 308f57b090ea9586c92905d9113b79063ef2105c3e7ed4c90c7406d15400d376
output_sha256: d7dae9cd7e8f0305824b30a2b8683abe6d7d8c9bd0509a0d84026322fd65344c
semantic_sha256: b6f83cd9cc405b4e5006ecfbd8ddfb25c115d7a55c270abc55d8d3b2a97b9a18
hash_basis: LF-normalized bytes
audit: >
  An independent hostile audit rederived the complete layer census, unique
  gcd-18 row, screen and Farkas directions, all 30 top-divisor masks, strict
  two-high gap, complete one-high terminal closure, and exact ledger change.
  It independently recomputed the row and byte-matched the frozen normal,
  optimized, and stored evidence, with no assertions, floats, or randomness.
  The divisor-product boundary was accepted as arithmetic only.
---

# THM-3261 -- projected-k3 z216 unique gcd18 terminal descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The first occupied layer below THM-3251

THM-3251 leaves the projected `k=3` necessary atlas at cap `216` and ledger
`373284`.  Its complete first occupied layer is

```text
z1=216: 480 rows = 447 wall + 33 order,                    (1)
row SHA256 53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649.
```

For a row with ruler `L=14 lcm(E)`, stratify by `g=gcd(216,L)`.  The exact
census is

```text
(g,count)=(8,19),(18,1),(24,135),(36,15),(72,310).        (2)
```

Thus there is a unique `g=18` row.  It is atlas row `185`,

```text
E=(1,5,7,9,11,13),   L=630630,   high=62108,
first_d=L/18=35035,  wall=true.                                (3)
```

No other row is classified or closed by the theorem.

## 2. Exact screen and the top divisor point

Running the inherited exact ray/status screen on `(3)` gives

```text
143 states = 0 crude + 113 exact-status + 30 residual,     (4)
direct Farkas=0, inherited exact Farkas=113.
```

The screen, row-semantic, row-failure, mask, and residual digests are

```text
5e2721f49b5f1d9f7dadd156fa7277cccb3a0693b648c2375d714bf47765c4c5,
33ba863d05950923a9a10719f4c17ab4c4b6df77cf8ca3f355b61618b3250d25,
e45b0879387d0df2907a0d17e021659a1f4eafd63959d040cae3d8ca78d7a1e6,
2d753b2dc75f3ba79d3ae0b3f79999669807f6906562ff5072c71bb0ed2fd45e,
652e914b0e9e11ec997c8f168286e6684864770d42b04b94676c4299ddecb812.
```

For every one of the `30` residual denominator masks, if `D=lcm(ds)`, then

```text
D/first_d=18, hence D=L.                                   (5)
```

The full pointed divisor interval in this row is the divisor lattice of
`18=2*3^2`, abstractly the product of a two-chain and a three-chain, and all
residuals occupy its top point.  Equation `(5)` is only denominator metadata:
it constructs no `C2` or `C3` action, common carrier, free product, owner map,
or phase transport.

## 3. Complete one-high terminal exclusion

Because `(3)` is a wall row, every actual residual assignment must contain a
high slot.  The duplicate-permitting two-high union bound has the strict gap

```text
2739297996619 / 899643669274350 > 0,                       (6)
```

so at most one high slot can occur.  Every actual residual therefore lies in
the complete one-high bank.  The exact terminal census is

```text
1 row, 30 masks,
30 zero-high hostile relaxations,
30 one-high cases over 1 low-label set,
30 coarse-cardinality closures,
0 exact-cell fallbacks, 0 max-gap fallbacks, 0 failures,
minimum cardinality slack 8.                               (7)
```

The zero-high cases in `(7)` are not asserted physical; the wall condition
removes them.  In every one-high case the high ray is replaced by its exact
supremum and the low slots range over a complete enlarged state space.  Each
coarse cardinality certificate already contradicts the necessary projected
support size.  Thus all `30` residual masks, and hence the unique row `(3)`,
close.

The terminal semantic, tagged closure, and case-vector digests are

```text
607e7cd3166f078ace15835d71dd3f7f0e732624c957a9273c3f7a277b475e0f,
88af81e3326b00bceb14c280d2557c22210ba9c33050074814baeaaec9e340af,
4dd6438ccc4109f81388a4eb8ee50015e76b6b9050e1521c00763167408b0148.
```

## 4. Exact consequence and boundary

Removing one row gives

```text
373284-1=373283.                                           (8)
```

The projected cap remains `z1<=216`, and the residual occupied layer is

```text
479 rows = 446 wall + 33 order.                            (9)
```

This is a one-row decrement in a necessary projected atlas.  It is not a
physical-cover classification, a `k<=1` result, a final-rung theorem, or an
LRC(14) proof.  In particular the simultaneous appearance of the primes two
and three in `(5)` is not the modular co-occurrence object: the missing datum
is still an action on one physical carrier.

## 5. Exact evidence

The companion pins the THM-3139 engine and the promoted THM-3251 source,
output, and semantic hashes.  It freshly derives all `480` rows, checks `(1)`
and `(2)`, isolates `(3)` uniquely, reruns the complete screen, binds all
residual masks and digests, and recomputes the complete terminal probe.
Normal and optimized runs LF-normalize byte-for-byte to the stored transcript
with empty standard error.  All truth-bearing checks use explicit exceptions,
not optimization-sensitive Python assertions.
