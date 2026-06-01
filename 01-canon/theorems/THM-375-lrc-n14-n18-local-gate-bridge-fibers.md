---
id: THM-375
name: lrc-n14-n18-local-gate-bridge-fibers
status: PROVED (exact finite computation)
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on:
  - THM-357
  - THM-360
---

# THM-375: The n=14 and n=18 local gate invoices have exact bridge fibers

## Statement

Consider the local LRC endpoint-cover problem at denominator `n`: include a
single `n`-gate speed, look only at the `2n` endpoints owned by that gate, and
ask which lower columns `1,...,n-1` strictly protect all of those endpoints.

For `n=14`, the exact minimum lower cover size is `8`.  Every minimum cover
contains the forced odd fan

```text
(1,3,5,7,9,11,13)
```

and exactly one even bridge from

```text
{2,4,6,8,10,12}.
```

For `n=18`, the exact minimum lower cover size is also `8`.  Every minimum
cover contains the forced fan

```text
(1,5,7,9,11,13,17)
```

and exactly one bridge from

```text
{6,12}.
```

Thus the local `18`-gate invoice is not larger than the local `14`-gate
invoice; it is more rigid, with a two-element bridge fiber.

## Verification Record

The statement is certified by exact rational endpoint enumeration in

```text
04-computation/lrc_n14_n18_tournament_pingpong_s481.py
```

with stored output:

```text
05-knowledge/results/lrc_n14_n18_tournament_pingpong_s481.out
```

The relevant output rows are:

```text
n=14 forced=(1,3,5,7,9,11,13), exact=8, covers=6
n=18 forced=(1,5,7,9,11,13,17), exact=8, covers=2
```

The computation enumerates the finite cover family using the same strict
endpoint-protection predicate as THM-357/THM-360.

## Proof

For a fixed `n`, there are finitely many candidate lower columns and finitely
many endpoints owned by the `n`-gate.  The script constructs the exact
incidence relation

```text
column v protects endpoint e
```

using rational arithmetic and strict containment in the forbidden interval.
It then exhaustively searches cover subsets by size.

For `n=14`, no subset of size `7` covers all gate-owned endpoints, and the
six size-`8` covers are exactly the forced fan plus one even bridge.  For
`n=18`, no subset of size `7` covers all gate-owned endpoints, and the two
size-`8` covers are exactly the forced fan plus one of `6,12`.

Since the search space and incidence relation are finite and exact, this is a
complete proof of the scoped local statements.

## Significance

This upgrades the bridge-fiber part of HYP-1942 from heuristic to exact local
fact.  The first-even cases `14=2*7` and `18=2*9` share the same cover size,
but `n=18` has a square-torsion bridge fiber rather than the six-choice even
fiber of `n=14`.

## Related

- HYP-1910
- HYP-1942
- THM-357
- THM-360
- `07-reflections/lrc-n14-n18-tournament-pingpong-s481.md`
