---
id: THM-836
title: Thin two-sheet owner shells force a modular packing obstruction
status: RESERVED — candidate exact reduction; proof and replay are in progress and no branch deduction may use this file
source: codex-2026-07-15-S10 continuation
depends_on: [THM-797, THM-824]
related: [THM-774, THM-803, HYP-6820]
---

# THM-836 — thin two-sheet owner-shell packing obstruction

This namespace is reserved for an exact candidate consequence of THM-797's
signed-wall exit obligations, with THM-824's centre-membership guard retained.
For the common folded pair `(13d,5d)`, write the maximum speed as `B` and put

```text
s=2B-13d,
delta(d)=min({6d}_13,13-{6d}_13).
```

The candidate reduction is

```text
delta(d) <= 11s(13d+s)/(2(117d+11s)).
```

It would imply that the shell `s=1` is empty, and that the shell `s=3` can
survive only for `d=2,11 mod 13`, with the two distinguished exit owners
consecutive.  The intended proof passes through the exact directional
owner-exit coefficient multiset

```text
(9,20,31,42,53,64,75,86,97,108).
```

This is a reservation, not a theorem.  Promotion requires the complete
signed-wall derivation, strict endpoint audit, modular-packing proof,
fraction-exact replay over a nontrivial finite bank, and Tournament Analysis
or an explicit explanation of why the theorem-bearing object is not binary.
