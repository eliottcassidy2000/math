---
id: THM-837
title: One order-three Hamming-five context closes by exact five-comb recursion
status: RESERVED — exact replay is being finalized; no all-context or full Hamming-five conclusion is claimed
source: codex-2026-07-15-S10 continuation
depends_on: [THM-810, THM-815, THM-820, THM-823]
related: [HYP-6820]
---

# THM-837 — one-context five-comb recursive closure

This namespace is reserved for one exact branch of the arbitrary-height
order-three Hamming-five interface.  The context is

```text
C={1,5,8,12},  b=10,  R=C union {b},
pair bits=(1,1,1),
least CRT bases={1:16,5:28,8:37,10:4,12:10} mod 39.
```

The candidate proof recursively intersects the strict-safe residual with the
five labelled danger combs.  Its state retains active endpoint runs, the
remaining CRT progressions, and the last speed; a raw global cell mask is not
sound because unused endpoints create false component cuts.  Preliminary
exact counts are 75,371 states, zero covering prefixes, and 57 nonempty
terminal rows.

This reservation is not a theorem.  Promotion requires a complete proof of
the recursive least-speed bound at every node, exact endpoint conventions,
stored deterministic replay and certificate, and a clear statement that the
other 95 directed-flag/parity contexts remain open.
