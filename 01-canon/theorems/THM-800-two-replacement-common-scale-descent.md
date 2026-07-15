---
id: THM-800
title: Two-replacement common-scale descent in the full-residue chart
status: CLAIMED / PROOF AUDIT IN PROGRESS — deck-capacity reduction derived; normalized scale-one closure not yet proved
source: codex-2026-07-15-S10 (Hamming-two continuation)
depends_on: []
related: [THM-769, THM-770, THM-795, HYP-6775, HYP-6800, HYP-6820]
verification: forthcoming exact Fraction replay
---

# THM-800 — Two-replacement common-scale descent in the full-residue chart

This file reserves THM-800 while the proof and exact replay are audited.  The
claimed theorem number is intentionally backed by a narrow, explicit current
claim rather than by the stronger still-open Hamming-two rigidity statement.

Let `13 not divide c`, let `r,s` be distinct members of `[12]`, and put

```text
A=(c[12]\{cr,cs}) union {w_r,w_s},
w_r=cr (mod 13),             w_s=cs (mod 13).
```

For each replacement define its lifted-splice deck order

```text
D_i=c/gcd(c,w_i).
```

## Claimed scale-descent statement

If `M(A)=1/13`, then

```text
c divides w_r  and  c divides w_s.                       (1)
```

Consequently `A=cB`, where `B` is a labelled Hamming-two perturbation of
`[12]`.  Thus uniform Hamming-two rigidity at every AP scale reduces exactly
to the scale-one chart; no mixed deck order and no genuine higher deck order
can support a tight packet.

## Derived proof spine awaiting final audit

1. Over any missing-owner splice family, a replacement of deck order `D`
   occupies at most

   ```text
   gcd(c,w) * (floor(2D/13)+1)
   ```

   sheets.  If two replacements cover all `c` sheets, at least one deck order
   is `1` or `2`.  If neither is `1`, restriction to the uncovered parity
   class forces both orders to be `2`.

2. In the order-two/order-two chart, the own replacement covers the parity on
   which the half-step is zero.  A cross replacement can cover the other
   parity at the `r` splice only when

   ```text
   s/r in {6,7}  (mod 13).
   ```

   Reversing `r,s` would also require its inverse to lie in `{6,7}`, but
   `{6^(-1),7^(-1)}={11,2}`.  Hence this chart cannot cover both missing-owner
   splice families.  The complementary-owner case is to be checked separately
   in the replay because the core-safe germ becomes two-sided.

3. Once one replacement has order one, its cross phase at a separated missing
   owner is constant and outside the closed `1/13` danger arc.  The other
   replacement must cover that entire sheet family and hence also has order
   one.  For complementary missing owners, use the core's two-sided safe germ:
   the descended replacement leaves one side of its boundary tooth uncovered,
   forcing the second replacement to be dangerous at every lifted splice and
   hence to have order one.

The remaining work before changing the status to `PROVED` is to formalize the
closed-germ argument without orientation ambiguity, replay every residue/parity
case exactly, and attach Tournament Analysis to the sheet--replacement
incidence carrier.  The stronger assertion that the resulting scale-one
Hamming-two packet must equal `[12]` is **not** yet part of this theorem.

