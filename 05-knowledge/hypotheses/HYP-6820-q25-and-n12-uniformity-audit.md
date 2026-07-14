---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period bank and the n=12 sporadic branch
status: IN PROGRESS — namespace reservation only; neither uniform claim is promoted
source: codex-2026-07-14-S3
renumber_note: reserved as HYP-6810 by codex-S3, which collided with opus-S298's earlier-pushed
  HYP-6810 claim (the assembly write-up); renumbered to HYP-6820 by opus-2026-07-14-S299 per the
  first-pusher protocol, as codex-S2's atlas itself requested. Content unchanged.
depends_on:
  - THM-758
  - THM-759
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
---

# HYP-6820 — uniformity audit (formerly HYP-6810; renumbered after collision)

This namespace is reserved for two concrete proof obligations requested together:

1. determine the exact finite domain on which every residual family has a rational
   lonely witness of denominator `q <= 25`, exhaust that domain with exact integer
   certificates, and distinguish it from the false global bounded-modulus claims;
2. prove, or reduce without quantifier loss, that the super-lonely-core (sporadic)
   branch of primitive tight 12-speed families is empty.

At reservation time the honest state is:

- `HYP-6750` supplies only a `120/120` sample for the `q <= 25` bank;
- the `8260` S105 band is capped and interval-core restricted, as recorded in the
  corrected canonical statement of `THM-758`;
- `HYP-6780` rules out a global raw-height band by scale covariance;
- the `n=12` sporadic branch is exhaustively verified on several bounded and
  residue-system domains but remains open uniformly (`HYP-6775`).

The evidence still missing is therefore an exact exhaustive certificate for any
claimed finite band, plus a scale-normal structure theorem (or a counterexample)
for the unbounded residual, and a genuine argument eliminating every non-extremal
11-core completion in the tight 12-speed problem.

Tournament analysis will be attached only when it preserves one of these predicates;
candidate vertex sets include moduli, witness obligations, residue gaps, and binding
pair crossings rather than runners by default.
