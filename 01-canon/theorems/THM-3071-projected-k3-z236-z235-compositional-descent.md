---
id: THM-3071
title: "Projected k3 z236 z235 compositional descent"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-lrc-z236-z235-composition-2026-08-01
depends_on:
  - THM-3061-projected-k3-z239-gap238-z237-compositional-descent
script: 04-computation/lrc14_j7_k3_z236_z235_compositional_descent_thm3071.py
output: 05-knowledge/results/lrc14_j7_k3_z236_z235_compositional_descent_thm3071.out
script_sha256: c8735a0e1328b08e98e9f27b86f901d2b0c491f832381d98f7a8a14f11e4f345
output_sha256: c04e6f57ac7025100645f1a0f546e3e5f79f6444fa2269469c09b38e746e772f
semantic_sha256: d4cfb4c7e497007f044041e9f7e56d3b109107d033421d592f7f3f74e0995f02
hash_basis: LF-normalized bytes
---

# THM-3071 -- projected k3 z236 z235 compositional descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Statement

Inside the same pinned lossless projected `k=3` necessary atlas as THM-3061,
the inherited compositor closes exactly

```text
the one occupied z_1=236 row,
all eleven occupied z_1=235 rows.                     (1)
```

These twelve rows are disjoint from THM-3061's already-closed `4+44` rows at
`z_1=239,237`. The proved consequence is therefore

```text
374780-12=374768,
projected k=3 cap: z_1<=234.                          (2)
```

The next layer `z_1=234` is occupied by `381` atlas rows and is not inspected
or closed here. This candidate has no physical-cover conclusion, no `k<=1`
conclusion, and no LRC(14) conclusion.

## 1. Pinned universe and additive set difference

The exact companion pins LF-normalized source, output, and semantic hashes of
THM-3061 as well as the projected-body atlas hash. It parses all `6,060`
structured `row=` records and independently scans their literal `;z1=...;`
tokens. Both routes give the neighboring census

```text
z_1=239:   4 rows =   2 wall +  2 order;
z_1=238:   0 rows;
z_1=237:  44 rows =  33 wall + 11 order;
z_1=236:   1 row  =   1 wall +  0 order;
z_1=235:  11 rows =   5 wall +  6 order;
z_1=234: 381 rows = 330 wall + 51 order.              (3)
```

Here `wall` means below the strict-order floor in the inherited atlas grammar,
and `order` means at or above it. The positive layers on both sides guard
against a parser miss masquerading as a gap. The twelve selected rows have
row-order SHA-256

```text
560fd7ba77b5170e234e2cde486438583791b1ab7696a516ffe23204e980349e.
```

THM-3061 closed the disjoint set `{z_1=239} union {z_1=237}` of size `48`.
The present selected set `{z_1=236} union {z_1=235}` has size `1+11=12`.
The levels are pairwise distinct and the atlas rows are keyed by their exact
body, so the ledger subtraction in `(2)` is an additive set difference, not a
recount of THM-3061's decrement.

## 2. Necessary screen and its first failure

On the twelve rows, the inherited screen gives

```text
579 states = 101 crude + 338 status + 140 residual;
 27 order states = 9 crude + 18 status + 0 residual.  (4)
```

The single `z_1=236` row has only two states and both close crudely. At
`z_1=235`, the screen is

```text
577 states = 99 crude + 338 status + 140 residual,    (5)
```

with the `140` residuals supported on exactly four wall bodies. Thus the
first failed implication is explicit: crude plus status closes `439` of
`579` states, not all `579`.

## 3. Literal terminal repair

For each of the four residual bodies the inherited THM-3061 terminal
compositor rebuilds the literal safe tables. Every body has a strictly
positive duplicate-two-high gap. The exact terminal split is

```text
zero-high scalar cases:      138;
one-high translated cases:   140;
complete-cardinality closes: 140.                    (6)
```

All `140` one-high cases pass strict complete-cell cardinality. The minimum
printed strict slack is `1`. No max-gap fallback and no unit-phase fallback
is used. The four rational gaps and per-body stage, residual, and terminal
hashes are frozen in the transcript, and no terminal survivor remains.

The logical direction is the inherited necessary one: every literal cover
maps into the relaxation and each terminal cell retains actual safe residues.
No null-set, representative-root, or measure-zero seam inference is made.

## 4. Exact evidence

The companion binds its twelve screen records and four terminal records to
checkpoint fingerprint

```text
ab43b8e07f31e2b37adcf445d2a1d6b76828dd7439e86b4f00beb7155dd473e9.
```

Fresh normal and optimized runs used disjoint empty checkpoint directories.
All `16/16` serialized envelopes were byte-identical, and both LF-normalized
transcripts were byte-identical to the stored output.  The independent hostile
audit also rederived the atlas set difference, the order/wall split, the
positive two-high gaps, the exact one-high exhaustiveness, and the minimum
strict cardinality slack `1`:

```text
python 04-computation/lrc14_j7_k3_z236_z235_compositional_descent_thm3071.py --processes 4 --checkpoint-dir <fresh-normal> --output <normal-output>
python -O 04-computation/lrc14_j7_k3_z236_z235_compositional_descent_thm3071.py --processes 4 --checkpoint-dir <fresh-optimized> --output <optimized-output>
```

The output SHA-256 is
`c04e6f57ac7025100645f1a0f546e3e5f79f6444fa2269469c09b38e746e772f`,
and the semantic SHA-256 is
`d4cfb4c7e497007f044041e9f7e56d3b109107d033421d592f7f3f74e0995f02`.
Every truth-bearing check uses explicit `require` and therefore remains active
under `python -O`.

## 5. Scope

This is a proved exact theorem inside an already-defined projected necessary
atlas.  Its only conclusion is the twelve-row projected decrement and cap
`z_1<=234`. The `381` rows at `z_1=234`, all
physical-carrier reconstruction, the `k<=1` debt, and LRC(14) remain open.
