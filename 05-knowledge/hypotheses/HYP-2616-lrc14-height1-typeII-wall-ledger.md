---
id: HYP-2616
title: LRC(14) height-1 type-II support-six wall ledger
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S12
depends_on:
  - THM-538
  - HYP-2608
  - HYP-2612
  - HYP-2613
  - HYP-2614
related:
  - HYP-2610
  - HYP-2611
  - HYP-2615
  - OPEN-Q-108
---

# HYP-2616 - LRC(14) Height-1 Type-II Wall Ledger

## Claim

The height-1 large-offset resonances that refuted the naive span-only
Minkowski count are finite and harmless in the one-large-offset wide regime.

For the binding rows `k=8,9,10`, let `B(k)=16,15,13` be the bounded-spread
finite-check ceiling.  Consider primitive sets

```text
E = {0} union C union {M},
|C| = k-2, C subset {1,...,B(k)}, M > B(k),
```

for which there is a height-1 support-six relation touching `M`:

```text
+/- M + sum_{e in S} +/- e = 0,     |S| >= 5.
```

Every such row in the exact ledger has `meas(S7(E)) <= cap_k`, with large
margin.

## Evidence

Script:

- `04-computation/lrc14_height1_typeII_wall_ledger_codex_s12.py`
- output: `05-knowledge/results/lrc14_height1_typeII_wall_ledger_codex_s12.out`

Exact exhaustive counts:

| `k` | `B(k)` | primitive height-1 one-large rows | `M` range | over cap | worst `meas(S7)` | cap | margin | worst row |
|---:|---:|---:|---|---:|---:|---:|---:|---|
| 8 | 16 | 226046 | `17..81` | 0 | `213159/968240 ~= 0.220151` | `2243/5880 ~= 0.381463` | `93713/580944 ~= 0.161312` | `{0,1,4,7,10,13,16,19}` |
| 9 | 15 | 250264 | `16..84` | 0 | `547/1470 ~= 0.372109` | `1979/4004 ~= 0.494256` | `51353/420420 ~= 0.122147` | `{0,1,2,3,4,5,6,7,21}` |
| 10 | 13 | 54173 | `14..76` | 0 | `403/840 ~= 0.479762` | `55/91 ~= 0.604396` | `1361/10920 ~= 0.124634` | `{0,2,4,5,6,7,8,10,12,14}` |

This directly includes the HYP-2612 resonances
`21=1+2+3+4+5+6` and the k=10 height-1 signed wall around `22`; both are
resonant, but not proof-threatening.

## What This Changes

HYP-2612 showed that the theorem "large span implies large first coefficient"
is false.  HYP-2616 says the first obstruction to that theorem can be moved
into a finite ledger:

```text
bounded finite check
+ finite height-1 one-large wall ledger
+ height>=2 / multi-large / signed theta tail and signed-mass sequence spine
```

This is a real improvement to the HYP-2608 residual, because the most visible
anti-coset walls are no longer part of the analytic tail.

## What Remains Open

This does not prove LRC(14).  The ledger covers only:

- one large offset;
- coefficient height exactly `1`;
- support-six or larger relations touching that large offset.

The live proof obligations are:

1. height-2 and higher one-large walls;
2. multi-large low-height walls;
3. the relative signed support-six permanent/theta tail from HYP-2613/HYP-2614;
4. the HYP-2615 signed-mass sequence bookkeeping spine;
5. no-scale cluster-collapse rows.

## Tournament Analysis

Vertices in the script are offsets, not runners or arcs.  The observable is
height-1 wall incidence: `p(v)` counts rows in which offset `v` participates.
The switch orients `v -> w` when `p(v) >= p(w)`, with numeric order breaking
ties.  The resulting tournaments are transitive in all three rows
(`0` directed 3-cycles).

This quotient preserves the LRC predicate relevant to the support-six tail:
which offsets carry finite low-height proof obligations.  It destroys
witness-time geometry and phase location.  The challenged assumption is that
LRC tournament vertices must be runners; here the useful vertices are finite
wall obligations/offset incidences.

## Status

Partially confirmed by exhaustive exact computation of the stated finite
ledger.  Not a theorem about all wide-spread rows.
