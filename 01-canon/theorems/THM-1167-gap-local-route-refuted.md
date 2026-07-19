---
id: THM-1167
title: Exact finite witnesses refute the uniform minimum-over-all-gaps four-comb certificate
status: The displayed sufficient condition is valid and its universal premise is REFUTED by exact-rational witnesses in the stated finite family. THM-1161 strengthens the negative to an exact infinite legal family and sharp factor one. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.72; corrected codex-S75)
depends_on:
  - THM-1160
  - THM-1142
related: [THM-1137, THM-1097, THM-1161, THM-1162, MISTAKE-169]
script: 04-computation/min_config_kps_S128c72.py (+ .out)
---

# THM-1167 — the uniform all-gap certificate is refuted

THM-1160 proposed a three-tooth spacing statement inside one `k1`-gap,
with a measured margin of `3.05` against a required `1.295`.  The obvious
next step was to minimize that margin.  It does not survive.

## 1. The sufficient condition that was tested

Define

```text
W(k1,k2,k3,k4)
  = min over all k1-gaps in [0,1]
      (longest piece surviving the k2,k3,k4 teeth).
```

If `7*k4*W>1`, the four-comb theorem holds for that quadruple regardless of
which full `k1`-gap the core leaves available.  The 495-core atlas supplies a
component of length at least `1/70`, and every legal `k1>13*max(P)>=104`, so
such a component contains a full `k1`-gap.  This is a valid sufficient
condition; it would have removed core position from the problem.

## 2. Exact finite witnesses show that it fails

The exact-rational search over the stated consecutive-type family gives:

| `k1` | minimizing quadruple | `7*k4*W` |
|---:|---|---:|
| 157 | `(157,158,159,160)` | 0.79013 |
| 197 | `(197,198,199,200)` | 0.78193 |
| 237 | `(237,238,239,240)` | 0.77651 |
| 277 | `(277,278,279,280)` | 0.77267 |
| 317 | `(317,318,319,320)` | 0.76980 |

Thus the proposed universal certificate does not exceed `1.295`; it does
not even reach `1` on these witnessed rows.

## 3. Mechanism

The worst gap sits near `j=k1/4`.  THM-1142 gives the raw adjacent-tooth gap

```text
(a-jd)/(ab).
```

For consecutive killers, `d=k4-k1=3`, and at `j=k1/4` this is about
`1/(4k4)`.  Subtracting the three tooth widths leaves about
`0.77/(7k4)`.  The same linear descent that creates useful nonuniformity
also pushes the worst gap below threshold.

## 4. Correct scope

The **uniform minimum-over-all-gaps certificate** cannot prove the four-comb
theorem.  This does not exclude every argument one might call gap-local.
A successful proof must use which gaps the whole core-safe set actually
makes available.  THM-1161 strengthens this bounded negative to an exact
infinite legal family and proves sharp local factor one.  THM-1162 records
positive whole-safe-set telemetry, but only on finite killer banks.

Uniform `r=5` therefore remains open.  Do not retry a universal all-gap
spacing bonus; the live object is global selection and phase coupling among
wall-event stalks.
