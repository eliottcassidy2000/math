---
id: HYP-2949
title: LRC14 Baire-Haar any-angle carrier
status: PROOF-INTERFACE / measurable taut-wave guardrail, not a proof
source: codex-2026-06-24-S147
related:
  - HYP-2946
  - HYP-2945
  - HYP-2944
  - HYP-2943
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2937
  - HYP-2932
  - HYP-2832
  - THM-572
  - OPEN-Q-108
---

# HYP-2949: LRC14 Baire-Haar Any-Angle Carrier

The prompt's Borel/Baire/Haar and any-angle path-planning material should be
read as an event-algebra guardrail for LRC14.

On the compact metric group `R/Z`, Haar's theorem gives a unique normalized
translation-invariant measure.  For a finite LRC14 row `S`, each danger event

```text
D_v = {t in R/Z : ||v t|| < 1/14}
```

is a finite union of arcs.  Therefore finite Boolean combinations of the
`D_v` are Borel and Baire, with finite boundary.  The finite boundary is
Haar-null and meagre.  Endpoint choices are real proof debt, but they are not
bulk measure.

## Computation

The script
`04-computation/lrc14_baire_haar_anyangle_codex_s147.py` stores output at
`05-knowledge/results/lrc14_baire_haar_anyangle_codex_s147.out`.

It computes the regular-open danger union at threshold `1/14` for the core
rows:

```text
AP                 danger_mu=1        safe_mu=0        endpoint-only residual
GW 12->24          danger_mu=1        safe_mu=0        endpoint-only residual
near/K33 12->36    danger_mu=1259/1260 safe_mu=1/1260 positive open witness
petal 10->20       danger_mu=979/980   safe_mu=1/980  positive open witness
petal 13->26       danger_mu=181/182   safe_mu=1/182  positive open witness
```

So the Baire/Haar event split sees the same local story as the exact Farey
frontier:

```text
tight rows: endpoint-only, no open witness interval
loose rows: positive regular-open witness interval
```

The point is not that measure solves LRC14.  The point is that the event
algebra says which endpoint walls and owner/carry labels must be retained
before a scalar estimate is legal.

## Any-Angle Path-Planning Import

The existing any-angle families contribute proof-carrier roles:

```text
Field D*: interpolate continuous cost, but keep labels before scalarizing.
Theta*: direct line-of-sight shortcut; LRC analogue is direct witness shortcut.
Lazy Theta*: delay expensive visibility/owner checks until expansion.
Block A*: finite local database; LRC analogue is the AP/GW/petal/K33 atlas.
ANYA: interval nodes and taut paths; LRC analogue is safe-time intervals.
CWave: geometric wavefront primitives; LRC analogue is Haar wavefront arcs.
```

The creative sixth carrier is:

```text
Haar-Baire Taut Wave*
```

with state

```text
(regular-open Baire set U,
 Haar mass mu(U),
 finite boundary debt,
 C27/K33 owner label,
 parent taut wall).
```

Its line-of-sight check is not `is this segment visible?` but:

```text
does a positive regular-open Haar/Baire witness survive after deleting
the known C27/K33 obstacle walls?
```

This matches the current LRC14 proof shape: AP/GW are endpoint-only residuals,
while near/K33 and petal rows expose positive open witness intervals.  The
state-lift branch must preserve the owner of the boundary wall that kills the
open component.

## Haar/Baire Guardrails

1. **Haar theorem.**  On `R/Z` and finite products of `R/Z`, normalized Haar
   measure is the invariant mass behind `p0`, `GOOD`, `denseSet`, and
   `safeSet` events.
2. **Borel sets.**  Finite arc unions are Borel, so exact event definitions can
   be formalized without measurability ambiguity.
3. **Baire sets.**  In the metric circle, these finite arc events also have
   the Baire property; finite boundaries are meagre.
4. **Anti-collapse.**  Positive Haar measure, comeagreness, endpoint existence,
   and exact LRC loneliness are different predicates.  They agree only after
   regular-open and endpoint labels are retained.

## Tournament Analysis

Tournament vertices are carriers:

```text
Field D*
Theta* / Lazy Theta*
Block A*
ANYA
CWave
Haar-Baire Taut Wave*
```

The pairwise observable is:

```text
measurability retention,
interval-node strength,
taut-wall retention,
finite-atlas fit,
dynamic-update fit,
scalar-decoy resistance.
```

The conservative order is transitive:

```text
Haar-Baire Taut Wave*
> CWave
> ANYA
> Block A*
> Theta* / Lazy Theta*
> Field D*.
```

## Proof Targets

1. **Regular-open endpoint lemma.**  Replace finite LRC14 event sets by
   regular-open interiors plus finite boundary debt; prove Haar mass and Baire
   open components are unchanged.
2. **Taut-wall ANYA lemma.**  If a loose row has a positive witness interval,
   its endpoints should wrap tautly around named C27 shell or K33 owner walls.
3. **Lazy visibility lemma.**  Delay owner/carry checks until a candidate
   witness interval is expanded, as Lazy Theta* delays line-of-sight checks.
4. **Block database lemma.**  Cache the finite AP/GW/petal/K33 local atlas as
   a proof database.
5. **Haar wavefront lemma.**  Formalize `p0`, `GOOD`, `denseSet`, and `safeSet`
   as Haar-measurable regular-open events modulo finite boundary debt.
6. **Baire obstruction lemma.**  Separate endpoint-only tight residuals from
   positive open witnesses by Baire category before scalar measure estimates.

## Proof-Route Use

HYP-2949 suggests the next LRC14 route should not be another scalar count.  It
should be a measurable taut-wave proof search:

```text
exact M/Farey branch
-> regular-open Haar/Baire event carrier
-> C27/K33 wall owner
-> finite endpoint debt
-> HYP-2908/THM-572 state lift if the wall is nonunit/K33.
```

This is not a proof of LRC14.  It is a carrier discipline for deciding what
must be preserved before measure, category, or any-angle path analogies can
enter the proof.
