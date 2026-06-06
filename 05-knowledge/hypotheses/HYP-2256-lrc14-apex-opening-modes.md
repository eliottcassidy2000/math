---
id: HYP-2256
status: OPEN proof-route refinement for LRC14
source: codex-2026-06-06-S678b
related:
  - HYP-2253
  - HYP-2255
  - HYP-2254
  - HYP-2252
  - HYP-2241
  - HYP-2168
  - HYP-2167
  - HYP-2164
---

# HYP-2256: LRC14 Primitive Apex Debt Splits Into Local Opening Modes

## Claim

HYP-2253's primitive apex-debt target

```text
primitive apex debt => p_0(V,1/14)>0
```

should be proved by a local opening-mode package rather than by a single
global covering estimate.

For an apex speed `a` divisible by `14`, four local certificate modes appear:

1. **Clock shutter:** delete the apex speed; at some unit clock `t=j/14`, the
   non-apex row has a one-sided safe cone wider than the apex danger shutter
   radius `1/(14a)`.
2. **Apex-free side door:** the full row already has a positive safe interval
   whose two endpoint owners are non-apex speeds.
3. **One-apex hinge:** the full row has a positive safe interval with exactly
   one apex-owned endpoint and one non-apex endpoint.
4. **Apex-period aperture:** both endpoints are apex-owned, but the non-apex row
   is safe through a gap between neighbouring apex danger slits.

So the proof route becomes:

```text
apex debt
  -> deleted-clock cone OR endpoint-owner interval
  -> one of four local opening lemmas
  -> p_0 > 0.
```

## Exact Evidence

S678b adds `04-computation/lrc14_apex_opening_modes_s678b.py` with stored output
in `05-knowledge/results/lrc14_apex_opening_modes_s678b.out`.

The script reuses the exact S677 coherent carry atlas and audits all primitive
apex-debt rows:

```text
primitive apex-debt rows: 414
p_0-positive rows:       414
p_0-wall rows:             0
```

Using the priority order `clock_shutter`, `apex_free_side_door`,
`one_apex_hinge`, `apex_period_aperture`, the certificate histogram is:

```text
clock_shutter             235
apex_free_side_door       106
one_apex_hinge             73
apex_period_aperture        0
```

The independent endpoint-owner classification, before clock priority, is:

```text
apex_free                 230
one_apex                  182
apex_both_only              2
```

The two `apex_both_only` rows are not failures; they are already caught by the
clock-shutter lemma.  In the raw safe-interval population the endpoint-owner
histogram is:

```text
one_apex      5410
apex_free     3502
apex_both     2494
```

The smallest clock surplus in the audit is `1/2352`, at
`Vstar/interval_block/[0,4]:h1`, with apex `28`, clock `11/14`, shutter
`1/392`, and right cone `1/336`.

## Interpretation

This turns HYP-2241's owner-private deletion bit into a boundary-owner question:
do not ask whether the apex is globally harmless.  Ask which boundary of a
specific safe interval it owns.

The local modes also explain why the naive `n`-clock proof is incomplete.  The
clock-shutter lemma covers more than half the primitive apex-debt probes, but
side carries can erase the literal unit-clock cone.  The endpoint-owner
side-door and hinge modes are the repair.

## Tournament Analysis

Vertices are proof-certificate modes, not runners.  Candidate vertex sets
considered included runners, unit clocks, apex-deletion rows, safe intervals,
endpoint owners, carry congruence sites, and proof obligations.

The certificate-mode tournament uses observable

```text
(locality, exactness, coverage, deletion-owner, quotient-stability, action)
```

and is transitive:

```text
clock_shutter
> apex_free_side_door
> one_apex_hinge
> apex_period_aperture
> raw_p0_scalar
> raw_res27_shadow
```

with `directed_3cycles=0`, singleton SCCs, and one Hamiltonian path.

## Next Lemma Target

Prove the four local lemmas abstractly:

```text
L1 clock shutter:
  if a deleted-apex unit-clock cone width exceeds 1/(14a), then p_0>0.

L2 side door:
  if a safe interval has non-apex owners on both endpoints, it survives apex debt.

L3 hinge:
  if one endpoint is apex-owned and the other is non-apex-owned, the crossing
  is transversal and opens an interval.

L4 aperture:
  if both endpoints are apex-owned and non-apex depth is zero between them,
  the apex period itself supplies a positive aperture.
```

Then prove that every normalized primitive apex-debt row falls into at least
one of these modes, using the `Res_27` carry-owner quotient plus endpoint-owner
labels.
