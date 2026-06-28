---
id: HYP-3247
title: LRC14 shell-lag commutator and contact-support sidecar
status: EVIDENCE / exact bounded-bank scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1345
technique: LTI-345
tournament_technique: LTT-245
script: 04-computation/lrc14_shell_lag_contact_sidecar_codex_20260628.py
result: 05-knowledge/results/lrc14_shell_lag_contact_sidecar_codex_20260628.out
reflection: 07-reflections/lrc14-shell-lag-contact-sidecar-codex-20260628.md
related:
  - HYP-3245
  - HYP-3246
  - HYP-3244
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3239
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3228
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3138
  - OPEN-Q-108
---

# HYP-3247: LRC14 Shell-Lag Commutator and Contact-Support Sidecar

## Claim

HYP-3245 treats AP as the triangular ordinary-autocorrelation law and records
outward lag transport on the HYP-3202 trap boundary.  HYP-3228 gives the exact
shell functional `10q0+q3+10q6`.  HYP-3247 asks whether those two projections
commute on the bounded k=8 bank.

They do not.  Ordinary support autocorrelation is a lossy quotient for the
shell packet, and even the residue histogram modulo `7` is still not enough.
The missing bounded-bank repair is an ordered contact-support sidecar:

```text
contact_support(E) = positions of the non-unit gaps in the anchored gap word.
```

The new bounded-bank lesson is:

```text
ordinary lag profile              forgets too much shell data
ordinary lag + residue histogram  still forgets shell data
ordinary lag + residue histogram + ordered contact support
                                 repairs the shell packet on the bounded bank
```

This is the controlled-forgetting form of the HYP-3245/HYP-3228 interface.
Shell magic does not only see how much pair mass moved outward; it also sees
where the long gaps sit in the ordered endpoint/contact word.

## Exact Scout

The scout reuses the exact `3432`-row anchored bounded k=8 bank
`E={0} ∪ A, A⊂{1,...,14}, |A|=7` and imports `row_moments` from HYP-3200's
bank script.  For each row it records:

```text
ordinary support autocorrelation
residue histogram mod 7
residue word mod 7
gap multiset
contact support = positions of non-unit gaps
shell magic = 10q0+q3+10q6
```

Exact fiber counts:

```text
ac fibers                  = 1747
mixed magic on ac fibers   = 1677

ac + hist fibers                 = 3370
mixed magic on ac+hist fibers    = 62

ac + hist + contact_support fibers = 3432
mixed magic there                  = 0
```

So residue data repairs most of the loss, but not all of it.  There remains a
finite family of genuinely shell-visible commutators.

The structural comparison inside those `62` mixed `ac+hist` fibers is the key
new fact:

```text
adding gap_multiset      leaves all 62 fibers unresolved
adding contact_support   resolves all 62 fibers
```

Thus the hidden coordinate is positional, not purely metric.  Gap sizes alone
do not repair the shell packet; the positions of the non-unit gaps do.

## Sharp Collision Pair

The scout exhibits a concrete bounded-bank collision:

```text
E_A = (0,1,2,3,4,12,13,14)
E_B = (0,1,2,10,11,12,13,14)
```

They have:

```text
same ordinary lag profile
same residue histogram mod 7
same residue word mod 7 = (0,1,2,3,4,5,6,0)
same q6 = 1/98
```

but different shell data:

```text
magic(E_A) = 211/98
magic(E_B) = 71394/35035
q0+q6(E_A) = 1507/7644
q0+q6(E_B) = 77309/420420
q3(E_A)    = 347/1911
q3(E_B)    = 41819/210210
```

The visible difference is only the location of the long gap `8` in the
anchored gap word:

```text
gaps(E_A) = (1,1,1,1,8,1,1), contact_support(E_A) = (4,)
gaps(E_B) = (1,1,8,1,1,1,1), contact_support(E_B) = (2,)
```

This is the exact shell-lag commutator.  The lag profile and the residue
profile are unchanged, but the endpoint/contact placement moves and the shell
functional responds.

## Proof-Frontier Use

This bounded-bank result sharpens several nearby packets:

- HYP-3245: outward lag transport is useful, but not terminal; it needs a
  contact-sidecar.
- HYP-3228: shell magic is not a pure lag functional.
- HYP-3204: the ordered-tail exchange route was already positional; HYP-3247
  identifies one exact bounded-bank positional payload.
- HYP-3243: the repair wants a circle endpoint arrangement cell rather than a
  scalar lag count.
- HYP-3244: the repair is a controlled-forgetting sidecar, not an analogy.

The next theorem-facing statement should therefore be:

```text
shell magic = lag transport + ordered contact-support sidecar + named repair
```

with the repair routed to endpoint cells, Toeplitz/Green trap chambers,
tiling/half-tiling descent, or named residual debt.

## Tournament Analysis

The proof-carrier tournament is transitive:

```text
half_tiling_descent_sidecar
-> circle_endpoint_arrangement_cell
-> lag_plus_contact_word
-> lag_plus_contact_support
-> lag_plus_residue_histogram
-> lag_plus_gap_multiset
-> lag_profile_only
```

The guiding rule is that the better carrier is the one that preserves more of
the shell/contact/endpoint payload while still staying finite and theorem-facing.
