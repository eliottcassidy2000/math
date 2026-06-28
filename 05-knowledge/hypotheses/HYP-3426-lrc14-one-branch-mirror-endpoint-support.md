---
id: HYP-3426
title: LRC14 one-branch mirror reduction and endpoint-owner support certificates
status: EVIDENCE / exact finite-ruler mirror and endpoint-support audit; not an LRC14 proof
source: codex-2026-06-28
tangent: T1387
technique: LTI-387
tournament_technique: LTT-287
script: 04-computation/lrc14_one_branch_mirror_endpoint_support_codex_20260628.py
result: 05-knowledge/results/lrc14_one_branch_mirror_endpoint_support_codex_20260628.out
reflection: 07-reflections/lrc14-one-branch-mirror-endpoint-support-codex-20260628.md
related:
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3419
  - HYP-3417
  - HYP-3415
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3426: LRC14 One-Branch Mirror / Endpoint-Support Certificates

## Claim

HYP-3425's two-color obstruction has a mirror reduction.  For

```text
S = O union 2E,     u = 2t,
```

branch `0` uses `t=u/2`, while branch `1` uses `t=(u+1)/2`.  The involution

```text
u -> 1-u
```

preserves `E_safe` and sends branch-0 survivors to branch-1 survivors.  Thus
the two-branch target can be attacked through the one-color interval statement

```text
E_safe is not contained in B0_odd,
```

where `B0_odd` is the union of odd near-integer bad intervals in `o*u/2`.
Once branch `0` survives, branch `1` survives by mirror symmetry.

## Exact Readout

Script:

```text
04-computation/lrc14_one_branch_mirror_endpoint_support_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_one_branch_mirror_endpoint_support_codex_20260628.out
```

The audit covers `162` rows: the HYP-3425 bank, the canonical `84m` tower
through `m=40`, extended two-tail rows, and `60` deterministic primitive
covering rows.

```text
mirror identity branch1=mirror(branch0): 162/162
branch0 measure equals branch1:         162/162
positive one-branch survivor:           162/162
selected branch0 score >= 1/14:         162/162
endpoint-labelled survivors:            162/162
endpoint support histogram:             {1: 353, 2: 13103, 3: 72}
max endpoint support size:              3
```

The tight canonical row

```text
S = {1,2,3,4,5,6,7,8,9,10,11,13,84}
```

has branch-0 and branch-1 measures both equal to `563/105105`.  Its four
branch-0 survivor intervals have endpoint supports

```text
{13,42}, {11,42}, {5,42}, {7,42}.
```

The smallest branch-0 measure remains the tight canonical row.  The largest
observed endpoint support size is `3`, and only `72` survivor intervals in the
whole audit need support size `3`.

## Proof Target

The next finite lemma can be sharpened from HYP-3425's two-color statement to
the one-branch cover statement:

```text
For every primitive covering 13-row S=O union 2E,
E_safe is not contained in B0_odd.
```

Equivalently, at least one component of `E_safe` has a positive subinterval
whose endpoints are owned by the finite arrangement of even-safe walls and
odd near-integer walls.  The endpoint audit suggests a local certificate
theorem:

```text
Every surviving one-branch interval has an endpoint-owner certificate using
at most three speed owners.
```

This is not yet proved uniformly, but it gives a smaller proof search space:
owner triples, not full rows.

## Guardrail

Do not confuse the one-branch statement with an immediate proof of LRC14.  It
is still the hard interval-piercing lemma.  The gain is structural: branch
selection is removed by the mirror involution, and the surviving intervals are
controlled by small endpoint-owner certificates.  Owner-current labels from
HYP-3417/HYP-3419 remain exception routers, not substitutes for the interval
proof.

## Tournament Analysis

Tournament vertices are proof obligations and retained information channels,
not runners or raw scalar masses.

```text
score_hist={32:1, 55:1, 57:1, 60:1, 61:1, 62:1, 65:1}
directed_3cycles=0
hamiltonian_path=
  one_branch_interval_piercing_lemma
  -> endpoint_owner_triple_certificate
  -> mirror_involution_branch_equivalence
  -> component_local_cover_decomposition
  -> two_color_bad_core_identity
  -> owner_current_exception_router
  -> raw_branch_mass_lower_bound
```

The leading proof route is one-branch interval piercing, with endpoint-owner
triples as the finite certificate layer.
