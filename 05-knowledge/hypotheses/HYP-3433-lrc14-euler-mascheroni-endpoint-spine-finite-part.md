---
id: HYP-3433
title: LRC14 Euler-Mascheroni endpoint-spine finite-part ledger
status: EVIDENCE / exact canonical-tail audit; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3432/HYP-3431/HYP-3430/HYP-3429 after Euler-Mascheroni prompt
tangent: T1394
technique: LTI-394
tournament_technique: LTT-294
script: 04-computation/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.py
result: 05-knowledge/results/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.out
reflection: 07-reflections/lrc14-euler-mascheroni-endpoint-spine-finite-part-codex-20260628.md
related:
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3421
  - HYP-3419
  - HYP-3417
  - HYP-3129
  - HYP-2982
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3433: LRC14 Euler-Mascheroni Endpoint-Spine Finite Part

## Claim

HYP-3430 blocks the tempting scalar shortcut: a harmonic intercept can
calibrate scale but cannot replace endpoint certificates.  HYP-3431 gives an
all-`m` corridor-fence certificate for the canonical `84m` tower, and HYP-3432
shows that harmonic wall budgets are useful only as sidecars until exact
branch/wall/interval labels return.  HYP-3429 showed that the HYP-3425
two-adic relocation survivor is often controlled by a low-rank endpoint spine.
HYP-3433 isolates the canonical tail family

```text
S_m = {1,2,...,11,13,84m}
```

and tests whether Euler-Mascheroni belongs to the LRC14 proof as a raw scalar
or as a labelled tail normalizer.

The exact audit through `m=180` says the latter.  From `m=5` onward, the best
endpoint spine is the one-wall even label

```text
E:84m
```

with exact survivor-window length

```text
1/(49m).
```

Even better, the selected component has a simple modular address.  Put

```text
a_m = ceil(48m/7).
```

Then the best checked window is exactly

```text
[(14a_m+1)/(588m), (14a_m+13)/(588m)].
```

This is one safe component of the tail half-speed `42m`: its length is
`(12/14)/(42m)=1/(49m)`.  Thus the next proof target is an address lemma:
prove that this tail component survives all fixed small-speed odd/even walls
for every `m >= 5`.

## Exact Readout

Script:

```text
04-computation/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.out
```

Aggregate readout:

```text
rows audited:                       180
canonical family:                   {1,...,11,13,84m}
eventual law begins at m:           5
failures on checked tail:           []
tail label law:                     best label = E:84m
tail length law:                    best_len = 1/(49m)
tail address law:                   a_m=ceil(48m/7)
tail address failures:              []
exact scaled law residual:          0
endpoint-weight exchange residuals: []
```

The selected best branch is both-branch through `m=20`, branch `1` from
`m=21` to `m=26`, both-branch again at `m=27`, and branch `1` for all checked
`m>=28`.

## Euler-Mascheroni Role

Summing the selected endpoint-spine lengths gives

```text
sum_{m=5}^M 1/(49m) = (H_M - H_4)/49.
```

Therefore

```text
sum_{m=5}^M 1/(49m) - log(M)/49 -> (gamma - H_4)/49.
```

At `M=180`, the finite-part value is `-0.030680458422`, while the limiting
constant is `-0.030737095274`.

This is the right way for Euler-Mascheroni to enter this lane: as the finite
part of a labelled harmonic endpoint tail.  It is not a proof of the LRC14
floor and not a legal replacement for endpoint labels, branch masks,
component geometry, owner-current sidecars, or signed-SPEC/Rprime data.

## Candidate Lemma

For the canonical tower `S_m={1,...,11,13,84m}`, prove for every `m>=5` that
the component

```text
I_m = [(14ceil(48m/7)+1)/(588m),
       (14ceil(48m/7)+13)/(588m)]
```

satisfies:

1. both endpoints are active `E:84m` walls;
2. the interior lies in the even-safe set for the even half;
3. at least one branch, usually branch `1` after the transient, avoids the
   odd two-color bad core;
4. no fixed-speed wall from `{1,...,11,13}` cuts the component.

The likely proof is a finite residue check modulo `7` plus monotone interval
separation: the address `ceil(48m/7)` keeps the component near the limiting
coordinate `u=8/49`, while the tail width shrinks like `1/m`.

## Proof Use

This packet sharpens the HYP-3431/HYP-3432/HYP-3429 lane in two ways.

First, it turns the canonical E-only exception into an explicit address
family.  A proof can try to certify `I_m` directly rather than search over all
components.

Second, it explains how harmonic and Mertens-style constants should be used
in the LRC14 program.  Infinite packet tails may need renormalized finite
parts, but those finite parts are sidecars attached to labelled endpoint
families.  This complements HYP-3432's collision ledger: the budget can rank
debt, while the labelled finite part can summarize a certified tail.  Forgetting
the label creates an illegal scalar quotient.

The next useful generalization is to seek the same pattern in other endpoint
towers:

```text
tail speed V(m)
-> tail half-speed e(m)
-> address a_m
-> one-wall E component length c/m
-> finite-part constant after subtracting c log M.
```

## Tournament Analysis

Tournament vertices are proof carriers and tail normalizers, not runners.

```text
score_hist={17:1,30:1,43:1,53:1,57:2,64:2}
directed_3cycles=0
hamiltonian_path=
  eventual_E_endpoint_harmonic_law
  -> H3429_endpoint_spine_certificate
  -> gamma_finite_part_sidecar
  -> endpoint_weight_exchange_rate
  -> Mertens_loglog_tail_normalizer
  -> total_good_measure_density_probe
  -> raw_window_count_growth
  -> raw_gamma_constant_scalar
```

## Assumption Challenge

Alternate vertices considered:

```text
runners, the m-parameter, survivor components, endpoint labels, branch masks,
harmonic tail atoms, total measures, Mertens/gamma constants, proof obligations.
```

The chosen quotient keeps endpoint labels plus the actual survivor window.  It
destroys raw runner identity and most component geometry, so it is legal only
for tail normalization and address-lemma routing, not as the final relocation
predicate.
