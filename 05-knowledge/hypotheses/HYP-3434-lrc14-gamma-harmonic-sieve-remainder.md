---
id: HYP-3434
title: LRC14 Euler-Mascheroni harmonic-sieve overlap remainder
status: EVIDENCE / exact harmonic-sieve remainder audit; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3433/HYP-3432/HYP-3431/HYP-3430/HYP-3429/HYP-3426
tangent: T1395
technique: LTI-395
tournament_technique: LTT-295
script: 04-computation/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.py
result: 05-knowledge/results/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.out
reflection: 07-reflections/lrc14-overlap-tax-harmonic-sieve-remainder-codex-20260628.md
related:
  - HYP-3433
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3418
  - HYP-3417
  - HYP-3415
  - HYP-3401
  - HYP-3236
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3434: LRC14 Gamma Harmonic-Sieve Remainder

## Claim

HYP-3430 blocks the scalar shortcut: Euler-Mascheroni / harmonic intercepts
calibrate denominator tails but do not determine endpoint-spine class.
HYP-3431 proves the canonical `84m` corridor-fence case, HYP-3432 shows
reciprocal endpoint budgets are only wall-ranking sidecars, and HYP-3433
extracts the labelled endpoint-tail finite part for the canonical tower.
HYP-3426 reduces the two-branch covering-floor target to the one-branch
piercing statement

```text
branch0 = E_safe minus B0_odd is nonempty.
```

HYP-3434 identifies the exact first-sieve remainder that survives those
predecessors:

```text
|branch0|
  = |E_safe| - sum_o |E_safe cap B0_o|
    + (sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|).
```

The first two terms are the naive union-bound slack.  The last term is the
overlap tax.  The finite theorem target is now sharper:

```text
either naive_slack >= 0,
or endpoint-spine structure certifies overlap_tax > -naive_slack.
```

This is the concrete compression failure.  Set union is associative, but the
scalar compression that replaces a union by a sum of masses is not a legal
homomorphism unless the overlap sidecar is retained.

## Euler-Mascheroni Role

Euler-Mascheroni is useful here only as a denominator-prefix calibration:

```text
H_N = log N + gamma + error_N
sum_{odd n <= N} 1/n
  = 1/2 (log N + gamma + log 2) + odd_error_N.
```

Thus gamma is a reminder that harmonic endpoint budgets have a bounded
discrete-to-log remainder.  It is not the proof.  The proof object is the
exact rational overlap tax on `E_safe`.

## Exact Readout

Script:

```text
04-computation/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.out
```

The script audits the same `150`-row endpoint-spine bank used by HYP-3429.

Aggregate readout:

```text
exact sieve identity rows:            150/150
positive one-branch survivor rows:    150/150
naive slack positive rows:            91/150
naive slack negative rows:            59/150
positive overlap-tax rows:            150/150
negative-slack rows rescued by tax:   59/59
tightest tax/deficit rescue:          1.090875 (covering_AP_with_84)
odd-prefix gamma residual range:      [0.002106737, 0.044081619]
```

The tightest branch and tightest rescue is still the canonical covering row
`{1..11,13,84}`:

```text
|E| = 107/245
single_sum = 31259/63063
union = 9068/21021
branch0 = 563/105105
naive_slack = -18586/315315
overlap_tax = 4055/63063
tax/deficit = 1.090875
```

The worst naive slack row is `random_covering_12`; it is not a counterexample,
because the overlap tax is much larger than the deficit:

```text
naive_slack = -941318888891/12839701214340
overlap_tax = 105709968060469/526427749787940
tax/deficit = 2.739021
```

## Function And Information Reading

View every proposed invariant as a function.

```text
row
  -> exact interval family (E_safe, B0_o)
  -> compressed scalar packet.
```

The LRC predicate is not a function of the scalar packet `|E|` and
`sum_o |E cap B0_o|` alone.  The missing coordinate is the overlap tax.  In
information-theory language, the scalar packet has destroyed mutual
intersection information.  In algebra language, the compression fails to
respect the union operation after replacing sets by additive masses.

This also records a positive/negative duality: the odd bad intervals are the
negative cover, while the overlap tax is the negative-negative correction that
returns positive survivor mass.

## Graph And Topology Route

The next structural model should build a conductance or cut graph whose
vertices are `E_safe` components and odd bad intervals, with incidence weights
given by exact intersections.  Then:

```text
overlap_tax = collision / shared-incidence mass inside E_safe.
```

HYP-3429 says the relevant survivor windows have endpoint-spine rank at most
`2`; that suggests the graph has a low-rank local Menger cut, bounded
treewidth, or a Green-current certificate that lower-bounds the overlap tax.
Algebraic connectivity should be tested only after retaining negative leakage
and endpoint labels, following the HYP-3236 conductance guardrail.

## Candidate Lemma

For every primitive covering row `S = O union 2E`, exactly one of the following
proof exits should hold:

1. `|E_safe| - sum_o |E_safe cap B0_o| >= 0`;
2. a rank-`<=2` endpoint spine certifies
   `sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|`
   exceeds the naive deficit;
3. a named owner-current, exact-period, state-lift, or two-adic loss debt is
   emitted.

This lemma would join HYP-3433's endpoint-tail finite part, HYP-3432's
wall-budget sidecar, HYP-3431's corridor fence, HYP-3430's scalar firewall,
HYP-3426's one-branch mirror reduction, HYP-3429's endpoint-spine certificate,
and HYP-3428's loss ledger.

## Tournament Analysis

Tournament vertices are proof carriers and compression residues, not runners.

```text
score_hist={13:1, 29:1, 47:1, 61:1, 62:1, 65:1, 66:1}
directed_3cycles=0
hamiltonian_path=
  exact_overlap_tax_sieve
  -> endpoint_spine_tax_certificate
  -> conductance_graph_overlap_flow
  -> two_adic_even_child_loss_ledger
  -> odd_prefix_gamma_cap
  -> raw_harmonic_union_bound
  -> scalar_gamma_slogan
```

## Guardrail

HYP-3434 does not prove LRC14 and does not promote gamma, Mertens tails, or
harmonic sums to theorem vertices.  The constants are calibration shadows.
The load-bearing data are exact interval intersections, endpoint-spine labels,
and the overlap-tax inequality.
