---
id: HYP-2789
title: Consecutive one-far binding rows have an exact endpoint-period certificate
status: VERIFIED finite certificate/addendum to THM-563 for consec-base single-far branch; bridge target narrowed
source: codex-2026-06-21-S75
depends_on:
  - THM-563
  - HYP-2786
  - HYP-2788
  - HYP-2785
  - HYP-2784
  - HYP-2782
  - HYP-2779
  - THM-546
related:
  - HYP-2684
  - HYP-2694
  - HYP-2745
  - HYP-2772
  - HYP-2720
  - THM-557
  - OPEN-Q-108
---

# HYP-2789: One-Far Endpoint-Period Certificate

## Claim

For the consecutive binding bases

```text
B = {0,1,...,k-2},   k=8..12,
```

the single-far deviation

```text
Delta_w(B) = p0(B union {w}) - Phi(B)
```

has an exact endpoint-period certificate for every wide `w>=15`.  After
incoming `THM-563`, the periodicity lemma itself is canonical there; this
HYP is an independent exact table/check for all consecutive rows k=8..12 and
a proof-obligation refinement for the remaining non-consecutive bounded bases.

The Abel endpoint identity writes

```text
Delta_w(B) = S_B(w) / w,
```

where `S_B(w)=w*Delta_w(B)` is a signed sum of the periodized primitive `G0`
evaluated at rational one-miss arc endpoints.  Therefore `S_B(w)` is periodic
in `w` with period equal to the lcm of the endpoint denominators.  Scanning one
complete endpoint period gives:

```text
k  period  max positive S_B(w)   maxS/(15*margin)   true max Delta/margin
8     420  1                         0.360626              0.114075
9    2940  43/49                     0.442682              0.081192
10   5880  1007/980                  0.437695              0.104288
11  17640  1289/980                  0.451942              0.083918
12  17640  2034/1715                 0.310189              0.042799
```

Hence for every consecutive binding row and all `w>=15`,

```text
Delta_w(B) <= maxS_B / 15 < cap_k - Q(k-1).
```

This closes the consecutive-base one-far branch exactly.  It does not by itself
close arbitrary bounded bases; the remaining bridge is to combine `THM-563`
with `HYP-2788` and a finite period-max/slack certificate for the bounded bases
that arise after removing one wide perturbation.

After incoming `THM-563`, HYP-2789 should be read as a subcertificate rather
than the main theorem: it extends the stored consecutive table through k=11,12,
cross-checks the endpoint implementation, and records the row-pressure
tournament.  The proof route should now pass through `THM-563` plus the
`HYP-2788` near-cap dichotomy, not through an infinite signed-tail estimate.
Incoming mac-mini S6 period-max work further verifies the top dangerous k=8
and k=9 bounded bases, and the continuous period-max addendum closes dilated
consecutive bases produced by scale reduction.  Any residual finite work is
therefore non-consecutive bounded-base slack/certification, not rechecking
consecutive or dilated consecutive bases.

## Why This Does Not Contradict HYP-2785

HYP-2785 refutes a residue-only `w mod 7` table.  HYP-2789 uses a much larger
endpoint period:

```text
420, 2940, 5880, 17640, 17640.
```

The finite object is not a small residue table.  It is the exact periodic
endpoint numerator of the Abel identity.  This is why the result is compatible
with the Dedekind/equidistribution warning in HYP-2785 and the low-head
localization in HYP-2786.

## Exact Scout

Script:

```text
04-computation/lrc14_onefar_endpoint_period_codex_s75.py
```

Stored output:

```text
05-knowledge/results/lrc14_onefar_endpoint_period_codex_s75.out
```

The all-`w>=15` maxima match the short HYP-2786 low-head scout:

```text
k=8  w=21
k=9  w=68
k=10 w=22
k=11 w=16
k=12 w=71
```

The larger `maxS` rows occur at much larger `w`, but division by `w` makes them
safe.  The proof certificate therefore separates two quantities:

```text
1. endpoint numerator maxS_B, used for the rigorous all-w bound;
2. actual Delta_w maximum, used to locate the pressure row.
```

## Proof Obligation

The remaining OPEN-Q-108 bridge after `THM-563` can be phrased as:

```text
P1. treat the Abel/Dedekind endpoint-period lemma as THM-563;
P2. formalize the consecutive-base endpoint-period tables above as an exact
    k=8..12 regression/addendum to THM-563;
P3. finish any remaining finite period-max/slack ledger for non-consecutive
    bounded bases that can appear in the HYP-2788 single-perturbation regime;
P4. formalize HYP-2788's near-cap/genuine-wide dichotomy so the wide proof
    routes to the single-perturbation certificate or to the slack floor.
```

## Tournament Analysis

Vertices are the consecutive one-far binding rows `k=8..12`.  Pairwise
observable:

```text
larger maxS/(15*margin), then larger true max Delta/margin.
```

Switch/gauge: compare endpoint-period numerators before scalarizing by `1/w`.
The tournament is transitive, with pressure path:

```text
k=11 > k=9 > k=10 > k=8 > k=12.
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1}`, no directed 3-cycles, one
Hamiltonian pressure path.

## Challenged Assumption

Rejected assumption: the one-far signed bound must be an infinite analytic
Dedekind tail even for a fixed consecutive base.

Preserved predicate: the exact `Delta_w`, the cap-margin comparison, and the
HYP-2785 warning that `w mod 7` is not enough.

Destroyed information: the full Fourier phase distribution is compressed to
the Abel endpoint numerator `S_B(w)`.

What remains: extend the `THM-563` period-max certificate from the consecutive
binding bases to the finite non-consecutive bounded-base list left by HYP-2788.
