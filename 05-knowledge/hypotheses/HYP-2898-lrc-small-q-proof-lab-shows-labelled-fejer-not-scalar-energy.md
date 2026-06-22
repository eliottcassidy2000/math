---
id: HYP-2898
status: COMPUTATIONAL SIGNAL / proof-route discriminator
source: codex-2026-06-22-S111
tags: [lrc14, small-q, witness-floor, p0, additive-energy, fejer, tournament-analysis]
related:
  - HYP-2895
  - HYP-2896
  - HYP-2890
  - HYP-2889
  - HYP-2866
  - HYP-2832
  - THM-527
  - THM-534
  - OPEN-Q-108
results:
  - 04-computation/lrc_small_q_proof_lab_codex_s111.py
  - 05-knowledge/results/lrc_small_q_proof_lab_codex_s111.out
---

# HYP-2898: the smaller even-q cases select the labelled Fejer route

The exact small-q laboratory applies the current LRC(14) proof atoms to the
easier even thresholds

```text
q = 8, 10, 12, 14
threshold = 1/q
cluster GOOD_q = {maxgap(frac(e*x) : e in E) > 2/q}.
```

This is not a separate proof of the known smaller LRC cases.  It is a controlled
test of which LRC14 proof carriers are stable below q=14.

## Exact audit

Script:

- `04-computation/lrc_small_q_proof_lab_codex_s111.py`
- output: `05-knowledge/results/lrc_small_q_proof_lab_codex_s111.out`

For each even `q`, the script scans all primitive anchored bounded shapes
`E={0,...}` with `|E|=k`, `max(E)<=q`, and hard cluster sizes
`k=q/2+1,...,q-1`.  It computes exactly:

- `nu(E)=meas(GOOD_q(E))`;
- `D(E)=1-nu(E)`;
- `p0_q(E)`, the q-sector cover atom;
- `cap_q,k = min meas(G_P)` over small parts of size `q-1-k`;
- additive energy and difference-profile data.

## Findings

The Bonferroni floor is uniformly comfortable in the small-q lab:

```text
min_{q,k} [nu(AP_k)+cap_q,k-1] = 11/36 ~= 0.305556 at (q,k)=(12,7).
```

The bounded p0-cap margin is also always positive:

```text
min_{q,k} [cap_q,k - max_bounded p0_q(E)] = 1/25 at (q,k)=(10,6).
```

For q=14 itself, the smallest bounded p0 margin is the known k=8 pressure:

```text
cap_8 - p0(AP_8) = 319/5880 ~= 0.054252.
```

Consecutive/AP is the exact bounded-bank extremizer for the witness side:

```text
min nu(E) = nu(AP_k)
max D(E)  = D(AP_k)
AP difference-profile failures = 0.
```

This matches HYP-2866's denominator-prefix reading: consecutive owns the
aggregate low-denominator dense mass, even though individual buckets and local
compressions are not monotone.

The p0 side is more subtle.  AP remains cap-safe everywhere, but is not always
the p0 maximizer in the bounded bank:

```text
p0_AP_not_unique_or_not_max = 8.
```

Those non-AP p0 leaders still stay below `cap_q,k`, so they are not proof
threats.  They are warnings against replacing the cap theorem by a literal
"AP maximizes every p0 bank" statement.

The scalar additive-energy route fails immediately below q=14:

```text
scalar p0 inversions = 12706
scalar D inversions  = 12139
```

The q=14 worst rows reproduce the HYP-2890 residual-leak signal.  For example,
at `(q,k)=(14,9)` the worst scalar p0 inversion is

```text
low-energy row  (0,2,4,6,7,8,10,12,14), p0=99/245
higher-energy row (0,2,3,4,5,7,8,9,10), p0=689/5880.
```

This is the even-AP plus midpoint bridge already identified as the worst
residual-leak pattern in HYP-2890, not random noise.

## Proof-route implication

The smaller cases say:

```text
good route:
  Bonferroni floor
  + consecutive/AP nu-D extremality
  + bounded p0 cap
  + AP-facing difference-profile/Fejer majorization
  + labelled residual leak

bad route:
  scalar additive-energy monotonicity.
```

This explains why applying stronger techniques to smaller cases does not
produce a new shortcut.  It instead validates the current LRC14 architecture:
the proof must keep sector/Fourier/residue labels after extracting the positive
same-frequency additive-energy packet.

## Concurrent signal

Post-fetch incoming work HYP-2895/HYP-2896 sharpens what this lab is and is not
claiming.  HYP-2895 treats smaller `N` as a covering/tiler training atlas:
AP/Goddyn-Wong and apex-denominator sporadics are boundary exact-tilers, while
bounded covering rows have compactness margins and one-large rows go to
equidistribution.  HYP-2896 closes the one-tail zeta `-1/12` disproof branch by
explicit q-witness/binding-pair arithmetic.

The present HYP-2898 is complementary: it does not classify exact tilers or
one-tail rows.  It tests the analytic cap/floor atoms after the boundary cases
are routed away, and it says the remaining proof branch must be labelled
Fejer/residual control rather than scalar additive-energy monotonicity.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
Bonferroni_floor,
bounded_nu_minimizer,
bounded_p0_cap,
AP_difference_profile,
scalar_additive_energy.
```

Pairwise observable: exact failure count and smallest margin across the
`q=8,10,12,14` bounded banks.

Hamiltonian path:

```text
Bonferroni_floor
  > bounded_nu_minimizer
  > bounded_p0_cap
  > AP_difference_profile
  > scalar_additive_energy.
```

Assumption challenged: proving a scalar additive-energy monotonicity theorem on
smaller cases should close the q=14 cap theorem.  It does not; scalar energy is
already non-monotone in the easier q=8 and q=10 labs.  The preserved predicate is
the labelled AP-facing Fejer packet plus residual leak, not the scalar energy
value.

## Next proof target

Promote the small-q pattern into a theorem schema:

1. Prove consecutive maximizes the aggregate dense-set denominator prefixes
   for `GOOD_q` (HYP-2866 style).
2. Prove bounded p0-cap with the p0 leaders allowed to be non-AP but cap-safe.
3. Keep HYP-2890's same-frequency packet and prove the residual-leak inequality
   over labelled hidden-fold/support-cycle terms.

The small-q lab says these are not q=14 artifacts; they are the stable proof
objects.
