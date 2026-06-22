---
id: HYP-2889
status: EVIDENCE / proof-target; exact finite scout, no LRC14 proof claimed
source: codex-2026-06-22-S103
tags: [lrc14, additive-energy, fejer, majorization, p0, Ly, interval-extremality, signed-tail, tournament-analysis]
related:
  - HYP-2885
  - HYP-2886
  - HYP-2890
  - HYP-+2888
  - HYP-2873
  - THM-534
  - OPEN-Q-108
results:
  - 04-computation/lrc_additive_energy_majorization_codex_s103.py
  - 05-knowledge/results/lrc_additive_energy_majorization_codex_s103.out
  - 04-computation/lrc_additive_moment_coefficients_kps.py
  - 05-knowledge/results/lrc_additive_moment_coefficients_kps.out
  - 04-computation/lrc14_gamma_frequency_tail_codex_s104.py
  - 05-knowledge/results/lrc14_gamma_frequency_tail_codex_s104.out
---

# HYP-2889: additive-energy majorization is anchored AP-facing, not monotone

HYP-2885 correctly identifies additive energy / Fejer concentration as the
right extremal carrier for the LRC14 cap branch, but the carrier must be used in
an AP-facing way.  The exact S103 scout refutes three tempting stronger
statements:

```text
scalar A(E) monotone => p0(E) monotone;
pairwise difference-profile majorization => p0(E) monotone;
one-step compression toward an AP => p0(E) monotone.
```

What survives is sharper:

```text
the interval/AP difference profile majorizes every tested row, and the AP is
still the p0 and L_y maximizer in the same exact banks.
```

This also fits the concurrent HYP-+2888 boundary/rational-witness result and
the S39 refinement: maximal AP-like energy reaches exact unsafe coverage `1`,
but exact tiling is scaling-invariant and anchored, not translation-invariant.
Length-13 AP translates have the same additive energy as `{1,...,13}` but
positive safe measure, while `d*{1,...,13}` has the explicit boundary witness
`t=1/(14d)`.  Therefore the cap branch must prevent non-AP rows from exceeding
the anchored AP boundary value; it should not try to prove a universal
positive safe-measure floor at the strict threshold.

Thus the proof route should not be `L_y <= G(A(E))` for a scalar monotone
function `G`.  It should be:

```text
anchored AP-facing Fejer majorization
  + signed sector/Fourier remainder lemma
  => L_y(E) <= L_y(AP_k)
  => p0(E) <= cap_k by THM-534.
```

## Exact S103 evidence

Script:

- `04-computation/lrc_additive_energy_majorization_codex_s103.py`
- output: `05-knowledge/results/lrc_additive_energy_majorization_codex_s103.out`

Exact primitive banks:

```text
k=8:  E=(0 plus 7 from 1..13), rows=1716
k=9:  E=(0 plus 8 from 1..13), rows=1287
k=10: E=(0 plus 9 from 1..12), rows=220
```

For all three banks:

```text
AP diff-profile majorization failures = 0
rows with p0 > p0(AP)                 = 0
rows with L_y > L_y(AP)               = 0
```

The AP profiles are the expected Fejer profiles:

```text
k=8:  (7,6,5,4,3,2,1)
k=9:  (8,7,6,5,4,3,2,1)
k=10: (9,8,7,6,5,4,3,2,1)
```

Here the profile is the sorted list of positive-difference multiplicities

```text
d_E(h) = #{ unordered pairs {a,b}: |a-b| = h }.
```

This is the finite Fejer object because

```text
|Ehat(x)|^2 = k + 2 sum_h d_E(h) cos(2*pi*h*x).
```

The word "anchored" matters.  Difference profiles and additive energy are
translation-invariant, but LRC coverage at threshold is not: the exact tilers
are the scaled consecutive multiples `d*{1,...,13}`, not arbitrary AP
translates with the same energy.  The Fejer layer can control the interval
profile, but the sector/origin labels must remain available for the final
coverage inequality.

## Refuted shortcuts

Scalar additive energy is a useful correlation heuristic but not a monotone
proof order.  The exact banks have many strict inversions:

```text
scalar A monotonicity inversions for p0:
  k=8:  1686
  k=9:  1254
  k=10: 197
  total: 3137
```

Worst displayed examples:

```text
k=8:
lower-A  (0,2,4,6,7,8,10,12) A=272 p0=188/735
higher-A (0,2,3,4,5,7,8,9)  A=280 p0=341/5880
gap = 1163/5880

k=9:
lower-A  (0,1,2,3,4,5,6,7,12) A=401 p0=1093/2940
higher-A (0,2,3,4,5,7,8,9,10) A=405 p0=689/5880
gap = 499/1960
```

Full pairwise difference-profile majorization is also not a valid monotone
order for `p0`:

```text
pairwise profile-majorization => p0 monotone violations:
  k=8:  351261
  k=9:  199693
  k=10: 5798
```

This is the important correction.  Difference-profile majorization is useful
only with the interval/AP as the top profile, not as an arbitrary comparison
principle between two non-AP rows.

Local compression is also too naive.  Among one-step moves `e -> e-1`, S103
finds many cases where the difference profile improves but `p0` decreases:

```text
profile-up but p0-down moves:
  k=8:  572
  k=9:  520
  k=10: 85
```

Example:

```text
(0,1,2,3,4,5,6,7,12) -> (0,1,2,3,4,5,6,7,11)
A: 401 -> 417
p0: 1093/2940 -> 493/1470
```

So a proof by local compression must carry sector labels or Fourier signs; the
unlabelled Fejer profile alone loses the information that controls `p0`.

The incoming KPS S31l moment-coefficient check explains why scalarization fails
without contradicting the additive-energy mechanism.  The frequency-1
`s=2` coefficient is positive and dominant, but higher additive-moment
coefficients are signed: for example k=9 has `s=2:+8.910e-04`,
`s=3:-7.138e-05`, `s=4:-4.448e-05`, and k=12 has negative coefficients from
`s=4` onward.  Thus additive energy is the leading Fejer term, while the tail is
a signed cancellation problem.  The right analogy is the H-max Jensen/convexity
proof, not a term-by-term positive Savchenko-style comparison.

## Proof target

The next lemma should be stated as an anchored AP-facing signed-majorization
theorem.

Candidate form:

```text
For k=8,9,10 (and then the slack k>=11 cases), write

  L_y(E) - L_y(AP_k) = FejerDef(E) + SectorRem(E).

FejerDef(E) <= 0 follows from AP difference-profile majorization.
SectorRem(E) <= -FejerDef(E) or SectorRem(E) <= 0 follows from the labelled
sector/Fourier packet structure in the THM-534 certificate.
```

This is deliberately weaker than pairwise monotonicity and local compression.
It asks only that the AP top profile remains the top profile after the
specific signed sector kernel defining `L_y` is applied.

A plausible proof split:

1. Prove the Pollard/Karamata-style interval layer:

```text
d_AP^down majorizes d_E^down for every k-set E.
```

2. Express the THM-534 `L_y` certificate as a low-support sector Fourier
functional while preserving the anchor/origin labels.  The positive Fejer part
is controlled by the previous layer, but AP translates are not equivalent at
the strict LRC threshold.

3. Route the signed remainder by labelled sector packets:

```text
low relation-depth packets -> finite AP/Freiman atlas;
high relation-depth packets -> HYP-2636/HYP-2884 Abel or L2 cancellation;
exact-period witness packets -> HYP-2886 if the finite denominator route is used.
```

## Tournament Analysis

Vertices were proof carriers, not runners:

```text
AP_diff_majorization,
Ly_certificate,
pairwise_diff_majorization,
local_profile_compression,
scalar_additive_energy.
```

The pairwise observable was implication reliability on exact banks.  Lower
counterexample count wins.  S103 fingerprint:

```text
AP_diff_majorization      0
Ly_certificate            0
local_profile_compression 1177
scalar_additive_energy    3137
pairwise_diff_majorization 556752
```

The tournament is effectively transitive:

```text
AP_diff_majorization
  > Ly_certificate
  > local_profile_compression
  > scalar_additive_energy
  > pairwise_diff_majorization.
```

Assumption challenged: additive energy is not a scalar monotone sufficient
statistic for LRC coverage.  It is a Fejer/interval extremal carrier whose
labels must be retained until after the THM-534 sector functional is applied;
the higher-moment tail is signed, and exact threshold tiling is anchored rather
than affine-translation invariant.

## S104 addendum: exact tail target for the labelled remainder

HYP-2890 gives a concrete analytic version of the "labelled signed sector
remainder" above.  The full same-frequency additive-energy packet is positive
and absolutely convergent:

```text
Gamma_k(m)=C_{k,r}/m^4,  all C_{k,r}>0 for k=8..13.
```

But it overpredicts AP, so the signed sector/Fourier remainder is not small
noise.  With

```text
R_sf(E)=p0(E)-p0_decorr(k)-Gamma_k^sf A*(E),
```

Anchored AP-facing extremality becomes the residual-leak bound

```text
R_sf(E)-R_sf(AP) <= Gamma_k^sf(A*(AP)-A*(E)).
```

S104 exact anchored scans found no violations for k=8 (`3432` rows) or k=9
(`3003` rows); the worst k=9 ratio `0.933` occurs at the scaling-like row
`(0,2,4,6,7,8,10,12,14)`.  This is the first finite target for the labelled
remainder lemma proposed here.
