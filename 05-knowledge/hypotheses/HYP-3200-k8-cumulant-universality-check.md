# HYP-3200: k=8 Cumulant Universality Check

**Status:** PARTIALLY-TRUE / exact bounded-bank census; not an LRC14 proof.

## Claim Tested

The HYP-3160/S31ai thread left two live claims around the k=8 hard node:

1. Consecutive speeds should maximize the degree-two total covariance
   `Sigma kappa_2` of the six inner-sector emptiness indicators.
2. The associativity-compression defect
   `Sigma kappa_3 / S3` for consecutive speeds looked numerically close to
   `1/7`, suggesting a possible universal apex-prime law.

HYP-3200 tests both claims exactly on the anchored bounded bank

```text
E = {0} union A,  A subset {1,...,14}, |A|=7.
```

The bank has `3432` rows and `3431` primitive rows.  The nonprimitive row
`(0,2,4,6,8,10,12,14)` is the dilation twin of consecutive speeds.

Incoming HYP-3161/S31aj is the broader companion: it checks
`Sigma kappa_3/S3` for consecutive rows across `k=4..10`, verifies that the
near `1/7` value is just a k=8 crossing, and records the stronger all-bank
readout that consecutive speeds have no `Sigma kappa_2` beaters among all
`3432` bounded k=8 clusters.  In that same pass, `2002` is identified as
`C(14,5)=2*7*11*13`, the k=10 binding-row configuration count.  Thus `2002`
belongs to Pascal/binomial configuration currency, not to an entropy or
minimum-description proof route.

## Exact Readout

For `E=(0,1,2,3,4,5,6,7)`:

```text
Sigma kappa_2      = 6237419/8643600
Sigma kappa_3      = 407891843/2117682000
S3                 = 991/735
Sigma kappa_3 / S3 = 407891843/2855269200
ratio - 1/7        = -3757/2855269200
```

Therefore the `1/7` claim is false exactly:

```text
consec_ratio_equals_1/7 = False
rows_with_ratio_exactly_1/7_all = 0
rows_with_ratio_exactly_1/7_primitive = 0
```

The positive result is the covariance statement:

```text
primitive: consec_rank_Sigma_kappa2_MAX = 0/3431
all:       consec_rank_Sigma_kappa2_MAX = 1/3432
```

The all-bank rank is only displaced by the nonprimitive dilation twin.  In
primitive normal form, consecutive speeds are the exact top row for the tested
bounded bank.

The entropy route is also refuted exactly in this bank: consecutive speeds
have high entropy, with `consec_rank_entropy_MIN = 3428/3431` and
`consec_rank_entropy_MAX = 2/3431` among primitive rows.  Extremality is spread
as covariance, not concentrated as low entropy.

The `kappa_4` sign remains useful as a phi4 stabilizer sidecar, but it is not
itself an exact minimum theorem in this bank: consecutive speeds have
`consec_rank_Sigma_kappa4_MIN = 9/3431`.

## Tournament Analysis

Vertices are proof signals rather than runners or sectors.  The pairwise
observable is whether a signal survives exact bounded-bank verification.  The
gauge directs `A -> B` when verification score of `A` is higher.  The resulting
priority path is:

```text
total_covariance_sigma_k2
-> bounded_bank_exactness
-> kappa4_phi4_stabilizer
-> associator_sigma_k3_sidecar
-> random_3002_probe
-> raw_scalar_defect_ratio
-> entropy_min_description_route
-> exact_one_seventh_ratio_claim
```

This explicitly challenges the default vertex choice: the proof-facing
tournament is over retained information channels, not over runners, gaps, or
raw sectors.

## Proof-Frontier Use

HYP-3200 sharpens HYP-3160 into a cleaner split:

- Even/commutative target: prove a primitive-normal-form inequality for
  `Sigma kappa_2`, ideally by reflection-fold, Perron/covariance, or
  Ferromagnetic comparison arguments.
- Odd/non-associative target: keep `Sigma kappa_3`, `S3`, Worpitzky
  `-9S3`, and minority-edge data as sidecars.  Do not compress them into a
  scalar `1/7` law.
- Stabilizer sidecar: use `kappa_4` as a phi4 sign/stability diagnostic, not
  as a standalone minimization theorem.

The next theorem-facing task is not "prove `1/7`."  It is:

```text
For primitive k=8 bounded-core rows, prove consecutive speeds maximize
total empty-sector covariance Sigma kappa_2, then attach the odd Worpitzky /
associator residual without forgetting its sidecar coordinates.
```

## Artifacts

- Script: `04-computation/lrc_k8_cumulant_universality_codex_20260628.py`
- Result: `05-knowledge/results/lrc_k8_cumulant_universality_codex_20260628.out`

## Links

HYP-3200 refines HYP-3160 and should be read with HYP-3161, HYP-3154,
HYP-3153, HYP-3152, HYP-3151, HYP-3150, HYP-3147, HYP-3144, HYP-3142,
HYP-3139, HYP-3138, HYP-3132, HYP-3122, THM-577, and the Pascal-pair-mass
`2002=C(14,5)` thread.
