---
id: HYP-2614
title: LRC(14) support-6 residue-cusp tail
status: OPEN
source: codex-2026-06-19-S12
depends_on:
  - THM-538
  - HYP-2608
  - HYP-2612
  - HYP-2613
related:
  - HYP-2606
  - HYP-2611
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2614 - LRC(14) Support-6 Residue-Cusp Tail

## Claim

The HYP-2608(a) wide-spread tail should be proved as a signed reciprocal
theta-sum with boundary/cusp control, not as a bare absolute Minkowski volume
count.

For every exact six-support vector with nonzero coefficients not divisible by
`7`,

```text
K(n_1,...,n_6,0,...,0) = C_d(n_1,...,n_6 mod 7)/(n_1...n_6),
```

where `d` is the ambient nonzero-offset dimension and `C_d` depends only on the
six residue addresses modulo `7`.  Thus each support-6 relation hyperplane
`sum e_i n_i=0` is a finite combination of residue-addressed reciprocal sums.

The likely finishing lemma is:

```text
finite low-height wall ledger
+ residue-addressed signed reciprocal hyperplane tail
+ no-scale cluster quotient
< cap margin.
```

The new part is the middle term: use summation by parts / cotangent-Dedekind
sum technology on the relation hyperplane, with proper boundary faces controlled
by the support-6 floor and the already separated low-height wall ledger.

## Evidence From S12

Script:

- `04-computation/lrc14_support6_residue_cusp_codex_s12.py`
- output: `05-knowledge/results/lrc14_support6_residue_cusp_codex_s12.out`

The exact reciprocal residue identity was verified on `320` random
six-support vectors across ambient dimensions `6,7,8,9`; worst numerical error
was `2.948e-23`.

The shell ledgers show that the absolute mass is mostly a boundary/cusp shadow,
while signed sums are much smaller:

- AP support `(1,2,3,4,5,6)` in ambient `d=7`: through `H=28`, `absK=0.920275`
  but signed sum is `0.0316608`; on last nonempty shell `h=27`, the one-face
  boundary has `5,291,542` relations, `absK=0.0163888`, signed only
  `8.10e-5` (`abs/signed≈202`).
- Resonant support `(1,2,3,4,5,21)`: through `H=28`, `absK=0.508434` but
  signed sum is `-0.00234285`; on shell `h=27`, the one-face boundary has
  `abs/signed≈3485`.
- Wide sampled support `(2,3,4,5,6,68)` in ambient `d=8`: through `H=30`,
  `absK=0.100030` but signed sum is `8.94e-5` (`abs/signed≈1119`).
- `k=10` height-one wall support `(1,2,4,7,8,22)` in ambient `d=9`: through
  `H=24`, `absK=0.117563` but signed sum is `0.000262635`
  (`abs/signed≈448`).

These are exact support-6 relation sums using the THM-538 kernel, not the blunt
`64*c1^6/prod |n_i|` envelope.

## Guardrail

The tempting stronger shortcut is false: the raw residue coefficient does **not**
have zero one-coordinate marginals.  S12 measured worst one-coordinate
marginals:

```text
d=6: 0.384780076212
d=7: 0.193862704055
d=8: 0.107446514548
d=9: 0.0637696095737
```

So the needed cancellation is not "sum over one residue coordinate and vanish."
It must use the integer relation hyperplane, plus finite deletion of low-height
anti-coset walls.  This prevents the next proof attempt from overclaiming a
nonexistent pointwise marginal identity.

## Reframing

S11 found the signed permanent layer; S12 attaches the missing address
coordinate.  The correction term has the shape

```text
scalar/product collision:    sum e_i n_i = 0
address coordinate:          n_i mod 7 and boundary-face type
derivative sum:              reciprocal hyperplane sum 1/(n_1...n_6)
```

This is the same projection-repair grammar recurring elsewhere in the repo:
the scalar collision alone is lossy, the residue/fraction address prevents the
false absolute envelope, and the derivative/reciprocal sum is where cancellation
lives.

## Tournament Analysis Note

The useful tournament vertices are proof-obligation quotients, not runners or
arcs:

- free product envelope;
- coupled absolute hyperplane;
- residue permanent constant;
- signed reciprocal shell;
- boundary-face cancellation;
- low-height wall ledger;
- cluster translation quotient.

The S12 route tournament is transitive with Hamiltonian proof path

```text
boundary_face_cancellation
> signed_reciprocal_shell
> low_height_wall_ledger
> cluster_translation_quotient
> residue_permanent_constant
> coupled_absolute_hyperplane
> free_product_envelope.
```

This preserves the LRC predicate "wide support-6 correction below cap margin"
but destroys witness-time geometry.  That is acceptable here because the
obligation is analytic tail control, not construction of a witness time.

## Status

Open.  LRC(14) is not proved here.  HYP-2614 sharpens the single residual:
prove a signed reciprocal hyperplane tail theorem after finite low-height wall
deletion, instead of trying to force a uniform absolute Minkowski count.
