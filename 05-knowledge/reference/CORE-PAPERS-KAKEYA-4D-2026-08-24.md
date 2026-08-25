# Four-dimensional Kakeya: primary pins and transfer boundaries

> **Freshness/status audit: 2026-08-24.**  `OPEN` below refers to the
> Euclidean linear Kakeya set/maximal-function conjecture in dimension four.
> Multilinear, sticky, plany, finite-field geometric, and Arithmetic-Kakeya
> statements have different hypotheses and conclusions and must not be merged.

## Current Euclidean status

**OPEN.**  A Kakeya (Besicovitch) set in `R^4` contains a unit line segment in
every direction.  The conjecture says every such set has Hausdorff and
Minkowski dimension four; the corresponding sharp linear maximal-function
estimate is also open.  Joshua Zahl's
[2025--2026 survey](https://arxiv.org/abs/2512.09397) is the current overview
used here.

### General four-dimensional bounds: published and preprint

Katz--Zahl,
[*A Kakeya maximal function estimate in four dimensions using planebrushes*](https://arxiv.org/abs/1902.00989v3),
is **PUBLISHED, with the corrected 2025 v3 used here**.  It proves the general
Hausdorff-dimension lower bound

```text
3+(sqrt(17665)-97)/600 = 3.059849573...
```

and a distinct maximal-function estimate at dimension

```text
3+(2195737-13925 sqrt(17665))/6959096 = 3.049570924....
```

The v3 correction repairs Lemma 7.1.  Do not attach the Hausdorff number
`3.059...` to this maximal-function theorem.  These estimates do not prove
dimension four.

Borges--Chan--Chen--Liu--Xi--Zhan,
[*Restriction and Kakeya maximal estimates in R^4*](https://arxiv.org/abs/2511.22824),
is a **PREPRINT v1 CLAIM, submitted 2025-11-28**.  Its Theorem 1.2 claims the
improved general maximal estimate at

```text
d_0=(159+sqrt(145))/56 = 3.0543....
```

This supersedes `3.049570924...` as the current claimed maximal-function
record, while Katz--Zahl remains the corrected published baseline.  The
preprint estimate implies only the smaller Hausdorff lower bound `d_0`, so it
does not displace Katz--Zahl's general Hausdorff record `3.059849573...` and
does not prove dimension four.

Guth--Zahl,
[*Polynomial Wolff axioms and Kakeya-type estimates in R^4*](https://arxiv.org/abs/1701.07045),
provides the polynomial-Wolff and trilinear broad input.  Its corrected
general conditional exponent is `3+1/40`, not the superseded `3+1/28`; the
paper does not prove that every Kakeya family satisfies its polynomial Wolff
axioms.

### Broad/multilinear and narrow/plany regimes

Guth,
[*The endpoint case of the Bennett--Carbery--Tao multilinear Kakeya conjecture*](https://arxiv.org/abs/0811.2251),
proves the endpoint `n`-linear theorem, hence the four-linear theorem in
`R^4`.  It controls genuinely transverse direction families.  It does not
imply the linear conjecture because narrow or plany configurations can have
small wedge products.

The planebrush is the complementary narrow mechanism in the corrected
Katz--Zahl argument.  Łaba--Rai Choudhuri--Zahl's
[*A bound for plany Kakeya sets in F_q^4 using the planebrush method*](https://arxiv.org/abs/2507.09605)
is a finite-field plany theorem, not a Euclidean dimension-four solution.

Rai Choudhuri,
[*An improved bound on the Hausdorff dimension of sticky Kakeya sets in R^4*](https://doi.org/10.4171/RMI/1621),
proves `dim_H K >= 13/4` under the extra sticky hypothesis; see also the
[arXiv version](https://arxiv.org/abs/2410.23579).  This is **PUBLISHED** in
*Revista Matemática Iberoamericana* 42 (2026).  It neither treats arbitrary
Kakeya sets nor proves the predicted sticky dimension four.

## Finite-field geometry is a separate solved cardinality problem

Dvir,
[*On the size of Kakeya sets in finite fields*](https://arxiv.org/abs/0803.2336),
proves that a subset of `F_q^n` containing a line in every direction has size
at least `C_n q^n`.  Bukh--Chao,
[*Sharp density bounds on the finite field Kakeya problem*](https://arxiv.org/abs/2108.00074),
gives the asymptotically sharp density bound

```text
|K| >= (2-1/q)^(-(n-1)) q^n.
```

These are geometric finite-field Kakeya theorems.  They do not transfer
finite cardinality directly to Euclidean Hausdorff dimension, and a list of
directions without chosen affine lines is not a finite-field Kakeya set.

## Arithmetic Kakeya is a third problem

The Katz--Tao additive-projection or Arithmetic-Kakeya conjecture is distinct
from Dvir's geometric finite-field theorem.  Green--Ruzsa,
[*On the arithmetic Kakeya conjecture of Katz and Tao*](https://arxiv.org/abs/1712.02108),
gives equivalent formulations and explains its Euclidean implications.
Katz--Tao,
[*New bounds on Kakeya problems*](https://arxiv.org/abs/math/0102135),
proves the quantitative implication

```text
SD(alpha)  ==>  upper Minkowski dimension in R^n >= 1+(n-1)/alpha.
```

Tao's 2025 status paper
[*Sum-difference exponents for boundedly many slopes, and rational complexity*](https://arxiv.org/abs/2511.15135)
pins the best upper bound on the relevant exponent at
`alpha=1.6751308705...`, the root in `(1,2)` of `alpha^3-4alpha+2=0`.  In dimension
four, the Katz--Tao implication then gives only

```text
1+3/alpha = 2.790904850...,
```

below the corrected general Hausdorff bound `3.059849573...`.  To improve the
latter through this implication alone would require

```text
alpha < 3/(3.059849573...-1) = 1.456417031....
```

The repository's certificate cells prove none of these external Arithmetic-
Kakeya exponents; they remain a search language and exact finite audit.

## Repo import: THM-4035

[THM-4035](../../01-canon/theorems/THM-4035-sixty-clock-separation-and-finite-kakeya-spine.md)
imports only the algebraic transversality predicate.  Over `F_61`, one common
`C_60` clock admits both a full-spark twisted-cubic representation and a
two-plane representation with every four-minor zero.  Three independent
clocks chart the nonzero torus of `P^3(F_61)`.

The import preserves cyclic phase and determinant rank.  It destroys or omits
affine line basepoints, direction completion, incidence multiplicity,
Euclidean angular size, shadings, two-ends, scale, and Wolff nonconcentration.  Those
coordinates are mandatory sidecars before any Kakeya theorem can consume the
finite model.

## Explicit non-implications

None of the cited results or THM-4035 proves:

- the Euclidean Kakeya conjecture in `R^4`;
- a new Euclidean dimension or maximal-function bound;
- a finite-field-to-Euclidean transfer theorem;
- that one periodic orbit supplies all projective directions;
- Arithmetic Kakeya at any new exponent; or
- LRC(14).
