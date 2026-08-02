---
id: THM-3112
title: "Cycle-weighted Young-subgroup gap and uniform-octopus boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every finite
  nonnegative alphabet has a cycle-weighted symmetric-group gap operator
  which is a nonnegative mixture of Young-subgroup projection gaps.  The
  normalized one-row THM-3110 target is exactly a signed combination of these
  operators.  At degree six it is outside the conjugacy-averaged span of all
  uniform alpha_A octopus templates, so a cycle-weighted or rooted refinement
  is necessary.  This is a reduction and obstruction, not product-Gamma or
  Gaussian-moment closure.
source: root/multiscale-newton-flag-2026-08-02
audit: >
  An independent hostile audit rederived the coefficientwise cycle-colouring
  identity, colour-weight normalization, subgroup projection order,
  Frobenius/Kostka spectra, singleton alpha normalization, signed beta
  identity, direct-versus-Cesi sign convention, conjugacy-sum formula,
  generalized-octopus span obstruction, both degree-six rank/spectrum tables,
  and the degree-five monomial hostile.  Fresh normal and optimized runs
  byte-match the stored output; both LF hashes and documentation checks pass.
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
script: 04-computation/gmc_cycle_weighted_young_subgroup_gap_thm3112.py
output: 05-knowledge/results/gmc_cycle_weighted_young_subgroup_gap_thm3112.out
script_sha256: 5f285c7f4a8dcecd0f40e29a12abc76237538d42dfd17fa72aa8886f3e19f3bc
output_sha256: 45f91150d23d465eca9571eedd4bee6fdd7f541d1faf6fb971837b0b94def616
hash_basis: LF-normalized bytes
---

# THM-3112 -- cycle-weighted Young-subgroup gap and uniform-octopus boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem isolates the exact operator hidden behind the one-row extremum
suggested by the post-THM-3110 Young-cover census.  Its positive result is a
universal Young-subgroup projection decomposition.  Its negative result is
equally important: the published uniform-subgroup octopus cone is already too
small in degree six.  The missing inequality must retain cycle weights,
colour-fibre refinement, or a rooted insertion path.

## 1. Cycle-weighted central operators

Let `S=(s_1,...,s_d)` be a finite labelled alphabet of nonnegative real
numbers, with zeros and repeated entries allowed.  For `sigma in Sym_N` put

```text
p_type(sigma)(S)=product_(C cycle of sigma) sum_i s_i^|C|.       (1)
```

In the real group algebra define

```text
Omega_N(S)=1/N! sum_(sigma in Sym_N) p_type(sigma)(S) sigma,
beta_N(S)=h_N(S) identity-Omega_N(S).                            (2)
```

For a colouring `f:[N]->[d]`, let

```text
H_f=product_i Sym(f^(-1)(i)),
P_f=|H_f|^(-1) sum_(tau in H_f) tau,
a_f=N!^(-1)|H_f| product_(j=1)^N s_(f(j)).                      (3)
```

Then

```text
Omega_N(S)=sum_f a_f P_f,             sum_f a_f=h_N(S),
beta_N(S)=sum_f a_f(identity-P_f) >= 0.                         (4)
```

Here `>=0` is the direct positive-semidefinite order in every finite
dimensional representation of `Sym_N`.

### Proof

The coefficient of `sigma` on the right of the first identity in `(4)` is

```text
1/N! sum_(f:sigma f=f) product_j s_(f(j)).                       (5)
```

A fixed colouring is constant on every cycle.  Choosing its colour
independently on each cycle turns `(5)` into `(1)/N!`, proving the first
identity coefficient by coefficient.

For a fixed occupancy vector `(n_1,...,n_d)`, there are
`N!/product_i n_i!` colourings, while each has

```text
a_f=(product_i n_i!)/N! product_i s_i^n_i.                       (6)
```

Their total is the corresponding monomial in `h_N(S)`.  Summing all weak
compositions proves the second identity.  Finally `P_f` is a self-adjoint
idempotent in every unitary representation: it is the orthogonal projection
onto the `H_f`-invariant subspace.  Thus every `identity-P_f` is positive
semidefinite and the coefficients `a_f` are nonnegative.  This proves `(4)`,
including zero and repeated alphabet entries.

If `H subset K`, then `P_H P_K=P_K P_H=P_K`, so

```text
identity-P_H <= identity-P_K.                                  (7)
```

Consequently a stochastic coupling which coarsens colour fibres gives an
operator comparison.  Equation `(7)` is the precise star--mesh/refinement
order left available by the theorem.

## 2. Frobenius and tableau shadows

For a partition `lambda` of `N`, Frobenius character inversion gives the
scalar of `Omega_N(S)` on the Specht module `V_lambda`:

```text
Omega_N(S)|V_lambda=s_lambda(S)/f^lambda identity.               (8)
```

Hence `(4)` also says

```text
h_N(S)-s_lambda(S)/f^lambda >=0.                                 (9)
```

There is a second direct proof of `(9)`.  In
`s_lambda=sum_mu K_(lambda,mu)m_mu`, standardizing equal entries from left
to right injects the semistandard tableaux of any fixed content `mu` into
the `f^lambda` standard tableaux.  Thus `K_(lambda,mu)<=f^lambda`, whereas
`h_N=sum_mu m_mu`.  The inequality is coefficientwise before specializing
to the nonnegative alphabet.  Zero-padding the alphabet simply kills
monomials and changes nothing.

For the one-letter alphabet `(1)`, the only colouring has `H_f=Sym_N`, so

```text
N beta_N((1))
 =N identity-1/(N-1)! sum_(sigma in Sym_N)sigma
 =alpha_[N].                                                      (10)
```

Thus `(4)` is literally a cycle-weighted extension of the `alpha_A`
operator used in the octopus inequality, not merely a spectral analogy.

## 3. Exact row-normalized product-Gamma operator

Let a finite signed alphabet bank be

```text
Phi(f)=sum_R c_R f(S_R),                                         (11)
```

and distinguish an alphabet `Q` with positive coefficient `c_Q`.  Whenever
`Phi(h_N)>0`, put

```text
r_N=Phi(h_N)/(c_Q h_N(Q)),
D_N=sum_R c_R Omega_N(S_R)-r_N c_Q Omega_N(Q).                   (12)
```

The trivial scalar of `D_N` is zero by definition.  Combining `(4)` and
`(12)` cancels the identity terms exactly and gives

```text
D_N=-sum_R c_R beta_N(S_R)+r_N c_Q beta_N(Q).                    (13)
```

By `(8)`, its scalar on `V_lambda` is

```text
[Phi(s_lambda)-r_N c_Q s_lambda(Q)]/f^lambda.                    (14)
```

Therefore `D_N>=0` is exactly the assertion that the normalized response is
minimized by the row partition `(N)`.  In the THM-3110 banks every residual
alphabet is coordinatewise dominated by `Q`.  Because Schur functions have
nonnegative monomial coefficients, `s_lambda(Q)=0` then forces
`s_lambda(S_R)=0` for every `R`; hence all terms in `(14)` vanish.  Otherwise
`(14)>=0` is equivalent to

```text
Phi(s_lambda)/(c_Q s_lambda(Q)) >= Phi(h_N)/(c_Q h_N(Q)).        (15)
```

In Alon--Kozma--Puder's direct operator order, `(15)` is `D_N>=0`.  In
Cesi's Laplacian convention
`Gamma(G)={w: I(w)-rho(w)>=0 for every rho}`, augmentation `I(D_N)=0`
makes the equivalent statement `-D_N in Gamma(Sym_N)`.  The sign conventions
must not be identified.

Equations `(13)--(15)` are an exact all-degree reformulation.  They do not
prove that the signed combination in `(13)` is positive.

## 4. Why the published uniform octopus cone cannot suffice

For `2<=k<=N`, conjugacy-sum `alpha_A` over all `k`-subsets `A` of
`[N]`.  (Dividing by `binom(N,k)` would only rescale the generator.)  At a
nonidentity permutation moving exactly `s` points, its
coefficient is

```text
-binom(N-s,k-s)/(k-1)!                                           (16)
```

when `k>=s`, and zero otherwise.  It depends only on support size, not on
cycle type.  Hence every conjugacy-symmetrized linear combination of all
`alpha_A` lies in the `(N-1)`-dimensional support-size span.

Every operator difference in the ordinary octopus inequality and in the
disjoint, common-`(k-1)`, and co-size-one generalized octopus inequalities
of Alon--Kozma--Puder is a linear combination of `alpha_A`.  If a central
`D_N` were a nonnegative combination of any such possibly uncentred
templates, averaging the claimed identity over conjugation would put `D_N`
in the span described by `(16)`.

This is false already for both THM-3110 banks at `N=6,(a,b)=(1,2)`.  The ten
nonidentity conjugacy classes and the five averaged generators have exact
rank five.  Appending `D_6` raises the rank to six.  More transparently, the
two cycle types `(4,1,1)` and `(2,2,1,1)` both move four points but have

```text
                         q_D(4,1,1)          q_D(2,2,1,1)
I1                  -9741852/74353       -32006680/74353
I2                -22944341/185384      -73957487/185384.        (17)
```

The coefficients differ in each bank.  Thus `D_6` is outside even the
linear span, a fortiori outside the conic hull, of every published
`alpha_A` template.  The unique positive nonidentity class coefficient is
the transposition class `(2,1^4)`:

```text
I1: 146677056/74353,             I2: 368938435/185384.           (18)
```

Nevertheless all nontrivial Specht scalars are positive; the exact minima
are `737590/74353` and `112775/23173`.  The obstruction is to the uniform
template decomposition, not to the desired operator inequality at degree
six.

## 5. A sharp false shortcut and the remaining coupling problem

The positive decomposition of each individual `beta_N(S)` does **not** make
the signed combination `(13)` coefficientwise positive.  Already at
`N=5,(a,b)=(1,2)` and monomial/occupancy type `(2,1,1,1)`, the response ratios
are

```text
             row rate       monomial rate          difference
I1             6/8447          10/3007       -66428/25400129
I2            16/62019           5/4484      -238351/278093196.  (19)
```

For any fixed set partition of that block type, the coefficient in `(13)`
is a positive shape-dependent factor times the displayed difference.
Therefore it is negative.  A termwise Young-subgroup-gap proof is impossible.

The surviving target is narrower and more structured:

1. couple several signed bank atoms simultaneously in the colour-fibre
   partition space;
2. use subgroup inclusion `(7)` or a nonlinear star--mesh elimination before
   centralizing; and
3. retain a rooted Jucys--Murphy/tableau insertion path, because the degree-six
   transposition correction and the exact content-order hostiles show that
   neither support size nor scalar content is sufficient.

This is an Ewens-deformed generalized-octopus problem.  It is not solved by
the present theorem.

## 6. Literature boundary

Alon--Kozma--Puder, *On the Aldous--Caputo Spectral Gap Conjecture for
Hypergraphs* ([arXiv:2311.02505](https://arxiv.org/abs/2311.02505)),
Definition 1.5 and Theorems 1.6--1.9, supply the `alpha_A` operators and the
three generalized octopus templates used in Section 4.  Cesi,
*A few remarks on the octopus inequality and Aldous' spectral gap conjecture*
([arXiv:1310.6156](https://arxiv.org/abs/1310.6156)), Section 3, supplies the
Laplacian cone convention distinguished after `(15)`.  Neither source states
the cycle-weighted identity `(4)` or the THM-3110 signed inequality.

## 7. Exact companion and scope

Run

```text
python 04-computation/gmc_cycle_weighted_young_subgroup_gap_thm3112.py
python -O 04-computation/gmc_cycle_weighted_young_subgroup_gap_thm3112.py
```

and compare with

```text
05-knowledge/results/gmc_cycle_weighted_young_subgroup_gap_thm3112.out.
```

The companion checks `(4)` coefficient by coefficient on three finite
group-algebra controls, derives both degree-six operators from the live
THM-3110 banks, verifies their augmentation and every Specht scalar, computes
the rank jump and `(17)--(19)`, and uses only integer or rational arithmetic.

This theorem does not prove the all-degree row extremum, arbitrary anchored
product-Gamma goodness, SFC in width three, GMC(2), NC2, LRC(14), JC(2), or
DC(2).
