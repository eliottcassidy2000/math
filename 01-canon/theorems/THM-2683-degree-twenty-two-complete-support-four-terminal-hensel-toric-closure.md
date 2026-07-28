---
id: THM-2683
title: "Degree-twenty-two complete support-four terminal Hensel and toric closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the inherited
  genuine nonsplit degree-twenty-two branch, all five support-four strata are
  empty.  All scaled eliminants have the same irreducible fixed quintic.  At
  the first order beyond t-degree ten, BDEW and CDEW have full-rank
  three-monomial systems.  BCDW's unique kernel line misses a necessary toric
  quadric; BCEW's projective kernel line misses two toric quadrics by a
  nonzero Sylvester resultant; and BCDE's projective kernel plane has three
  toric quadrics whose degree-four Macaulay map has full row rank fifteen.
  Every certificate is uniform in both the degree-five root and degree-ten
  unordered-root-pair fields.  THM-2671's fixed-place Kummer and y=0
  mechanisms then empty the five charts.  Only full support remains in this
  branch; split/even descent, integral raising, JC(2), and DC(2) remain open.
source: root-long-frontiers-2026-07-28-support-four-toric-hensel
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2636-degree-twenty-two-BCD-triple-spectral-square-Kummer-closure
  - THM-2671-degree-twenty-two-complete-support-three-weighted-hensel-closure
script: 04-computation/jc2_degree22_support_four_toric_hensel_thm2683.py
output: 05-knowledge/results/jc2_degree22_support_four_toric_hensel_thm2683.out
script_sha256: 6f4636725542a7a3804be978bf9c4f3c79cd961d2c10f15c6af918ab0e353d1e
output_sha256: 3efbadaffde08f978ef8b80df5b476f224c4e21a9b1c14884b46d816f82a2e59
hash_basis: LF-normalized bytes
---

# THM-2683 -- every support-four stratum is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2671 closes every degree-twenty-two coefficient stratum of support at
most three.  The five support-four charts are

```text
BCDE, BCDW, BCEW, BDEW, CDEW.                            (1)
```

The linear terminal-rank method used for support three is insufficient on
the first, second, and third charts in `(1)`: their terminal monomial vectors
have more columns than the matrix rank.  The correct object is the monomial
image itself, not its ambient vector space.  Its toric binomials close BCDE,
BCDW, and BCEW; ordinary terminal rank closes BDEW and CDEW.

## 1. Fixed section and factor atlas

Use THM-2671's normalized variables

```text
y=11s,                v=u/y^2,                zeta=Z/y^3,
(b_2,...,b_6)=(B/y^2,C/y^3,D/y^4,E/y^5,W/y^6).           (2)
```

Eliminating `zeta` from the two normalized fluxes gives the same primitive
sixty-term pre-scale eliminant `R_gen` of raw content `28,344,976` used in
THM-2671.  On the five charts choose a coefficient root `rho` by

```text
BCDE, BCDW, BCEW, BDEW: rho^2=B,
CDEW:                   rho^3=C,                         (3)
```

and put `t=rho/y`.  The active normalized coefficients are their nonzero
parameter times `t` to the corresponding weights.  In each case the
specialized eliminant

```text
P(t,v) in C[t,v]                                         (4)
```

is monic of `v`-degree five and has `t`-degree ten.  Its fixed section is
the same irreducible squarefree quintic

```text
L_5(v)=
 155624547606v^5+3215383215v^4-1700698560v^3
 +58124770v^2-855470v+2583.                              (5)
```

As in THM-2671, if `(4)` factors, monic Gauss integrality lets the factors be
chosen in `C[t,v]`.  The smaller fixed factor has degree one or two.  Its
choices are exhausted uniformly by the degree-five root field `A_1` and the
degree-ten unordered-root-pair field `A_2` of `(5)`.

## 2. The common terminal Hensel system

Write a proposed factorization as

```text
P=QS,
Q=sum q_n(v)t^n,             S=sum s_n(v)t^n.             (6)
```

The fixed factors are coprime in `A_1` and `A_2`, so Hensel uniqueness gives

```text
f_n=r_n-sum_(i=1)^(n-1)q_i s_(n-i),
q_n=rem_(q_0)(f_n (s_0 mod q_0)^(-1)),
s_n=(f_n-q_n s_0)/q_0.                                  (7)
```

Since `deg_t P=10`, degree additivity in an actual polynomial factorization
forces

```text
q_11=s_11=0.                                             (8)
```

The companion implements `(7)` over the genuine sparse three-parameter rings
`A_i[p,q,r]`, checks every division and reconstructed coefficient, and
extracts the five scalar equations in `(8)`.  No parameter is encoded into
another exponent and no finite-field specialization is used.

## 3. Two full-rank closures

On `BDEW`, in parameter order `(D,E,W)`, the complete terminal monomial
support is

```text
(E, EW, DE).                                              (9)
```

The five-by-three coefficient matrix has rank three in both `A_1` and `A_2`.
It has seven nonzero maximal minors in either field.  Their first primitive
numerators have respectively degree/term pairs `(4,5)` and `(9,10)`, and
every nonzero minor is coprime to the corresponding field modulus.  Thus
`(8)` forces `(E,EW,DE)=0`, impossible on the active torus.

On `CDEW`, after the scale choice `C=1`, the complete support is

```text
(E, EW, D^2).                                             (10)
```

The same rank, seven-minor, degree/term, and modulus-coprimality statements
hold.  Equation `(8)` would force `(E,EW,D^2)=0`, again impossible.

## 4. BCDE: three kernel-plane quadrics generate every quartic

In parameter order `(C,D,E)`, use the terminal support

```text
m=(E,DE,C,CD,CD^2,C^2E,C^3).                            (11)
```

The five-by-seven terminal matrix has rank four in both fields, with pivot
columns `(0,1,2,3)`.  Its projectivized kernel is therefore a plane with the
last three columns as homogeneous coordinates `[lambda:mu:nu]`.

Every monomial point `(11)` satisfies the three toric quadrics

```text
Q_1=m_0m_3-m_1m_2,
Q_2=m_2m_4-m_3^2,
Q_3=m_2m_5-m_0m_6.                                      (12)
```

Substitute the exact three-vector kernel basis into `(12)`.  Multiply each
resulting ternary quadratic by all six degree-two monomials.  The resulting
degree-four Macaulay coefficient map has shape

```text
15 quartic monomials by 18 quadratic multiples.          (13)
```

It has full row rank fifteen in both fields, with pivot columns

```text
(0,1,2,3,4,5,6,7,8,9,10,12,13,14,15).                  (14)
```

Thus every ternary quartic belongs to `(Q_1,Q_2,Q_3)`.  A projective common
zero would in particular force `lambda^4=mu^4=nu^4=0`, which is impossible.
The selected fifteen-by-fifteen minor is nonzero and modulus-coprime.  Its
primitive numerator has degree four and five terms over `A_1`, with digest

```text
3df995446b6ce494aa6fcd37da780361aac48174e812eaf4dbe3b080e620ff72,
```

and degree nine and ten terms over `A_2`, with digest

```text
2b04c535fd5d8428843d85681ef7df95756adcecb389130e868b9cb9bd949b38.
```

Hence the kernel plane contains no monomial point and BCDE is closed.

## 5. BCDW: the kernel line misses one toric quadric

In parameter order `(C,D,W)`, the BCDW terminal vector is

```text
m=(C,CD,CD^2,C^3,CW)=C(1,D,D^2,C^2,W).                  (15)
```

Its five-by-five coefficient matrix has rank exactly four in both factor
fields.  Rows `(0,1,2,3)` give a cofactor kernel vector

```text
k=(k_0,...,k_4),                                         (16)
```

whose five coordinates are nonzero and modulus-coprime.  The full fifth row
annihilates `(16)`, so `(16)` spans the complete kernel at every embedding.

Every monomial vector `(15)` obeys the toric identity

```text
m_0 m_2-m_1^2=0.                                        (17)
```

If `(8)` held, then `m=lambda k` for some nonzero `lambda`, and `(17)` would
force

```text
k_0 k_2-k_1^2=0.                                        (18)
```

Exact reduction gives a nonzero modulus-coprime numerator in both fields.
For `A_1` it has degree four, five terms, and digest

```text
32c3befea4674e4ff2900c5210753260fdc8f0a84f2d8c06e29894d7c35718bd;
```

for `A_2` it has degree nine, ten terms, and digest

```text
8de0da3bd0aed7bdf62b3b706240dbeadaf5c7d250341a58398ea90c33a9e8d7.
```

Thus `(18)` is impossible uniformly.  The rank defect was an artifact of
linearizing the monomial image; its first toric equation closes the chart.

## 6. BCEW: the kernel line misses two toric quadrics

In parameter order `(C,E,W)`, use the terminal support order

```text
m=(E,EW,C,CW,C^2E,C^3).                                 (19)
```

The five-by-six terminal matrix has rank four in both fields.  Exact RREF
has pivot columns `(0,1,2,3)`, so its projectivized kernel is a line with
homogeneous coordinates `[lambda:mu]`, using the last two columns as the
free coordinates.

Every monomial vector `(19)` satisfies

```text
F=m_0m_3-m_1m_2=0,
G=m_2m_4-m_0m_5=0.                                      (20)
```

After substituting the two exact kernel basis vectors, `F` and `G` become
binary homogeneous quadratics in `[lambda:mu]`.  Every one of their six
coefficients is nonzero and modulus-coprime.  Their four-by-four Sylvester
resultant is also nonzero and modulus-coprime in both fields.  Its primitive
numerator has degree four and five terms over `A_1`, with digest

```text
a7846c6cd848e619fe62d19dd482171ea4a9f1f57ee0b2e0ce2b344e39876833,
```

and degree nine and ten terms over `A_2`, with digest

```text
d8d17d914d6de0d97a15acec27fa579e511e454787492f62e3a3e4b6efd653b1.
```

Hence the two quadratics have no common point in the projective kernel line.
No nonzero monomial vector `(19)` can satisfy `(8)`, closing BCEW.

## 7. Irreducibility, Kummer closure, and consequence

Sections 3--6 exclude every line and quadratic fixed factor.  Section 1
exhausts the smaller side of every possible factorization of the monic
quintic `(4)`, so each of the five signed-`t` eliminants is absolutely
irreducible.

THM-2671's remaining argument is uniform in the scale choice `(3)`.  At all
five roots of `(5)`, `t` is a uniformizer and the first flux makes `zeta` a
unit.  Since

```text
H=rho^3 zeta/t^3=Z=T^2,                                 (21)
```

the retained square has odd order `-3` at five fixed places.  Its connected
double cover has genus at least two.  A nonconstant `t` would embed this
cover into `C(x)`, impossible; a constant `t` makes successively `y,v,zeta`,
and `T,q` constant, contradicting the genuine nonsplit deck.  THM-2671's
two uniform quintics separately close `y=0` before any scale is chosen.

Therefore

```text
BCDE=BCDW=BCEW=BDEW=CDEW=empty.                          (22)
```

Together with the support-at-most-three closure, any survivor in the
inherited genuine nonsplit degree-twenty-two branch now has full coefficient
support `BCDEW`.  This does not close that chart, the split/even short-edge
branch, integral `2`-adic raising, another degree, `JC(2)`, or `DC(2)`.

## 8. Reproduction

Run

```bash
python3 04-computation/jc2_degree22_support_four_toric_hensel_thm2683.py
python3 -O 04-computation/jc2_degree22_support_four_toric_hensel_thm2683.py
```

Both executions byte-match the stored output and declared hashes.  An
independent hostile audit reconstructed the five terminal supports and both
factor fields, verified the native three-parameter ring and Hensel recursion,
checked the BCDW cofactor ordering and toric defect, and independently derived
the BCEW two-quadric Sylvester and BCDE degree-four Macaulay obstructions.

QED.
