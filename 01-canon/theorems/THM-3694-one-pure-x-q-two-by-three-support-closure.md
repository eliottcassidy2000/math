---
id: THM-3694
title: "One-pure-x-Q normalized two-by-three support closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In
  the normalized planar chart P=x plus exactly two distinct nonlinear
  monomials and Q=y plus exactly three distinct nonlinear monomials, suppose
  both P supports have positive x exponent and exactly one Q support is pure
  x.  Then Jac(P,Q) is not constant, with no hypothesis on the six cross
  determinants.  The exact saturated support atlas forces the pure-x exponent
  to be 1/2 or 3/2 on every nonempty algebraic stratum, contradicting its
  nonnegative-integral nonlinear support.  Other axis branches and JC(2)
  remain open; no counterexample is claimed.
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS -- root independently checked the eleven-label expansion, exact
  activity, partition completeness and affine exits, the four principal-open
  cover, complex-superset saturation, specialized ideal membership on every
  positive-dimensional stratum, the characteristic-zero lattice conclusion,
  the 180,500-point independent census, and both controls.  No typed dual is
  inferred.
depends_on: []
related:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3689-fully-transverse-two-by-two-sparse-support-gate
  - THM-3690-complete-normalized-two-by-two-sparse-support-closure
  - THM-3692-ordinary-two-by-three-shift-root-peeling
script: 04-computation/jc2_one_pure_x_q_two_by_three_support_closure_thm3694.py
output: 05-knowledge/results/jc2_one_pure_x_q_two_by_three_support_closure_thm3694.out
script_sha256: 22164ce7cab89d2617861a2a18649865089372a4e8e0d04b5cf6fca32f5854dd
output_sha256: 33f3036832b828a021ba9d53d87b1b3e06e2dc7e9cc6c96484b36ad2d8095ad5
semantic_sha256: 1148a4e75ef85db80d89cdcdf552d97f0321b9a330a30f08a3dc97015aca8adf
hash_basis: raw LF bytes for files; affine-exit ledger, saturated Groebner bases, exponent obstruction, and bounded census for semantic hash
---

# THM-3694 -- one-pure-x-Q normalized two-by-three support closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem closes the first axis cell left by THM-3692.  Its terminal
obstruction is not a bounded-degree absence: every nonempty complex support
stratum forces the nominal pure-`x` exponent off the integer lattice.

Work over a characteristic-zero field `K`.  Write `X^(r,s)=x^r y^s` and
`|(r,s)|=r+s`.

## 1. Statement

Let

```text
A,C,B,D in Z_{>=0}^2,                  k in Z_{>=2},
A != C,                                B,D,(k,0) pairwise distinct,
|A|,|C|,|B|,|D| >=2,
A_1 C_1 B_2 D_2 !=0.                                      (1)
```

For arbitrary nonzero coefficients `a,c,b,d,e in K*`, put

```text
P=x+a X^A+c X^C,
Q=y+b X^B+d X^D+e x^k.                                    (2)
```

Then

```text
Jac(P,Q) is not constant.                                  (3)
```

No cross-transversality is assumed.  Any subset of

```text
det(A,B), det(A,D), det(A,(k,0)),
det(C,B), det(C,D), det(C,(k,0))                            (4)
```

may vanish.  The theorem does not cover two pure-`x` terms in `Q`, a pure-`y`
term in `P`, or the typed dual obtained by interchanging the two components.

## 2. Eleven labels and the finite equality core

The eleven potential terms of `Jac(P,Q)-1` have buckets

```text
A-e_1, C-e_1, B-e_2, D-e_2, (k,-1),
A+B-e_1-e_2, A+D-e_1-e_2, A+(k,0)-e_1-e_2,
C+B-e_1-e_2, C+D-e_1-e_2, C+(k,0)-e_1-e_2.                (5)
```

The fifth divergence label is inactive because the pure-`x` term has zero
`y` derivative.  The other four divergence labels are active by `(1)`.  Each
cross label is active if and only if its determinant in `(4)` is nonzero,
because all five displayed source coefficients are nonzero and the field has
characteristic zero.

Every active term in `(5)` has positive total degree.  Hence a constant
Jacobian would be exactly `1`, and every active bucket would contain at least
two labels.  For each of the `64` cross-activity masks, enumerate every
no-singleton partition of the four fixed active divergence labels and the
declared-active cross labels.  This gives exactly

```text
52,922 activity/partition systems.                          (6)
```

Impose the axis equation for the fifth label and all two-coordinate bucket
equalities within each partition block.  Exact affine-linear solving removes
`52,710` systems with the ordered ledger

```text
inconsistent                                      42,447
A=C                                                 7,782
two Q supports equal                                1,813
|A| forced below 2                                    226
|C| forced below 2                                    226
|B| forced below 2                                    108
|D| forced below 2                                    108.               (7)
```

The ordering only assigns one diagnostic to a system that may have several
defects.  Every exit contradicts `(1)`.  Exactly `212` affine residual systems
remain.

## 3. Exact principal-open saturation

On each residual system, every declared-inactive determinant in `(4)` is
imposed as an equation.  Let `S_0` be the product of

```text
- the four active divergence coordinates A_1,C_1,B_2,D_2;
- all declared-active determinants in (4);
- |U|(|U|-1) for U in {A,C,B,D,(k,0)}.                     (8)
```

For nonnegative integer vectors, the last factors are nonzero exactly when
the displayed supports have total degree at least two.  This is a lossless
relaxation of `(1)` to complex exponent coordinates.

Distinctness is a union of principal opens.  Choose one of the two coordinate
differences witnessing `A!=C`, and independently one witnessing `B!=D`.  The
active `y` coordinates of `B,D` already distinguish them from `(k,0)`, so no
additional witness is needed.  The four choices cover the full distinctness
condition in `(1)`.

For each choice, multiply its two witnesses into `S_0`, call the result `S`,
and introduce an inverse variable `s` with

```text
s S - 1 = 0.                                               (9)
```

Equation `(9)` is the exact principal-open encoding; it never divides by an
unchecked factor.  The ledger is

```text
212 residual systems x 4 witness charts = 848 saturations,
0 identically empty witness products,
800 unit ideals,
48 nonunit saturated ideals.                               (10)
```

The `48` nonunit bases have exact profiles

```text
32: one exponent parameter, two affine-linear basis equations;
16: three exponent parameters, three basis equations of degrees 18,1,1.    (11)
```

Thus the latter `16` charts are positive-dimensional; they are not discarded
or replaced by selected sample points.

## 4. Uniform half-integral obstruction

Let `q` denote the nominal pure-`x` exponent `k`, now temporarily treated as a
complex support coordinate.  For **every one** of the `48` nonunit reduced
Groebner bases, exact polynomial reduction gives

```text
(2q-1)(2q-3) = 4q^2-8q+3 = 0.                              (12)
```

Equivalently, the polynomial in `(12)` belongs to every surviving saturated
ideal.  This ideal-membership statement applies to the full
positive-dimensional components in `(11)`, not only to extracted points.
Every algebraic point therefore has

```text
q in {1/2,3/2}.                                             (13)
```

But `(1)` requires `q=k` to be an integer at least two.  Hence none of the
`48` algebraic strata contains an admissible exponent packet.  Together with
`(7)`, this proves `(3)`.

## 5. Independent census and controls

An independent direct census through total degree six constructs the active
buckets from the exponent vectors, without activity partitions, affine
solving, Groebner bases, or saturation.  Its exact universe is

```text
25 nonlinear monomials,
180,500 support pairs/triples satisfying the typed axis conditions,
0 no-singleton survivors.                                  (14)
```

The closest small hostile illustrates the terminal debt.  Take

```text
P=x+a x y+c x^3,
Q=y+b y^2+d x^2y+e x^4.                                   (15)
```

Its debt is

```text
Jac(P,Q)-1
 = (a+2b)y +(3c+d)x^2 +2ab y^2
   +(-ad+6bc)x^2y +(-4ae+3cd)x^4.                          (16)
```

The choices

```text
b=-a/2,                  d=-3c,                  e=-9c^2/(4a)
```

cancel all four nonterminal buckets, but leave the debt `-a^2 y^2`.  Equivalently,
the `y^2` coefficient `2ab` is forced nonzero under the chart hypotheses.  This
is the smallest visible form of the algebraic obstruction, not a counterexample.

Support drop remains a real positive control:

```text
P=x+alpha y^2,                    Q=y+delta P^2
```

has one nonlinear term in `P`, three in `Q`, and `Jac(P,Q)=1`.

## 6. Reproduction and consequence

Run

```bash
python3 -B 04-computation/jc2_one_pure_x_q_two_by_three_support_closure_thm3694.py
python3 -O -B 04-computation/jc2_one_pure_x_q_two_by_three_support_closure_thm3694.py
```

Both streams must agree byte-for-byte with the stored output.  The companion
uses exact rational polynomial arithmetic, checks every ledger in `(6)--(14)`,
reduces `(12)` against every nonunit basis, verifies the hostile and tame
controls, and rejects inactive Python `assert` statements.

Together with THM-3692, the normalized two-by-three chart with both `P`
divergence labels active is now closed when zero or one `Q` divergence label is
inactive.  The next typed cells have two pure-`x` terms in `Q`, or at least one
pure-`y` term in `P`.  No symmetry claim identifies those cells, and no Keller
counterexample is constructed here.  **QED.**
