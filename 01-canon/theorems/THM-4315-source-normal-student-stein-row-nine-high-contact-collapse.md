---
id: THM-4315
title: "Source-normal Student--Stein row-nine and high-contact collapse"
status: >
  PROVED RELATIVE TO THM-4308 AND THM-4312 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED.  The source-normal bracket tangent cokernel at row
  m is stationary expectation for an explicit Student/Pearson diffusion.
  At row nine this gives one exact scalar obstruction E9.  On the cubic-
  corner k=1 locus, E9 cuts out ten distinct finite source points, the exact
  row-nine depth fibre above each point is affine seven-space, and all ten
  points have L1 nonzero; the j=1728 high-contact locus cannot reach row
  nine.  Although the row-eight and row-nine fibres both have dimension
  seven, truncation from the latter has singleton image.  An exact F19
  proposal therefore supplies a hostile to probabilistic-extinction
  arguments.  This is a finite row-nine theorem, not a row-ten or all-row
  lift, seam-entry, JC(2), or DC(2) result.
source: root / planar-Jacobian stochastic-process bridge session, 2026-09-01
audit: >
  PASS.  The primary SymPy certificate reconstructs the Student--Stein
  operator, invariant density, posterior filter, literal G9, bracket
  obstruction, row-nine depth modules, corner quintic, high-contact Bezout
  witness, fibre truncation, and a good-prime stochastic hostile.  A
  dependency-free standard-library Fraction/sparse implementation derives
  the same data, including all characteristic-19 ranks and augmented ranks.
  Normal and optimized runs byte-match both frozen transcripts.
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4312-source-normal-cubic-corner-repeated-face-collapse
related:
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
  - THM-3499-regular-shortlex-languages-have-logarithmic-density
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-4316-source-normal-row-ten-cubic-corner-extinction
primary_script: 04-computation/jc2_source_normal_student_stein_row9_thm4315.py
primary_output: 05-knowledge/results/jc2_source_normal_student_stein_row9_thm4315.out
primary_script_sha256: c9625d0b974c4c579388ec14bcafaacfac5bf4a525d3d02c649d1bbf510b29fa
primary_output_sha256: e67eb0aa3435a3e36fe7b0a4fd2a9081a3ad9def759b5fdcee646477aabdeca6
independent_audit_script: 04-computation/jc2_source_normal_student_stein_row9_thm4315_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_student_stein_row9_thm4315_independent_audit.out
independent_audit_script_sha256: a4012c1037084e088953ab758b1e69ee8e66030088d201ffd0fa9cef7d9e2f56
independent_audit_output_sha256: ad2e59e9a0910eb1495db049430b2c6cebd12c5d09ddd283c7a83ceba4bfe74d
hash_basis: raw LF bytes
---

# THM-4315 -- Source-normal Student--Stein row-nine and high-contact collapse

**PROVED RELATIVE TO THM-4308 AND THM-4312 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.  THIS IS AN EXACT FINITE ROW-NINE PROJECTION THEOREM.
IT EXCLUDES THE `L_1=0` / `j=1728` CUBIC-CORNER BRANCH AT ROW NINE, BUT IT
DOES NOT ASSERT A ROW-TEN OR ALL-ROW LIFT, SEAM ENTRY, `JC(2)`, OR `DC(2)`.**

## 1. Statement, inheritance, and the stochastic bridge

Work over an algebraically closed field `k` of characteristic zero in the
fixed source-normal, residual-weight-at-most-twelve gauge of THM-4308.  Let
`v_n=(A_n,C_n)`, retain its degree caps, and write

```text
q=-(x^2+6)/2.
```

THM-4308 proves that changing the terminal tangent before source row `m` by
`theta v_0'`, with `deg(theta)<=m-1`, changes the next compatibility by

```text
D_m theta=q' theta-(q/m)theta'
         =[(x^2+6)theta'-2mx theta]/(2m).                 (1)
```

The inheritance pass was:

- closest proved mechanism: THM-4308's one-dimensional row cokernel and
  exact row-eight depth fibre;
- canonical hostile: THM-3163 and MISTAKE-354, which show that an abstract
  posterior Markov realization is automatic and carries no source algebra;
- corrected near miss: equal affine dimensions at two horizons do not imply
  that randomness transports between them;
- least-used sidecar: the actual cokernel functional, retained with its
  degree cap and Poisson primitive rather than only its scalar value;
- live concepts: bracket tangent, stationary law, posterior killing,
  row-eight cloud, unique spine, corner quintic, and high-contact cubic.

The precise connection is not merely an analogy.  It sends a row-defect
polynomial `f` to the same polynomial viewed as an observable under the
probability law `mu_m` below.  It preserves exactly the predicate
"`f` lies in `im(D_m)`", destroys the gauge primitive and all source/depth
labels, and requires `theta`, the degree cap, and the depth modules as
sidecars.  Its cheapest decisive test is the finite moment row in `(5)`.

## 2. The Student--Stein row theorem

For every `m>=1`, put

```text
mu_m(dx)=c_m (x^2+6)^(-(m+1)) dx                  on R.   (2)
```

Then `mu_m` is stationary for the Pearson diffusion generator

```text
L_m h=(x^2+6)h''-2mxh'.                                  (3)
```

Indeed, if `rho_m` is the density in `(2)`, then

```text
d[(x^2+6)rho_m]/dx+2mx rho_m=0.                          (4)
```

Equivalently, `2mD_m` is a Stein operator.  Integration by parts gives
`E_mu_m[D_m theta]=0` for `deg(theta)<=m-1`; the boundary term vanishes in
that entire range.  The map `D_m:k[x]_(<=m-1)->k[x]_(<=m)` is injective by
THM-4308.  Both its image and the kernel of the nonzero expectation
functional have dimension `m`, so

```text
f in im(D_m)                 iff E_mu_m[f(X)]=0.          (5)
```

For the polynomial degrees used here, the rational moments are

```text
E_mu_m[X^(2r)]
 =6^r(2r-1)!! / product_(j=1)^r(2m-2j+1),
E_mu_m[X^(2r+1)]=0,

with 0<=2r<=m in the first line and 0<=2r+1<=m in the second. (6)
```

Thus `(5)` is a rational coefficient identity and transports to every
characteristic-zero field.  It reproduces THM-4308's rows `m=5,...,8`; at
row nine its primitive integer generator is

```text
(12155,0,4290,0,5148,0,11880,0,45360).                  (7)
```

There is also an exact row-index transition.  Given `X~mu_m`, retain `X`
with probability

```text
6/(X^2+6).                                               (8)
```

The survival probability is `(2m+1)/(2m+2)`, and the conditional survivor
has law `mu_(m+1)`.  The density ratio canonically yields `(8)` as the
maximal same-state keep filter, normalized to have retention one at `x=0`.
Unlike THM-3163's universal posterior construction, it preserves the row
functional; it still does not remember the source or depth sidecars.

## 3. The exact row-nine bracket gate

Direct expansion of the fixed source gives

```text
G_9=(20U+10W+4Z)x^6
    +(10alpha_11+6beta_11)x^7
    +(5upsilon_5+4xi_10)x^8
    +(eta+zeta_3)x^9.                                    (9)
```

Apply `(7)` to `G_9` minus the bracket prediction from rows zero through
eight.  Before the Hasse gate, the resulting numerator is

```text
E9_full=
 68294026800 Delta^2+3653910000 Delta Theta
-5288166000 Delta xi_10+176911616000 Delta
+1547488800 Phi^2+3987447750 Phi alpha_11
+24602292000 Phi eta+4222003500 Phi zeta_3
+2258685000 Theta^2+225494104640 Theta
+1993723875 eta^2+263331993600 xi_10
+105193437167616,                                         (10)
```

with Student obstruction `E9_full/328050`.  On THM-4308's exact gates

```text
Delta=896/15,        Theta=512/75,       zeta_3=-3Phi/2,
```

equation `(10)` is `-39/5` times

```text
E9=
 613527750 Phi^2-511211250 Phi alpha_11
-3154140000 Phi eta-255605625 eta^2
+6736896000 xi_10-46483785515008.                        (11)
```

The seven-dimensional THM-4308 terminal depth fibre is mapped with rank
seven by the restricted `G_9` equation.  Hence:

> **Row-nine bracket iff.**  At a fixed source point of the THM-4308 gate,
> a degree-capped bracket extension through `G_9` exists if and only if
> `E9=0`.  When it exists, its truncation to row eight is one uniquely
> selected point of the old affine seven-space.

This statement is bracket-theoretic.  Row-nine depth membership is audited
on the cubic-corner intersection in the next section.

## 4. Cubic-corner intersection and exact depth fibre

On THM-4312's surviving `k=1` corner, one has

```text
xi_10=(4343625Phi^2+124805668864)/12798000,

eta=(2091705253888/107983125-(2839/1185)Phi^2)/Phi,

U=-13(820125Phi^2+13056802816)/57591000,
c_2=11(14625Phi^2+404652032)/474000,

beta_11=-alpha_11,              alpha_11^2=4Uc_2,        (12)
```

with `Phi`, `U`, and `c_2` nonzero.  Equation `(11)` fixes

```text
alpha_11=
 [6971519208442078125 Phi^4
 -14082869793796263936000 Phi^2
 -74378924775425263164981248]
 /[396452079682031250 Phi^3].                            (13)
```

Put `X=Phi^2`.  Squaring `(13)` and using `(12)` gives the primitive
quintic

```text
Q(X)=
316016952601619726458584136962890625 X^5
+14163685391496771581808513584548828125000 X^4
+137633556412285429978153325875719168000000000 X^3
-6709927871370175861935855782936495259648000000 X^2
+16759499408238096044037088463875607198378754048000 X
+44257795605986960276324945338517826145242081533100032.
                                                               (14)
```

Exact Euclidean algorithms give

```text
gcd(Q,Q')=1,
gcd(Q,X(820125X+13056802816)(14625X+404652032))=1.        (15)
```

Therefore `(14)` has five distinct allowed `X`-roots.  Each has two
distinct square roots `Phi`, and `(12)--(13)` then fix every remaining
source coordinate.  Thus the row-nine source gate consists of exactly ten
geometric points over the algebraic closure.

At each of those points, the exact row-nine projected depth universes are

```text
pi_9(P_2): 75 ambient rows, 160 columns, rank 59,
            left nullity 16;

pi_9(P_3): 85 ambient rows, 251 columns, rank 73,
            left nullity 12.                             (16)
```

The combined 28 residuals have rank three on the ten row-nine tangent
coordinates, with the augmented rank also three and pivots in coordinates
`7,8,9`.  They impose no further source equation and leave a new affine
seven-space `J_9`.

## 5. High-contact collapse

THM-4312's genuine high-contact equation `L_1=0` is

```text
R(X)=
7547170421607067494140625 X^3
+164114458618573873612800000000 X^2
+2284603892441775363795663716352000 X
+2970579390109346274816679296272171008.                  (17)
```

The exact gcd is `gcd(Q,R)=1`.  A compact witness reduces the two
polynomials modulo 19 to

```text
Q=-5X^5+5X^4-6X^3+2X^2-7X-4,
R= 3X^3+3X+8,

(-X^2+9X+1)Q
+(-8X^4+4X^3-X^2+9X+3)R=1              in F_19[X].       (18)
```

Consequently no `L_1=0` point reaches the row-nine bracket gate.  All ten
points in `(14)` lie on the `L_1!=0`, first-carrier `j=0` branch.  This does
not say that THM-4312's row-eight high-contact locus was empty; it dies only
when the row-nine condition is imposed.

## 6. Terminal clouds, a unique spine, and the finite-field hostile

At each fixed row-nine source point,

```text
J_8 isomorphic A^7,          J_9 isomorphic A^7,
dim image(J_9 -> J_8)=0.                                  (19)
```

The last equality is the rank-seven `G_9` selection in Section 3.  Hence the
pushforward of every probability law on `J_9` is a Dirac mass on `J_8`.
Uniform finite-field laws, Gaussian tangent laws, or any other nondegenerate
full-fibre marginals on both spaces cannot be projectively consistent.
Equal affine dimensions merely mean that seven old directions die and seven
new terminal directions appear.

More generally, fixed source rows and the inherited boundary data admit at
most one all-row degree-capped formal bracket lift.  If two such lifts first
differ in row `n`, their difference is `theta_n v_0'`.  Equality of their
next prescribed source row gives `D_(n+1)theta_n=0`; injectivity of `(1)`
forces `theta_n=0`, a contradiction.  Thus a finite terminal cloud can have
at most one point on an infinite spine.  Almost-sure death of a randomly
sampled leaf cannot prove that the spine is absent.

The independent audit makes that warning exact over `F_19`.  The source
point

```text
(Delta,Theta,K,Phi,eta,zeta_3,upsilon_5,xi_10,
 U,W,Z,alpha_11,beta_11)
= (4,1,5,8,6,7,9,4,12,14,12,15,4)                       (20)
```

satisfies all row-five through row-nine gates, the corner square relation,
and every projected-rank condition; `X=Phi^2=7`, the three forbidden factors
are `(7,10,10)`, and `L_1=5`.

For one explicitly declared proposal, choose each determinant tangent
uniformly and independently, then apply the depth filters.  Its successive
codimensions are

```text
(5,6,7,8,2,7,3).                                       (21)
```

Hence its acceptance fractions are `19^(-28)` through row eight and
`19^(-38)` through row nine.  If `I_i` is survival through stage `i` and
`s_i` is the partial sum of `(21)`, then

```text
Z_i=19^(s_i) I_i                                        (22)
```

is an exact mean-one likelihood martingale.  Point `(20)` survives despite
the tiny proposal probability.  These fractions are coordinate- and
proposal-dependent rank checks, not intrinsic probabilities on source data
and not an extinction theorem.

This is the sharp boundary with the older stochastic canon.  THM-3163 says
that an abstract posterior chain is always available; Section 2 is useful
because it intertwines the bracket operator.  THM-3499 shows how genuine
finite-state recurrence can control density, but here a density-zero unique
spine is still possible.  THM-2370's clone hostile likewise warns that an
average or squared energy can lose the unique terminal orientation.

## 7. Reproduction and scope

Run

```text
python3 -B 04-computation/jc2_source_normal_student_stein_row9_thm4315.py
python3 -B -O 04-computation/jc2_source_normal_student_stein_row9_thm4315.py
python3 -B 04-computation/jc2_source_normal_student_stein_row9_thm4315_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_student_stein_row9_thm4315_independent_audit.py
```

The normal and optimized streams must byte-match the declared frozen
outputs.  The independent path imports neither the primary certificate nor
a computer-algebra package.

The theorem proves the exact row-nine bracket gate, and the exact row-nine
depth projection only on its cubic-corner `k=1` intersection.  It does not
prove a row-ten extension, membership in the infinite depth modules,
polynomial termination, existence of a Keller pair, entry into the reduced
seam, `JC(2)`, or `DC(2)`.  THM-4316 supplies the separate row-ten test.

**QED.**
