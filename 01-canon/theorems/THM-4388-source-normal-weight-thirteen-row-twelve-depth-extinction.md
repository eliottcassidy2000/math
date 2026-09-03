---
id: THM-4388
title: "Source-normal weight-thirteen row-twelve depth extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4389 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. In the complete fixed THM-4308 source-normal
  residual-weight-at-most-thirteen family, the Phi=0 branch dies by
  row-twelve bracket compatibility. The Phi!=0 branch has exactly twenty-six
  reduced geometric row-twelve bracket source-projection points, all
  disjoint from the required row-twelve joint projected P_2/P_3 depth
  equation. Higher residual weights, chart or seam entry, all-row lifting,
  polynomial termination, Keller pairs, JC(2), and DC(2) remain OPEN.
source: root + weight13_anchor + independent referee / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4389-source-normal-weight-thirteen-row-ten-nonisotrivial-elliptic-pencil
related:
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-4380-source-normal-weight-twelve-row-twelve-extinction
  - THM-4385-source-normal-row-ten-elliptic-sign-quotient
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
mistake_firewall:
  - MISTAKE-287
  - MISTAKE-540
  - MISTAKE-541
primary_script: 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388.py
primary_output: 05-knowledge/results/jc2_source_normal_weight13_face_row12_extinction_thm4388.out
primary_script_sha256: 4dd8ef97659205aa5a8dd9e430d003b9c573108033f020e3cd1095d97b6572be
primary_output_sha256: b48bfc5f18e8ad5ab32a148cbc2f789068869e79391eb207615c2efba5861d14
independent_audit_script: 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.out
independent_audit_script_sha256: 57a3bf13ef33f7314514de728e6f249147c0a8e1dcbd0134c65764d49c87e375
independent_audit_output_sha256: 1df8b3e33eb164c83810597813d30cd53c59198404a6235ec1cca887d2a7eeca
hash_basis: raw LF bytes
audit: >
  PASS. The primary imports only the audited THM-4308/4315 row operators,
  reconstructs the literal enlarged source, and checks both Phi strata,
  every new bracket selector, the complete projected P_2/P_3 universes,
  localization boundaries, the affine quotient scheme, and the final
  Bezout extinction. The independent audit does not import the primary; it
  reorganizes the linear eliminations and executes 172 explicit checks.
  It does reuse the audited THM-4308/4315 row operators and is therefore not
  called import-free. Normal, optimized, and hash-seeded executions of both
  certificates byte-match their frozen LF outputs.
---

# THM-4388 -- source-normal weight-thirteen row-twelve depth extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4389 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS CLOSES THE COMPLETE FIXED THM-4308 SOURCE-NORMAL
RESIDUAL-WEIGHT-AT-MOST-THIRTEEN FAMILY THROUGH ROW TWELVE. CHART OR SEAM
ENTRY, HIGHER RESIDUAL WEIGHTS, AN ALL-ROW LIFT, POLYNOMIAL TERMINATION, A
KELLER PAIR, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Exact universe and theorem

Work over an algebraically closed field `k` of characteristic zero. Retain
the fixed source-normal chart, degree caps, bracket recursion, and projected
depth modules of THM-4308/4315. Extend its residual polynomial by

```text
H_13=H_12+c51*p^5*y+c23*p^2*y^3.                         (1)
```

Since `wt(p)=2` and `wt(y)=3`, the nonnegative solutions of `2a+3b=13`
are exactly `(a,b)=(5,1),(2,3)`. Thus `(1)` is the complete
residual-weight-thirteen face in this fixed chart. The conclusion is

```text
no point of the fixed residual-weight-at-most-thirteen family satisfies
all bracket equations and the required projected P_2/P_3 depth equations
through row twelve.                                      (2)
```

This is not a monotone closure in the weight cap: a term of weight at least
fourteen can enter a row already used below. Statement `(2)` therefore gives
neither higher-weight extinction nor entry into this chart from a putative
global Keller pair.

## 2. Visibility of the odd face

The two new monomials have exact source expansions

```text
p^5*y   =x*t^7*(1+x^2*t)^6,
p^2*y^3 =x^3*t^8*(1+x^2*t)^5.                            (3)
```

With out-of-range binomial coefficients interpreted as zero, their complete
row increments are

```text
Delta G_n=
 [binom(6,n-7)c51+binom(5,n-8)c23]*x^(2n-13),
                                                        7<=n<=13. (4)
```

In particular,

```text
Delta G_7 =c51*x,
Delta G_8 =(6c51+c23)*x^3,
Delta G_9 =(15c51+5c23)*x^5,
Delta G_10=(20c51+10c23)*x^7,
Delta G_11=(15c51+10c23)*x^9,
Delta G_12=(6c51+5c23)*x^11,
Delta G_13=(c51+c23)*x^13.                               (5)
```

Every displayed power is odd. THM-4328's algebraic Student--Stein
functionals kill the direct odd-face contributions. The primary and
independent certificates separately reconstruct the exact recursion and show
that the scalar bracket/depth gates through row nine are indeed unchanged;
the stochastic parity is an explanation of direct invisibility, not the
extinction proof.

The face itself has not been forgotten. Its first labelled channel map is

```text
(h_7,h_8)=(c51,6c51+c23),
c51=h_7,                 c23=h_8-6h_7,                   (6)
```

so the sidecar is lossless and the parameters reappear in later tangent and
depth equations. In inherited notation, the unchanged row-eight gate is

```text
Delta=896/15,       Theta=512/75,       zeta_3=-3Phi/2, (7)
```

and the unchanged row-nine equation is

```text
613527750P^2-511211250P*a-3154140000P*e
-255605625e^2+6736896000X-46483785515008=0,             (8)
```

where `P=Phi`, `e=eta`, `X=xi_10`, and `a=alpha_11`.

## 3. The `P=0` branch

On `P=0`, row ten again leaves the two reduced nonzero roots of

```text
F(e)=18612736875e^2-4820239249178624.                    (9)
```

The inherited graphs include `beta_11=-6e/5`. Unlike the weight-twelve
cap, row eleven does not kill these points. Its remaining bracket equation
has literal `c51` coefficient `323e/432`, a unit on `V(F)`, and gives

```text
c51=
-3016054591963389705285537747767*e
 /6517871654473709985808572150                 in k[e]/(F). (10)
```

Every row-eleven projected-depth residual then vanishes modulo `F`.

At row twelve the bracket selector has rank eight and leaves two equations
linear in `c23`. Their literal `c23` coefficients are

```text
323e/360,                    -5/7.                       (11)
```

Their compatibility determinant reduces modulo `F` to the nonzero scalar

```text
6281860880237041075158211534217247003921830342656
------------------------------------------------------- .
25491226342562384812197955789004609878125               (12)
```

Thus the two equations are incompatible at both roots of `F`: the complete
`P=0` branch dies by row-twelve bracket compatibility.

## 4. The `P!=0` coefficient graphs

Put

```text
D(P,e)=7231154026500P^3+50541940696500P^2e
      +6793915500000Pe^2+353642000625e^3
      -631918028977864704P-91584545734393856e.           (13)
```

This is the cubic defining THM-4385's weight-twelve elliptic carrier. In the
complete weight-thirteen family, row-ten joint depth is instead

```text
D(P,e)-707284001250P^2*c51=0.                            (14)
```

Because `P` is already localized, `(14)` uniquely solves

```text
c51=D(P,e)/(707284001250P^2).                            (15)
```

Hence the old elliptic curve is exactly the transverse slice `c51=0`, not a
higher-weight invariant. THM-4389 refines this correction: fixed-`c51`
fibres form a non-isotrivial genus-one pencil, whereas the total row-ten
carrier is rational.

At row eleven there is one bracket residual, with literal `c23` coefficient

```text
323P/504.                                                (16)
```

It therefore solves `c23` uniquely with only powers of the already inverted
`P` in the denominator. Row-eleven projected depth is automatic. This leaves
a two-dimensional `(P,e)` base for the row-twelve test, with both new face
coordinates fixed by response graphs.

## 5. The saturated row-twelve bracket scheme

The two row-twelve residuals have simultaneous-sign characters, one odd and
one even. Put

```text
r=e/P,                         Y=P^2.                    (17)
```

After dividing the odd residual by `P`, let their primitive quotient
polynomials be `Q_o(Y,r)` and `Q_e(Y,r)`. Their complete coefficients are
frozen in the primary transcript, and

```text
bideg(Q_o)=(3,5),              bideg(Q_e)=(3,4),        (18)
```

where the first degree is in `Y`. The raw resultant factors as

```text
Res_Y(Q_o,Q_e)=u*(9r+8)*K_13(r),             u in Q^*,  (19)
```

with `K_13` irreducible and squarefree of degree thirteen. The factor
`9r+8` is not an affine solution: direct specialization gives

```text
gcd_Y(Q_o(Y,-8/9),Q_e(Y,-8/9))=1.                       (20)
```

It is the common-leading-coefficient degree drop at `Y=infinity`, so it is
not discarded without saturation evidence. The affine lexicographic
Groebner basis has bidegrees

```text
(1,12),                       (0,13),                   (21)
```

and the second member is `K_13`. The first is linear in `Y`; its leading
coefficient, constant term, and `K_13` are pairwise coprime in the ways
needed for a unique nonzero `Y` over each root. Consequently the affine
quotient scheme has exactly thirteen reduced geometric points. Each has two
distinct lifts `P=+-sqrt(Y)`, `e=rP`, after which `(15)--(16)` determine
`c51,c23`. Thus:

```text
the row-twelve P!=0 bracket source projection consists of exactly
twenty-six reduced geometric points.                    (22)
```

They are not asserted rational or real, and they are not generally on
THM-4385's old elliptic curve. The point

```text
F_23: (r,Y,P,e,c51,c23)=(6,9,3,18,6,14)                (23)
```

is a positive bracket-carrier control. It is not used to infer the
characteristic-zero count.

## 6. Row-twelve depth kills all 26 bracket points

The exact row-twelve projection universes are

```text
pi_12(P_2): 117 x 267, rank 87;
pi_12(P_3): 130 x 424, rank 105.                         (24)
```

The joint terminal system has shape `55 x 13`, rank four, and pivots
`(9,10,11,12)`. Its only remaining mismatch is, up to a nonzero scalar,

```text
N(P,e)=272008125P^2-43740000Pe+10213932924928.           (25)
```

The literal joint residual is `-N/288360000`. Independently selecting the
`P_2` and `P_3` solutions gives their mismatch as `N/108135000`.

In quotient coordinates,

```text
N_q(Y,r)=(272008125-43740000r)Y+10213932924928.          (26)
```

If its `Y` coefficient is zero, its constant is nonzero. Otherwise solve
`(26)` for `Y`, substitute into `Q_o,Q_e`, and clear the third power of the
slope. The resulting primitive polynomials have degrees six and five. Their
degree-preserving, unit-normalized reductions modulo `11` are

```text
f_6=3r^6-4r^5-2r^2-2r-2,
f_5=5r^5-5r^3+5r^2+2r+5.                                (27)
```

They satisfy the literal identity

```text
(2r^4+4r^3+3r^2-5r+2)f_6
 +(r^5-3r^4-2r^3-3r^2+5r+1)f_5=1       in F_11[r].    (28)
```

The reductions preserve degree, so `(28)` proves that the characteristic-zero
eliminants are coprime. No one of the twenty-six bracket points satisfies
`N=0`. At the finite control `(23)`, `N_q=4 mod 23`, providing the hostile
depth check. Together with Section 3, this proves `(2)`.

## 7. Quotient typing and stopping boundary

The two-sheet lift in `(22)` is a coefficient-normalization deck, not a map
from the physical source curve. On the live-seam coefficient atlas write

```text
A5=a^5,              Rtilde=R/gamma,
f21=[p^2*y]Rtilde,    f31=[p^3*y]Rtilde,
alpha^2=A5.                                               (29)
```

The normalization bookkeeping gives

```text
P=-alpha*A5^3*f21/2,       e=-alpha*A5^4*f31/2,
r=A5*f31/f21,              Y=A5^7*f21^2/4.              (30)
```

Changing the normalization root `alpha` changes the signs of `(P,e)` and
fixes `(r,Y)`. The coefficients in `(29)` are constant along the physical
source line. Therefore `(30)` supplies no nonconstant
`P^1_X -> Ebar` or `P^1_X` to the thirteen-point quotient, and no
Riemann--Hurwitz argument is available from it. THM-4012's genuine
`P^1`-to-elliptic mechanism has a differently typed target and map; importing
that argument here would lose the required coordinate.

The theorem ends at the fixed finite chart and at row twelve. The next
uncut source operation is the complete residual-weight-fourteen face, whose
channels can enter earlier rows. Separately, the source-normal coefficient
atlas still lacks a proved chart/seam entry map from arbitrary hypothetical
Keller pairs.

## 8. Reproduction

Artifacts:

```text
04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388.py
05-knowledge/results/jc2_source_normal_weight13_face_row12_extinction_thm4388.out
04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.py
05-knowledge/results/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.out
```

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388.py
python3 -B -O 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388.py
python3 -B 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_weight13_face_row12_extinction_thm4388_independent_audit.py
```

Normal, optimized, and hash-seeded streams of both artifacts byte-match the
frozen LF outputs. The independent certificate imports no THM-4388 code, but
it does reuse the audited THM-4308/4315 row operators; its independence is in
the theorem-specific reconstruction and elimination, not the inherited
compiler.
