---
id: THM-3651
title: "Russell-cylinder degree-seven double-zero sixth-order closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the two
  conjugate degree-seven folds on both the principal zero-second-debt and
  zero-fourth-debt walls, exact target-ring coefficients
  give an actual lift J_0=1, J_1=J_2=0.  Separately, the complete arbitrary-
  target-two-form local systems survive total source order five but become
  inconsistent at total order six.  A 25-support exact left certificate over
  Q(sqrt(94)) proves the latter obstruction.  This is a fixed-compiler,
  fixed-fold local closure; it constructs no Keller pair and proves no case
  of JC(2).
source: root/rebuild-thm3651/2026-08-24
depends_on:
  - THM-3650-russell-cylinder-principal-zero-second-debt-fourth-stable-hyperplane
related:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
  - THM-3649-russell-cylinder-qstar-actual-lift-fourth-stable-closure
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
script: 04-computation/jc2_russell_cylinder_degree_seven_double_zero_sixth_order_closure_thm3651.py
output: 05-knowledge/results/jc2_russell_cylinder_degree_seven_double_zero_sixth_order_closure_thm3651.out
script_sha256: 782d770fcb890f2459943ef66eed03903f0f7d3d92b9ce56a249647c5cfe36c5
output_sha256: 276db219037b6929554fc30da9924f274428d787c1d8a9afc567c33453aea381
historical_lift_script_sha256: acc022c0715d262d55cfd2dc23d8bae382d1058fcfdc90dc36a06b2ba143f5f4
historical_lift_transcript_sha256: 70e5e2876e1175e0ee18328d07e86451ac06db8feab0c7575cd85b6f52773946
certificate_plus_sha256: acd75c71e93c95a163fdc91b5d241f3041a7ea98ef6118fee237e974a224a753
certificate_minus_sha256: 599a55ad0e15634dfde922c17c8e98a397a7f2b141a51b8b5ee36707730cb4dc
hash_basis: >
  Raw LF bytes for files.  K-coefficient hashes use a_num/a_den+s*b_num/b_den;
  target-vector hashes bind nonzero coefficients to nested (b,c,e) exponents;
  certificate hashes bind nonzero weights to canonical local row labels.
audit: >
  PASS.  Static hostile audit passed after adding field, compiler,
  differentiated-sidecar, common-target, truncation, normalization, and
  fold-parameterization gates.  The frozen combined companion passed its
  final normal and optimized replays at 1300 active gates; both transcripts
  are byte-identical to the 3514-byte stored LF output, with SHA-256
  276db219037b6929554fc30da9924f274428d787c1d8a9afc567c33453aea381.
  Independent hostile dynamic controls passed.  The recovered historical
  lower-lift companion also replayed all 970 gates with its pinned byte-
  identical transcript.
---

# THM-3651 -- degree-seven double-zero sixth-order closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The theorem has
two deliberately separate assertions.  First, each of two conjugate
algebraic folds has an actual target-ring lift only
through `J_2`.  Second, the larger space of arbitrary local target two-forms
survives through total source order five and fails at order six.  The second
assertion obstructs a constant Jacobian for the fixed folds, but it does not
promote the first lift to orders three, four, or five.

All algebra below is exact.  Let

```text
K=Q(s),                         s^2=94.                (1)
```

Every statement for the plus sign has a minus-sign mate obtained by the
field automorphism `s -> -s`.

## 1. Compiler and the two double-zero folds

Use the Russell-cylinder target ring and compiler

```text
R_K=K[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),       (2)
```

followed by the quadratic fold

```text
q=Q(x)+t^2,                         w=t.                (3)
```

Put

```text
Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,
P=x^2(x^2-1)^2,

beta_+/-=-211/130 +/- (99/260)s,
alpha_+/-=-1683/260 -/+ (11/65)s,

Q_+/-=Q_1+P(alpha_+/- + beta_+/- x).                  (4)
```

Thus, written without paired-sign notation,

```text
beta_+=-211/130+(99/260)s,     alpha_+=-1683/260-(11/65)s,
beta_-=-211/130-(99/260)s,     alpha_-=-1683/260+(11/65)s. (5)
```

The companion verifies directly that both degree-seven polynomials have

```text
Q_+/-(-1,0,1)=(-3,-3/4,-3),
(Q_+/-)'(-1,0,1)=(-9/2,1,9/2).                        (6)
```

Moreover,

```text
alpha=-(259+16beta)/36,
520beta^2+1688beta-5717=0.                            (7)
```

Thus THM-3650 applies: both folds lie on the principal zero-second-debt
plane and on its zero-fourth-debt hyperplane.  The triple is ordinary because
the retained tangent directions are `(1,-9),(1,4),(1,9)`.

## 2. Actual target-ring lift through `J_2`

For the plus fold, there are coefficients `F_i,G_i in R_K` such that

```text
F^#=c+wF_1+w^2F_2+w^3F_3,
G^#=e+wG_1+w^2G_2+w^3G_3,                             (8)
```

and, on pulling back through `(3)`,

```text
Jac_(x,t)(F^#,G^#)=sum_(n>=0)t^nJ_n(x),

J_0=1,                         J_1=J_2=0              (9)
```

as identities in `K[x]`.  Conjugating all coefficients proves the same
statement for `Q_-`.

The certificate uses raw target monomials `b^i c^j e^k` in nested increasing
`(i,j,k)` order.  Their restriction degrees are

```text
(deg gamma(b),deg gamma(c),deg gamma(e))=(27,19,16).  (10)
```

The cutoffs below select one reproducible solution; no minimality is claimed.

### 2.1 The `J_0` square

At raw weighted cutoff `160`, there are `141` monomials.  The ascending
coefficient system has `179` rows and `282` columns, ordered as

```text
[c' gamma(m)]_m,                 [-e' gamma(m)]_m,    (11)
```

for the `G_1,F_1` blocks.  Reduction modulo

```text
p=1000033,                    s=20857 mod p,
20857^2=94 mod p                                      (12)
```

has operator and augmented rank `178`.  Its deterministic skeleton is

```text
pivot-column hash
7784507c6b1b77d78bfbd33fa2385eaa0c5a0c4e1a457c50faa6c27737a07768,

pivot-row hash
2113e21fc204fdc9fc3be762bc1d022a2388bb14a8390ce4b8603ac661de7454. (13)
```

The selected `178 x 178` square is rebuilt over `K`.  Replacing each entry
`a+sb` by the faithful rational block

```text
[ a   94b ]
[ b    a  ]                                             (14)
```

gives a `356 x 356` rational square.  Its exact solution and the full,
unselected polynomial residual prove

```text
c' gamma(G_1)-gamma(F_1)e'=1.                          (15)
```

The representative-bound data are

```text
F_1: support 89,
 target  0ee493528fa3b1bf9ca6a824af0a6b7cdfd1588c2a99f358b171a7cea0f219d8,
 gamma degree 160,
 gamma   89178122fa5660177bc2ec8f2c097232fee112efa5ccba496f5ec1a2574e13bc,
 delta degree 153,
 delta   ce48af3830f332279a6f9c6e5ad26883baf59d9b3c33ce488901f6d1e247f019;

G_1: support 86,
 target  6b02b7b4af6d9a50ec526e0367802a811ebfbc03ab263aa660cd67e772346a38,
 gamma degree 157,
 gamma   1e21a0e832129ef01b121061019923d04df1d7faae2f7440601dc516e5049d3f,
 delta degree 150,
 delta   a53ae91d4966f864882344580697213481687c81da1ee13c96f840a1c4e75d26. (16)
```

Here `delta` is the `q`-derivative of the chosen raw target representative
after restriction.  Freezing it with the target coefficient vector prevents
an illicit replacement by an equal restriction with a different lift jet.

### 2.2 The coupled `J_1,J_2` square

At cutoff `301`, there are `744` monomials.  Order the unknown blocks as

```text
[F_2,G_2,F_3,G_3]                                    (17)
```

and put all `J_1` rows before all `J_2` rows, each in ascending `x` degree.
The row envelopes are `320+461=781`, and the operator has `2976` columns.
Modulo `(12)`, operator and augmented rank are both `776`.  The pivot hashes
are

```text
columns a9697a7a900ca6ddb98a14981d01847456e04cfb48932e8c300d07065c7084d6,
rows    449e1c433c85fe550e6dda35664a1584b70fe5c4da4846003173ade228fbbfef. (18)
```

The modular computation only selects the skeleton.  The selected `776 x 776`
square over `K` is solved through its faithful `1552 x 1552` rational block,
then substituted into every row of the full `K[x]` equations.  Both complete
residuals vanish.  The chosen higher coefficients are bound by

```text
F_2: support 230,
 target 7e38ac68a6cbd4702e1c50ba566849fa9d0f9198b058523111f1a6c0070d863c,
 gamma degree 301,
 gamma  7e7c7072644f5870d9798a1b157759f92bfb4406a4d8deaf560c76e3f53b1483;

G_2: support 228,
 target 0acaf6b4655dc7cea05fae327c4cb89de446c3c618b6768c1f290be22ce4f2b9,
 gamma degree 298,
 gamma  884f6b51f22cd3f22a1834f9663a9e7146823b24fd04cbeb8b28e08fe684d830;

F_3: support 230,
 target c4007f9663bbd876cb663afbb796fd7e7904f346e8c6746e60cec8175fe5987e,
 gamma degree 301,
 gamma  d066e6a37d2d878dc648f937c2467863d6561985d6d1038323263045d79ad62c;

G_3: support 86,
 target d7b163e2ec6949316a0dcccdbf105302dc69e36bbcc05ec67d61a2979b044ccb,
 gamma degree 300,
 gamma  757b39f18070c1bdb46a9eb199edc42dd0efaf7b4f1d0aba1866f2c1b1dd3a92. (19)
```

This proves `(9)`.  Since the folds satisfy `(7)`, substituting `J_0=1` and
`J_2=0` in PROVED THM-3650 also gives the retained identity

```text
lambda(J_4)=0.                                        (20)
```

Equation `(20)` is one three-branch quotient value.  It is **not** the
polynomial identity `J_4=0`, and this construction proves neither `J_3=0`
nor any continuation past `J_2`.

## 3. Complete arbitrary-two-form local systems

This section starts a different problem.  At the common target point put

```text
y=c/3,                         z=e+3.                 (21)
```

Then `(y,z,w)` are regular local target coordinates.  Write one arbitrary
target two-form germ

```text
Omega=A(y,z,w) dy^dz+B(y,z,w) dy^dw+C(y,z,w) dz^dw.  (22)
```

Indeed the common point is `(b,c,e,w)=(0,0,-3,0)`, and the derivative with
respect to `b` of `c^2e-b(b+4)` is `-4` there.  Thus eliminating `b` makes
`(c,e+3,w)`, equivalently `(y,z,w)`, a regular coordinate system.

The coefficients of `A,B,C` are **common target coefficients at the one
shared target point**.  They are evaluated along all three source branches;
the matrix does not grant three independent coefficient jets branch by
branch.

At a retained branch `a in {-1,0,1}`, let `h=x-a`.  If the pullback density is
`j_a(h,t)`, define `M_N` to map the common coefficient monomials of `(22)` of
total degree at most `N` to all Taylor coefficients

```text
[h^i t^j]j_a,              i+j<=N.                  (23)
```

These are Taylor coefficients, not raw mixed derivatives; a raw derivative
would differ by `i!j!`.  The row order is branches `-1,0,1`; then total degree
ascending; then `t^d,h t^(d-1),...,h^d`.  The column order is the three slots
`A,B,C`, followed by increasing `(deg_y,deg_z,deg_w)`.

The exact universe sizes are

```text
rows    =3 binom(N+2,2),
columns =3 binom(N+3,3).                              (24)
```

Every coefficient monomial beyond degree `N` is invisible because `y,z,w`
vanish at the target point.  There is one essential truncation precaution:
`y` and `z` must be retained through total degree `N+1` **before** their
derivatives are taken.  For `N=6`, the companion retains degree seven and
checks the regression coefficient

```text
[h^6] y_h at a=-1
 =1253248661/4225+s(477739493/25350).                 (25)
```

Premature degree-six truncation would erase this term.

Use the denominator-clearing target vector

```text
j_a(h,t)=12 mod (h,t)^(N+1)             for all a.   (26)
```

Scaling by `1/12` gives the equivalent unit-density problem.  Faithful
rational block rank computations give

```text
        shape       rank(M_N)   rank(M_N|12)
N=4     45 x 105        40            40
N=5     63 x 168        57            57
N=6     84 x 252        77            78.            (27)
```

Therefore one common arbitrary target two-form jet matches constant density
on all three branches through total order five, whereas none does through
total order six.  The `N=5` solution is not asserted to be closed,
decomposable, or of the form `dF^# wedge dG^#`; it must not be combined with
the actual lift of Section 2.

## 4. The 25-support sixth-order certificate

The following row `ell_+` annihilates every one of the `252` columns of
`M_6`.  An entry `(a,i,j):u+sv` is the weight on the Taylor coefficient
`[h^i t^j]j_a`, not on a raw derivative.  Omitted entries are zero.

```text
(-1,0,0) -416265780017/23562825 +s*129267792/54925
(-1,1,0)    2464348103/966680   +s*75921/10985
(-1,0,2)  -19437583987/5800080  +s*163017/21970
(-1,2,0)             6975/143    +s*20/3
(-1,1,2)           -47295/572    -s*95/13
(-1,3,0)                25/11
(-1,0,4)             39555/286    +s*90/13
(-1,2,2)              -75/22
(-1,1,4)              225/44
(-1,0,6)             -675/88

( 0,0,0) -160054491789/7854275
( 0,1,0)        1342208/1859      +s*10944/169
( 0,0,2)    46861679809/4833400   -s*163017/21970
( 0,1,2)             1944/143      -s*108/13
( 0,0,4)           -69093/143      -s*90/13
( 0,0,6)             1215/44

( 1,1,0)     -1521051607/371800
( 1,0,2)   -14152473763/2230800
( 1,2,0)             1467/11       +s*4/3
( 1,1,2)             9459/44       +s
( 1,3,0)              -65/11
( 1,0,4)             7587/22
( 1,2,2)             -195/22
( 1,1,4)             -585/44
( 1,0,6)            -1755/88.                         (28)
```

The table uses the right-hand-side normalization `12` from `(26)`.  Its
pairing is

```text
<ell_+,12>=
 -275824386272/604175+s(1551213504/54925),            (29)

Norm_(K/Q)<ell_+,12>
 =177370060592178176/1329185 !=0.                    (30)
```

For unit right-hand side, the pairing is one twelfth of `(29)`:

```text
<ell_+,1>=
 -68956096568/1812525+s(129267792/54925) !=0.         (31)
```

Thus `(28)` proves the strict rank jump in the last row of `(27)` without
depending on rank software.  Conjugating every coefficient gives the minus-
fold certificate.  With the hash convention in the front matter,

```text
ell_+ acd75c71e93c95a163fdc91b5d241f3041a7ea98ef6118fee237e974a224a753,
ell_- 599a55ad0e15634dfde922c17c8e98a397a7f2b141a51b8b5ee36707730cb4dc. (32)
```

The nonzero norm in `(30)` stays nonzero under both embeddings `K -> C`, and
the exact nonzero rank minors used in `(27)` stay nonzero under either
embedding as well.  Hence the obstruction holds for both resulting complex
folds, not merely for a formal presentation over `K`.

The hostile deletion is active.  Removing only the last weight
`(1,0,6)=-1755/88` leaves nonzero pairings on `47` of the `252` columns; on
the named column `(dy^dz)w` the residual is

```text
-5265/22 !=0.                                         (33)
```

## 5. THM-3683 specialization cross-control

This section is a cross-control, not a proof dependency on independently
hostile-audited THM-3683.  First, the companion types the fold into the
parameterization used by that theorem.  With

```text
Q_6=Q_1-(259/36)P,              R_1=P(1-x^2),
R_2=P(4-9x),

Q_beta=Q_6+beta P(x-4/9),
Q_(p,r)=Q_6+pR_1+rR_2,                                  (34)

p(r)=(520/9)r^2-(1688/81)r-5717/729,                  (35)
```

equation `(7)` gives

```text
p(-beta_+/9)=0,
Q_beta_+=Q_(0,-beta_+/9).                              (36)
```

Thus the cross-control really specializes at the same degree-seven fold,
not merely at a numerically related parameter.  Inheriting THM-3683's
displayed quartic formula, put `r=-beta_+/9` and evaluate

```text
F(r)=72783360r^4-77822208r^3-28419741r^2
                         +7849770r-1276420.           (37)
```

The present companion independently performs the specialization arithmetic
and checks

```text
F(-beta_+/9)
 =2187(21544632beta_+-97639283)/33800.                (38)
```

Consequently the cross-control value

```text
D_6=-256F(r)/3^13
 =275824386272/200201625-s(210658624/2471625)         (39)
```

satisfies

```text
<ell_+,12>=-(3645/11)D_6.                             (40)
```

This exact agreement connects the complete total-degree obstruction to the
retained even-order debt calculation.  It is not used to prove `(27)--(31)`.

## 6. Mechanism, boundary, and non-consequences

The mechanism is a genuine extension failure.  At total orders four and
five, constant density lies in the image of the complete common-target-two-
form jet map.  At order six, one new source layer creates a left-cokernel
class whose constant response has nonzero field norm.  Retaining degree seven
before differentiation is the missing coordinate that makes the class
visible.

The result implies for each fixed fold `Q_+/-` that no
arbitrary target two-form, hence no form `dF^# wedge dG^#` from an actual
target pair, can pull back to a nonzero constant density.  It does not imply
any of the following:

- that the Section 2 actual pair has `J_3=0`, polynomial `J_4=0`, or `J_5=0`;
- that the arbitrary-form survivor through order five is an actual pair;
- a classification of other points on the THM-3650 zero-fourth hyperplane;
- a statement for another compiler, a nonquadratic or noncritical fold, or
  another retained tangent chart;
- a polynomial Keller pair or a proof or refutation of `JC(2)`.

The planar Jacobian conjecture remains **OPEN**.

## 7. Reproduction and audit state

Reproduce the combined package with

```bash
python3 -B 04-computation/jc2_russell_cylinder_degree_seven_double_zero_sixth_order_closure_thm3651.py
python3 -O -B 04-computation/jc2_russell_cylinder_degree_seven_double_zero_sixth_order_closure_thm3651.py
```

The expensive exact step is the `1552 x 1552` rational block solve.  The
companion reconstructs all target restrictions and `q`-derivative sidecars,
checks the pinned supports and hashes, verifies the full unselected `K[x]`
residuals, builds all three complete local matrices from the literal compiler,
checks every one of the `252` certificate annihilations, runs the hostile
deletion, evaluates the quartic cross-control, and checks an assertion-free
AST.

As an independent historical control, the recovered standalone lower-lift
script with SHA-256

```text
acc022c0715d262d55cfd2dc23d8bae382d1058fcfdc90dc36a06b2ba143f5f4
```

was freshly replayed under concurrent load.  All `970` gates passed and its
stdout SHA-256 was exactly

```text
70e5e2876e1175e0ee18328d07e86451ac06db8feab0c7575cd85b6f52773946.
```

That historical control confirms Section 2 only.  The frozen combined package
passed independent static and dynamic hostile audit, and its normal,
optimized, and stored transcripts are byte-identical.  **QED.**
