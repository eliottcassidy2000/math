---
id: THM-3650
title: "Russell-cylinder principal zero-second-debt fourth-stable hyperplane"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the ordinary
  principal zero-second-debt plane, an exact
  endpoint identity computes the retained fourth-stable cokernel for every
  target two-form.  Under J_0=1 and J_2=0, its zero set is a fiberwise affine
  hyperplane in the endpoint third jets.  Exact Q6 and Q* controls miss that
  hyperplane, whereas Qh and Q9 lie on it; Q9 also has an actual finite local
  target-ring control through retained order four.  The scope is only the
  retained three-branch quotient, with no global Keller-pair or JC(2) claim.
source: root/zero-debt-next-order/2026-08-21
depends_on:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
related:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3639-russell-cylinder-all-retained-cells-universal-second-stable-debt
script: 04-computation/jc2_russell_cylinder_principal_zero_second_debt_fourth_stable_hyperplane_thm3650.py
output: 05-knowledge/results/jc2_russell_cylinder_principal_zero_second_debt_fourth_stable_hyperplane_thm3650.out
script_sha256: 40743857d18b08bcccec87de63637ef0582b5a27670e81ec23a666a295151a38
output_sha256: ae5db24875bb66d13132d0d1cc1d2965202505dadf4fbbcb29e5705a35474f8d
hash_basis: raw LF bytes
audit: >
  PASS -- independent hostile reconstruction solved the full endpoint ansatz,
  recovered every displayed coefficient, and checked all 105 target-two-form
  jets, the affine zero hyperplane, beta discriminant, Q6/Q*/Qh/Q9 and uniform
  controls, and the actual local Q9 F4 cancellation.  Normal, optimized, and
  stored transcripts are byte-identical at 212 active gates; hashes, AST0,
  docs, diff, line-ending, and whitespace checks pass.
---

# THM-3650 -- principal zero-second-debt fourth-stable hyperplane

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
generalizes the two fourth-stable identities of
THM-3642 across the ordinary principal zero-second-debt plane of THM-3641.
The exact companion and an independent symbolic reconstruction verify every
displayed identity.

All rings, germs, and closed points are over `C`.  The symbolic certificate is
over a rational function field over `Q` and then base-changes.

## 0. Compiler, fold, and principal jet plane

Use the Russell-cylinder compiler

```text
R=C[b,c,e]/(c^2 e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3).        (1)
```

Put `y=c/3`, `z=e+3`, and pull a target two-form back through

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                             (2)
```

Assume the principal retained values and slopes

```text
Q(-1,0,1)=(-3,-3/4,-3),
Q'(-1,0,1)=(-9/2,u,9/2).                                (3)
```

The retained tangent vectors in `(y,z)` are

```text
(1,-9),                 (1,4u),                 (1,9).  (4)
```

Write

```text
A=9-4u,             B=9+4u,
r_i=Q''(i),          k_i=Q'''(i),        i=-1,0,1.      (5)
```

The triple is ordinary exactly when

```text
u!=+-9/4,                     equivalently AB!=0.       (6)
```

This theorem is restricted to the zero-second-debt plane

```text
A r_-+B r_+=-243,

r_+=(-243-A r_-)/B.                                    (7)
```

Its normalized retained cokernel is

```text
lambda_u(P)=[AP(-1)-18P(0)+BP(1)]/18.                  (8)
```

For a pulled-back target two-form, write the source density as

```text
sum_(n>=0) t^n J_n(x).                                  (9)
```

Primes below are source `x` derivatives.  No decomposability, closedness, or
actual-pair hypothesis is imposed on the two-form identity.

## 1. Full universal fourth-stable endpoint identity

On `(7)`, define

```text
P_-=
 5184r_0r_-+69984r_0
+1280r_-^2u^2+2304r_-^2u-11664r_-^2
+62208r_-u^2+132192r_-u-389286r_-
-2304k_-u^2+11664k_-
-1152k_+u^2-5184k_+u-5832k_+
+622080u^2+1819584u-2184813,

P_+=
 5184r_0r_-+69984r_0
+1280r_-^2u^2-8064r_-^2u+11664r_-^2
+62208r_-u^2-287712r_-u+240570r_-
+1152k_-u^2-5832k_-
+2304k_+u^2+10368k_+u+11664k_+
+622080u^2-2799360u+177147.                            (10)
```

Set

```text
alpha_- = P_-/(26244B),
alpha_+ =-P_+/(26244B),

beta_- =-2(80r_-u^2-72r_-u-243r_-+1728u^2-4374)/(729B),
beta_+ =-2(20r_-u-27r_-+432u-243)/729,

chi_0  = 4(2r_-+27)/(9B),
delta  = (8r_-u-18r_-+216u-243)/243.                  (11)
```

Then every target two-form satisfies the exact retained identity

```text
lambda_u(J_4)=
 alpha_- J_0(-1)+beta_- J_0'(-1)-(A/81)J_0''(-1)
                   +chi_0 J_0'(0)
+alpha_+ J_0(1) +beta_+ J_0'(1) -(B/81)J_0''(1)
+delta[J_2(-1)-J_2(1)]
                   +(A/27)J_2'(-1)-(B/27)J_2'(1).     (12)
```

This chosen representative has no `J_0(0)`, `J_0''(0)`, `J_2(0)`, or
`J_2'(0)` term.  More importantly, `(12)` uses neither `J_1` nor `J_3`.
All coefficients, including the full endpoint coefficients `(10)`, are part
of the statement.

## 2. Constant specialization and the zero hyperplane

Under the polynomial identities

```text
J_0=1,                           J_2=0,                 (13)
```

equation `(12)` reduces to

```text
lambda_u(J_4)=-2N_4/(243B),

N_4=12A r_-^2+486(A-3)r_--AB k_-+B^2k_+
                                      +243(22A-153).   (14)
```

Thus, on the ordinary locus, the retained fourth-stable debt vanishes if
and only if

```text
k_+=[AB k_--12A r_-^2-486(A-3)r_-
                         -243(22A-153)]/B^2.           (15)
```

For fixed `u` and curvature packet, `(15)` is a genuine affine hyperplane in
the endpoint third jets: its `k_+` coefficient before division is `B^2!=0`.
Globally in `(r_-,k_-,k_+)` it is the displayed polynomial hypersurface.
The central jets `r_0,k_0` cancel from `(14)`.

There is therefore no universal fourth-stable nonvanishing theorem on the
zero-second-debt plane.  The order-three Hermite evaluation matrix at
`-1,0,1` has rank `12`; a polynomial of degree at most `11` realizes every
prescribed value, slope, curvature, and third-jet packet.  In particular,
both sides of `(15)` occur for genuine collision polynomials.

## 3. Why 105 monomials prove the arbitrary-two-form identity

In the regular local target coordinates `(y,z,w)`, write

```text
omega=CalA dy^dz+CalB dy^dw+CalC dz^dw.                (16)
```

If `(Y,Z,w)` denotes its pullback, the source density is

```text
CalA(Y_xZ_t-Y_tZ_x)+CalB Y_x+CalC Z_x.                 (17)
```

Both sides of `(12)` are linear in the three coefficient slots.  Only target
Taylor monomials of total degree at most four can enter the `J_4` values,
`J_2` first derivatives, or `J_0` second derivatives.  There are `35` such
monomials in each slot, hence `105` basis two-forms total.  The companion
checks `(12)` separately on all `105` over

```text
Q(u,r_-,r_0,k_-,k_0,k_+),                              (18)
```

after the substitution `(7)`.  It also checks all `63` exact-degree-five
monomial slots are invisible; higher degree is then excluded by the target
maximal-ideal order.  No finite sampling or numeric interpolation is used.

Consequently `(12)` holds for arbitrary polynomial or formal local target
two-forms, and therefore for `dF wedge dG` from every actual target-ring
pair.  This implication does not assert that any such pair satisfies the
global polynomial identities `(13)`.

## 4. The beta line and the two nonzero controls

Let

```text
Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,
L=x(x^2-1),

Q_beta=Q_1+L^2[beta x-(259+16beta)/36].                (19)
```

Every `Q_beta` has `u=1` and satisfies `(7)`.  Its exact fourth debt is

```text
lambda(J_4)=-64(520beta^2+1688beta-5717)/6561.         (20)
```

The zero parameters are exactly

```text
beta=-211/130 +- (99/260)sqrt(94).                     (21)
```

Indeed the quadratic discriminant is `396^2*94`.  The rational points
`beta=0,-1` give the hostile nonzero controls

```text
Q_6=Q_1-(259/36)L^2
   =-(259x^6-36x^5-680x^4+72x^3+502x^2-36x+27)/36,

Q_*=-x^7-(27/4)x^6+3x^5+18x^4-3x^3
                      -(27/2)x^2+x-3/4.                (22)
```

They recover the THM-3642 invoices

```text
Q_6: lambda(J_4)=365888/6561,
Q_*: lambda(J_4)=5440/81.                              (23)
```

Thus zero second debt alone does not imply zero fourth debt.  Conversely,
the irrational roots `(21)` already give explicit algebraic collision
polynomials on the fourth-zero hyperplane.

## 5. Rational zero controls `Q_h` and `Q_9`

The THM-3624 polynomial

```text
Q_h=1/5408(
 44069x^11+112059x^10-154749x^9-406377x^8
+188107x^7+513081x^6-82835x^5-230931x^4
+5408x-4056)                                            (24)
```

has

```text
u=1,       (r_-,r_0,r_+)=(0,0,-243/13),
           (k_-,k_0,k_+)=(0,0,10449/169).              (25)
```

At `u=1,r_-=0`, equation `(15)` is exactly

```text
65k_--169k_++10449=0,                                  (26)
```

the fourth-order THM-3624 equality.  Hence `Q_h` has both zero second and
zero fourth retained debt.

A lower-degree rational control is

```text
Q_9=Q_6+(5717/324)L^3

   =(5717/324)x^9-(5717/108)x^7-(259/36)x^6
     +(5825/108)x^5+(170/9)x^4-(6365/324)x^3
     -(251/18)x^2+x-3/4.                               (27)
```

It has the same values, slopes, and curvatures as `Q_6`, while

```text
(k_-,k_0,k_+)=(35234/27,-6365/54,13094/27).            (28)
```

Direct substitution gives `N_4=0`.  Thus `Q_9` is a rational simultaneous
zero-second/zero-fourth control.

For completeness, every ordinary `u` has a uniform rational-function
control.  Put

```text
Q_5,u=ux^5-2ux^3+ux+(9/2)x^4-(27/4)x^2-3/4,

H_u=-[64u^2x+756ux^2+144ux-432u
                         +1944x^2+243x-972]/(16B),

K_u=[2048u^3x^2-1536u^3+9216u^2x^2+8784u^2x-6912u^2
     +23328ux^2+51516ux-9720u
                         +9477x^2+51759x-4374]/(32B^2),

Q_u=Q_5,u+L^2H_u+L^3K_u.                              (29)
```

Then `deg Q_u<=11` and its nonfixed jets are

```text
(r_-,r_0,r_+)=(0,0,-243/B),
(k_-,k_0,k_+)=(0,0,243(88u-45)/B^2).                  (30)
```

Equations `(7)` and `(15)` follow immediately.  At `u=1`, `(29)` specializes
exactly to `Q_h`.

## 6. Actual local `Q_9` fourth-value control

The hyperplane is nonvacuous for actual target-ring data at the retained
triple.  For `Q_9`, take `G=w` and

```text
F=F_0(y,z)+w^2F_2(y,z)+w^4F_4(y,z),                    (31)

F_0=-1/118272960(
 6725752777y^3-1223385280y^2z-147386304y^2
-23380969yz^2+45143280yz-118272960y-2074176z^2),

F_2=-1/1995856200(
 28674065232y^2+23094463735yz+106774200y
-4336511632z^2-508472640z),

F_4=-(779580313/39917124)y-(428886352/249482025)z.     (32)
```

These are finite actual target-ring polynomials because `y=c/3` and
`z=e+3`.  At each of the three retained branches their pulled-back source
Jacobian satisfies

```text
(J_0,J_0',J_0'')=(1,0,0),
J_1=J_3=0,
(J_2,J_2')=(0,0),
J_4=0.                                                  (33)
```

Before adding `w^4F_4`, the fourth-value vector is

```text
(4049599153,26351689457,34929416497)/997928100,         (34)
```

which is nonzero but has `lambda_1` equal to zero.  The linear `F_4` term
kills the entire vector.  This is an actual local retained-jet control, not
a claim that `J_0=1` or `J_2=0` holds as a polynomial in `x`.

## 7. Strict stopping boundary

The exact positive conclusion is only the retained implication

```text
principal ordinary zero-second-debt plane,
polynomial J_0=1 and J_2=0
       ==> lambda_u(J_4) is given by (14),

lambda_u(J_4)=0 iff the third jets satisfy (15).        (35)
```

Equation `(12)` is universal for arbitrary local target two-forms, but it
only controls one three-branch cokernel value.  It does not construct an
actual global pair satisfying `(13)`, does not show `J_4=0` as a polynomial,
does not address later stable orders or points away from `-1,0,1`, and does
not classify nonordinary `u=+-9/4`.

In particular, neither the zero hyperplane nor the local `Q_9` control gives
a Keller map on `A^2`, and no proof or refutation of the two-dimensional
Jacobian conjecture follows.  `JC(2)` remains **OPEN**.

## 8. Exact reproduction

Run

```bash
python3 04-computation/jc2_russell_cylinder_principal_zero_second_debt_fourth_stable_hyperplane_thm3650.py
python3 -O 04-computation/jc2_russell_cylinder_principal_zero_second_debt_fourth_stable_hyperplane_thm3650.py
```

The frozen normal and optimized transcripts match the stored output
byte-for-byte.  The companion also gates zero Python `assert` AST nodes, so
optimization removes no mathematical check.
