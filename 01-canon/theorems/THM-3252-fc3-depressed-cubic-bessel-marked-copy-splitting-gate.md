---
id: THM-3252
title: "FC(3) affine-coordinate cubics: Bessel-order functional nonsplitting and the ordinary-basis gate"
status: >
  PROVED STRUCTURAL + VERIFIED-EXACT / PERIOD-VALUE EXCLUSION CONDITIONAL.
  After exact cubic depression, the `B=0` branch is THM-3250/3251.  For
  `A*B!=0`, the two
  endpoint primitives form one rank-two regular-singular connection with
  residue exponents `-1/3,-2/3`, and its coordinates satisfy Bessel equations
  of orders `1/3,2/3`.  A collision-safe theorem classifies every rational
  endpoint-exponential splitting of one weighted packet: noncritical grouped
  sources must vanish, while a critical source must lie on
  `span((2B/3,2eta))`.  Triangle moment identities force a nonzero source
  transverse to every such line.  For noncollinear and three-distinct
  collinear images, a hypothetical endpoint-elementary mixed period would satisfy
  both unequal Bessel-order equations and hence vanish packetwise, a
  contradiction.  For a doubled knot, one cyclic row recovers the whole
  nonsplit packet.  These conclusions prove functional nonsplitting, but do
  not by themselves produce a linearly independent derivative-closed
  E-function basis ordinary at `s=1`.  The former Beukers specialization was
  invalid; see MISTAKE-356.  Excluding period values zero and `1/2` in the
  non-pure cubic branch remains conditional on that ordinary-basis gate.
audit: >
  An independent audit rederived both Euclidean divisions, the collision-safe
  rational splitting classification (including its critical top-degree
  obstruction), the turn and divided-difference moment packets, the critical
  fibre directions, the 1/3-versus-2/3 quotient defect, doubled-knot cyclicity,
  and the functional nonsplitting statements.  A later hostile audit found
  that the selected independent scalar family was not shown derivative-closed;
  the full ordinary augmented system is not known independent.  THM-3250/3251
  survive this audit because their exact independent block vectors are already
  derivative-closed and ordinary at one.  Fresh normal and optimized runs
  byte-match the stored structural transcript and both declared hashes.
source: codex-2026-08-02-fc3-cubic-frontier
depends_on:
  - THM-3250-fc3-noncollinear-pure-power-turn-current-exclusion
  - THM-3251-fc3-collinear-pure-power-spline-residue-exclusion
related:
  - THM-3039-the-FC-n-exponential-period-bridge-forced-level-is-the-simplex-volume
  - THM-3203-fc3-complex-affine-coordinate-quadratic-cauchy-green-cycle-exclusion
  - THM-3230-marked-c3-trace-centered-norm-and-terminal-prefactor-recovery
external:
  - "Cauchy--Green formula."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc3_depressed_cubic_bessel_marked_extension_thm3252.py
output: 05-knowledge/results/fc3_depressed_cubic_bessel_marked_extension_thm3252.out
script_sha256: 90936f7ac2e4e096bc07a825252759120b8c4a97308acd6208fcf09107ea96ce
output_sha256: 1d93829c3b5a6678984426521b998d901fdda54ad91dd92dc46a274e04547e1a
---

# THM-3252 -- affine-coordinate cubics: exact nonsplitting and an open value gate

**PROVED STRUCTURAL + VERIFIED-EXACT / PERIOD-VALUE EXCLUSION CONDITIONAL.**

## 1. Result and sharp scope

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},          area(Delta)=1/2,
ell in Qbar[u,v] affine and nonconstant,
Q in Qbar[z],                              deg Q=3.          (T0)
```

The desired period conclusion is

```text
int_Delta exp(Q(ell(u,v)))du dv !=0,
int_Delta exp(Q(ell(u,v)))du dv !=1/2.                      (T1)
```

For the non-pure depressed branch, (T1) is **not yet proved**.  What is proved
below is the exact primitive connection, its collision-safe rational splitting
classification, and functional nonsplitting of the scalar period packets from
the endpoint exponentials.  The missing step is to realize the relevant
independent functions inside a derivative-closed E-function basis whose system
is ordinary at `s=1`.  Conditional on that gate, (T1) follows.  Pure cubics
remain covered unconditionally by THM-3250/3251.  Nothing here covers a
genuinely bivariate cubic, nonalgebraic coefficients, the projective-leading-
form reduction, or `FC(3)` itself.

The coupled mechanism is stated first.  Let

```text
p(t)=At^3+Bt,                     A,B in Qbar, A*B!=0,       (1)
E_x(s)=exp(sp(x)),
J_(0,x)(s)=integral_0^x exp(sp(t))dt,
J_(1,x)(s)=integral_0^x t exp(sp(t))dt.                     (2)
```

Every endpoint `x` is algebraic.  Put

```text
Y_x=(J_(0,x),J_(1,x))^T,
kappa=2B/(3A),

Cmat=[ 0             2B/3       ],
     [ -2B^2/(9A)    0          ]

D=diag(-1/3,-2/3).                                        (3)
```

Then the primitive packet has the exact connection

```text
Y_x'=(Cmat+D/s)Y_x
     +1/(3s){(x,x^2+kappa)^T E_x-(0,kappa)^T}.              (4)
```

Consequently the two scalar homogeneous equations are

```text
s^2 f''+s f'+(lambda^2 s^2-1/9)f=0,
s^2 g''+s g'+(lambda^2 s^2-4/9)g=0,
lambda^2=4B^3/(27A).                                       (5)
```

Thus the general depressed cubic does not create many unrelated
transcendental primitives: it couples the two nontrivial cubic residue
classes into one rank-two Bessel connection.

The structural core below is an exact **one-copy splitting classification**.
It identifies every way a weighted packet `Y_w=sum_x w_xY_x` can become a
rational combination of its endpoint exponentials.  Sections 8--10 then use
the geometry retained by the triangle weights to rule out the remaining
mixed-coordinate cancellation and prove functional nonsplitting.  Section 10
states the additional specialization hypothesis that would imply (T1).

## 2. Exact depression of an arbitrary affine-coordinate cubic

Start with

```text
Q(z)=a z^3+b z^2+c z+d,             a!=0.                  (6)
```

Set

```text
t=z+b/(3a),
A=a,
B=c-b^2/(3a),
C_0=d-bc/(3a)+2b^3/(27a^2).                                (7)
```

Direct expansion gives

```text
Q(z)=At^3+Bt+C_0.                                           (8)
```

For a triangle phase `Q(ell)`, translation only changes the three endpoint
marks from `ell(v_j)` to `ell(v_j)+b/(3a)`, while the constant contributes
the harmless external factor `exp(C_0s)`.  There are three sharp boundaries:

* `A=0` is not cubic and returns to the quadratic theory;
* `A!=0,B=0` is the translated pure-cubic branch already closed by
  THM-3250 and THM-3251; and
* the present coupled calculation assumes `A*B!=0`.

No division by `B` below is silently used at the pure-power boundary.

## 3. Derivation of the rank-two primitive system

The two Euclidean divisions carrying all of the calculation are

```text
p=(t/3)p'+(2B/3)t,                                          (9)

tp=(t^2/3+2B/(9A))p'-(2B^2/(9A)).                         (10)
```

Integrating the `p' exp(sp)` terms by parts on the path from `0` to `x`
gives

```text
3sJ_(0,x)'+J_(0,x)-2BsJ_(1,x)=xE_x,                       (11)

3sJ_(1,x)'+2J_(1,x)+(2B^2/(3A))sJ_(0,x)
  =(x^2+kappa)E_x-kappa.                                  (12)
```

Equations (11)--(12) are precisely (4).  The constant `-kappa` in (12) is
the lower-boundary term from the constant part of the quotient in (10).
It is load-bearing when `p(x)=0`: the endpoint exponential then belongs to
the same class as this constant source and must be grouped with it.

## 4. Scalar Bessel equations and the pure-cubic degeneration

Discard the inhomogeneous endpoint source temporarily and write
`Y=(f,g)^T`.  The homogeneous part of (4) is

```text
f'=(2B/3)g-f/(3s),
g'=-(2B^2/(9A))f-2g/(3s).                                  (13)
```

Because `B!=0`, either coordinate can be eliminated.  Eliminating `g` gives
the first equation of (5), and eliminating `f` gives the second.  The local
exponents at `s=0` are therefore the two nonintegral cubic classes

```text
-1/3,                         -2/3.                         (14)
```

At `B=0`, both off-diagonal entries of `Cmat` vanish before any elimination,
and (4) splits into

```text
sJ_(0,x)'=-(1/3)J_(0,x)+(x/3)exp(Ax^3s),
sJ_(1,x)'=-(2/3)J_(1,x)+(x^2/3)exp(Ax^3s).                 (15)
```

These are exactly the residue `1/3` and `2/3` blocks used in
THM-3250/3251.  Thus pure cubics are not a disconnected special trick:
they are the diagonal degeneration of the general cubic connection.

## 5. Collision-safe source data

Let `X` be a finite multiset of algebraic endpoints and let `w_x in Qbar`
be weights, with repeated endpoints combined if desired.  Define

```text
Y_w=sum_(x in X) w_xY_x.                                   (16)
```

Group endpoints by their phase value `eta=p(x)`.  The source vector on the
distinct exponential `exp(eta s)` is

```text
b_eta=sum_(x:p(x)=eta) w_x(x,x^2+kappa)^T
      -1_(eta=0)(sum_x w_x)(0,kappa)^T.                    (17)
```

The second term is global because every endpoint copy in (4) contributes
the same lower-boundary constant.  Formula (17), rather than the individual
endpoint vectors, is the invariant input when two or all three cubic
endpoint values collide.

Let `mathcal E_X` denote the `Qbar(s)`-span of the distinct exponentials

```text
{exp(p(x)s):x in X} union {1}.                              (18)
```

The eigenvalues of `Cmat` are exactly the two critical values of `p`.
Indeed, if `p'(r)=3Ar^2+B=0`, then

```text
p(r)=2Br/3,
p(r)^2=-4B^3/(27A),                                        (19)
```

and `det(Cmat-eta I)=eta^2+4B^3/(27A)`, so the same quadratic
is the characteristic equation of `Cmat`.

## 6. Exact one-copy splitting theorem

The following are equivalent:

```text
(i)  both coordinates of Y_w belong to mathcal E_X;

(ii) for every phase class eta,
       b_eta=0,                              eta not critical,
       b_eta in L_eta,                       eta critical,    (20)

     where
       L_eta=span_Qbar((2B/3,2eta)^T).
```

This criterion is exact, including `eta=0`, endpoint collisions, and a
phase class containing three roots.

### 6.1 Necessity

Distinct exponentials are linearly independent over `Qbar(s)`.  If (i)
holds, write

```text
Y_w=sum_eta R_eta(s)exp(eta s),             R_eta in Qbar(s)^2. (21)
```

Coefficient comparison in (4) gives, with `M_eta=Cmat-eta I`,

```text
R_eta'=M_eta R_eta+(D/s)R_eta+b_eta/(3s).                  (22)
```

First, every rational solution of (22) is polynomial.  At a finite nonzero
pole, differentiation creates an unmatched pole of one higher order.  At a
pole of order `m>=1` at zero, the leading coefficient would lie in the
kernel of

```text
-mI-D=diag(-m+1/3,-m+2/3),                                 (23)
```

which is invertible.  Hence no finite pole exists.

Write a polynomial solution as `R=sum_(k=0)^N r_ks^k`.  Multiplying (22)
by `s` and comparing coefficients gives

```text
D r_0+b_eta/3=0,
M_eta r_(k-1)+(D-kI)r_k=0,                 1<=k<=N,
M_eta r_N=0.                                                (24)
```

If `eta` is not critical, `M_eta` is invertible.  The last equation and
downward induction force every `r_k=0`, and then `b_eta=0`.

Suppose `eta` is critical.  Then

```text
ker M_eta=span((2B/3,eta)^T),
im M_eta =span((2B/3,-eta)^T).                              (25)
```

If `N>=1`, the top two equations in (24) would require

```text
(NI-D)(2B/3,eta)^T in im M_eta.                             (26)
```

The two coordinate ratios in (26) say

```text
N+1/3=-(N+2/3),                                            (27)
```

or `2N+1=0`, impossible for an integer `N>=1`.  Thus `R_eta`
is constant, lies in `ker M_eta`, and the first equation of (24) gives

```text
b_eta=-3D r_0 in span((2B/3,2eta)^T)=L_eta.                (28)
```

No polynomial of positive degree provides an unrecorded exceptional split.

### 6.2 Sufficiency

For a noncritical class with `b_eta=0`, take `R_eta=0`.  For a critical
class with `b_eta` in `L_eta`, equation (28) uniquely supplies a constant
`R_eta in ker M_eta`.  Subtract the resulting rational exponential
particular solution from `Y_w`.  The difference is an entire solution `H`
of the homogeneous connection.

There is no nonzero entire homogeneous solution.  If
`H=sum_(n>=N)h_ns^n`, `h_N!=0`, has its first nonzero Taylor coefficient,
the lowest term of

```text
sH'=sCmat H+DH                                                (29)
```

would give `(NI-D)h_N=0`; its diagonal entries are `N+1/3` and
`N+2/3`.  This contradiction proves `H=0`, hence (i).

## 7. Hostile boundary tests and what the theorem excludes

The exceptional critical line is real structure, not an artifact of the
proof.  It is nevertheless unavailable to one isolated endpoint.  If `r`
is a critical point, its natural endpoint source is

```text
v_r=(r,r^2+kappa)^T=(r,B/(3A))^T.                          (30)
```

Using `eta=p(r)=2Br/3`, one obtains

```text
det[v_r,(2B/3,2eta)^T]=-2B^2/(3A)!=0.                      (31)
```

Thus even a critical endpoint is transverse to its only splitting line.
A split can arise only after sources with a common phase exponent have been
combined, or after the zero-phase class has been combined with the constant
lower-boundary source.

The exact verifier includes three deliberately hostile controls for
`p=t^3+3t`:

* all three roots of `p(t)=4` are grouped before testing their source;
* the critical endpoints `t=+/-i` at values `+/-2i` are each transverse to
  their exceptional line; and
* the nonzero root `i sqrt(3)` of `p(t)=0` is grouped with the constant
  source in (17).

It also solves the rational-polynomial coefficient system through degree
five at both critical values and a noncritical value.  This finite audit is
a regression control for the general pole-and-leading-degree proof above,
not a replacement for it.

The theorem therefore excludes exactly this failure mode:

```text
one weighted rank-two primitive packet becomes endpoint-elementary
without its grouped source satisfying the critical-line criterion.       (32)
```

It does **not** say that either coordinate of a packet is separately
non-elementary, nor that a scalar expression using coordinates from two
different packets cannot cancel.

## 8. Three distinct endpoints: the mixed copies cannot cancel

This is where the apparent marked-copy obstruction closes.  Both relevant
three-distinct geometries have one common moment packet.

### 8.1 Geometry and moments

In the noncollinear Cauchy--Green geometry of THM-3250, let `omega_j=c_j`
be the vertex-turn weights, let `x_j` be the translated endpoint values,
and put `Omega=W`, the nonzero oriented double-area factor.  In the
three-distinct collinear geometry of THM-3251, let `omega_j=alpha_j` be the
B-spline divided-difference weights and put `Omega=1`.  In both cases,

```text
sum_j omega_j=0,
sum_j omega_jx_j=0,
sum_j omega_jx_j^2=-Omega!=0.                              (33)
```

For the turn weights, the first two identities telescope around the edge
cycle and the third is the degree-zero Cauchy--Green identity.  For the
B-spline weights they are the first three divided-difference moments.

Define two differently marked sections of the same primitive connection:

```text
U=sum_j omega_jY_(x_j),
V=sum_j omega_jx_jY_(x_j),
F=U_1-V_0.                                                  (34)
```

The area-period formula is uniformly

```text
Omega K(s)=exp(C_0s)F(s).                                  (35)
```

### 8.2 The unmarked packet `U` cannot split

Because `sum omega_j=0`, the constant correction in (17) cancels for `U`.
Its grouped source at phase value `eta` is

```text
b_eta=sum_(j:p(x_j)=eta) omega_j(x_j,x_j^2+kappa)^T.       (36)
```

The sum over all phase classes is `(0,-Omega)^T` by (33), so at least one
`b_eta` is nonzero.  A nonzero source at a noncritical value violates the
first alternative of (20).

It remains to audit a critical collision, not assume generic endpoints.
If `r` is critical and `eta=p(r)`, then

```text
p(x)-eta=A(x-r)^2(x+2r),
B=-3Ar^2,                  kappa=-2r^2.                    (37)
```

Thus every endpoint in this phase class is `r` or `-2r`, and its source is
respectively

```text
(r,r^2+kappa)=r(1,-r),
(-2r,4r^2+kappa)=-2r(1,-r).                               (38)
```

Every nonzero grouped critical source therefore lies on `span((1,-r))`.
But the exceptional splitting line in (20) is

```text
L_eta=span((1,2r)),                                        (39)
```

and these two lines are transverse because `r!=0`.  Hence the nonzero
source guaranteed by (33) fails the splitting criterion in every possible
phase class.  Therefore

```text
U notin mathcal E_X^2.                                     (40)
```

This proof includes a double critical root, the third root with the same
critical value, and all three-endpoint phase collisions.

### 8.3 The Bessel-order mismatch kills cross-copy cancellation

Suppose for contradiction that `F` belonged to `mathcal E_X`.  Work modulo
that derivative-stable space and write lower-case letters for classes.
Both `u` and `v` obey the homogeneous connection (13), while (34) gives

```text
u_1=v_0=:f.                                                 (41)
```

By (5), the left side of (41), as a second coordinate, satisfies the
order-`2/3` scalar equation, while the right side, as a first coordinate,
satisfies the order-`1/3` equation.  Subtracting those two equations gives

```text
-(1/3)f=0.                                                  (42)
```

Thus `u_1=0`.  The second row of (13), and `A*B!=0`, then forces `u_0=0`.
This says that the whole packet `U` splits, contradicting (40).  Hence

```text
F notin mathcal E_X.                                       (43)
```

No simplicity theorem or abstract Schur argument is needed.  The two marked
copies cannot glue because the scalar coordinate they would have to share
has two incompatible cubic local exponents.

## 9. One doubled knot: a cyclic row recovers the packet

Suppose two vertex values equal `a` and the remaining value is `c!=a`.
Put `L=c-a` and

```text
T=Y_a-Y_c,
G=(-c,1)T=-cT_0+T_1.                                      (44)
```

The triangular marginal formula of THM-3251 is

```text
L^2K(s)=exp(C_0s)G(s).                                     (45)
```

The packet `T` has weights `1,-1`, so its constant lower-boundary source
cancels.  It cannot split:

* if `p(a)!=p(c)`, each phase group is a singleton with nonzero source;
  a noncritical singleton violates (20), while a critical singleton is
  transverse by (31);
* if `p(a)=p(c)`, its single grouped source is
  `(a,a^2+kappa)-(c,c^2+kappa)`, which is nonzero because its first
  coordinate is `a-c`; if the common value is critical, (37)--(39) again
  put this source on `span((1,-r))`, transverse to the splitting line.

Thus

```text
T notin mathcal E_{a,c}^2.                                 (46)
```

Now suppose `G` were endpoint-elementary.  Modulo the exponential space,
`g=0` gives `t_1=ct_0`.  Put

```text
Gamma=2B^2/(3A)!=0.                                        (47)
```

The two rows of (11)--(12), modulo their sources, are

```text
3st_0'=-t_0+2Bst_1,
3st_1'=-Gamma st_0-2t_1.                                  (48)
```

Differentiate `t_1=ct_0` and substitute (48).  The result is

```text
[c+s(2Bc^2+Gamma)]t_0=0.                                  (49)
```

The bracket is never the zero polynomial: if `c!=0` it has nonzero
constant term, and if `c=0` its coefficient is `Gamma!=0`.  Hence
`t_0=t_1=0`, contradicting (46).  Therefore

```text
G notin mathcal E_{a,c}.                                   (50)
```

The doubled-knot scalar row is cyclic for the rank-two connection.  This is
the multiplicity boundary that the three-distinct mixed-coordinate proof
cannot cover directly.

## 10. E-functions and the remaining ordinary-basis gate

The splitting theorem remains unchanged if `mathcal E_X` is enlarged by
finitely many algebraic exponentials having zero source.  In the proof of
section 6, every such coefficient satisfies (22) with `b_eta=0`; the same
pole and leading-degree audit forces it to vanish, including at a critical
value because the constant equation has invertible `D`.  We may therefore
adjoin `exp(-C_0s)` to the exponential spaces in (43) and (50).

For `k=0,1`, the expansion

```text
J_(k,x)(s)=sum_(n>=0) s^n/n!
             integral_0^x t^k(At^3+Bt)^n dt               (51)
```

has algebraic coefficients with the required conjugate and denominator
growth, and (11)--(12) are holonomic.  Hence these primitives, their
weighted combinations, and the algebraic exponentials are E-functions.
The **full augmented primitive vector** satisfies a first-order system with
denominator only `s`, so that system is ordinary at `s=1`.  Equations (43)
and (50), after adjoining `exp(-C_0s)`, prove that the selected scalar `F` or
`G` is functionally independent from the distinct endpoint exponentials.

These two facts cannot be spliced without another argument.  The full
augmented vector is not known to be linearly independent, while the selected
independent list containing `F` or `G` is not shown closed under
differentiation.  Replacing a dependent ordinary system by an independent
basis can introduce an apparent singularity at `s=1`.  The elementary pair
`(s-1,exp(s))` is the minimal warning: it is functionally independent, but its
first value vanishes and the minimal scalar equation for `s-1` is singular at
one.  Thus Beukers Corollary 1.4 does not presently apply to the selected
family.

The precise remaining gate is any proof of one of the following equivalent-
purpose certificates:

1. a linearly independent, derivative-closed E-function vector ordinary at
   `s=1` that contains `F` (respectively `G`) and the needed distinct
   exponentials;
2. a basis of the generated differential module, regular at `s=1`, in which
   those functions have coefficients regular at one and nonzero specialized
   coordinate rows; or
3. a direct specialization theorem controlling the full relation module at
   `s=1`.

Conditional on such a certificate, Beukers transfers the proved functional
independence to `Qbar`-linear independence of the relevant values.  Indeed, in a
three-distinct geometry, (35) shows that either forbidden value would give
respectively

```text
F(1)=0,
F(1)-(Omega/2)exp(-C_0)=0.                                 (52)
```

In the doubled-knot geometry, (45) gives respectively

```text
G(1)=0,
G(1)-(L^2/2)exp(-C_0)=0.                                   (53)
```

Each would contradict the conditional ordinary-basis certificate.  Hence
(T1) follows **conditionally** for `A*B!=0`.

If the depressed coefficient `B` is zero, THM-3250 and THM-3251 prove the
same conclusion using independent derivative-closed split residue blocks
ordinary at one; their Beukers step is unaffected by MISTAKE-356.
Noncollinear, three-distinct collinear, and doubled-knot images exhaust every
nonconstant affine coordinate on a triangle.  Thus the exact structural
reduction is complete for every cubic `Q` in (T0), while the non-pure cubic
period-value conclusion remains at the stated gate.

## 11. Precise C3 dictionary with the Jacobian sidecar

There is a typed structural parallel with THM-3230, not a theorem transfer.
The two local Mellin exponents in (14) are the two nontrivial characters of
the cubic monodromy packet.  THM-3230 writes a trace-centered cubic Kummer
element as

```text
delta=pi A_1+pi^2 B_1,                                     (54)
```

whose two terms have character degrees `1` and `2` modulo three.  The map is

```text
source:       local exponent classes {-1/3,-2/3} in (4);
target:       tame C3 character degrees {1,2} in (54);
map:          local monodromy exponent modulo Z;
preserves:    the two-component nontrivial C3 packet;
destroys:     relative amplitude, cyclic ordering, and geometric marking;
sidecar here: Bessel order of the marked scalar coordinate.               (55)
```

In THM-3230 a norm keeps a cubeclass but loses the relative character
amplitude, so a marked sheet is required.  Here the scalar period initially
loses the relative amplitude between `U` and `V`, but their coordinate marks
carry different Bessel orders; the defect (42) recovers enough information
to forbid gluing.  This connection explains the recurring two-channel
`C3` packet.  It does not import the Jacobian theorem or prove `FC(3)`.

## 12. Reproduction contract

Run

```bash
python3 04-computation/fc3_depressed_cubic_bessel_marked_extension_thm3252.py
python3 -O 04-computation/fc3_depressed_cubic_bessel_marked_extension_thm3252.py
```

The verifier checks, over exact SymPy arithmetic:

1. the affine cubic depression (7)--(8);
2. both Euclidean divisions and 60 Taylor-coefficient instances of
   (11)--(12) over three complex/algebraic controls;
3. both scalar eliminations in (5);
4. the `B=0` residue degeneration;
5. the critical and noncritical splitting-source ranks through rational
   polynomial degree five for `p=t^3+3t`;
6. all three collision/boundary controls listed in section 7;
7. the uniform critical-group direction, splitting-line transversality, and
   cross-coordinate Bessel-order defect; and
8. exact turn/B-spline moment controls and the doubled-knot cyclic row; and
9. 18 direct Taylor-coefficient comparisons for the noncollinear,
   three-distinct collinear, and doubled-knot period formulas under
   `p=t^3+3t`.

Normal and optimized runs must be byte-identical to the archived output.  The
verifier certifies the structural identities and nonsplitting mechanisms; it
does not certify specialization at `s=1`.
