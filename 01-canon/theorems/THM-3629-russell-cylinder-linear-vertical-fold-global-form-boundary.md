---
id: THM-3629
title: "Russell-cylinder linear vertical-fold global-form boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For an even
  retained-collision polynomial Q
  and H in t C[t] with H'(0) nonzero, every hypothetical global regular
  target pair with nonzero constant source Jacobian is already a
  noninjective polynomial Keller map.  If H is nonlinear, neither output
  can be surface-only or stable-only: both must genuinely mix the
  Danielewski and stable coordinates.  If H=ht is affine, a global exact
  polynomial target two-form pulls back to a unit, but decomposing that
  exact form as dF wedge dG is precisely the still-open surface Darboux
  problem.  The complete affine-linear cell in the global inverse-cylinder
  coordinates is empty.  Arbitrary higher mixed pairs remain OPEN; no
  JC(2) counterexample is claimed.
source: root / audit_thm3623 linear-boundary continuation, 2026-08-21
audit: >
  PASS -- an independent hostile derivation checked the shifted D=0
  collision, nonlinear critical-value unit gate, global exact form and Russell
  transport, decomposition determinant, complete affine-linear q-degree
  staircase, and local/nonclosed boundaries; normal, optimized, and stored
  73-gate transcripts are byte-identical.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3623-russell-cylinder-even-general-vertical-fold-all-order-closure
related:
  - THM-3618-compiler-one-graph-observable-fibre-separator-no-embedding
  - THM-3622-compiler-one-observable-graph-closure-normalization-arm-debt
  - THM-4054-exceptional-affine-simple-zero-retained-packet
script: 04-computation/jc2_russell_cylinder_linear_vertical_global_form_boundary_thm3629.py
output: 05-knowledge/results/jc2_russell_cylinder_linear_vertical_global_form_boundary_thm3629.out
script_sha256: fa7c8c0edaf2ee31140f810923001126bf935ee1b3da3a21dc8cfd2203e54613
output_sha256: 182dcf40daa4eec4c8a02118e0dbedeb72b517fba4dc5e926e80336c03a77e6e
hash_basis: raw LF bytes
---

# THM-3629 -- Russell-cylinder linear vertical-fold global-form boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
identifies the exact globality
debt left by the `H'(0)!=0` local escape in THM-3623.  It proves several
unrestricted no-entry lemmas, but it deliberately does **not** close every
global target pair in this boundary.

All rings, differentials, and closed points are over `C`.

## 0. Setup and statement

Use the exponent-two Danielewski surface and its cylinder

```text
Y_2=Spec R,
R=C[b,c,e]/(c^2e-b(b+4)),                 A=R[w].       (1)
```

The THM-3561 compiler is

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3).                                             (2)
```

Let `Theta:Y_2 x A1_w -> Y_1 x A1_S` be the polynomial Russell
isomorphism of THM-3605 and put

```text
Psi=Theta o (Phi x id).                                (3)
```

Assume

```text
Q in C[x] even,       Q(0)=-3/4,       Q(1)=Q(-1)=-3,
H in t C[t],          h=H'(0)!=0,                         (4)
```

and define

```text
E_(Q,H)(x,t)=(x,Q(x)+H(t),t).                          (5)
```

For regular functions `F,G` on the exponent-one Russell target, transport
them through `Theta` to `F^#,G^# in A`.  Write

```text
f=F o Psi o E_(Q,H),          g=G o Psi o E_(Q,H).     (6)
```

The boundary theorem has five parts.

1. If `Jac_(x,t)(f,g)` is a nonzero constant, then `(f,g)` is already a
   noninjective polynomial Keller map and hence a counterexample to `JC(2)`.
2. If `deg(H)>=2`, neither `F^#` nor `G^#` lies in `R`, and neither lies in
   `C[w]`.  Both outputs must genuinely mix surface and stable coordinates.
3. If `H=ht`, there is an explicit global exact polynomial two-form on the
   Russell target whose source pullback is a nonzero constant.
4. Decomposing that particular form as `dF wedge dG`, even with arbitrary
   stable dependence allowed, is equivalent to the still-open polynomial
   Darboux system on `Y_2`.
5. For `H=ht`, no pair affine-linear in the global inverse-cylinder
   coordinates `(b,c,e,w)` has nonzero constant source Jacobian.  This last
   statement does not require evenness or the collision values of `Q`.

Arbitrary higher-degree genuinely mixed pairs remain **OPEN**.

## 1. Every survivor would already refute `JC(2)`

At `t=0`, equations `(2),(4)` give

```text
(x,q)=(-1,-3),(0,-3/4),(1,-3)
          |-> (b,c,e,w)=(0,0,-3,0).                    (7)
```

Therefore the same three distinct points have one common image under
`Psi o E_(Q,H)`.  If the two polynomials `(f,g)` in `(6)` had nonzero
constant Jacobian, they would form a polynomial Keller map identifying the
three points in `(7)`.  It would be a genuine counterexample to `JC(2)`.

This implication is not used as a proof of nonexistence.  It specifies the
actual consequence of any positive construction.

## 2. A shifted `D=0` sign collision on every even graph

Write

```text
Q(x)=Qbar(x^2).                                         (8)
```

The two values in `(4)` make `Qbar` nonconstant.  Fix any `s in C` and put

```text
p_s(y)=1+y(Qbar(y)+s).                                  (9)
```

This is a nonconstant polynomial with constant term one.  It therefore has
a root `y_s!=0`.  Choose `r_s` with `r_s^2=y_s`.  At

```text
x=+-r_s,          q=Qbar(y_s)+s=-1/y_s,                (10)
```

one has `D=0`, and direct substitution into `(2)` gives

```text
(b,c,e)=(-4,0,-3/y_s)                                  (11)
```

at both signs.  Thus every shifted even graph `q=Q(x)+s` has a pair of
distinct compiler points with the same target.  Fixing the same stable value
at the two points preserves the collision under the Russell isomorphism.

The nonzero root in this argument is load-bearing: it makes `r_s` and
`-r_s` distinct.

## 3. Nonlinear vertical terms force both outputs to mix

Assume now that `deg(H)>=2`.  The polynomial `H'` is nonconstant, so choose

```text
tau in C with H'(tau)=0,             s=H(tau).          (12)
```

Suppose first that `F^# in R`, while `G^# in A` is arbitrary.  Along the
source line `t=tau`, the vertical derivative of `f` vanishes and the
vertical derivative of `g` is its explicit `w` derivative:

```text
f_t(x,tau)=0,
g_t(x,tau)=(partial_w G^#) o (Phi(x,Q(x)+s),tau).       (13)
```

If the source Jacobian were `kappa in C*`, equations `(13)` would give an
identity in `C[x]`:

```text
f_x(x,tau) g_t(x,tau)=kappa.                           (14)
```

The only units of `C[x]` are nonzero constants.  Hence

```text
f_x(x,tau)=alpha in C*,
f(x,tau)=alpha x+beta.                                 (15)
```

But Section 2 supplies distinct points `+-r_s` with the same surface target,
so the left side of `(15)` has equal values there.  The right side differs by
`2alpha r_s!=0`, a contradiction.

If `G^# in R`, the symmetric identity at `t=tau` is

```text
-f_t(x,tau)g_x(x,tau)=kappa,                           (16)
```

which forces `g_x` to be a nonzero constant and gives the same contradiction.
Therefore

```text
deg(H)>=2 and Jac(f,g) in C*
        imply F^#,G^# notin R.                         (17)
```

There is a separate stable-only obstruction valid for every `H` in `(4)`.
If, for example, `F^#=p(w)`, then

```text
-p'(t)g_x(x,t)=kappa.                                  (18)
```

Both factors in `(18)` are polynomial units, so `p' in C*` and

```text
g(x,t)=gamma x+r(t),                 gamma in C*.       (19)
```

Equation `(19)` cannot take one common value at the three points `(7)`.
The same argument applies with the two outputs reversed.  Combining this
with `(17)` proves the precise nonlinear boundary:

```text
every possible nonlinear-H survivor must genuinely mix R and w
in each of its two outputs.                             (20)
```

## 4. The affine boundary passes the global two-form gate

On `Y_2`, put

```text
eta=db wedge dc/c^2.                                   (21)
```

Although `(21)` is written on `c!=0`, it has a global polynomial
representative.  Define

```text
alpha=-(ce/4)db+(e(b+2)/4)dc+(c(b+2)/8)de.             (22)
```

Differentiation in the Kähler module of `R` gives

```text
eta=d alpha
   =(e/2) db wedge dc
    +(3c/8) db wedge de
    -(b+2)/8 dc wedge de.                              (23)
```

Thus `eta` is global, polynomial, and exact.  The exact compiler minors are

```text
Jac_(x,q)(b,c)=-3c^2,
Jac_(x,q)(b,e)= 6ce,
Jac_(x,q)(c,e)= 6(b+2).                                (24)
```

Substitution of `(24)` into `(23)` and use of `c^2e=b(b+4)` gives

```text
Phi^*eta=-3 dx wedge dq.                               (25)
```

The form can be written completely polynomially on the exponent-one Russell
target.  Put

```text
g(B)=-B^2(B+6)/8,             K_0(B)=(B-2)(B+6)/64,

b_*=g(B)+C^2S,
e_*=Y^2K_0(B)+2(g(B)+2)S+C^2S^2.                       (26)
```

These are the inverse-cylinder coordinates of THM-3605, with `c_*=C`.
Define

```text
eta_tilde
 =(e_*/2) db_* wedge dC
  +(3C/8) db_* wedge de_*
  -(b_*+2)/8 dC wedge de_*.                            (27)
```

Every coefficient in `(27)` is polynomial.  By construction,

```text
eta_tilde=(Theta^(-1))^*eta,
Psi^*eta_tilde=-3 dx wedge dq.                         (28)
```

If `H=ht`, equations `(5),(28)` give

```text
E_(Q,ht)^*Psi^*eta_tilde=-3h dx wedge dt.              (29)
```

In particular, the global target form

```text
                     -(4/h)eta_tilde                   (30)
```

pulls back to the normalized constant `12 dx wedge dt`.

Thus the affine boundary is not obstructed by regularity, exactness, or
existence of a global unit two-form.  Its remaining debt is decomposition by
two global regular functions.

## 5. Decomposing the global form is exactly the surface Darboux problem

Suppose that regular target functions satisfy the identity

```text
dF wedge dG=lambda eta_tilde,             lambda in C*. (31)
```

Transport through `Theta` and use the relative differential `d_Y` on
`A=R[w]`:

```text
dF^#=d_YF^#+F^#_w dw,
dG^#=d_YG^#+G^#_w dw.                                  (32)
```

Their wedge is

```text
d_YF^# wedge d_YG^#
 +(G^#_w d_YF^#-F^#_w d_YG^#) wedge dw.                (33)
```

The surface part of `(31)` says

```text
d_YF^# wedge d_YG^#=lambda eta!=0.                     (34)
```

Hence `d_YF^#,d_YG^#` are linearly independent over `Frac(R[w])`.  The
vanishing `dw` coefficient in `(33)` then forces

```text
F^#_w=G^#_w=0.                                         (35)
```

Characteristic zero makes `(35)` equivalent to `F^#,G^# in R`.  Conversely,
surface functions with

```text
{F^#,G^#}=lambda                                       (36)
```

give `(31)`.  Therefore exact global decomposition of `(27)` is equivalent
to the residual polynomial Darboux system `(36)` of THM-3561.  That system
is **OPEN** in arbitrary support.

For `H=ht`, restrict explicitly to surface-only pairs `F^#,G^# in R`.  The
same calculation gives the exact equivalence

```text
Jac_(x,t)(f,g)=kappa in C*
 iff {F^#,G^#}=-kappa/(3h).                             (37)
```

The implication uses that `(x,t)->(x,q=Q(x)+ht)` is a polynomial
automorphism and that the compiler image is dense in `Y_2`.

An arbitrary mixed target pair is not required to satisfy `(31)`: its
two-form may differ from a scalar multiple of `eta_tilde` by a form killed
on the chosen source.  This is the exact remaining conormal exit.

## 6. Complete affine-linear global cell

Remain in the affine case `H=ht`, but now allow `Q` to be any polynomial.
Suppose `F^#,G^#` are affine-linear in `(b,c,e,w)`.  Their wedge has the
form

```text
A db wedge dc+B db wedge de+C dc wedge de
+D_0 db wedge dw+E_0 dc wedge dw+L de wedge dw,         (38)
```

with constant coefficients.  Ignoring the decomposability relation among
these six constants only enlarges the tested universe.

Use the polynomial coordinate change `q=Q(x)+ht` and put

```text
delta_Q=partial_x+Q'(x)partial_q.                       (39)
```

The six source coefficients in the order displayed in `(38)` have the
following leading `q` terms:

| target wedge | `deg_q` | leading coefficient |
|---|---:|---|
| `db wedge dc` | 4 | `-3h x^10` |
| `db wedge de` | 4 | `6h x^7` |
| `dc wedge de` | 3 | `6h x^6` |
| `db wedge dw` | 3 | `6x^5` |
| `dc wedge dw` | 2 | `5x^4` |
| `de wedge dw` | 2 | `2x` |

Indeed, the first three rows are `(24)` multiplied by `h`, while the last
three are `delta_Q(b),delta_Q(c),delta_Q(e)`.

If `(38)` pulled back to a constant, its `q^4` coefficient would give

```text
-3h A x^10+6h B x^7=0,                                 (40)
```

so `A=B=0`.  The `q^3` coefficient then gives

```text
6h C x^6+6D_0 x^5=0,                                   (41)
```

so `C=D_0=0`.  Finally, the `q^2` coefficient gives

```text
5E_0 x^4+2Lx=0,                                        (42)
```

so `E_0=L=0`.  Thus the only constant source coefficient in this enlarged
affine-linear cell is zero.

This no-go is all-degree in `Q` and makes no parity assumption.  It is not an
extrapolation from a bounded target search.

## 7. Sharp local controls and OPEN boundary

On the local chart `b+4!=0`, put `a=e/(b+4)`.  Since

```text
eta=da wedge dc,                                        (43)
```

the affine fold has the exact local Darboux pair

```text
F=a,                  G=-4c/h,                          (44)
```

with source Jacobian `12`.  Equations `(23),(30)` show that the two-form of
this pair extends globally even though `a` does not.  Polynomial
decomposition, not two-form regularity, is the missing global coordinate.

For nonlinear `H` with `H'(0)!=0`, the enlarged-form survivor from THM-3623
is

```text
-4/H'(w) eta.                                           (45)
```

It is not closed unless `H` is affine:

```text
d(-4eta/H')=4H''/(H')^2 dw wedge eta.                  (46)
```

The tempting local functions `F=a,G=-4c/H'(w)` do not silently repair
closedness.  Their source coefficient is

```text
12+4 delta_Q(a) c H''/(H')^2,                          (47)
```

not the constant produced by the enlarged form `(45)`.

What is proved by this theorem:

1. every positive global pair would already be a `JC(2)` counterexample;
2. every nonlinear survivor must genuinely mix `R` and `w` in both outputs;
3. the affine boundary has a global exact unit two-form;
4. exact decomposition of that form is the residual surface Darboux problem;
5. the complete affine-linear inverse-cylinder cell is empty.

What remains **OPEN**:

- affine `H=ht`: arbitrary higher-degree mixed pairs and the surface Darboux
  system `(36)`;
- nonlinear `H` with `H'(0)!=0`: pairs in which both outputs genuinely mix
  surface and stable coordinates; and
- target two-forms differing from `(27)` by a form killed on the chosen
  source.

No global target pair, Keller example, or `JC(2)` counterexample is
constructed.

## 8. Exact verification

The companion verifies the compiler and retained triple, the exact `D=0`
sign collision, the polynomial representative `(23)`, all compiler minors
and the sign in `(25)`, the Russell inverse transport, affine normalization,
the decomposition determinant, nonlinear and stable-only hostile controls,
the complete affine-linear `q`-degree staircase, and the local/nonclosed
boundaries.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_linear_vertical_global_form_boundary_thm3629.py
python3 -O 04-computation/jc2_russell_cylinder_linear_vertical_global_form_boundary_thm3629.py
```

Both transcripts must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_linear_vertical_global_form_boundary_thm3629.out
```

The frozen companion reports

```text
CHECKS=73
RESULT=PASS
```

with zero Python assertion statements.  They are finite controls for the
written proof, not a replacement for it.
