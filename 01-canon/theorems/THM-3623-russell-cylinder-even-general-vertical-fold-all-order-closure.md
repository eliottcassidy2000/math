---
id: THM-3623
title: "Russell-cylinder even general-vertical-fold all-order closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let Q be even
  with Q(0)=-3/4 and
  Q(1)=Q(-1)=-3, and let 0!=H in t^2 C[t].  On the closed source plane
  q=Q(x)+H(t), w=t, no pair of regular functions on the exponent-one
  Russell target cylinder has a nonzero constant source Jacobian.  If
  k=ord_0(H) and lambda=[t^k]H, the nth comparison invoice occurs at
  source order k(n-1) and has coefficient
  -16(2lambda/3)^(n-1)/(n-1)!.  The affine-linear and static vertical
  controls are method/local boundaries, not global counterexamples.
  No JC(2) counterexample is claimed.
source: root / general_vertical_fold comparison-curve extension, 2026-08-21
audit: >
  PASS -- an independent hostile derivation checked the exact chart,
  epsilon-divisible error ledger, generalized coefficient and shift orders,
  all-order induction, and k=1/H=0 hostile boundaries; normal, optimized,
  and stored 276-gate transcripts are byte-identical.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3619-russell-cylinder-even-fold-higher-jet-staircase
related:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
  - THM-3614-russell-cylinder-collision-free-full-linear-projection-rigidity
script: 04-computation/jc2_russell_cylinder_general_vertical_fold_all_order_thm3623.py
output: 05-knowledge/results/jc2_russell_cylinder_general_vertical_fold_all_order_thm3623.out
script_sha256: a808ab75fd6d7893d1a010116b4e9202e6cf2e2a35b72522f6f3106f3a4984a4
output_sha256: 6c981576bac2e02691061d97e24cb967c5422d2b5595b5b80dc74c952e0f742e
hash_basis: raw LF bytes
---

# THM-3623 -- Russell-cylinder even general-vertical-fold all-order closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
extends the exact moving-triple
mechanism of THM-3619 from the quadratic vertical fold to every nonzero
polynomial vertical displacement having a critical point at the retained
triple collision.  Its quantifier over target pairs is unrestricted, but the
statement remains local to the specified family of closed source planes.

All rings, formal germs, derivatives, and differential forms are over `C`.

## 0. Statement, geometry, and scope

Let

```text
Q in C[x] even,        Q(0)=-3/4,        Q(1)=Q(-1)=-3,
0!=H in t^2 C[t].                                           (1)
```

Define

```text
E_(Q,H):A2_(x,t) -> A3_(x,q,w),
E_(Q,H)(x,t)=(x,Q(x)+H(t),t).                              (2)
```

The image is the principal hypersurface

```text
q-Q(x)-H(w)=0,                                             (3)
```

so `(2)` is a closed polynomial copy of `A2`.  Since `(1)` forces
`deg(H)>=2`, the projection to `(x,q)` has image `C[x,H(t)]`, which does not
contain `t`; hence the plane is genuinely not a polynomial graph
`w=h(x,q)`.

Let `Psi` be the stabilized Russell map of THM-3605.  For arbitrary regular
functions `F,G` on its exponent-one Russell target cylinder, put

```text
J_(Q,H)(F,G)
 =Jac_(x,t)(F o Psi o E_(Q,H),G o Psi o E_(Q,H)).          (4)
```

The theorem is

```text
                         J_(Q,H)(F,G) notin C*.            (5)
```

Write

```text
k=ord_0(H)>=2,              lambda=[t^k]H!=0,
q_n=Q^(n)(1).                                               (6)
```

The comparison forces

```text
q_1=9/2,
q_n=-(n+1)q_(n-1),                         n>=2,           (7)
```

and the `n`th constraint appears at source order `k(n-1)`.  Thus the
quadratic THM-3619 staircase is the special case `k=2,lambda=1`; no parity
or monomial assumption on `H` is needed.

The result does **not** cover non-even `Q`, `H'(0)!=0`, `H=0`, implicit
source planes outside `(2)`, or arbitrary stabilized surfaces.  The boundary
controls in Section 6 show failure of this local enlarged-form method, not
global target counterexamples.  In particular, `(5)` is not a proof of
`JC(2)` and constructs no Jacobian counterexample.

The information ledger is

| item | retained or lost |
|---|---|
| source | the closed plane `(3)` with one stable triple collision |
| target | the full exponent-one Russell cylinder |
| retained | even side parity, stable coordinate `w=t`, and the complete moving collision germ |
| varied | the vertical ramification order `k` and all higher coefficients of `H` |
| tested predicate | a nonzero constant ordinary `(x,t)` Jacobian for any regular target pair |
| decisive sidecar | the order `k(n-1)` moving-side comparison invoice |
| lost boundary | the `P`-direction when `H'(0)!=0`, or all side motion when `H=0` |

## 1. Local Darboux chart and exact general pullback

Before the polynomial Russell-cylinder isomorphism, retain the THM-3561
coordinates

```text
D=1+x^2q,
a=q/D^2,                   c=xD(D+2),
b=ac^2=(D-1)(D+2)^2,       e=a(b+4)=q(D+3).              (8)
```

At the common target point

```text
(b,c,e,w)=(0,0,-3,0),                                      (9)
```

the function `b+4` is a unit.  Thus

```text
a=e/(b+4),                     c,                     w    (10)
```

are regular local target parameters.  Direct differentiation gives the
global rational identity

```text
Jac_(x,q)(a,c)=-3.                                         (11)
```

Write the actual target two-form `dF wedge dG`, or more generally an
arbitrary formal target two-form, as

```text
Omega=P da wedge dc+K da wedge dw+R dc wedge dw.           (12)
```

On `(2)`, total `x` derivatives at fixed `t` use `q_x=Q'(x)`, while
`q_t=H'(t)` and `w_t=1`.  Equations `(11),(12)` give the exact pullback

```text
E_(Q,H)^*Omega=j(x,t) dx wedge dt,
j=-3H'(t)P+a_xK+c_xR.                                     (13)
```

At the middle collision branch `x=t=0`, evenness of `Q` and `(1)` give

```text
a=-3/4,             c=0,             a_x=0,             c_x=3. (14)
```

Suppose `(4)` is a nonzero constant and scale one output so that this
constant is `12`.  Since `H'(0)=0`, `(13),(14)` force

```text
R(-3/4,0,0)=4.                                            (15)
```

This is the first place where the criticality of `H` is load-bearing.  The
proof below permits every formal `P,K,R`, a strict enlargement of forms
`dF wedge dG`.  An obstruction in this larger universe therefore obstructs
all regular target pairs.

## 2. Exact comparison triple and epsilon-divisible errors

Put

```text
s=H(t),
g=(1-4s/3)^(-1/2),                X=g-1,
Q_infinity(x)=-3/4-9/(4x^2).                              (16)
```

The square-root branch with constant term `1` is unique in `C[[t]]`.  From
`(6)` one has

```text
X=(2lambda/3)t^k+O(t^(k+1)).                              (17)
```

The three comparison points

```text
gamma_0=(0,-3/4+s,t),
gamma_-=(-g,-3/g^2,t),          gamma_+=(g,-3/g^2,t)      (18)
```

share the exact target

```text
(a,c,w)(gamma_0)=(a,c,w)(gamma_-)
                 =(a,c,w)(gamma_+)=(-3/4+s,0,t).          (19)
```

Indeed, `D=1` at the middle and `D=-2` at the sides, while

```text
Q_infinity(g)+s=-3/g^2.                                  (20)
```

The actual fold points over `x=+-g` instead have

```text
q=Q(g)+s=-3/g^2+epsilon,
epsilon=Q(g)-Q_infinity(g),                               (21)
```

where evenness of `Q` makes the same `epsilon` valid on both sides.  The
error is not merely heuristic.  At the positive actual side one has exactly

```text
D=-2+g^2 epsilon,
c=g^3 epsilon(-2+g^2 epsilon),

a-(-3/4+s)
 =epsilon(3g^2 epsilon-8)/(4(g^2 epsilon-2)^2).           (22)
```

The denominator in `(22)` is a unit, and the negative side has the same `a`
and the opposite `c`.  Hence every target-value displacement in `(12)` is
divisible by `epsilon`.  The corresponding errors in `a_x,c_x` are also
`epsilon`-divisible because their denominators remain units at `D=-2`.

Let `j_-,j_0,j_+` denote the source germs of `(13)` in the local coordinates

```text
xi=x+1,                     xi=x,                    xi=x-1
```

at the three collision preimages, and define

```text
C_(Q,H)(t)=j_-(-X,t)+j_+(X,t)-2j_0(0,t).                 (23)
```

At the exact comparison target, parity gives

```text
a_x(-g)=-a_x(g),                  c_x(-g)=c_x(g),         (24)
```

and direct differentiation, retaining the actual side slope `Q'(g)`, gives

```text
c_x=12-2g^3Q'(g)                 at each side,
c_x=3                            at the middle.           (25)
```

Thus the exact comparison `P` terms cancel among the three branches, the two
side `K` terms cancel, and the middle `K` term vanishes.  With

```text
rho(t)=R(-3/4+s,0,t),                      rho(0)=4,       (26)
```

the divisibility ledger `(22)` yields formal series `U,V in C[[t]]` such
that

```text
C_(Q,H)(t)
 =rho(t)[18-4g^3Q'(g)]+epsilon U(t)+H'(t)epsilon V(t)

 =-4rho(t)g^3(Q'-Q_infinity')(g)
    +epsilon U(t)+H'(t)epsilon V(t).                     (27)
```

The terms have the following exact roles:

| coefficient | exact comparison | actual-side error |
|---|---|---|
| `P` | cancels with weights `1,1,-2` | `H' O(epsilon)` |
| `K` | cancels by `(24)` | `O(epsilon)` |
| `R` | gives the main term in `(27)` | `O(epsilon)` |

This is the complete nuisance-term invoice.

## 3. Generalized coefficient and vertical-source invoice

Let

```text
q_n^infinity=Q_infinity^(n)(1)
            =(-1)^(n-1)9(n+1)!/4.                        (28)
```

For `n>=1`, suppose the jets of `Q` and `Q_infinity` agree through order
`n-1` at `x=1`; when `n=1`, this means only their common value `-3`.  Put

```text
delta_n=q_n-q_n^infinity.                                (29)
```

Taylor expansion at `g=1+X` gives

```text
epsilon=O(X^n)=O(t^(kn)),

(Q'-Q_infinity')(1+X)
 =delta_n X^(n-1)/(n-1)!+O(X^n).                         (30)
```

Set

```text
N=k(n-1).                                                (31)
```

The `K,R` errors in `(27)` begin at order `kn=N+k`, and the `P` error begins
at order

```text
kn+k-1=N+2k-1.                                           (32)
```

Both lie strictly above `N`.  Using `(15),(17),(27),(30)` gives the exact
leading coefficient

```text
[t^N]C_(Q,H)(t)
 =-16(2lambda/3)^(n-1) delta_n/(n-1)!.                   (33)
```

For the intrinsic vertical version, define

```text
Delta_N(j)=[t^N]j_-(0,t)-2[t^N]j_0(0,t)+[t^N]j_+(0,t).  (34)
```

Suppose the source coefficients of total degree below `N` equal those of the
constant `12`.  A shifted source monomial `xi^i t^ell`, with `i>=1`, can
contribute to `[t^N]C_(Q,H)` only if

```text
N-ell>=ki.                                               (35)
```

Consequently

```text
i+ell<=N-(k-1)i<N,                                       (36)
```

because `k>=2`.  Every shifted contribution therefore belongs to an already
vanishing lower total-degree source row.  Only `i=0` survives, so

```text
[t^N]C_(Q,H)=Delta_N(j)                                  (37)
```

modulo those lower rows.

Under the previous side-jet matches, `(28)` obeys

```text
q_n^infinity=-(n+1)q_(n-1)^infinity,                     (38)
```

and hence, for `n>=2`,

```text
delta_n=q_n+(n+1)q_(n-1).                               (39)
```

Combining `(33),(34),(37),(39)` gives the generalized invoice

```text
Delta_(k(n-1))
 =-16(2lambda/3)^(n-1)/(n-1)!
    (q_n+(n+1)q_(n-1)),                       n>=2.       (40)
```

At `n=1`, the order-zero invoice is

```text
Delta_0=-16(q_1-9/2).                                    (41)
```

For `H=t^2`, `(40)` reduces exactly to the THM-3619 coefficient
`-2^(n+3)/(3^(n-1)(n-1)!)`.

## 4. All-order induction and polynomial contradiction

Assume for contradiction that `(4)` is the normalized constant `12`.
Then `C_(Q,H)=0` identically.  Equation `(41)` first forces

```text
q_1=9/2=q_1^infinity.                                    (42)
```

If the side jets match through order `n-1`, `(33)` has a nonzero coefficient
because `lambda!=0` and characteristic is zero.  It therefore forces

```text
q_n=q_n^infinity,                                        (43)
```

closing the induction for every `n`.

But every jet in `(28)` is nonzero, whereas the derivatives of a polynomial
vanish above its degree.  Thus no polynomial `Q` can have all the forced
jets.  Equivalently, the forced germ satisfies

```text
xQ'(x)+2(Q(x)+3/4)=0,                                    (44)
```

whose rational solution with `Q(1)=-3` is `Q_infinity`; a nonconstant
polynomial term `a_d x^d` would enter the left side of `(44)` with nonzero
coefficient `(d+2)a_d`.  This contradiction proves the provisional claim
`(5)`.

## 5. Ramification controls

The exact companion checks four deliberately different vertical
polynomials:

```text
H_2=t^2+5t^5,                 k=2, lambda=1,
H_3=-2t^3+7t^4,               k=3, lambda=-2,
H_4=3t^4/5-t^7,               k=4, lambda=3/5,
H_7=11t^7+t^11,               k=7, lambda=11.             (45)
```

For each it checks `(17),(32),(33),(36)` through `n=6`.  These controls show
that neither evenness of `H`, purity of its leading monomial, nor small
ramification order is being silently used.  The proof itself is all-order;
the finite gates are controls, not its replacement.

## 6. Sharp excluded boundaries

### 6.1 The `H'(0)!=0` enlarged-form escape

Let `h=H'(0)!=0`.  At the middle collision branch, `(13)` now gives

```text
12=-3hP_0+3R_0,                                          (46)
```

while the order-zero comparison gives only

```text
R_0(18-4q_1)=0.                                          (47)
```

Thus `R_0=0` is a genuine branch, rather than the forced normalization
`R_0=4` of `(15)`.  In the enlarged formal two-form universe it is realized
for every `Q` by

```text
Omega=-4/H'(w) da wedge dc,                              (48)
```

because `H'(w)` is a unit and `(13)` pulls `(48)` back identically to `12`.
For nonlinear `H`, `(48)` need not be closed and is only an enlarged-form
survivor.

The affine case is sharper.  If `H(t)=ht`, then

```text
F=a,                         G=-4c/h                     (49)
```

is an exact local Darboux pair with source Jacobian `12`, for every `Q`.
For example,

```text
H=t,                 Q=-3/4-9x^2/4                      (50)
```

has the collision values but `Q'(1)=-9/2`, contradicting the first forced
jet of `(7)` while `(49)` survives exactly.

This does **not** contradict the global graph-slice theorem: `a=e/(b+4)` is
regular only on the local chart `b+4!=0`, not necessarily on the whole affine
target cylinder.  Equations `(48)--(50)` are exact method/local boundaries,
not global target counterexamples.  If one separately assumes `R_0!=0`, the
direct moving comparison does force the recurrence even when `k=1`; however,
the vertical quotient `(37)` ceases to isolate because `(36)` can become an
equality at total degree `N`.

### 6.2 The static `H=0` boundary

If `H=0`, then

```text
s=0,                         g=1,                         X=0. (51)
```

The comparison triple never moves.  Central constancy gives
`R(-3/4,0,t)=4`, and `(27)` reduces exactly to

```text
C_(Q,0)(t)=4(18-4Q'(1)).                                 (52)
```

Thus it forces only `Q'(1)=9/2`; no higher derivative of `Q` occurs.  The
even polynomial

```text
Q_static=-3/4-27x^2/4+9x^4/2                            (53)
```

has all collision values and `Q_static'(1)=9/2`, but

```text
Q_static''(1)=81/2!=-27/2.                               (54)
```

The static comparison cannot distinguish `(53)` from the formal rational
germ at second order.  Again, this is a precise stopping witness for the
comparison method, not an asserted Darboux pair or global counterexample.

## 7. Exact verification

The companion verifies the chart identities, pullback `(13)`, comparison
geometry `(18)--(22)`, epsilon divisibility, four ramification scales through
six recurrence rows, every `P/K/R` order inequality, the shifted-source-row
gate, the rational germ, and both excluded boundaries.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_general_vertical_fold_all_order_thm3623.py
python3 -O 04-computation/jc2_russell_cylinder_general_vertical_fold_all_order_thm3623.py
```

The normal and optimized transcripts must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_general_vertical_fold_all_order_thm3623.out
```

The frozen companion reports

```text
CHECKS=276
RESULT=PASS
```

with zero Python assertion statements.  They are finite controls for the
all-order proof, not a replacement for it.
