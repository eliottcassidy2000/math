---
id: THM-3612
title: "Russell-cylinder even-fold non-graph collision-jet rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The closed planes q=Q(x)+t^2, w=t are genuine non-graphs over (x,q)
  and retain the full THM-3561 triple collision at t=0.  For every even
  Q with Q(0)=-3/4 and Q(±1)=-3, and for every pair of regular functions
  on the exponent-one Russell target cylinder, a nonzero constant source
  Jacobian is impossible if Q'(1)!=9/2.  At the unique first-jet escape
  Q'(1)=9/2, it remains impossible if Q''(1)!=-27/2, by an exact
  second-jet discrete-difference invoice.  The doubly tuned locus
  Q'(1)=9/2, Q''(1)=-27/2 remains OPEN and has an explicit second-order
  formal survivor; no JC(2) counterexample is claimed.
source: root / implicit_cylinder_slice non-graph wildcard, 2026-08-21
audit: >
  PASS -- independent hostile audit reconstructed the fold geometry, complete
  first- and second-jet systems, cubic-coefficient cancellation, sharp
  survivors, and OPEN doubly tuned boundary; normal, optimized, and stored
  241-gate transcripts are byte-identical.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
related:
  - THM-3607-russell-cylinder-mixed-projection-degree-seven-gate
  - THM-3610-russell-cylinder-full-linear-projection-collision-rigidity
  - THM-3611-russell-cylinder-arm-separating-nonlinear-first-coordinate-rigidity
script: 04-computation/jc2_russell_cylinder_even_fold_jet_rigidity_thm3612.py
output: 05-knowledge/results/jc2_russell_cylinder_even_fold_jet_rigidity_thm3612.out
script_sha256: fa210b98b493f5426f77567769c8f3356440de0b2dad00fca3c029b71efa0f43
output_sha256: abb71fdac67e97e0c759cf970e3d787b0cbe60c802638303e1d237eab4f59c4b
hash_basis: raw LF bytes
---

# THM-3612 -- Russell-cylinder even-fold non-graph collision-jet rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem closes an infinite, genuinely non-graph family of polynomial
planes using only the common collision point and its branch jets.  The
quantifier over target pairs is unrestricted: the outputs may be arbitrary
regular functions on the Russell target cylinder, not merely linear
projections of its four displayed coordinates.

All rings, germs, derivatives, and differential forms are over `C`.

## 0. Statement and dimension ledger

Retain the THM-3561 compiler functions

```text
D=1+x^2q,
c=xD(D+2),
b=(D-1)(D+2)^2,
e=q(D+3),                         c^2e=b(b+4).            (1)
```

The Russell map from THM-3605 has polynomial target coordinates

```text
B=b+cw,
C=c,
Y=ce+(2b+4)w+cw^2,
S=((b+2)(e+3w^2)+cw(3e+w^2))/8,                         (2)
```

which satisfy `CY=B(B+4)`.  Let `Q in C[x]` be even and satisfy

```text
Q(0)=-3/4,                       Q(1)=Q(-1)=-3.           (3)
```

Define the polynomial embedding

```text
E_Q:A2_(x,t) -> A3_(x,q,w),
E_Q(x,t)=(x,Q(x)+t^2,t).                                  (4)
```

Its image is the principal hypersurface

```text
H=q-Q(x)-w^2=0,                                           (5)
```

so `C[x,q,w]/(H)=C[x,w]` and `(4)` is a closed copy of `A2`.
The projection to `(x,q)` corresponds to

```text
C[x,q] -> C[x,t],             q |-> Q(x)+t^2,             (6)
```

whose image is `C[x,t^2]`.  It is generically two-to-one under
`t <-> -t`, and `t` is not in that image.  Thus `(4)` is genuinely not a
polynomial graph `w=h(x,q)`.

Put

```text
p_- =(-1,0),                 p_0=(0,0),
p_+ =(1,0)                                                        (7)
```

in the `(x,t)`-plane.  Equations `(1)--(4)` send all three to

```text
y_0=(B,C,Y,S)=(0,0,0,-3/4).                               (8)
```

For arbitrary

```text
F,G in O(Y_1 x A1_S),
J_Q(F,G)=Jac_(x,t)(F o Psi o E_Q,G o Psi o E_Q),          (9)
```

where `Psi` is the stabilized map of THM-3605.  Write

```text
v=Q'(1),                         r=Q''(1).                (10)
```

Then the theorem is

```text
v!=9/2                       ==> J_Q(F,G) notin C*,
v=9/2 and r!=-27/2           ==> J_Q(F,G) notin C*.       (11)
```

The conclusion holds for **every regular target pair**.  The doubly tuned
case

```text
v=9/2,                           r=-27/2                  (12)
```

is not covered and remains **OPEN**.

The information ledger is

| item | retained or lost |
|---|---|
| source | a closed polynomial `A2` in the stabilized source `A3` |
| target | the full smooth Russell target threefold near `(8)` |
| retained | polynomiality, one full stable triple collision, and a unit residue two-form |
| lost | graph structure over `(x,q)` and the rest of the stable collision curve |
| tested predicate | nonzero constant ordinary source Jacobian for any regular target pair |
| decisive sidecar | the three branch tangent planes, then one discrete vertical second jet |
| open boundary | the doubly tuned jet locus `(12)` and non-even folds on its broader tangent-match locus |

## 1. Stable-volume residue on the fold

THM-3605 proves

```text
Psi^*(omega_1 wedge dS)=-3 dx wedge dq wedge dw.          (13)
```

For `H` in `(5)`, direct exterior multiplication gives

```text
dH wedge (3 dx wedge dw)=-3 dx wedge dq wedge dw.         (14)
```

Hence the fold carries the nowhere-vanishing residue two-form

```text
eta_Q=3 dx wedge dt.                                      (15)
```

This does not itself produce a target Darboux pair.  Selecting the plane
also selects a normal direction; the three pushed normal sidecars at the
common target point need not agree.  Sections 2--3 measure that mismatch.

## 2. The universal first-jet obstruction

Let

```text
Z=S+3/4.                                                  (16)
```

At `y_0`, the differential of the target equation is

```text
d(CY-B(B+4))=-4dB.                                       (17)
```

Thus the target is smooth there, `dB=0`, and `(dC,dY,dZ)` is a cotangent
basis.  Since `Q` is even,

```text
Q'(0)=0,                        Q'(-1)=-v.                (18)
```

Differentiating `(1),(2),(4)` gives the exact branch table

```text
         dC                         dY                         dZ
p_0      3 dx                       -9 dx+4 dt                 0
p_+      (12-2v) dx                 (-36+6v) dx+4 dt           (9-v)dx/2
p_-      (12-2v) dx                 (-36+6v) dx+4 dt          -(9-v)dx/2.
                                                                  (19)
```

For any regular target pair, write its common cotangent two-form as

```text
(dF wedge dG)_(y_0)
 =A dC wedge dY+K dC wedge dZ+R dY wedge dZ.             (20)
```

The three pullback coefficients on `dx wedge dt` are

```text
j_0=12A,
j_+=(48-8v)A+2(v-9)R,
j_-=(48-8v)A+2(9-v)R.                                   (21)
```

Suppose these are one common nonzero constant.  If `v!=9`, subtracting the
last two equations gives `R=0`; comparison with `j_0` then forces
`v=9/2`.  If `v=9`, the side values are `-24A` while the middle value is
`12A`, so equality forces `A=0` and the common value is zero.  Therefore
`v=9/2` is the unique first-jet escape, proving the first implication in
`(11)`.

At this escape, `(21)` forces `R=0` and `A!=0`; the coefficient `K` remains
free.  After adding constants, scaling one output, and applying a constant
determinant-one target basis change, any hypothetical pair can therefore
be normalized to

```text
J_Q(F,G)=12,
F=C+O(m^2),                       G=Y+beta Z+O(m^2),       (22)
```

in the completed local ring at `y_0`.

## 3. The exceptional second-jet invoice

Write the quadratic Taylor parts in the convention

```text
F=C+f_CC C^2/2+f_CY CY+f_CZ CZ
    +f_YY Y^2/2+f_YZ YZ+f_ZZ Z^2/2+O(m^3),

G=Y+beta Z+g_CC C^2/2+g_CY CY+g_CZ CZ
    +g_YY Y^2/2+g_YZ YZ+g_ZZ Z^2/2+O(m^3).              (23)
```

Let `r=Q''(1)=Q''(-1)`.  The six equations obtained by setting the two
first source derivatives of the Jacobian to zero at each of the three
preimages have rank five.  Their complete solution is

```text
f_CC=-g_CY-3 beta/8,
f_CY=-g_YY-7 beta/32,
f_CZ=1+4r/27-7 beta^2/64-beta g_YY/2-g_YZ/2,
(f_YY,f_YZ,f_ZZ)=ell(1,beta,beta^2),                    (24)
```

with all six `g`-Hessian entries and `ell` free.

Now allow arbitrary cubic target jets in both `F` and `G`.  For a source
Jacobian `J`, let `[t^2]_(p_i)J` denote the coefficient of `t^2` in its
restriction to the vertical germ through `p_i`.  Exact substitution of
the three source germs into `(23),(24)` gives

```text
[t^2]_(p_-)J-2[t^2]_(p_0)J+[t^2]_(p_+)J
                  =-16(2r+27)/3.                       (25)
```

Every cubic target coefficient cancels from `(25)`; target terms of degree
at least four cannot affect a second source jet.  If `J` were the constant
`12`, the left side would vanish.  Therefore `(25)` is impossible whenever
`r!=-27/2`, proving the second implication in `(11)`.

This is an all-regular-function statement: `(C,Y,Z)` are formal parameters
at the smooth point `y_0`, so `(23)` is the unrestricted target Taylor
expansion.  The proof does not assume that `F,G` are linear or polynomial
in the four displayed ambient coordinates.

## 4. Sharp controls and the honest stopping boundary

### 4.1 Cheapest first-jet no-go

The quadratic interpolation

```text
Q_0(x)=-3/4-9x^2/4                                      (26)
```

has `Q_0'(1)=-9/2`, so Section 2 kills every regular target pair already
at the common cotangent space.  Explicitly, the three values in `(21)` are

```text
12A,                         84A-27R,       84A+27R.     (27)
```

### 4.2 The `-576` first-order survivor

Put

```text
Q_*(x)=-3/4-27x^2/4+9x^4/2.                            (28)
```

Then

```text
Q_*'(1)=9/2,                      Q_*''(1)=81/2.         (29)
```

The global target polynomials

```text
F_*=C+7C(S+3/4),                  G_*=Y                  (30)
```

give source Jacobian `12` with both first source derivatives zero at all
three collision preimages.  Their vertical `t^2` coefficients are

```text
(-141,147,-141)              at (p_-,p_0,p_+),          (31)
```

whose discrete second difference is `-576`, exactly `(25)` at `r=81/2`.
Thus `(30)` is a sharp first-order survivor, not a Keller pair.

### 4.3 Doubly tuned second-order survivor -- OPEN

The even sextic

```text
Q_dag(x)=-3/4-27x^2/2+18x^4-27x^6/4                    (32)
```

satisfies `(3)` and the doubly tuned invoices

```text
Q_dag'(1)=9/2,                    Q_dag''(1)=-27/2.      (33)
```

Hence `(25)` vanishes.  This is not merely a vacuous equality.  Put
`Z=S+3/4` and take the global target polynomials

```text
F_dag=C-CZ-209C^3/144+C^2Y/16+7CY^2/64-755CZ^2/108,
G_dag=Y.                                                   (34)
```

For the `Q_dag` fold, the source Jacobian of `(34)` is `12` modulo the
third power of the maximal ideal at each of `p_-,p_0,p_+`; equivalently,
all source terms of total degree one and two vanish at all three branches.
This exact second-order survivor demonstrates that the boundary `(12)` is
real.  No global constant-Jacobian claim is made: higher jets of `(34)` are
nonzero.  The next honest test is the simultaneous third- and higher-jet
recursion on `(12)`.

## 5. Affine closure and collision-losing positive controls

The three collision points in the original `(x,q)`-plane are noncollinear:

```text
det [[1,-9/4],[-1,-9/4]]=-9/2.                          (35)
```

Therefore any affine two-plane in `A3_(x,q,w)` containing one common-`w`
point over all three collision branches is the unique affine span
`w=w_0`.  It is a constant graph, already routed to the graph analysis of
THM-3605, rather than an affine non-graph exit.  This observation does not
claim the still-open arbitrary nonlinear-`Y` graph case.

Two coordinate planes show why non-graphness by itself is not an
obstruction.

On `x=0`, one has

```text
(B,C,Y,S)=(0,0,4w,q+3w^2/4),
Jac_(q,w)(Y,S)=-4,                                      (36)
```

with inverse `w=Y/4`, `q=S-3Y^2/64`.  This plane contains only the middle
collision branch.  On `q=0`, the polynomial inverse-cylinder coordinate

```text
W=Y(B+2)/8-CS                                           (37)
```

restricts to `w`, and

```text
(C,W)=(3x,w),                     Jac_(x,w)(C,W)=3.      (38)
```

This plane contains none of the three collision points.  The fold theorem
therefore isolates collision-compatible branch geometry, not a blanket
prohibition on polynomial planes in the cylinder.

## 6. Scope and open exits

The result closes arbitrary regular target pairs on every even fold `(4)`
outside the doubly tuned locus `(12)`.  Its conclusions are global no-go
statements proved by local contradiction at the retained collision.

It does not close:

- the doubly tuned even folds `(12)`, for which `(34)` survives through
  second source order;
- non-even folds on the broader first-jet tangent-match locus;
- implicit polynomial planes not equivalent to a quadratic fold over the
  stable coordinate; or
- source embeddings that do not present one coordinate as `w=t`.

All four exits are **OPEN**.  The first should be tested by a simultaneous
higher-jet recursion before any global degree search.

## 7. Exact companion and reproduction

The deterministic exact companion checks:

- the compiler and Russell target relations;
- closedness, generic two-sheetedness, collision retention, and the residue
  sign for the fold;
- every entry and sign in the common-cotangent branch table `(19)--(21)`;
- the complete exceptional Hessian parametrization `(24)`;
- cancellation of every arbitrary cubic target coefficient in `(25)`;
- the `Q_0`, `Q_*`, and `Q_dag` value and jet invoices;
- the sharp first-order and second-order survivors `(30),(34)`; and
- the affine-span and both collision-losing coordinate-plane controls.

Finite rows verify the displayed algebra.  The quantifier over all regular
target pairs is discharged by the unrestricted formal Taylor expansions
and the first nonzero collision-jet contradiction.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_even_fold_jet_rigidity_thm3612.py
python3 -O 04-computation/jc2_russell_cylinder_even_fold_jet_rigidity_thm3612.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_even_fold_jet_rigidity_thm3612.out`.
