---
id: THM-3608
title: "Russell-cylinder nonlinear target-shear formal rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  polynomial source graph, a
  first output L=u(C)B+F(C) with u(0) nonzero and a second output
  M=S+K(L,C) cannot have nonzero constant ordinary Jacobian.  This includes
  arbitrary separated shears K=H(C)+G(L), all fibre-triangular polynomial
  target automorphisms of this order, and even noninvertible u with u(0)
  nonzero.  The formal arm equation forces
  u(C)(ac+h)+(F(C)-F(0))/C to be constant; the resulting rational graph has
  the unavoidable D=0 pole.  An exact formal-algebraic hostile shows that
  u(0) nonzero is sharp for this proof.  Arbitrary nonlinear target
  coordinates that identify the two arms remain open.
source: root/nonlinear_cylinder_shears target-shear niche, 2026-08-21
audit: >
  PASS.  An independent hostile audit rederived the full u,F,K determinant
  with K_T cancellation and fixed-T K_c, the boundary degree argument, the
  first formal-arm recurrence and completion injection, the forced endpoint
  and D-pole, the collision, the u(0)=0 formal-algebraic hostile, and the
  arm-identifying coordinate boundary.  Normal and optimized replays are
  byte-identical to the stored 167-gate transcript; the AST has no assertion
  gates, and documentation and diff checks pass.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3607-russell-cylinder-mixed-projection-degree-seven-gate
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_russell_cylinder_nonlinear_target_shear_thm3608.py
output: 05-knowledge/results/jc2_russell_cylinder_nonlinear_target_shear_thm3608.out
script_sha256: ffe3c69d6e243f67778ade2df0fe7d18790e76c3db85b481cec1d0d40840efb2
output_sha256: 8dcf9e228548ac36968ce3e932171e6f5365086c24412a041b130ed70102acc7
hash_basis: raw LF bytes
---

# THM-3608 -- Russell-cylinder nonlinear target-shear formal rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
strictly extends the mixed linear face of THM-3607.  All rings and
derivatives below are over `C`.

## 0. Statement and preservation/loss ledger

Retain the THM-3561 functions

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
b=(D-1)(D+2)^2=ac^2,
e=q(D+3)=a(b+4),                    Jac_(x,q)(a,c)=-3.  (1)
```

For an arbitrary polynomial graph `h in C[x,q]`, put

```text
A=ac+h,
B_h=b+ch=cA,
S_h=a+3A^2/4+cA^3/8.                                  (2)
```

The last two expressions are polynomial in `(x,q)`, by the Russell-cylinder
formula in THM-3605, even though `A` itself can have a pole.  Let

```text
u,F in C[c],       u(0)!=0,       K(T,c) in C[T,c],
L_h=u(c)B_h+F(c),
M_h=S_h+K(L_h,c).                                      (3)
```

Then the theorem is

```text
Jac_(x,q)(L_h,M_h) notin C*.                           (4)
```

Taking `K(T,c)=H(c)+G(T)` gives exactly the nonlinear target shears in the
motivating question.  The stronger `K(T,c)` scope is free: derivatives of
`K` in the `L` direction cancel from a two-dimensional Jacobian.

The information ledger is

| item | retained or lost |
|---|---|
| source | a polynomial graph `w=h(x,q)` in the stabilized THM-3561 source |
| target pair | `(u(C)B+F(C), S+K(u(C)B+F(C),C))` |
| retained | polynomial outputs; the full stable collision whenever the graph values agree |
| tested predicate | a nonzero constant ordinary planar Jacobian |
| decisive sidecar | the explicit `x=0` formal-arm completion and the unit `u(0)` |
| discarded | `Y`, the output `C`, and for nonconstant `u` invertibility of the target operation |
| sharp punctured survivor | the unique formal solution `Z=t`, displayed in Section 4 |
| first unsafe exit | `u(0)=0`, where a nonconstant formal Keller solution exists |

No collision assumption is used in `(4)`.  Hence `(4)` also holds for every
polynomial graph satisfying the full-stable collision condition

```text
h(0,-3/4)=h(1,-3)=h(-1,-3).                           (5)
```

## 1. Exact nonlinear transport equation

Write

```text
F_0=F(0),            R(c)=(F(c)-F_0)/c in C[c],
Z=u(c)A+R(c),        L_h=F_0+cZ,
Q(A,c)=3A/2+3cA^2/8.                                  (6)
```

Because `u(0)!=0`, it is a unit in `C[[c]]`, and locally

```text
A=(Z-R)/u,
A_a=Z_a/u,
A_c=(Z_c-R'-u'A)/u.                                   (7)
```

In `(a,c)` coordinates, differentiate `(2),(3)`.  In `K(T,c)`, write `K_c`
for the derivative with respect to its second argument while keeping `L`
fixed.  The `K_L` terms cancel, and direct expansion gives

```text
det partial(L_h,M_h)/partial(a,c)
 =-(Z+cZ_c+V(Z,c)Z_a),                                (8)

V(Z,c)=Q(A,c)(Z+cR'+cu'A)/u
       -cA^3/8-c K_c(F_0+cZ,c).                       (9)
```

Combining `(8)` with `Jac_(x,q)(a,c)=-3` fixes the source sign:

```text
Jac_(x,q)(L_h,M_h)=3(Z+cZ_c+V(Z,c)Z_a).               (10)
```

For `u=1`, `F=cR`, and `K(T,c)=H(c)+G(T)`, formula `(9)` becomes

```text
V=Q(Z+cR')-cA^3/8-cH'(c),        A=Z-R,               (11)
```

because `G(L)` contributes no Jacobian.  Constant `R=rho` and linear
`H=tau*c` recover the THM-3607 transport equation exactly.

## 2. The boundary polynomial is constant

Suppose the left side of `(10)` is the nonzero constant `3t`.  On the arm
`x=0`, equations `(1)` give `a=q,c=0`.  Put

```text
u_0=u(0),       R_0=R(0),       f(a)=Z(a,0) in C[a].  (12)
```

At `c=0`, equations `(7),(9)` give

```text
A=(f-R_0)/u_0,
V(f,0)=3f(f-R_0)/(2u_0^2).                            (13)
```

Restricting `(10)` therefore gives

```text
f+3f(f-R_0)f'/(2u_0^2)=t.                            (14)
```

If `deg f=d>=1`, the derivative term in `(14)` has degree `3d-1`,
strictly greater than `d`, and has nonzero leading coefficient.  It cannot
cancel `f-t`.  Thus `f` is constant, and `(14)` forces

```text
f=t.                                                   (15)
```

This is the only place where polynomiality along the boundary, rather than
mere formal regularity, supplies a degree argument.

## 3. Formal uniqueness forces `Z=t`

THM-3607 gives an explicit formal-arm isomorphism

```text
C[a][[c]] ~= C[q][[x]],                               (16)
```

obtained from the unique `Delta=1+O(c^2)` satisfying

```text
ac^2=(Delta-1)(Delta+2)^2,
x=c/(Delta(Delta+2)),             q=a Delta^2.        (17)
```

Thus `Z` belongs to `C[a][[c]]`.  Set `W=Z-t`.  By `(15)`, `c` divides
`W`.  If `W` is nonzero, choose its first coefficient

```text
W=z_N(a)c^N+O(c^(N+1)),       N>=1,       z_N!=0.     (18)
```

The coefficient of `c^N` in `(10)/3-t=0` is

```text
(N+1)z_N + 3t(t-R_0)z_N'/(2u_0^2)=0.                 (19)
```

The first term has degree `deg z_N`, while the derivative term has smaller
degree or vanishes.  Hence `(19)` has no nonzero polynomial solution.  This
contradicts `(18)`, so

```text
Z=t                                                       (20)
```

in the completion.  The localized source domain embeds in its `x`-adic
completion, so `(20)` holds already in its function field.

## 4. The forced graph and its sharp punctured control

Equations `(2),(6),(20)` force

```text
h=(t-R(c))/u(c)-ac.                                   (21)
```

The first term in `(21)` is regular at the prime divisor `D=0`, because
`c=0` there and `u(0)!=0`.  By contrast, `ac=b/c` has order `-1` there:
`c=xD(D+2)` has a simple `D` factor while `b=-4` modulo `D`.  Therefore
`(21)` has a genuine simple `D=0` pole and cannot be a polynomial `h`.  This
proves the claim `(4)`.

The rejected solution is exact and sharp on the punctured source.  For any
`t in C*`, set

```text
A_t(c)=(t-R(c))/u(c),        h_t=A_t(c)-ac.            (22)
```

On the open set `D*u(c)!=0`, one has

```text
L_t=F_0+tc,
M_t=a+3A_t(c)^2/4+cA_t(c)^3/8+K(F_0+tc,c),
Jac_(x,q)(L_t,M_t)=3t.                                (23)
```

Indeed the determinant in `(a,c)` is `-t`.  At all three collision points,
`a=-3/4,c=0,ac=0`, so `h_t` has the common value
`(t-R_0)/u_0`; the rational endpoint retains the full stable collision.

For the motivating subclass `u=1,F=cR,K=H+G`, formulas `(21)--(23)` reduce
to

```text
h=t-R(c)-ac,       L=tc,
M=a+3(t-R(c))^2/4+c(t-R(c))^3/8+H(c)+G(tc).           (24)
```

Thus a nonlinear target shear changes only the one-variable tail of the
punctured Darboux pair; it cannot pay the missing `D=0` pole invoice.

## 5. Why `u(0)!=0` is a sharp formal boundary

The unit condition is not cosmetic.  Take

```text
u(c)=c,        F(c)=c,        K=0,
L=cB+c=c(1+cA),             M=S.                     (25)
```

Put `v=ac^2` and let `phi=phi(v) in v C[[v]]` be the unique formal solution
of

```text
v=-phi/2-27phi^2/32-phi^3/8.                         (26)
```

Existence and uniqueness follow because the right side has derivative
`-1/2` at zero.  Its first terms are

```text
phi=-2v-27v^2/4-697v^3/16-89775v^4/256+O(v^5).       (27)
```

Define

```text
A=phi(ac^2)/c in C[a][[c]],       Z=1+cA=1+phi.       (28)
```

This `Z` is nonconstant.  Nevertheless, regarding `(phi,c)` as parameters,

```text
det partial(L,S)/partial(phi,c)
   =(6phi^2+27phi+8)/(16c^2),
det partial(a,c)/partial(phi,c)
   =-(6phi^2+27phi+8)/(16c^2).                        (29)
```

Hence

```text
det partial(L,S)/partial(a,c)=-1,
Jac_(x,q)(L,S)=3.                                     (30)
```

This is a formal-algebraic hostile, not a polynomial graph and not a
counterexample to `(4)` without its hypothesis.  It proves exactly that the
single-arm implication `constant Jacobian => Z is constant` fails once
`u(0)` is allowed to vanish.  Any polynomial no-filling theorem in that
singular branch needs a different, multiarm argument.

## 6. Polynomial target automorphisms and the arbitrary-coordinate exit

The theorem contains every fibre-triangular polynomial target automorphism
of the displayed order:

```text
(B,C,S) |-> (alpha B+F(C), C, beta S+K(alpha B+F(C),C)),
alpha,beta in C*.                                     (31)
```

The constant `beta` is removed by rescaling the second output.  More
generally, `(3)` allows nonconstant `u` with `u(0)!=0`, even though such a
map is not a polynomial automorphism: the target three-Jacobian contains
the nonunit factor `u(C)`.  Constant invertible recombinations of the two
outputs also preserve the conclusion.

Arbitrary polynomial coordinate changes are not covered.  A cheapest exact
hostile is the polynomial coordinate

```text
L_*=C+B(B+4),                                          (32)
```

which is a coordinate because `(B,C) |-> (B,L_*)` has inverse
`C=L_*-B(B+4)`.  On a graph, the Danielewski relation gives

```text
B_h(B_h+4)=c Y_h,
Y_h=ce+(2b+4)h+ch^2,
L_*=c(1+Y_h).                                         (33)
```

At `c=0`, `(32)` identifies the two target arms `B=0` and `B=-4`.  Thus it
destroys exactly the arm label used by the pole contradiction: a formal
central-arm inverse no longer determines which global algebraic branch
contains `B_h`.  This does not produce a polynomial Keller graph, but it
blocks any inference from the present one-arm proof to arbitrary target
coordinates.

The following remain **OPEN** here:

- the singular coefficient branch `u(0)=0` for polynomial source graphs;
- first outputs nonlinear in `B`, including `(32)`;
- projections involving `Y`, implicit non-graph planes, other cylinder
  isomorphisms, or changes of the fixed source coordinates.

## 7. Exact companion and reproduction

The deterministic companion verifies the Russell arm normal form, the full
symbolic `u,F,K` determinant and sign, concrete polynomial-source rows, the
boundary and first-coefficient degree gates, the sharp punctured controls
and collision values, the `D=0` residue, the formal-algebraic hostile and its
series coefficients, and the arm-identifying coordinate `(32),(33)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_nonlinear_target_shear_thm3608.py
python3 -O 04-computation/jc2_russell_cylinder_nonlinear_target_shear_thm3608.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_nonlinear_target_shear_thm3608.out`.
