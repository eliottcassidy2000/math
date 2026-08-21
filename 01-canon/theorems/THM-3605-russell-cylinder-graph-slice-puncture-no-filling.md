---
id: THM-3605
title: "Russell cylinder graph-slice and puncture no-filling gate"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.
  For the two-arm exponent-two Danielewski surface behind THM-3561, the
  Russell cylinder formula gives an explicit polynomial isomorphism to the
  cylinder over the exponent-one surface and transports the residue volume
  form up to sign.  The stabilized THM-3561 map remains etale and retains
  its triple collision.  In the only arm plane containing that collision,
  no polynomial graph over the source A2 makes the induced two-form a unit.
  Polynomial graphs over the punctured Darboux plane are classified exactly;
  the unit-form graphs merely apply a triangular change to the original
  rational Keller pair and retain its double pole.  The fixed second
  coordinate has no polynomial Jacobian mate.  Non-graph slices and
  projections genuinely mixing the cylinder coordinate remain open.
source: root / arm_chart_compactification stable-cylinder wildcard, 2026-08-21
audit: >
  PENDING INDEPENDENT AUDIT.  The exact companion checks both polynomial
  cylinder maps and their compositions in the hypersurface quotient rings,
  residue-volume transport, the THM-3561 specialization and collision curve,
  all graph identities, weight-ODE divisibility obstructions, pole hostiles,
  and positive punctured controls.  Normal and optimized replays are
  byte-identical to the stored output.  These checks do not promote the
  theorem before an independent proof audit.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
citation:
  - "Moser-Jauslin and Poloni, Isomorphisms between cylinders over Danielewski surfaces, arXiv:2002.12202, Section 5.1."
script: 04-computation/jc2_russell_cylinder_graph_slice_no_filling_thm3605.py
output: 05-knowledge/results/jc2_russell_cylinder_graph_slice_no_filling_thm3605.out
script_sha256: e66938c4da6394dd20aafd2e5e074846176eadb94093c95ce5fad41a47a692e8
output_sha256: b984b8e876aecb310a33f118b2e73bfd736ffa3712f72bcbf0c9fc1724609c40
hash_basis: raw LF bytes
---

# THM-3605 -- Russell cylinder graph-slice and puncture no-filling gate

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT
AUDIT.**  Nothing below is a polynomial counterexample to `JC(2)`.  The
cylinder map is three-dimensional, and every successful planar graph found
below is defined only on the punctured source of THM-3561.

All rings, varieties, and differential forms are over `C`.

## 0. Truth and dimension ledger

There is one cited input:

```text
Moser-Jauslin--Poloni, Section 5.1:
the Russell cylinder formula for W_(2,P) x A1 ~= W_(1,P) x A1.       (1)
```

The specialization of that formula, its displayed polynomial inverse, both
quotient-ring compositions, the volume identity, and every graph obstruction
below are internal derivations.  The citation is not used for any planar
no-filling conclusion.

The dimension ledger is

```text
A3_(x,q,w) --etale, noninjective--> Y_2 x A1
                                      |
                                      | Russell isomorphism
                                      v
                                    Y_1 x A1.                         (2)
```

Neither cylinder in `(2)` is asserted to be `A3`.  A graph inside one arm
plane times `A1` is an `A2`, but its pullback need not extend across the
puncture in the original `(x,q)`-plane.

## 1. The explicit exponent-two to exponent-one cylinder

Put

```text
P(T)=T(T+4),

Y_2=Spec C[b,c,e]/(c^2 e-P(b)),
Y_1=Spec C[B,C,Y]/(C Y-P(B)).                           (3)
```

The Bezout data used in the Russell formula are

```text
U(T)=(T+2)/8,             V(T)=-1/4,
U P'+V P=1,
g(T)=T-P(T)U(T)=-T^2(T+6)/8.                           (4)
```

Define a point map

```text
Theta:Y_2 x A1_w -> Y_1 x A1_S                        (5)
```

by

```text
C=c,
B=b+cw,
Y=ce+(2b+4)w+cw^2,
S=((b+2)(e+3w^2)+cw(3e+w^2))/8.                       (6)
```

These are polynomials and satisfy `CY=P(B)`.  The inverse is also polynomial.
Put

```text
H(T)=(T-2)(T+6)/64,             P(g(T))=P(T)^2 H(T).  (7)
```

Then, on `Y_1 x A1_S`,

```text
b=g(B)+C^2S,
w=Y(B+2)/8-CS,
e=Y^2H(B)+2(g(B)+2)S+C^2S^2.                           (8)
```

Indeed, `(7)` makes the last expression regular without division by `C`.
Substitution of `(8)` into `(6)`, and of `(6)` into `(8)`, gives the identity
modulo the respective relations in `(3)`.  Thus `(5)` is a polynomial
isomorphism; only its Russell-form starting point is cited, while `(4)`--`(8)`
give a self-contained verification of this specialization.

### 1.1 Stable residue-volume transport

Use the global residue two-forms

```text
omega_2=dc wedge db/c^2,          omega_1=dC wedge dB/C.             (9)
```

On the dense set `c!=0`, eliminate `e=P(b)/c^2`.  Direct differentiation of
`(6)` gives

```text
det partial(C,B,S)/partial(c,b,w)=-1/c.               (10)
```

Consequently

```text
Theta^*(omega_1 wedge dS)=-omega_2 wedge dw.           (11)
```

Both sides are regular top forms, so equality on the dense set proves `(11)`
globally.  This is a unit **three-form** statement.  Restricting a three-form
to a graph does not by itself produce a planar Keller pair; the arm two-form
must be checked separately.

## 2. Stabilizing the THM-3561 collision

Retain the exact THM-3561 functions

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
b=ac^2=(D-1)(D+2)^2,
e=a(b+4)=q(D+3).                                      (12)
```

They satisfy

```text
c^2e=P(b),                         Jac_(x,q)(a,c)=-3. (13)
```

Although `a` is rational, `(b,c,e)` is polynomial and defines the etale map
`Phi:A2->Y_2` of THM-3561.  Stabilize and compose:

```text
Psi=Theta o (Phi x id):A3_(x,q,w)->Y_1 x A1_S.        (14)
```

This is a polynomial etale threefold map.  From `(11)` and `(13)`,

```text
Psi^*(omega_1 wedge dS)=-3 dx wedge dq wedge dw.       (15)
```

The three source points

```text
r_0=(0,-3/4),          r_+=(1,-3),          r_-=(-1,-3)              (16)
```

all have

```text
(a,c)=(-3/4,0),              (b,c,e)=(0,0,-3).         (17)
```

For every common stable value `w`, their common image under `(14)` is

```text
(C,B,Y,S)=(0,0,4w, 3(w^2-1)/4).                       (18)
```

Thus the original triple collision becomes three source lines mapping
pointwise to the same target curve.  This is stable noninjectivity in
dimension three, not a planar polynomial collision.

## 3. The only collision-bearing arm plane

The common image `(18)` lies on the `B=0` arm of `Y_1`.  Hence the only arm
plane containing it is

```text
U_0 ~= A2_(A,C),
A=B/C,
B=CA,
Y=A(CA+4).                                             (19)
```

The `B=-4` arm chart deletes the entire collision line and cannot retain
`(18)`.  On the root-zero chart, `(6)` and `b=ac^2` give

```text
A=ac+w.                                                (20)
```

After `(12)`, this is

```text
A=xq(D+2)/D+w.                                         (21)
```

If `w=h(x,q)` is polynomial, `(21)` has a genuine simple pole along `D=0`.
Indeed, on that curve `x` is a unit, `x^2q=-1`, and

```text
xq(D+2)=-2/x !=0.                                     (22)
```

No polynomial `h` changes this principal part.  Thus no polynomial graph
over the full source `A2_(x,q)` even makes the transported arm coordinate
regular.

## 4. No polynomial source graph has a unit arm form

There is a stronger differential obstruction.  Let

```text
Gamma_h={w=h(x,q)},                    h in C[x,q].     (23)
```

Since `omega_1=dC wedge dA` on `(19)`, equations `(13),(20)` give

```text
(Psi|Gamma_h)^*omega_1
  =(3c-Jac_(x,q)(h,c)) dx wedge dq.                    (24)
```

Therefore a unit coefficient `lambda in C*` would force

```text
Jac_(x,q)(h,c)=3c-lambda.                              (25)
```

Equation `(25)` has no polynomial solution.  Give `x,q` weights `1,-2`, put

```text
u=x^2q,
p(u)=(u+1)(u+3),
c=xp(u).                                               (26)
```

The Hamiltonian derivation `delta=Jac(-,c)` raises weight by two.  For a
homogeneous rational expression `x^mF(u)`, direct differentiation gives

```text
delta(x^mF(u))
 =x^(m+2)(mF p'-F'p).                                  (27)
```

The weight-one part of `(25)` must come from the weight-minus-one component

```text
h_(-1)=x^(-1)F(u),                   F in u C[u].       (28)
```

Equations `(25),(27)` require

```text
-x(Fp)'=3xp,
Fp=-u(u+3)^2+K                                        (29)
```

for a constant `K`.  Divisibility by `u+3` forces `K=0`, but the remaining
numerator has value `4` at `u=-1`; it is not divisible by `p`.  This
contradiction proves the no-unit statement for every polynomial graph.

### 4.1 Collision condition and the affine-linear hostile

Because `Theta` is an isomorphism, a source graph retains the threefold
collision exactly when

```text
h(r_0)=h(r_+)=h(r_-).                                  (30)
```

For

```text
h=alpha x+beta q+gamma,                                (31)
```

condition `(30)` forces `alpha=beta=0`.  Hence the only affine-linear
collision-preserving graphs are constant.  Their coefficient in `(24)` is
`3c`, not a unit.  Arbitrary polynomial graphs can satisfy `(30)`, but `(29)`
excludes all of them from having a unit arm form.

## 5. Complete classification over the punctured Darboux plane

There is a sharp positive result if the graph is allowed to be polynomial in
the **punctured** coordinates `(a,c)`.  Let

```text
w=H(a,c),                              H in C[a,c].     (32)
```

Then `(20)` and `(13)` give

```text
dA wedge dc=(c+partial_a H) da wedge dc,
Jac_(x,q)(A,c)=-3(c+partial_a H).                       (33)
```

The nonzero Jacobian in `(13)` makes `a,c` algebraically independent.  Thus
the planar coefficient is a nonzero constant if and only if

```text
H(a,c)=(mu-c)a+G(c),
mu in C*,                              G in C[c].       (34)
```

For `(34)`,

```text
A=mu a+G(c),                Jac_(x,q)(A,c)=-3mu.       (35)
```

The original collision is automatic because all three points have the same
`(a,c)`.  The exact-image statement of THM-3561 also says that the punctured
source covers the root-zero arm plane; since `(35)` is a triangular plane
automorphism, the full stable image of `(32)` is the polynomial graph

```text
B=cA,
Y=A(cA+4),
S=cA^3/8+3A^2/4+(A-G(c))/mu                            (36)
```

inside `U_0 x A1_S`.  Formula `(36)` follows by substituting `(34)` into
`(6)`; it makes the target slice visibly isomorphic to `A2_(A,c)`.

This does **not** fill the source puncture.  Substitution of `(12)` into
`(35)` gives

```text
A=mu q/D^2+G(xD(D+2)).                                 (37)
```

The first term has a double pole along `D=0`, while the second is regular.
Since `mu!=0`, cancellation is impossible.  An affine-linear polynomial
`H(a,c)` has constant `partial_a H`, so `c+partial_a H` cannot be constant;
the successful family necessarily contains the quadratic term `-ca`.

### 5.1 Positive punctured control

The cheapest member of `(34)` is

```text
mu=1,             G=0,
w=(1-c)a,
A=a,
S=ca^3/8+3a^2/4+a.                                    (38)
```

It preserves the determinant `-3` and the triple collision exactly.  At the
collision its stable value is `w=-3/4` and `S=-21/64`.  It is defined on
`D!=0`, and its double pole is precisely the old THM-3561 puncture.  The
cylinder therefore reorganizes the rational counterexample but does not
polynomialize it.

## 6. The fixed second coordinate has no polynomial mate

One might allow a rational graph to cancel `(21)` and hope that the resulting
arm coordinate `F(x,q)` is polynomial.  Even then the unchanged second
coordinate

```text
c=x(1+x^2q)(3+x^2q)                                   (39)
```

has no polynomial Jacobian mate.  If `F in C[x,q]` satisfied

```text
Jac_(x,q)(F,c)=lambda in C*,                           (40)
```

the weight-zero part of `(40)` would come from

```text
F_(-2)=x^(-2)K(u),                    K in u C[u].      (41)
```

Using `(27)` with `m=-2` gives

```text
pK'+2p'K=-lambda,
(p^2K)'=-lambda p.                                    (42)
```

With

```text
I(u)=integral p(u)du=u^3/3+2u^2+3u,                   (43)
```

equation `(42)` says

```text
p^2K=-lambda I+K_0.                                   (44)
```

At `u=-3`, divisibility forces `K_0=0`; at `u=-1`, the right side is
`4lambda/3`, nonzero.  This proves `(40)` impossible.  It is a theorem about
the particular polynomial `(39)`, not a use of `JC(2)`.

## 7. Sharp boundary and open exits

The result closes exactly the following routes:

1. polynomial source graphs `w=h(x,q)` followed by the collision-bearing arm
   projection `(A,C)`;
2. all polynomial graphs `w=H(a,c)` on the punctured Darboux plane, including
   the complete unit-form family `(34)`;
3. rational graph choices that make `A` polynomial while keeping the second
   component `(39)`.

It does not close:

- an implicitly defined `A2` in the cylinder that is not a polynomial graph
  over `(x,q)` or `(a,c)`;
- a planar projection mixing `S` essentially with both arm coordinates;
- a different cylinder isomorphism followed by a genuinely different
  filling; or
- any construction that changes the fixed second coordinate `c`.

Those possibilities are **OPEN**.  Any future claim must provide an actual
polynomial `A2` source, two polynomial output functions, and their constant
ordinary Jacobian.  Stable three-volume preservation alone is insufficient.

## 8. Literature boundary

Moser-Jauslin and Poloni,
[*Isomorphisms between cylinders over Danielewski
surfaces*](https://arxiv.org/abs/2002.12202), Section 5.1, is cited only for
the Russell cylinder formula from which `(5)` starts.  The direct inverse
`(7),(8)`, specialization checks, volume calculation, THM-3561 transport,
and graph no-filling results are internal proof candidates.  No priority
claim is made for the general cylinder isomorphism.

## 9. Exact companion and reproduction

The deterministic exact companion checks:

- the Bezout and `P(g)=P^2H` identities;
- both polynomial cylinder relations and both inverse compositions modulo
  the two hypersurface ideals;
- the residue-volume determinant and stable THM-3561 pullback;
- the triple collision curve and the unique relevant arm formula;
- polynomial pole hostiles and affine-linear collision constraints;
- the weight formula `(27)`, both ODE contradictions, and one-arm positive
  boundary controls;
- the complete punctured graph family for several exact `mu,G` controls,
  including the target graph `(36)` and unavoidable double pole; and
- the special positive punctured slice `(38)`.

Finite rows verify the displayed algebra; the all-polynomial conclusions are
proved by the grading and divisibility arguments above.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_graph_slice_no_filling_thm3605.py
python3 -O 04-computation/jc2_russell_cylinder_graph_slice_no_filling_thm3605.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_graph_slice_no_filling_thm3605.out`.
