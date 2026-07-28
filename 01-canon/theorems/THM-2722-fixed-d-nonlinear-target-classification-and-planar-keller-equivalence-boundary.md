---
id: THM-2722
title: "Fixed-d nonlinear target classification and planar Keller equivalence boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  On every polynomial graph, a nonzero-constant-Jacobian target
  pair (P(A,d),d) forces d=ax+R(y-x^2/2) and
  P=alpha A+K(d), and is a triangular polynomial automorphism.  For an
  unrestricted pair (P(A,d),Q(A,d)), the same unit argument first makes
  (A,d) polynomial coordinates and then identifies the surviving problem
  exactly with an ordinary planar Keller pair.  Thus the fixed-d face is
  closed, while the full C[A,d] face is equivalent to JC(2), not a proof of
  it.  Targets involving B, nongraph source surfaces, JC(2), and DC(2)
  remain open.
source: thm2694-full-lift-fibre-scout-2026-07-28
depends_on:
  - THM-2702-polynomial-graph-coordinate-projection-keller-classification
related:
  - THM-2705-linear-target-planes-containing-d-polynomial-graph-keller-classification
  - THM-2715-nonlinear-d-target-shear-polynomial-graph-keller-classification
script: 04-computation/jacobian_fixed_d_nonlinear_target_boundary_thm2722.py
output: 05-knowledge/results/jacobian_fixed_d_nonlinear_target_boundary_thm2722.out
script_sha256: a2eac3f6d40022c47d8a2ba8286114bf1693f53ae5367c4591398a5d2045adcd
output_sha256: c71107ac8778b3c8ec52eacf9d4a19909d9d3320d580eda58bd8cfda3c5ee0d9
hash_basis: LF-normalized bytes
---

# THM-2722 -- fixed-`d` nonlinear targets and the exact planar boundary

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2705 classifies linear target planes containing `d`.  Allowing an
arbitrary polynomial in `A,d` looks much larger, but the Jacobian factors
through the coordinate pair `(A,d)`.  A nonzero constant makes both factors
units.  This closes the whole fixed-`d` nonlinear face and also identifies
the precise stopping boundary: if both outputs vary in `C[A,d]`, the residue
is exactly the planar Jacobian conjecture in renamed coordinates.

## 1. Statement

Work over `C` on a polynomial graph

```text
z=f(x,y),                  A=x^2-2y,                  d=z.       (1)
```

Put

```text
t=y-x^2/2,                 F(x,t)=f(x,t+x^2/2).                  (2)
```

Let `P,Q in C[A,d]` and define

```text
D_(P,Q)=P_A Q_d-P_d Q_A.                              (3)
```

Then the following hold.

1. The exact graph factorization is

   ```text
   Jac(P(A,d),Q(A,d))
     =2 D_(P,Q)(-2t,F) F_x.                           (4)
   ```

2. If the left side of `(4)` is `kappa in C^*`, there are

   ```text
   a,delta in C^*,             R in C[t],             (5)
   ```

   such that

   ```text
   F=ax+R(t),                  D_(P,Q)=delta,
   kappa=2a delta.                                      (6)
   ```

   The map

   ```text
   Theta:(x,t) |-> (A,d)=(-2t,ax+R(t))                 (7)
   ```

   is a polynomial automorphism.  Moreover

   ```text
   (P(A,d),Q(A,d))=(P,Q) o Theta,                      (8)
   ```

   so the graph target is a polynomial automorphism if and only if the
   planar Keller pair `(P,Q)` in the independent variables `(A,d)` is.

3. In the fixed-`d` case `Q=d`, equation `(4)` has a nonzero constant value
   if and only if

   ```text
   F=ax+R(t),             P=alpha A+K(d),              (9)
   a,alpha in C^*,        R,K in C[T].
   ```

   Every such target is a polynomial automorphism, with Jacobian
   `kappa=2a alpha`.

Conversely, every planar Keller pair `(P,Q)` occurs in the full `C[A,d]`
graph family: take `F=x`, so `(A,d)=(-2t,x)` is already a polynomial
coordinate system.  Therefore a nonautomorphic member of the unrestricted
family exists if and only if `JC(2)` is false.  This equivalence is a scope
boundary, not a solution of `JC(2)`.

## 2. The chain-rule factorization

The change `(x,y)->(x,t)` has Jacobian one, and `(2)` gives

```text
A=-2t,                       Jac(A,d)=2F_x.            (10)
```

For any two polynomials in `A,d`, the two-by-two chain rule gives

```text
Jac(P(A,d),Q(A,d))
 =det [P_A P_d; Q_A Q_d] Jac(A,d),                    (11)
```

which is `(4)`.  No algebraic independence has been used yet.

Suppose `(4)` equals `kappa!=0`.  The two factors

```text
D_(P,Q)(-2t,F),                 F_x                   (12)
```

belong to `C[x,t]` and have a nonzero constant product.  Both are therefore
units.  Write `F_x=a in C^*`; integration in `x` gives

```text
F=ax+R(t).                                             (13)
```

Now `(7)` has the explicit inverse

```text
t=-A/2,                  x=(d-R(-A/2))/a.             (14)
```

Thus substitution `(A,d)=(-2t,F)` is an isomorphism

```text
C[A,d]  ->  C[x,t].                                   (15)
```

The first factor in `(12)` is a scalar `delta`; injectivity of `(15)` now
forces the polynomial identity `D_(P,Q)=delta`.  This proves `(5)`--`(8)`.

The argument also explains why one must not replace the unit conclusion by
a statement only about values along a dependent curve.  The nonzero
Jacobian first forces `(A,d)` to be independent polynomial coordinates;
only then may the constant pullback be descended to `C[A,d]`.

## 3. Complete fixed-`d` classification

Set `Q=d`.  Then `D_(P,d)=P_A`, so Section 2 gives

```text
F=ax+R(t),                  P_A=alpha in C^*.          (16)
```

Integrating in the independent variable `A` yields exactly

```text
P(A,d)=alpha A+K(d).                                  (17)
```

Conversely `(13)` and `(17)` give

```text
U=P(A,d)=-2alpha t+K(ax+R(t)),
V=d=ax+R(t),                                            (18)
```

and direct differentiation gives `Jac(U,V)=2a alpha`.
The polynomial inverse is

```text
t=(K(V)-U)/(2alpha),
x=(V-R(t))/a,
y=t+x^2/2,                                             (19)
```

with `z=V` on the graph.  This proves necessity, sufficiency, uniqueness of
the normal form, and triangularity.  An affine rescaling or translation of
the fixed output `d` is absorbed by a target affine automorphism and gives
the same classification.

The linear choices of `K` recover the relevant branch of THM-2705.  The new
content is that no nonlinear `K(d)` and no higher mixed monomial in `P`
creates another graph: all nonlinear freedom is a harmless triangular
target shear.

## 4. The exact `JC(2)` boundary

For unrestricted `P,Q`, Section 2 proves that every nonzero-response graph
first straightens to the coordinate pair `(A,d)`.  After that
straightening, `(8)` contains precisely an ordinary planar Keller pair.
Both directions are literal:

```text
graph pair nonautomorphic
  <=> planar Keller pair (P,Q) nonautomorphic.          (20)
```

For the reverse implication one may always choose the graph `F=x`; then
`Theta(x,t)=(-2t,x)` is an automorphism and multiplies the planar Jacobian
only by `2`.  Thus attempting to classify all of `C[A,d]^2` without another
invariant simply restates `JC(2)`.

This is also the failure boundary of the fixed-`d` proof.  Fixing one output
to `d` makes `D_(P,d)=P_A`, whose constant-Jacobian equation integrates.
When both outputs vary, `D_(P,Q)=delta` is the general planar Keller
equation and has no known triangular normal form.

## 5. Scope

The theorem concerns polynomial graphs in the fixed THM-2696 quotient and
targets whose two components lie in the subalgebra `C[A,d]`.  It completely
closes the nonlinear face with one output fixed to `d` (or an affine
coordinate in `d`).

It does not classify targets involving `B`, nonlinear changes mixing all
three quotient coordinates, source surfaces that are not graphs in this
chart, or arbitrary planar or Weyl-algebra endomorphisms.  The equivalence
`(20)` proves neither `JC(2)` nor `DC(2)`; it identifies the exact point at
which this quotient family ceases to simplify them.

At zero response the conclusion is false.  For example `F=R(t)` makes
`Jac(A,d)=0`, while a polynomial `P` independent of `A` makes
`D_(P,d)=0`.  These are the sharp unit-gate failures.

## 6. Exact companion

Run

```text
python 04-computation/jacobian_fixed_d_nonlinear_target_boundary_thm2722.py
python -O 04-computation/jacobian_fixed_d_nonlinear_target_boundary_thm2722.py
```

and compare both transcripts with

```text
05-knowledge/results/jacobian_fixed_d_nonlinear_target_boundary_thm2722.out.
```

The companion uses exact symbolic arithmetic and explicit failure raises.
It checks the formal two-by-two chain rule, the graph factor, a generic
nonlinear fixed-`d` normal form and its two-sided inverse, a nontrivial
planar Keller composition, a variable-Jacobian hostile, and both zero-response
boundaries.  The unit and all-degree integration arguments are Sections 2--3,
not finite computations.

QED (candidate; independent hostile audit pending).
