---
id: THM-3563
title: "Nongraph conductor node-cycle obstruction to every target projection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On THM-2696's
  sharp constant-different plane
  s_3=s_1s_2-c, c nonzero, no pair of polynomial functions of all three
  target coordinates (A,B,d) pulls back to a planar map with nonzero
  constant Jacobian.  An explicit Laurent double curve makes every target
  function even, while polynomial Poincare exactness forces a nonzero odd
  inverse-cube difference.  Six normalization points form three source
  nodes whose cycle makes that difference equal to both r and -2r.  The
  theorem is confined to this one normalization plane and proves no JC(2)
  conclusion.
source: root-2026-08-18-planar-jacobian-counterexample-constructions
audit: >
  An independent hostile audit rederived the restricted target map,
  polynomial primitive and sign, even Laurent packet, all three distinct
  source nodes, alternating inverse-cube values, and the node-cycle
  contradiction.  It also checked the nonclosed minor-Bezout hostile,
  theorem scope, algebraic universe, active truth gates, byte-identical
  ordinary and optimized replays, hashes, documentation, and diff surface.
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
related:
  - THM-2722-fixed-d-nonlinear-target-classification-and-planar-keller-equivalence-boundary
  - THM-3554-punctured-kummer-collision-surface-normal-form
companion: 04-computation/jc_nongraph_conductor_node_cycle_thm3563.py
output: 05-knowledge/results/jc_nongraph_conductor_node_cycle_thm3563.out
script_sha256: 815839c28279a7dd749b465959deb3c2adbfc66726e7bd44d067840d71b8d5f2
output_sha256: 7410ac021b72b92be30f8d379e4242ec162ab35f163cf4d6952f4c652d4c5d32
hash_basis: LF-normalized bytes
---

# THM-3563 -- nongraph conductor node-cycle obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
closes every polynomial target projection of one fixed normalization plane.
The node identifications, primitive sign, algebraic universe, and relationship
to the older normalization packet have all survived independent hostile
audit.

The field has characteristic zero.  It is enough to work after base change
to an algebraic closure.  Fix `c!=0`.

## 1. Statement and the sharp constant-different plane

Restrict THM-2696's quotient

```text
(s_1,s_2,s_3) ->
(A,B,d)=(s_1^2-2s_2, s_2^2-2s_1s_3, s_3)             (1)
```

to its sharp level `s_1s_2-s_3=c`.  With source-plane coordinates
`(u,v)=(s_1,s_2)`, this is

```text
nu_c:A^2_(u,v) -> A^3_(A,B,d),

A=u^2-2v,
B=v^2-2u^2v+2cu,
d=uv-c.                                               (2)
```

The claim is that for every pair

```text
P,Q in k[A,B,d],                                      (3)
```

the pullback pair `p=P o nu_c`, `q=Q o nu_c` never has

```text
Jac_(u,v)(p,q)=kappa in k*.                            (4)
```

Equivalently, no polynomial projection `A^3_(A,B,d)->A^2`--linear or
nonlinear--turns this fixed normalization of THM-2696's singular nongraph
target surface into a planar Keller map.

## 2. The Laurent double curve

For `h!=0`, put

```text
gamma(h)=(u(h),v(h)),

u(h)=h/2+4c/h^2,
v(h)=2c/h-16c^2/h^4.                                  (5)
```

Direct substitution gives

```text
A(gamma(h))=(h^6+192c^2)/(4h^4),
B(gamma(h))=4c^2(h^6+192c^2)/h^8,
d(gamma(h))=-64c^3/h^6.                               (6)
```

All three expressions are even in `h`.  Consequently

```text
nu_c(gamma(h))=nu_c(gamma(-h)),                        (7)
```

and every target polynomial in `(3)` becomes even on `gamma`.

The source curve in `(5)` lies on

```text
K=8c^2-cu^3-10cuv+u^4v+2u^2v^2+v^3=0,                (8)
```

where the exact target identity is

```text
K=((c-3d)^2-AB)/2.                                    (9)
```

This is the explicit double-curve packet inside the normalization plane.
The proof below needs only the contained Laurent curve and its six special
points; identifying the whole scheme-theoretic conductor is not being
silently substituted for those exact identities.

## 3. A hypothetical Keller pair forces an odd inverse cube

Assume `(4)`.  In polynomial differential forms on `A^2`,

```text
d(p dq-kappa u dv)=dp wedge dq-kappa du wedge dv=0.    (10)
```

Polynomial Poincare exactness in characteristic zero gives a polynomial
`phi in kbar[u,v]` such that

```text
p dq=kappa u dv+dphi.                                  (11)
```

No analytic or formal primitive is being used.  If the closed form in
`(10)` is `a du+b dv`, one explicit polynomial primitive is

```text
integral_(0)^u a(t,v)dt + integral_(0)^v b(0,s)ds;     (12)
```

closedness makes its two derivatives `a,b`, and characteristic zero permits
the coefficientwise integrations.

Put `f(h)=phi(gamma(h))`.  By `(6)`, both `p(gamma(h))` and `q(gamma(h))`
are even, so

```text
p(gamma(h)) (q o gamma)'(h)                            (13)
```

is odd.  On the other hand, exact differentiation of `(5)` gives

```text
u(h)v'(h)
 =256c^3/h^7+24c^2/h^4-c/h.                           (14)
```

Its even part is exactly `24c^2/h^4`.  Pulling `(11)` to `gamma` and taking
even parts therefore gives

```text
(f_odd)'(h)=-24 kappa c^2/h^4,
f_odd(h)=8 kappa c^2/h^3,                              (15)
```

where `f_odd=(f(h)-f(-h))/2`; its additive constant is zero because it is
odd.  Thus every hypothetical pair `(3)--(4)` forces

```text
f(h)-f(-h)=16 kappa c^2/h^3.                           (16)
```

The sign in `(16)` is fixed twice: differentiating `+8 kappa c^2 h^-3`
gives the negative term in `(15)`, and `(11)` places `kappa u dv` on the
same side as `dphi`.

## 4. Six normalization points and the three-node contradiction

Choose

```text
epsilon^2=-3,
zeta=(1+epsilon)/2,             zeta^3=-1,
h_0^3=8 epsilon c,
h_j=h_0 zeta^j,                 0<=j<=5.               (17)
```

Then every `h_j` satisfies `h_j^6=-192c^2`.  Exact substitution in `(5)`
gives three identifications of **source points**:

```text
gamma(h_0)=gamma(h_5),
gamma(h_1)=gamma(h_2),
gamma(h_3)=gamma(h_4).                                 (18)
```

These are ordinary nodes of the source curve `(8)`.  Indeed both first
derivatives of `K` on `(5)` contain the factor `h^6+192c^2`, while the
Hessian determinant at those six parameter values is `108c^2!=0`.  The
three node images are distinct: they are exactly

```text
(u,u^2/2),                  u^3=8c/3.                  (18a)
```

The cubic is squarefree because `c!=0`.

Write

```text
f_0=f_5=a,             f_1=f_2=b,             f_3=f_4=e. (19)
```

Since `h_j^-3=(-1)^j h_0^-3`, equation `(16)` for `j=0,1,2` gives

```text
a-e=r,
b-e=-r,
b-a=r,

r=16 kappa c^2 h_0^-3 !=0.                            (20)
```

But the first two rows of `(20)` give `b-a=-2r`, whereas the third gives
`b-a=r`.  Hence `3r=0`, contradicting characteristic zero and `r!=0`.
This proves the claimed all-degree no-go.

## 5. Hostile control: a unit minor certificate is not integrable

The three natural tangent minors are

```text
{A,B}=-4(u^3+uv-c),
{A,d}= 2(u^2+v),
{B,d}=-2(v^2+u^2v-cu).                                 (21)
```

They admit the deceptively positive Bezout identity

```text
(-d/(4c^2)){A,B}
 +(B/(3c^2)){A,d}
 -(A/(6c^2)){B,d}=1.                                   (22)
```

Thus arbitrary source-polynomial minor coefficients are not the missing
gate.  The ambient two-form suggested by `(22)` is not closed: its closure
coefficient is

```text
partial_d(-d/(4c^2))
 -partial_B(B/(3c^2))
 +partial_A(-A/(6c^2))=-3/(4c^2).                     (23)
```

For an actual pair `P,Q`, the coefficient two-form is `dP wedge dQ` and is
closed automatically.  Equations `(15)--(20)` identify the global failure
behind `(23)`: closing the form would require the missing inverse-cube
normalization sidecar, and the three-node cycle prevents it from being the
restriction of a source polynomial.

## 6. Scope and the next construction escape

THM-2722 closes `k[A,d]^2` on arbitrary polynomial graphs only up to its
exact planar-Jacobian boundary and explicitly leaves targets involving `B`
open.  The present theorem is orthogonal:

- it allows arbitrary joint nonlinear dependence on all of `(A,B,d)`;
- it treats only the single sharp graph `s_3=s_1s_2-c`; and
- it uses the singular nongraph target surface and its normalization cycle.

It does **not** close `B`-targets on other graphs, nongraph source surfaces,
other constant-different carriers, arbitrary three-dimensional Keller maps,
or `JC(2)`.

The exact escape is also visible.  On the double curve, `h^-3` is a square
root of `-d/(64c^3)`.  A new architecture must deform the three-node cycle
or retain this Kummer/normalization owner until its finite branch has been
exported to infinity.  Simply projecting away the owner recreates the
obstruction `(20)`; simply adjoining it risks the puncture/filling boundary
of THM-3554.  This is a typed construction target, not a counterexample.

## 7. Exact companion

The companion checks `(5)--(9)`, `(14)--(16)`, the derivative and Hessian
node factors, all six source-coordinate reductions in `(18)`, the
alternating inverse cubes, the exact three-row contradiction, and the
nonclosed Bezout hostile `(21)--(23)`.  It uses explicit failure raises, not
Python `assert`.

Reproduce with

```bash
python3 04-computation/jc_nongraph_conductor_node_cycle_thm3563.py
python3 -O 04-computation/jc_nongraph_conductor_node_cycle_thm3563.py
```

Both modes must byte-match

```text
05-knowledge/results/jc_nongraph_conductor_node_cycle_thm3563.out.
```

```text
script_sha256 = 815839c28279a7dd749b465959deb3c2adbfc66726e7bd44d067840d71b8d5f2
output_sha256 = 7410ac021b72b92be30f8d379e4242ec162ab35f163cf4d6952f4c652d4c5d32
hash_basis    = LF-normalized bytes
```
