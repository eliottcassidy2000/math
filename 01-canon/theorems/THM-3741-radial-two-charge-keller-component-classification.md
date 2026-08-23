---
id: THM-3741
title: "Radial two-charge Keller-component classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over an algebraically
  closed characteristic-zero field, every nonsingular polynomial
  Q=X phi(XT)+T psi(XT) is exactly a nonzero linear form, a squarefree
  single-charge profile nonvanishing at the origin, or one of the two mixed
  monomial flanks aX+bX^nT^(n+1), aX^(n+1)T^n+bT with n>=2.  Among this
  complete nonsingular list, Q has a polynomial Jacobian mate if and only if
  Q is linear.  Thus the entire arbitrary-degree two-charge radial ansatz is
  empty as a planar Jacobian-counterexample search space.
source: root + jc_sparse_direct_search / 2026-08-22
audit: >
  PASS.  An independent hostile audit rederived the mixed-axis order boundary,
  the root-multiplicity quotient and off-profile torus root, the pure-charge
  squarefree criterion, both charge projections, the transfer of both mixed
  flanks to THM-3716, and the explicit linear mates.  Normal, optimized, and
  frozen output agree; script/output/semantic hashes and CHECKS=1126 match.
depends_on:
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3738-opposite-charge-radial-profile-critical-point-obstruction
related:
  - THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers
script: 04-computation/jc2_radial_two_charge_keller_component_classification_thm3741.py
output: 05-knowledge/results/jc2_radial_two_charge_keller_component_classification_thm3741.out
script_sha256: dd01fb8c51e89d0a119d5f8f7bd3a333e03ea7b6c075c4eb6314516660f86127
output_sha256: 8594c7c43da5f78c71d21a272227b22f29df07a1622c34c132f272d8a5c4c1c1
semantic_sha256: 25482d6fad2f12deffc6f12c41449ab0804d10a73cd2088b7a8cd029414f0077
hash_basis: raw LF bytes
---

# THM-3741 -- the full radial two-charge ansatz has no Keller pair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3738 killed
the generic boundary where both radial profiles are nonzero at the origin.
The excluded axes are not merely technical exceptions: they contain infinite
families of nonsingular noncoordinates.  This theorem classifies all of them
and then shows that none has a polynomial Jacobian mate.

Let `k` be an algebraically closed field of characteristic zero, let
`phi,psi in k[z]`, and put

```text
Q(X,T)=X phi(XT)+T psi(XT).                            (1)
```

Then `Q` has no critical point in `k^2` if and only if it belongs to one of
the following classes.  The list is made disjoint only to expose its
boundaries.

```text
(L) phi=a, psi=b are constant, with (a,b)!=(0,0);

(P+) psi=0 and phi is nonconstant, squarefree, and phi(0)!=0;
(P-) phi=0 and psi is nonconstant, squarefree, and psi(0)!=0;

(M+) phi=a, psi=bz^n, with a,b!=0 and n>=2;
(M-) phi=az^n, psi=b, with a,b!=0 and n>=2.            (2)
```

Moreover,

```text
there is P in k[X,T] with J(P,Q) in k*
                 if and only if Q is in the linear class (L).             (3)
```

Thus `(1)` supplies smooth noncoordinate components but no planar Keller
pair and hence no Jacobian-conjecture counterexample.

## 1. Universal derivatives and the torus eliminant

Write `z=XT` and

```text
A(z)=(z phi(z))',                 B(z)=(z psi(z))'.    (4)
```

Direct differentiation gives

```text
Q_X=A(z)+T^2 psi'(z),
Q_T=B(z)+X^2 phi'(z).                                  (5)
```

On the torus, multiplication by `X,T` packages the critical equations as

```text
[ A       z psi'] [X] = [X Q_X],
[ z phi'  B     ] [T]   [T Q_T].                      (6)
```

The determinant is the one-variable derivative

```text
A B-z^2phi'psi'=(z phi psi)'.                          (7)
```

Whenever `(7)` vanishes at `z_0!=0` while
`phi(z_0)psi(z_0)!=0`, the matrix has a kernel vector with both coordinates
nonzero.  Scaling that vector by a common `lambda` satisfying
`lambda^2XT=z_0` produces a point on the fibre `XT=z_0`; equations `(6)`
then give a critical point.  This is the torus-kernel mechanism of THM-3738,
and it remains valid even when one profile vanishes at the origin.

## 2. Classification at the two axes

At the origin the gradient is

```text
(Q_X,Q_T)(0,0)=(phi(0),psi(0)),                        (8)
```

so at least one constant term must be nonzero.  If both are nonzero,
THM-3738 says that nonsingularity is equivalent to both profiles being
constant.  This gives `(L)` with `ab!=0`.

It remains, up to the typed duality `(X,phi)<->(T,psi)`, to suppose

```text
phi(0)!=0,                         psi(0)=0.            (9)
```

### The pure-charge boundary

If `psi=0`, equations `(5)` become

```text
Q_X=(z phi)',                     Q_T=X^2phi'.         (10)
```

There is no critical point on `X=0`, because `Q_X=phi(0)`.  For `X!=0`, a
critical point exists exactly when `phi(z)=phi'(z)=0`: every such `z` is
realized by taking `T=z/X`.  Hence

```text
Q=Xphi(XT) is nonsingular
 iff phi(0)!=0 and gcd(phi,phi')=1.                    (11)
```

The constant case belongs to `(L)` and the nonconstant cases are precisely
`(P+)`.  Swapping `X,T` proves `(P-)`.

### The mixed boundary

Now assume `psi!=0`.  On the whole axis `X=0`, equations `(5)` read

```text
Q_T=0,                            Q_X=phi(0)+T^2psi'(0). (12)
```

If `psi'(0)!=0`, algebraic closure supplies a root of the second polynomial,
so nonsingularity forces

```text
k=ord_0(psi)>=2.                                      (13)
```

Factor

```text
phi(z)psi(z)=z^k R(z),               R(0)!=0.          (14)
```

If `R` is nonconstant, then

```text
(z phi psi)'=z^k D(z),
D(z)=(k+1)R(z)+zR'(z).                                (15)
```

For a root `alpha` of `R` of multiplicity `m`, necessarily `alpha!=0`, and
`ord_alpha(D)=m-1`.  Equivalently,

```text
gcd(D,R)=gcd(R',R).                                   (16)
```

If `R` has `s>=1` distinct roots, `D` has degree `deg R` and
`D/gcd(D,R)` has degree `s`.  None of its roots lies on `R`, and zero is not
a root because `D(0)=(k+1)R(0)!=0`.  It therefore supplies a
`z_0!=0` with

```text
(z phi psi)'(z_0)=0,                phi(z_0)psi(z_0)!=0. (17)
```

The torus-kernel argument following `(7)` gives a critical point.  Thus
`R` must be constant.  Unique factorization in `k[z]`, together with
`phi(0)!=0`, then forces

```text
phi=a in k*,                       psi=bz^k, b in k*.  (18)
```

Conversely the associated polynomial is

```text
Q=aX+bX^kT^(k+1).                                      (19)
```

Its second derivative in `(5)` is `(k+1)b(XT)^k`.  At any possible critical
point `XT=0`; because `k>=2`, the first derivative there is `a`, not zero.
So `(19)` is nonsingular.  The excluded equality case `k=1` is genuinely
singular: on `X=0` choose `T` with `a+bT^2=0`.  This proves `(M+)`, including
its sharp exponent boundary.  The other axis gives `(M-)`.

Together with `(8)--(9)` these cases exhaust all profile pairs and prove the
nonsingular classification `(2)`.

## 3. A charge projection closes the pure families

It remains to classify Jacobian mates.  Give `X,T` Euler weights `+1,-1`.
The pure component

```text
Q=Xphi(XT)                                               (20)
```

has weight `+1`.  If `J(P,Q)=lambda in k*`, only the weight `-1` part of
`P` can contribute to the weight-zero right side.  Every polynomial of
weight `-1` is

```text
P_(-1)=T F(XT),                    F in k[z].           (21)
```

The projected Jacobian equation is the exact one-variable antiderivative

```text
J(TF(XT),Xphi(XT))=-(zF(z)phi(z))'=lambda.             (22)
```

Integration and evaluation at zero give `Fphi=-lambda`.  Since the units of
`k[z]` are the nonzero constants, `phi` must be constant.  Thus every
nonlinear `(P+)` member has nonzero Hamiltonian debt despite its unimodular
gradient.  Typed duality gives

```text
J(XF(XT),Tpsi(XT))=(zF(z)psi(z))',                     (23)
```

and kills every nonlinear `(P-)` member in the same way.

## 4. The mixed flanks are monomial Broughton components

After nonzero source and target scalings, `(M+)` is

```text
x+x^n y^(n+1),                    n>=2,                (24)
```

and `(M-)` is its variable-swapped dual.  Invertible linear source changes
multiply a Jacobian by a nonzero scalar, so the existence of a constant
Jacobian mate is invariant under these normalizations.  THM-3716 applies to
`(24)` and proves that its Hamiltonian repair chain never terminates.
Therefore neither mixed family has a polynomial mate.

Finally, for a nonzero linear form `Q=aX+bT`, explicit normalized mates are

```text
P=-T/a if a!=0,                    P=X/b if a=0.        (25)
```

They satisfy `J(P,Q)=1`, proving both directions of `(3)`.  **QED.**

## 5. Exact controls and the next search boundary

Reproduce the exact audit surface with

```bash
python3 -B 04-computation/jc2_radial_two_charge_keller_component_classification_thm3741.py
python3 -B -O 04-computation/jc2_radial_two_charge_keller_component_classification_thm3741.py
```

The assertion-free companion verifies `(7)`, both charge projections, the
complete degree-at-most-two coefficient cube over `{-1,0,1}`, squarefree and
repeated-root pure controls, both mixed orientations through `n=8`, the
singular `n=1` boundary, and bounded mate-system obstructions through total
degree seven.  These are hostile controls; equations `(4)--(25)` prove the
arbitrary-degree theorem.

The structural lesson is sharper than the earlier generic obstruction.  A
counterexample component cannot be obtained by any arbitrary polynomial
choice of the two radial profiles, even on their exceptional axes.  The next
honest search must add a third Euler-charge sector, introduce genuinely
nonradial coupling, or move from component design to interacting
factorizations where no single radial potential is exposed.
