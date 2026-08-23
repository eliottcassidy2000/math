---
id: THM-3738
title: "Opposite-charge radial-profile critical-point obstruction"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.  Over an algebraically
  closed characteristic-zero field, let phi(0)psi(0) be nonzero.  The
  polynomial Q=X phi(XT)+T psi(XT) has no critical point if and only if phi
  and psi are both constant.  Thus every genuinely nonlinear mixture of the
  two opposite radial charges already loses gradient unimodularity and
  cannot be one component of a Keller map.  In particular, every nonzero
  lower/upper binomial-tower superposition from THM-3734 has a critical
  point, including unequal depths.  The proof is uniform in both degrees.
source: root + jc_sparse_direct_search / 2026-08-22
audit: >
  PENDING.  The depth-one resultant was independently derived before the
  all-degree eliminant/gcd proof.  A final hostile audit of the written
  degenerate-root argument, field assumptions, scope boundary, hashes, and
  transcript is still required.
depends_on: []
related:
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
script: 04-computation/jc2_opposite_charge_radial_profile_critical_obstruction_thm3738.py
output: 05-knowledge/results/jc2_opposite_charge_radial_profile_critical_obstruction_thm3738.out
script_sha256: 9513b252d566e54fcf99025c96eb088d3feb17790b326eacf242d38b8c0fee4f
output_sha256: d6e051ae56b1ad914f45719d07738775083f1d60fc0b5099001334f8f5835579
semantic_sha256: 4dc770e42f6183fa8a571016c2b1f29cc428000f2b2d2dc7b4511d1e7d8e89e6
hash_basis: raw LF bytes
---

# THM-3738 -- opposite radial charges cannot be mixed smoothly

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.**  THM-3734 exposes
two large families of nonsingular noncoordinate polynomials carrying charges
`+1` and `-1`.  The most direct counterexample design is to superpose the two
charges, hoping that their Hamiltonian top edges cancel.  They do not reach
that stage: every nonlinear two-charge radial superposition has a critical
point.

Let `k` be an algebraically closed field of characteristic zero.  For
`phi,psi in k[z]` satisfying

```text
phi(0) psi(0) !=0,                                    (1)
```

put

```text
Q(X,T)=X phi(XT)+T psi(XT).                            (2)
```

Then the following are equivalent:

```text
(i)  Q has no critical point in k^2;
(ii) the ideal (Q_X,Q_T) is the unit ideal;
(iii) phi and psi are both constant.                  (3)
```

Consequently, if `J(P,Q)` is a nonzero constant for some polynomial `P`,
then a `Q` of the form `(2)` subject to `(1)` must be linear.
The equivalence of `(i)` and `(ii)` is the weak Nullstellensatz; the proof
below establishes the substantive equivalence with `(iii)`.

## 1. The torus eliminant is one derivative

Write

```text
z=XT,
A(z)=(z phi(z))',                 B(z)=(z psi(z))'.    (4)
```

Direct differentiation gives

```text
Q_X=A(z)+T^2 psi'(z),
Q_T=B(z)+X^2 phi'(z).                                  (5)
```

On the torus `XT=z!=0`, multiply the equations in `(5)` by `X` and `T`.
They become

```text
[ A       z psi'] [X] = [X Q_X],
[ z phi'  B     ] [T]   [T Q_T].                      (6)
```

The determinant of the matrix in `(6)` is

```text
A B-z^2 phi'psi'
 =phi psi+z(phi'psi+phi psi')
 =(z phi psi)'.                                       (7)
```

Thus a critical point is controlled by a one-variable derivative, with the
only possible loss being that a kernel vector of `(6)` lies on an axis.

## 2. A degree/gcd count forces a non-axis kernel

Set

```text
P(z)=phi(z)psi(z),                 C(z)=(zP(z))'.      (8)
```

Suppose `P` is nonconstant of degree `n`, and let `s>=1` be the number of its
distinct roots.  Condition `(1)` gives `gcd(z,P)=1`, so

```text
gcd(C,P)=gcd(P+zP',P)=gcd(P',P).                       (9)
```

The last gcd has degree `n-s`.  Meanwhile `C` has degree `n`, because its
leading coefficient is `(n+1)lead(P)`.  Therefore

```text
deg[C/gcd(C,P)]=s>=1.                                 (10)
```

Over the algebraically closed field, `(10)` has a root `z_0` for which

```text
C(z_0)=0,                         P(z_0)!=0.            (11)
```

Also `z_0!=0`, since `C(0)=P(0)!=0`.

At `z_0`, the matrix in `(6)` is singular by `(7)`.  Any nonzero kernel
vector `(x_0,t_0)` has both coordinates nonzero.  Indeed, an `X`-axis kernel
would force `A(z_0)=phi'(z_0)=0`, hence `phi(z_0)=0`; a `T`-axis kernel would
similarly force `psi(z_0)=0`.  Either contradicts `(11)`.

Choose `lambda in k` with

```text
lambda^2 x_0t_0=z_0,                                  (12)
```

which is possible because `k` is algebraically closed.  Set
`X=lambda x_0`, `T=lambda t_0`.  Then `XT=z_0`, the scaled vector remains in
the kernel of `(6)`, and its two coordinates are nonzero.  Equations `(6)`
therefore give

```text
Q_X(X,T)=Q_T(X,T)=0.                                  (13)
```

So every nonconstant product `phi psi` produces a critical point.  Since
`k[z]` is a domain, a constant nonzero product forces both factors to be
constant.  Conversely, nonzero constant `phi,psi` make `(2)` a linear
polynomial with constant nonzero gradient.  This proves `(3)`.  **QED.**

## 3. Consequences for the Cohn divided-power towers

For THM-3734, write the lower and upper normalized potentials as

```text
q_r^L=X Phi_r^L(XT),             Phi_r^L(0)=1,
q_s^U=T Phi_s^U(XT),             Phi_s^U(0)=-1.        (14)
```

For arbitrary depths `r,s>=1` and nonzero scalars `lambda,mu`, the mixture

```text
lambda q_r^L+mu q_s^U                                  (15)
```

satisfies `(1)` and is nonlinear.  It therefore has a critical point.  The
depths need not match.  At `r=s=1`, `(15)` includes the first attempted
mixture

```text
X+X^2T + lambda(XT^2-T),                               (16)
```

whose explicit resultant obstruction is now explained by `(7)--(10)` rather
than by a special quartic computation.

This isolates a stricter design requirement for a planar Jacobian
counterexample.  Merely coupling the two opposite charge-one towers cannot
preserve the Bezout gradient.  A viable mixed-charge component must either
use additional charge sectors, cross the excluded axis boundary
`phi(0)psi(0)=0`, or abandon the radial `XT` profile.

Reproduce the exact audit surface with

```bash
python3 -B 04-computation/jc2_opposite_charge_radial_profile_critical_obstruction_thm3738.py
python3 -B -O 04-computation/jc2_opposite_charge_radial_profile_critical_obstruction_thm3738.py
```

The companion verifies the universal determinant identity, hostile repeated
and shared-root profiles, all 2,496 normalized integer profile pairs of
degree at most two, 16 linear equality controls, and all 36 unequal-depth
binomial pairs through depth six.  Its finite checks are controls; equations
`(4)--(13)` prove the all-degree statement.
