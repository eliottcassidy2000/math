---
id: THM-2570
title: "Global plane normalization, cusp-cylinder conductor, and survivor section of the sporadic Jelonek surface"
status: >
  PROVED over C, with all ring identities over Q.  The global finite
  birational and subintegral normalization by A^2, the normalization-defect
  module, the nonzero-c cusp-cylinder isomorphism, exact conductor and singular
  support, smooth c=0 boundary, normalized finite-survivor section, and
  theta=0,2,-1 strata are proved.  The companion verifies the parametrization
  kernel, all identities, both boundary sections, and the Jacobian ideal
  exactly in ordinary and optimized Python.  This is a theorem only about the
  fixed sporadic three-variable Keller map; no general Keller-tower or JC(2)
  consequence is claimed.
source: codex-2026-07-27
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2566-two-chart-saturated-cusp-atlas-and-parasitic-plane-ledger
related:
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - MISTAKE-287 (raw cusp/discriminant pullbacks require coefficient sidecars)
script: 04-computation/keller_jelonek_normalization_thm2570.py
output: 05-knowledge/results/keller_jelonek_normalization_thm2570.out
script_sha256: be7374e1a932e6725723ca4e3676488644a3c8023d718877897ba99e56d42e8d
output_sha256: 7e71c67c070c8cca59e87ed8f9b0006c43c2a5cc4c797f3e08f8bcf6375fe82a
---

# THM-2570 -- the sporadic Jelonek surface is a normalized cusp cylinder

Let

```text
L = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
A = Q[a,b,c]/(L).
```

THM-2473 identifies `Spec A=V(L)` as the Jelonek non-properness surface of
the fixed sporadic Keller map `F`.  THM-2566 shows that it is obtained from
two cusp pullbacks only after saturation.  The present theorem identifies the
surface itself: its normalization is an affine plane, and away from `c=0` it
is literally a cylinder over the ordinary cusp.

## 1. Global normalization by an affine plane

Put `B=Q[c,lambda]` and define

```text
nu : A -> B,

a |-> lambda^2(3-c lambda)/27,
b |-> lambda(4-c lambda)/3,
c |-> c.                                             (1)
```

Substitution kills `L`.  Conversely, eliminating `lambda` from the two graph
equations in (1) gives exactly the principal ideal `(L)`; the exact lexicographic
Groebner basis has the sole `lambda`-free member `L`.  Thus (1) is injective.

Two further identities explain the map without relying on elimination:

```text
lambda^2 - 3b lambda + 27a = 0,                     (2)
(3bc-4)lambda = 27ac - 3b.                          (3)
```

Equation (2) is monic, so `B=A[lambda]` is finite over `A`.  Equation (3)
puts `lambda` in `Frac(A)` because `3bc-4` is not the zero element of the
domain `A`; hence `A` and `B` have the same fraction field.  Finally `B` is a
polynomial ring and therefore integrally closed.  It follows that

```text
normalization(A) = Q[c,lambda].                     (4)
```

Section 2 strengthens this: the finite normalization is an elementary
subintegral extension, hence a universal homeomorphism.  In particular it is
set-theoretically bijective over `C`.  Its failure to be an isomorphism is a
cusp-style ring defect, not multiple geometric branches.

## 2. The nonzero-`c` chart is exactly a cusp cylinder

Define

```text
theta = 2 - c lambda,
T = 4 - 3bc,
S = 27ac^2 - 9bc + 8.                               (5)
```

Under (1), direct simplification gives

```text
T = theta^2,                 S = theta^3.            (6)
```

When `c` is invertible, the converse formulas are

```text
b = (4-theta^2)/(3c),
a = (theta^3-3theta^2+4)/(27c^2)
  = (theta-2)^2(theta+1)/(27c^2).                   (7)
```

Therefore the coordinate rings are exactly

```text
A[c^-1] = Q[c,c^-1,theta^2,theta^3],
B[c^-1] = Q[c,c^-1,theta].                          (8)
```

In geometric language,

```text
V(L) intersect D(c)  is  G_m x Spec Q[theta^2,theta^3],
its normalization    is  G_m x A^1_theta.           (9)
```

Moreover the `c`-chart map of THM-2566 becomes precisely the cusp
normalization:

```text
Phi_c o nu(c,theta) = (S,T) = (theta^3,theta^2).    (10)
```

There is also a global statement hidden in (6).  No inversion of `c` is
needed for the recovery identity

```text
lambda = [3(b-9ac)+3b(2-theta)]/4.                 (10a)
```

Indeed `T lambda=3(b-9ac)`, `c lambda=2-theta`, and
`4=T+3bc`.  Consequently

```text
B=A[theta],          theta^2=T in A, theta^3=S in A.  (10b)
```

Thus `A subset B` is an elementary subintegral extension.  Directly, over
any prime of `A`, if the residue of `theta^2` is nonzero then the only possible
residue of `theta` is `theta^3/theta^2`; if it vanishes then `theta` must
vanish in the residue field.  Hence there is one prime above each prime and
no residue-field extension.  Since `B` is already normal, the seminormalization
and normalization of `A` coincide with `B`.  This proves the claimed universal
homeomorphism and identifies the defect as genuinely cuspidal at every base
change.

Thus “cusp-shaped” is not merely a discriminant analogy on this chart; the
Jelonek surface is a product with the cusp ring.  Formula (7) also explains
the two components of its `a=0` section and their multiplicities.

## 3. The exact conductor

For `R_0=Q[c,c^-1]`, the elementary semigroup identity

```text
R_0[theta^2,theta^3] = R_0 + theta^2 R_0[theta]     (11)
```

shows that the conductor of the localized cusp ring inside its normalization
is `theta^2 R_0[theta]`.  Indeed every multiple of `theta^2` carries the whole
normalization into (11), while a nonzero constant term would produce a
forbidden `theta^1` term after multiplication by `theta`.

This conductor globalizes without an extra boundary component.  In `B`, put

```text
R = b - 9ac.
```

The exact identities

```text
T = (2-c lambda)^2,
T lambda = 3R,
S = 2T - 3cR                                      (12)
```

show that `TB` is an ideal of `B` contained in `A`.  Since (2) gives
`B=A+A lambda`, its contraction is

```text
TB = (T,R) in A.                                   (13)
```

It is the full conductor.  At primes where `c` is invertible this is the
localized cusp calculation (11).  At a prime containing `c`, the element
`T=4-3bc` is a unit, so `TB` already localizes to all of `B`.  Hence (13)
agrees with the conductor at every prime and therefore globally.

There is a second useful presentation.  From (12), `(T,S)` is contained in
`(T,R)`.  Conversely,

```text
cR=(2T-S)/3,                 4=T+3bc,
```

so

```text
R = (TR+3bcR)/4
```

belongs to `(T,S)`.  Therefore

```text
conductor_B/A = TB = (T,R) = (T,S) in A.            (14)
```

The whole normalization defect is one copy of its support.  Since
`B=A+A lambda`, multiplication by the class of `lambda` gives an exact
isomorphism

```text
A/(T,R)  ->  B/A,             f |-> f lambda mod A. (14a)
```

Surjectivity is the quadratic reduction (2), while the kernel is exactly the
conductor definition: `f lambda` lies in `A` iff `fB` lies in `A`.  Moreover
`T=R=0` gives `bc=4/3` and `a=4/(27c^2)`, so

```text
B/A  isomorphic to  A/(T,R)  isomorphic to  Q[c,c^-1]. (14b)
```

Thus the normalization defect is a trivial rank-one module along `E`; every
transverse cusp has delta invariant one.  This is stronger than equality of
set-theoretic supports.

On the normalization, it is the nonradical ideal
`((2-c lambda)^2)=(theta^2)`; in the cusp ring it is
`(theta^2,theta^3)`.  Its support is the single line `theta=0`, whose image is

```text
E = {(4/(27t^2),4/(3t),t) : t in C^*}.             (15)
```

This is exactly THM-2546's empty-fibre curve.

## 4. The singular locus is precisely the conductor support

The first target derivative has the structural form

```text
L_a = 2S.                                           (16)
```

If a point is singular, then `L=L_a=0`.  The cusp identity
`S^2-T^3=27c^2L` forces `T=0`; in particular `c!=0`.  On `T=0`, the next
derivative becomes

```text
L_b = 2(b-9ac) = 2R.                               (17)
```

Thus every singular point lies on `V(T,R)=E`.  Conversely, `T=R=0` gives
`bc=4/3`, `a=4/(27c^2)` and kills `L` and all three derivatives.  Hence

```text
Sing(V(L))_red = E = Supp(conductor).               (18)
```

The exact Jacobian-ideal Groebner basis is

```text
[16a-b^3c, (3bc-4)^2].                              (19)
```

Its square records the nonreduced cusp thickness; (18) is a statement about
the reduced singular locus, not an identification of the Jacobian scheme with
the conductor scheme.

## 5. The `c=0` boundary is smooth and remembers the missing slope

The boundary omitted by the cusp-cylinder localization is

```text
V(L,c) = V(16a-b^2,c),                              (20)
```

a smooth parabola because `L_a=16` there.  Restricting (1) to `c=0` gives

```text
(a,b,c) = (lambda^2/9, 4lambda/3, 0),               (21)
lambda = 3b/4.                                      (22)
```

Thus normalization is already an isomorphism along the whole boundary.  In
the cusp coordinate, however,

```text
theta = 2                                             (23)
```

on all of (20), and `Phi_c` collapses it to `(S,T)=(8,4)`.  The normalization
coordinate `lambda` is exactly the first-order sidecar lost by that collapse:
writing `theta=2-c lambda`, formulas (7) remain finite as `c->0` and give
(21).  This identifies the former “parasitic plane” mechanism at its actual
Jelonek boundary rather than treating `c=0` as a compactification chart.

The unique affine `F`-survivor on (21) is also polynomial:

```text
(x,y,z) = (0, 4lambda/3, -7lambda^2),               (24)
```

and direct substitution gives `F(x,y,z)=(lambda^2/9,4lambda/3,0)`.

## 6. The normalized finite-survivor section

On the cusp cylinder, one more exact factorization is

```text
q=b^2-12a = theta^2(theta-2)^2/(9c^2).              (25)
```

For `c theta !=0`, the unique finite point in the fibre over (7) is

```text
x =  2c/theta^2,
y =  (2-theta)(3theta+2)/(6c),
z = -theta^2(theta-2)^2(theta^2+4theta+2)/(8c^2).  (26)
```

These are the THM-2546 survivor formulas rewritten in the normalization
coordinate.  The companion also verifies (26) directly against the three
original components of `F`, so no inference from a coordinate eliminant is
being substituted for the actual fibre statement.

Three sections expose the geometry:

```text
theta=0:
  target E; conductor and singular curve; fibre empty;

theta=2 (inside c!=0):
  target (0,0,c); survivor (c/2,0,0);

theta=-1:
  target (0,1/c,c); survivor (2c,-1/(2c),9/(8c^2)). (27)
```

At `theta=0`, (26)'s `x=2c/theta^2` pole is the exact total-escape signal.
At `theta=2` and `theta=-1`, the factorization

```text
a=(theta-2)^2(theta+1)/(27c^2)                     (28)
```

identifies respectively the double target `c`-axis component and the simple
`a=0,bc=1` component.  The survivor values in (27) agree exactly with the
exceptional cases of THM-2546.

## 7. Mechanism and non-consequences

- **Underlying object:** the nonzero-`c` Jelonek surface is the semigroup ring
  `Q[c,c^-1,theta^2,theta^3]`; the missing exponent `1` is the entire
  normalization and conductor defect.
- **Map:** `theta` normalizes the cusp; `lambda=(2-theta)/c` is the global
  coordinate that continues across `c=0`.
- **Preserved predicate:** `theta=0` is exactly the singular conductor and the
  empty-fibre curve; `theta!=0` carries the unique finite survivor.
- **Destroyed information:** the raw pair `(S,T)=(theta^3,theta^2)` forgets
  `theta` at the cusp and forgets the boundary slope `lambda` when `c=0`.
- **Boundary:** the coefficient plane itself is not an infinity boundary.
  Its intersection with `V(L)` is the smooth parabola (20), on which the
  normalization is already an isomorphism.
- **No tower promotion:** none of (4), (9), or (14) is asserted for a general
  Keller map or a composite.  A composition law would have to track the
  normalization and conductor of every Jelonek component, not merely multiply
  discriminant polynomials.  No `JC(2)` consequence follows.

## Reproduction

```bash
python3 04-computation/keller_jelonek_normalization_thm2570.py
python3 -O 04-computation/keller_jelonek_normalization_thm2570.py
diff -u 05-knowledge/results/keller_jelonek_normalization_thm2570.out \
  <(python3 04-computation/keller_jelonek_normalization_thm2570.py)
sha256sum 04-computation/keller_jelonek_normalization_thm2570.py \
  05-knowledge/results/keller_jelonek_normalization_thm2570.out
```

Ordinary and optimized output agree byte-for-byte.  The script uses no
`assert`; it verifies the graph kernel, integral and birational relations,
cusp-cylinder inverse, conductor identities, exact Jacobian ideal, smooth
boundary, full normalized survivor, boundary survivor, and all three special
`theta` sections over `QQ`.
