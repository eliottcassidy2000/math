---
id: THM-2241
title: "Monic transverse response depth and the resultant nonproper quotient"
status: >
  PROVED. Let P,Q in C[x,y] have nonzero constant Jacobian, and suppose P is
  monic of y-degree d. For D=Jac(P,-), the pair (P,Q) is a polynomial
  automorphism if and only if D^(d+1)(x)=0; in the automorphism case
  D^d(x) is nonzero, so the cutoff is sharp. The proof combines the exact
  response kernel of THM-2230, a finite slice-Taylor descent, and monic
  properness. The leading X-coefficient of
  Res_Y(P(X,Y)-U,Q(X,Y)-V) cuts out exactly the Jelonek/nonproper set and
  transforms covariantly under target shears. Monic division gives a unique
  finite-channel section of every target-shear orbit, although the full
  response quotient has infinite C[P]-rank. These are exact decision and
  quotient theorems, not a proof of the planar Jacobian conjecture.
source: codex-2026-07-25-monic-response-depth-resultant
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - THM-2084-cubic-fiber-low-complement-gauss-manin-gate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
  - THM-2206-quadratic-deck-augmentation-and-hamiltonian-characteristic-incompatibility
---

# THM-2241 -- a finite transverse response tail

Let

```text
P,Q in C[x,y],          Jac(P,Q)=kappa in C*,
D=Jac(P,-),             s=Q/kappa.                    (1)
```

Then

```text
D(P)=0,                  D(s)=1.                      (2)
```

Assume throughout Sections 1--4 that `P` is monic in `y` of positive
degree

```text
d=deg_y P.                                            (3)
```

The monic orientation is not a restriction on a prospective planar Keller
counterexample: Section 3 explains the generic affine source change.

## 1. Slice-Taylor termination

THM-2230 gives the exact response kernel

```text
ker D=C[P].                                           (4)
```

Suppose `f in C[x,y]` satisfies `D^N f=0`. Then
`D^(N-1)f` belongs to (4), and

```text
f_1=f-s^(N-1)D^(N-1)f/(N-1)!
```

is killed by `D^(N-1)`. Descending induction gives

```text
D^N f=0                    implies f in C[P,s].       (5)
```

Equivalently, the finite Taylor projector

```text
pi(f)=sum_(j=0)^(N-1) (-s)^j D^j f/j!                (6)
```

lies in `ker D=C[P]`, and repeated application to `D^j f` reconstructs
`f` as a polynomial in `s` with coefficients in `C[P]`.

This step uses only (2), characteristic zero, and the kernel theorem. It
does not assume that `D` is locally nilpotent on the entire ring.

## 2. The monic transverse-response dichotomy

The exact criterion is

```text
(P,Q) is a polynomial automorphism
             iff
D^(d+1)(x)=0.                                         (7)
```

Moreover, on the automorphism branch,

```text
D^d(x)!=0.                                            (8)
```

### Termination implies invertibility

If the right side of (7) holds, (5) gives

```text
x in C[P,s].                                         (9)
```

Thus any sequence on which `(P,s)` remains bounded also has bounded `x`.
Write

```text
P(x,y)=y^d+a_(d-1)(x)y^(d-1)+...+a_0(x).             (10)
```

If `x` and `P(x,y)` are bounded while `|y|` tends to infinity, divide
(10) by `y^d`; its right side tends to `1` and its left side tends to
zero, a contradiction. Hence bounded `(P,s)` also bounds `y`.
The map

```text
F=(P,s):C^2 -> C^2                                   (11)
```

is therefore proper. It is etale because `Jac(P,s)=1`. A proper etale map
of the complex affine plane is a finite topological covering of the simply
connected plane, hence has degree one. Equivalently, finite etale
triviality of `A^2_C` makes (11) an isomorphism. Rescaling the second target
coordinate proves that `(P,Q)` is an automorphism.

### Invertibility gives the sharp depth

Conversely suppose `(P,s)` is an automorphism and write its polynomial
inverse as

```text
x=X(P,s),                y=Y(P,s).                   (12)
```

Because `P(x,Y)-T` is monic and irreducible over `C(x,T)`,

```text
[C(x,y):C(P,x)]=d.                                   (13)
```

Transporting (13) through the coordinate system `(P,s)` gives

```text
[C(P,s):C(P,X(P,s))]=d.                              (14)
```

For a nonconstant polynomial in the one variable `s` over `C(P)`, the
field degree in (14) is exactly its `s`-degree. Therefore

```text
deg_s X=d.                                           (15)
```

On `C[P,s]`, equation (2) says that `D` is differentiation with respect
to `s`. Equations (7)--(8) follow from (15).

In particular, a hypothetical counterexample in this monic orientation
must satisfy

```text
D^n(x)!=0                  for every n>=1.            (16)
```

Indeed, termination at any depth invokes (5) and the same properness
argument.

## 3. Generic monicization and exact scope

Let `m=deg P` and let `P_m` be the top homogeneous part. Choose a linear
source direction `v` with `P_m(v)!=0` and use it as the new `y`-axis.
After this invertible linear source change and a nonzero rescaling of `P`,
the transformed first component is monic of `y`-degree `m`. Vanishing of a
response iterate is unaffected by the nonzero scalar normalization.

Consequently planar JC is equivalent to the following finite polynomial
test in every generic monic orientation:

```text
D_P^(deg P+1)(x)=0.                                  (17)
```

Equation (17) is a **decision invariant**, not a proof that the displayed
polynomial always vanishes. Its useful content is that an all-depth local
nilpotence claim is unnecessary: one transverse orbit and one degree-sharp
cutoff suffice.

The monic hypothesis and the global properness step are load-bearing. For

```text
h=x y^2,
D_h=-2xy partial_x+y^2 partial_y,

D_h(x)=-2xy,       D_h^2(x)=2xy^2,       D_h^3(x)=0, (18)
```

although `h` is not a power of a linear form and `D_h` is not locally
nilpotent (`D_h^n(y)` never terminates). Thus termination for one
coordinate of a leading face does not justify a linear-power face
classification; lower faces and properness cannot be discarded.

## 4. The resultant cuts out nonproperness

Define

```text
R_F(X;U,V)
 =Res_Y(P(X,Y)-U,Q(X,Y)-V).                           (19)
```

A Keller map is etale and hence quasi-finite. Since `P` is monic in `Y`,
the finite roots of the specialization `R_F(X;u,v)` are exactly the
`X`-coordinates of `F^(-1)(u,v)`. The order of a root is the sum of the
local intersection multiplicities over points with that `X`-coordinate.
The nonzero Jacobian makes every such local multiplicity one. Therefore

```text
deg_X R_F(X;u,v)=#F^(-1)(u,v).                       (20)
```

Let `N` be the generic degree of `F` and put

```text
c(U,V)=[X^N]R_F(X;U,V).                              (21)
```

Then `c` is nonzero and

```text
c(u,v)=0
 iff #F^(-1)(u,v)<N
 iff (u,v) belongs to A(F),                          (22)
```

where `A(F)` is the Jelonek/nonproper set.

For completeness, the last equivalence is the standard missing-sheet
criterion for an etale generically finite map. If a target point has fewer
than `N` finite preimages, the missing nearby sheet escapes to infinity,
so the map is not proper there. Conversely a bounded-target escaping
sequence loses a sheet. Monicity makes the escape visible in the chosen
coordinate: bounded `X` and bounded `P` would bound `Y` by (10), so every
such sequence has `|X|` tending to infinity. This is exactly the degree
drop detected by (21).

Thus

```text
A(F)=V(c).                                           (23)
```

In particular, the Keller map is proper, hence an automorphism, exactly
when `c` is a nonzero constant.

## 5. Exact target-shear covariance

For `H in C[T]`, put

```text
Q_H=Q+H(P).                                          (24)
```

Modulo `P(X,Y)-U`,

```text
Q_H(X,Y)-V
 =Q(X,Y)-(V-H(U))
   +H(P(X,Y))-H(U),                                  (25)
```

and the final difference in (25) is divisible by `P(X,Y)-U`. Because that
polynomial is monic in `Y`, the root-product definition of the resultant
gives the exact identity

```text
R_(P,Q_H)(X;U,V)=R_F(X;U,V-H(U)).                    (26)
```

Consequently

```text
c_H(U,V)=c(U,V-H(U)).                                (27)
```

Properness is therefore an invariant of THM-2230's exact response fiber,
while the embedded nonproper curve is defined only up to the corresponding
vertical target shear. This is the correct quotient statement: setting
the curve equal as an embedded subset would discard a real gauge action.

## 6. A finite-channel section, but no finite-rank quotient

Monic division gives a unique expansion

```text
Q=sum_(j=0)^(d-1) a_j(x,P)y^j,
                  a_j in C[x,P].                    (28)
```

A target shear changes only `a_0`, by adding an element of `C[P]`.
Subtracting `a_0(0,P)` produces the unique representative satisfying

```text
a_0(0,P)=0.                                          (29)
```

Thus every target-shear orbit has a canonical `d`-channel polynomial
representative. The word “channel” must not be confused with finite rank:
because `P` and `x` are algebraically independent when `d>=1`, the classes

```text
[x],[x^2],[x^3],...                                  (30)
```

are `C[P]`-linearly independent in `C[x,y]/C[P]`. Hence the lossless
response quotient has infinite `C[P]`-rank.

The sharp next consumer is therefore finite **decision** data--the tail in
(17), the coefficient in (21), or the finite monic channels in (28)--not a
finite-rank replacement for the entire target-shear orbit. QED.
