---
id: THM-3236
title: "Contact spectrum, primitive-element gate, and root reconstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every THM-3232 etale first-contact stratum, the characteristic
  polynomial chi_m(T)=Norm(T-c_m) is the exact root-free resultant quotient
  Res_x(S_m,Tf'g'+Omega_m)/Res_x(S_m,f'g').  Its squarefree radical is the
  minimal polynomial of c_m.  The contact coefficient generates the whole
  root algebra iff chi_m is squarefree iff its discriminant is nonzero; in
  that case the root coordinate is a unique polynomial in c_m.  The exact
  identity Disc(chi_m)=Disc(S_m)kappa_m^2 identifies the gate with a Krylov
  basis determinant.  A jointly separating second observer generically
  resplits a collapsed spectrum outside at most binom(r,2) pencil values.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion computes every spectrum both as
  a multiplication characteristic polynomial and by resultant interpolation.
  It verifies five quadratic power-family gates, the discriminant/Krylov
  identity, exact root reconstruction, quadratic and cubic delayed
  resplitting, the cubic three-slope exceptional polynomial, and a sharp F2
  small-field pencil boundary.  Normal and optimized runs byte-match the
  stored transcript and the LF-normalized hashes below.  An independent
  hostile audit rederived the resultant normalization and constant term, the
  etale radical/minimal-polynomial statement, all five primitive-element
  equivalences, the Vandermonde/Krylov discriminant identity, affine and
  swap covariance, reconstruction, the pencil bound, and every exact hostile,
  including a raw cubic confirming that the unnormalized resultant quotient
  is already monic and exact, and found no defect.
depends_on:
  - THM-3232-root-free-contact-stratum-norm-and-discriminant-power
related:
  - THM-3229-hasse-pluecker-simple-root-contact-gcd-flag-and-degree-termination
  - THM-3227-selected-root-residue-contact-trie-primitive-carry-and-delayed-resplitting
  - THM-3178-squarefree-resultant-tangent-cone-and-first-witt-norm
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/contact_spectrum_primitive_element_thm3236.py
output: 05-knowledge/results/contact_spectrum_primitive_element_thm3236.out
script_sha256: 20bbaa807c32779a1e57752ea75b8df4a0c68ea898313b0d1a2dafee4d0bc840
output_sha256: cc72963331ca37917fc33d160078d03501d2857ea45df94bec069a41b044ba97
hash_basis: LF-normalized bytes
---

# THM-3236 -- contact spectrum, primitive-element gate, and root reconstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3232 takes the product of all first-live coefficients on one contact
stratum.  That norm cannot say whether two roots carry the same coefficient.
The full characteristic polynomial is the correct lossless symmetric lift:
it retains the entire unordered multiset of contact values, and its
discriminant decides exactly when those values recover the root algebra.

This is a root-free reconstruction theorem, not a root ordering theorem.

## 1. The contact spectrum

Retain the notation and hypotheses of THM-3232.  Thus

```text
S=S_m=G_(m-1)/G_m,       r=deg S>=1,
A=K[x]/(S),              d=f'g' in A^*,
c=c_m=-Omega_m/d in A^*.                                  (1)
```

Define the monic contact-spectrum polynomial

```text
chi_c(T)=Norm_(A/K)(T-c) in K[T].                          (2)
```

Since `S` is monic, the same resultant convention as THM-3232 gives the
exact root-free formula

```text
chi_c(T)
 =Res_x(S,Td+Omega_m)/Res_x(S,d).                          (3)
```

Indeed, `T-c=(Td+Omega_m)/d`; apply multiplicativity of the norm.  The
coefficient of `T^r` in the numerator of `(3)` is `Res(S,d)`, so the quotient
is monic without any extra normalization.  Its constant term recovers the
THM-3232 norm:

```text
chi_c(0)=(-1)^r Norm_(A/K)(c).                             (4)
```

After a splitting extension,

```text
chi_c(T)=product_(S(a)=0)(T-c(a)).                         (5)
```

Thus `chi_c` is precisely the unordered contact-value spectrum.  No root
factorization or ordering occurs in `(3)`.

## 2. Minimal polynomial and the primitive-element gate

Let `mu_c` be the minimal polynomial of the element `c in A`.  Because `A`
is etale, multiplication by `c` is diagonalizable over an algebraic closure.
Its eigenvalues are the values in `(5)`, with repetitions exactly when two
geometric roots have the same contact coefficient.  Consequently

```text
mu_c=rad(chi_c),                                           (6)
```

where `rad` means the product of the distinct monic irreducible factors.  In
particular

```text
dim_K K[c]=deg mu_c
            =number of distinct geometric values c(a).    (7)
```

The following conditions are equivalent:

```text
(i)   c(a) are pairwise distinct on the r roots of S;
(ii)  chi_c is squarefree;
(iii) Disc_T(chi_c)!=0;
(iv)  deg mu_c=r;
(v)   K[c]=A.                                              (8)
```

The implications `(i)--(iv)` follow from `(5)--(7)`.  Condition `(iv)` says
that the subalgebra generated by `c` already has the full dimension of `A`,
which is `(v)`.  This is the exact primitive-element gate.

If the conditions fail, `Spec(A)->Spec(K[c])` is the canonical quotient that
glues precisely the roots with equal first-contact value.  The norm remains
nonzero, but it cannot see this gluing.

## 3. The discriminant/Krylov identity

Use the power basis

```text
1,x,...,x^(r-1)                                            (9)
```

of `A`, and let `B_c` be the matrix whose columns are the coordinates of

```text
1,c,...,c^(r-1).                                          (10)
```

Put

```text
kappa_c=det(B_c).                                         (11)
```

Then

```text
Disc_T(chi_c)=Disc_x(S) kappa_c^2.                        (12)
```

To prove `(12)`, evaluate the two power bases at all geometric roots of `S`.
If `V_x` and `V_c` are the resulting Vandermonde matrices, then

```text
V_c=V_x B_c.                                              (13)
```

Squaring determinants gives `(12)`.  Since `S` is squarefree,
`Disc(S)!=0`; hence

```text
kappa_c!=0 iff Disc(chi_c)!=0.                            (14)
```

This identifies the spectral collision wall both as a resultant
discriminant and as failure of one finite Krylov basis.  Clearing the
nonzero powers of `Res(S,d)` in `(3)` turns it into an ordinary polynomial
wall in any bounded-degree coefficient chart.

## 4. Exact root reconstruction

Assume the equivalent conditions `(8)`.  Evaluation at `c` gives an
isomorphism

```text
K[T]/(chi_c)  ->  A,
T             |-> c.                                     (15)
```

Therefore there is a unique polynomial `R_c(T)` of degree less than `r` such
that

```text
x=R_c(c) in A.                                            (16)
```

Equivalently, at every geometric root,

```text
a=R_c(c(a)).                                              (17)
```

The coefficients of `R_c` are obtained by solving

```text
B_c coeff(R_c)=coeff(x).                                  (18)
```

Thus a surviving contact value identifies its root whenever the spectral
gate is open.  Formula `(16)` is canonical in the chosen affine `x`-chart,
but it does not order the roots: Galois acts simultaneously on both sides of
`(17)`.

## 5. Covariance and invariant content

For the affine source change `x=rho y+sigma` used in THM-3232,

```text
c~=rho^(m-1)c,
chi_c~(T)=rho^[r(m-1)] chi_c(T/rho^(m-1)).                 (19)
```

Hence pairwise distinctness and the primitive-element gate are affine-source
invariant, while

```text
Disc(chi_c~)=rho^[r(r-1)(m-1)]Disc(chi_c).                (20)
```

Constant independent row rescalings leave `chi_c` unchanged.  A row swap
gives

```text
chi_(-c)(T)=(-1)^r chi_c(-T),                             (21)
```

so it also preserves the gate.  A general constant `GL_2` target frame has
the root-dependent derivative multiplier from THM-3232; it can create or
remove spectral collisions and is not silently an invariance.

## 6. A two-observer delayed-resplitting lemma

The spectral construction is not special to one element.  Let `h_1,h_2 in A`
and suppose their joint geometric values separate the roots:

```text
a |-> (h_1(a),h_2(a)) is injective.                        (22)
```

For each unordered root pair, equality of

```text
h_1+lambda h_2                                             (23)
```

at that pair excludes at most one `lambda in Kbar`; if the two `h_2` values
are equal, `(22)` means it excludes none.  Therefore all but at most

```text
binom(r,2)                                                 (24)
```

values of `lambda` make `(23)` a primitive element of `A`.  In particular a
good `lambda in K` exists when `K` is infinite, or when
`|K|>binom(r,2)`.  Conversely, if `(22)` fails, no pencil member can separate
the glued pair.

For contact work one may take

```text
h_1=c_m,
h_2=d_j=-Omega_j/(f'g') in A                              (25)
```

for a later raw normalized coefficient difference.  When `j>m`, `d_j` is
not by itself the order-`j` transition coefficient after composition; it is
only the displayed raw germ difference.  Joint generation in `(22)` is an
additional hypothesis, and mixing different affine weights makes the
chosen pencil parameter chart-dependent.  The generated subalgebra and the
existence of a good pencil member remain well-defined in the fixed chart.

The field-size condition cannot be deleted from the base-field conclusion.
In `F_2^3`, the three joint values

```text
(0,0), (0,1), (1,0)                                       (26)
```

are distinct, but both available pencil members have only two distinct
values.  This abstract hostile concerns the pencil lemma; a contact algebra
still comes with the external primitive coordinate `x`.

## 7. Exact power-family controls

Take

```text
P=x^2-2,                  f=P,        g=P+P^m.             (27)
```

THM-3232 gives `c_m=(2x)^(m-1)`.  For `m=2,...,6`, the primitive/collapsed
pattern is exactly

```text
primitive, collapsed, primitive, collapsed, primitive.    (28)
```

At the first two depths,

```text
m=2: chi(T)=T^2-8,
m=3: chi(T)=(T-8)^2.                                     (29)
```

Both norms are nonzero, so `(29)` is a sharp proof that THM-3232's norm does
not imply root reconstruction.  At the collapsed order three, the next raw
difference is

```text
d_4=6x.                                                   (30)
```

Thus `h=c_3+d_4=8+6x` has

```text
chi_h(T)=T^2-16T-8,      Disc(chi_h)=288,
kappa_h=6,               x=(h-8)/6.                       (31)
```

This is exact delayed spectral resplitting.

A degree-three control makes the collision packet visible over `Q`.  For

```text
P=x^3-3x,                f=P,        g=P+P^2,              (32)
```

one has

```text
c_2=P'=3x^2-3,
chi_(c_2)(T)=(T+3)(T-6)^2,
mu_(c_2)(T)=(T+3)(T-6).                                  (33)
```

The next raw difference is `d_3=6x`; the pair `(c_2,d_3)` separates all
three roots.  The Krylov determinant of `c_2+lambda d_3` is

```text
54 lambda(4lambda^2-3),                                  (34)
```

so the three geometric exceptional slopes in `(24)` are attained exactly.

## 8. Scope and frontier consequences

The information ladder is now exact:

```text
selected roots and values
   -> contact spectrum chi_c
   -> norm (constant term)

chi_c squarefree
   -> coefficient-labelled root-algebra reconstruction
   -/> an ordered root or a physical carrier.             (35)
```

For Jacobian and PRS work, `Disc(chi_c)` is a finite algebraic wall deciding
whether one higher-contact coefficient can replace the root coordinate.  It
may remove a branch selector on the open gate, but it supplies neither a
Keller trajectory nor a common physical row.

For the THM-3227 contact trie, `(31)--(34)` are the root-free analogue of
delayed resplitting.  They retain the partition of roots by coefficient
values but destroy the trie's affine edge labels and additive prime carry.

For Gaussian moments and LRC, the theorem isolates a precise sufficient
sidecar: a physical construction exposing `c_m`, together with
`Disc(chi_c)!=0`, would recover the root address algebraically.  Current
canon supplies neither such carrier nor a uniform open-gate theorem.  No
moment noncancellation, semantic owner alignment, row decrement, JC(2), or
LRC(14) conclusion is asserted.

## 9. Connection contract

```text
source:      one THM-3232 etale contact algebra and its unit c_m;
map:         characteristic polynomial / resultant spectrum;
target:      unordered contact values, primitive gate, reconstruction R_c;
preserved:   all symmetric value data, Galois action, affine gate status;
destroyed:   root ordering, physical provenance, additive carry;
sidecar:     physical survival of c_m and an open spectral discriminant.
```

## 10. Exact companion

The assertion-independent companion

```text
04-computation/contact_spectrum_primitive_element_thm3236.py
```

uses exact rational arithmetic.  It computes spectra independently as
multiplication characteristic polynomials and interpolated Sylvester
resultants, verifies `(12)` and reconstruction, exhausts the five quadratic
power depths, checks both delayed-resplitting controls and the sharp cubic
exceptional pencil, and records the `F_2` boundary.  Normal and optimized
runs byte-match

```text
05-knowledge/results/contact_spectrum_primitive_element_thm3236.out
```

and the LF-normalized hashes are pinned in the frontmatter.

QED.
