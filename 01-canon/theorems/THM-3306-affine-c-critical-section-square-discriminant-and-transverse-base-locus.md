---
id: THM-3306
title: "Affine-c critical-section square discriminant and transverse base locus"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the fixed THM-3289 slice C=c+x and E'=1, in either THM-3212 cubic
  accessory response field, the coefficient a of the primitive linear
  subresultant is irreducible of degree 36.  Its failure to be a unit on the
  degree-53 saturated residual is exactly one squarefree irreducible
  degree-36 parameter divisor D: Res_x(a,H)=unit*D^2.  Along D the linear
  coefficient pair has a transverse base locus, the two localized gradient
  cubics have gcd exactly two in y, and the distinguished H-root has exact
  multiplicity two in x.  D is disjoint from the finite-clutch and owner
  boundary walls.
source: root/creative-synthesis-recover/2026-08-03
audit: >
  The primary exact companion constructs H(c,x), the primitive linear and
  quadratic subresultants, and multiplication by c_* in K_i[x]/(a); it also
  supplies additional bounded rational and residual-discriminant evidence
  outside this theorem's statement.  The independent audit does not import,
  execute, or pin the primary.  It rebuilds both THM-3212 response pairs,
  derives the generic subresultants, interpolates H from c=0,1,2,3, uses a
  separate dense-domain multiplication matrix and a full-rank Krylov check,
  and independently recomputes the normal, quadratic-row, H-derivative, and
  wall identities.  Normal, optimized, and stored transcripts agree for both
  routes; both sources have zero assertion nodes and zero floating literals.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go
related:
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3279-affine-transverse-clutch-critical-no-go
script: 04-computation/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.py
output: 05-knowledge/results/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.out
script_sha256: 6448f0a8d8238358adca613610cede3fca72dac210ba0487d126b6d466849697
output_sha256: 28e0395a4e90be88ebd413f27ebdaf82ee4ba51233d72fca42a423448bd6faa5
audit_script: 04-computation/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.py
audit_output: 05-knowledge/results/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.out
audit_script_sha256: 4afadde9302723c5af1a9a209525733c453f928332f21d82e1125f3ddb5662f8
audit_output_sha256: d72ed2fb78bd3e237ea5893a75573d54f2c58fa9080799a8633e73547aaba81d
hash_basis: LF-normalized bytes
---

# THM-3306 -- affine-c critical-section square discriminant and transverse base locus

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Fixed-slice statement

Let `K_i` be either cubic accessory field of
[THM-3212](THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch.md),
embedded in an algebraically closed characteristic-zero field.  Retain its
response pair and owner divisor

```text
V=4SDT^2/Gamma^2,       A=2SET/Gamma,
g=ST=rad(V),            2VA'-AV'=2V.                    (1)
```

Fix the simultaneous affine lane of
[THM-3289](THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go.md)
by

```text
C=c+x,                  d=C'=1,                  k=E'=1. (2)
```

Put

```text
L=y^2+y+(c+x)V,
R_1=2L(2y+1)+VA,
R_2=V^3+V^2y+L(-V'y+2V^2).                              (3)
```

These are the localized gradient cubics of THM-3289.  Their resultant has
the exact form

```text
Res_y(R_1,R_2)=V^3 K(c,x).                               (4)
```

For the two response fields, define

```text
beta_4111=S^3T^8x^9,
beta_3211=S^3T^8x^6(x-1)^3.                             (5)
```

Then `beta_i` divides `K`, and

```text
H(c,x)=K(c,x)/beta_i
      =H_0(x)+cH_1(x)+c^2H_2(x)                         (6)
```

has bidegrees

```text
deg_x(H_0,H_1,H_2)=(53,52,36).                          (7)
```

After dividing the penultimate standard subresultant row by the same
boundary and normalizing by one common element of `K_i^*`, it is

```text
a(x)y+b_0(x)+cb_1(x),
deg_x(a,b_0,b_1)=(36,44,36).                            (8)
```

The following statements hold in both accessory fields.

1. `a` is monic, squarefree, and irreducible over `K_i`.  In
   `A_i=K_i[x]/(a)`, both `b_1` and `H_2` are units.

2. If

   ```text
   c_*=-b_0/b_1 in A_i,                                  (9)
   ```

   then its reduced representative has degree 24 and

   ```text
   b_0+cb_1=b_1(c-c_*),
   H_0+cH_1+c^2H_2=H_2(c-c_*)^2             in A_i[c]. (10)
   ```

3. Let `D_i(c)` be the monic characteristic polynomial of multiplication by
   `c_*` on the 36-dimensional `K_i`-space `A_i`.  Then `D_i` is squarefree
   and irreducible of degree 36, and for some units `u_i,v_i in K_i^*`,

   ```text
   Res_x(a,b_0+cb_1)=u_i D_i(c),
   Res_x(a,H)=v_i D_i(c)^2.                              (11)
   ```

   Consequently, for every `c_0` in an algebraic closure of `K_i`,

   ```text
   gcd_x(a,H(c_0,x)) != 1       iff       D_i(c_0)=0.    (12)
   ```

   In particular, no rational parameter lies on this degeneracy divisor.

4. At the incidence cut out by `a=0` and `c=c_*`, both coefficients of the
   linear subresultant vanish.  The preceding quadratic subresultant remains
   nonzero, both cubics remain cubic, and hence

   ```text
   deg_y gcd(R_1,R_2)=2.                                 (13)
   ```

   The base ideal `(a,b_0+cb_1)` is transverse: its Jacobian determinant in
   `(x,c)` is

   ```text
   a'(x)b_1(x),                                           (14)
   ```

   which is a unit modulo `a`.  Moreover, at the corresponding residual root,

   ```text
   H_x=0,                    H_xx is a unit modulo a.     (15)
   ```

   Thus the distinguished `H`-root has exact multiplicity two in `x`.

5. Let

   ```text
   Delta(c,x)=1-A'(x)(c+x),
   F_i(c)=Res_x(g,Delta),
   B_i(c)=Res_x(g,H).                                    (16)
   ```

   Then

   ```text
   deg F_i=5,       factor degrees of F_i=(1,1,1,2),
   deg B_i=6,       factor degrees of B_i=(1,1,1,1,2),
   B_i=unit * F_i L_i,             deg L_i=1,
   gcd(D_i,B_iF_i)=1.                                    (17)
   ```

   Here `L_i` is the live `q_3` wall at the simple root of `S`.  If
   `t=S`, `V=v_1t+v_2t^2+...`, and `C_S=c+x_S`, the finite, live, and apparent
   exceptional values are

   ```text
   C_S^(finite)=1/2,
   C_S^(live)=-1/2-2v_2/(3v_1^2),
   C_S^(exceptional)=-v_2/(3v_1^2).                      (18)
   ```

   The exceptional value is the exact midpoint of the first two.  It is not
   a root of `B_i` or `F_i`, and `D_i` avoids all three values.

## 2. Quotient-norm proof of the square

The exact symbolic subresultant calculation gives `(6)--(8)`.  Factoring in
`K_i[x]` gives

```text
a irreducible of degree 36,
gcd(a,a')=gcd(a,b_1)=gcd(a,H_2)=1.                       (19)
```

Therefore `A_i` is a field and `(9)` is defined.  Exact reduction modulo `a`
gives

```text
b_0+b_1c_*=0,
H_1+2H_2c_*=0,
H_0-H_2c_*^2=0,                                         (20)
```

which is equivalent to `(10)`.

Let `M_*` denote multiplication by `c_*` on `A_i`.  Then

```text
D_i(c)=det(cI-M_*)=Norm_(A_i/K_i)(c-c_*).                (21)
```

The exact 36-by-36 characteristic polynomial is squarefree and has one
irreducible factor of degree 36.  Since `b_1` and `H_2` are units in `A_i`,
taking norms in `(10)` gives `(11)`.  The resultant criterion gives `(12)`.
This route proves the invariant zero locus and exponent two without declaring
any raw PRS representative canonical; all row normalizations are used only up
to a common field unit, as required by MISTAKE-360.

## 3. Fibre and normal structure on `D_i`

On `a=0,c=c_*`, the linear subresultant is zero.  Direct reduction of all
three coefficients of the preceding quadratic row gives three nonzero
elements of `A_i`.  Also `V'` is a unit modulo `a`, so the second polynomial
in `(3)` does not lose cubic degree.  The standard subresultant criterion now
gives `(13)`.

The Jacobian matrix of the coefficient pair is

```text
[ a'             0  ]
[ b_0'+cb_1'     b_1],                                  (22)
```

whose determinant is `(14)`; equation `(19)` makes it a unit along the base
locus.  Thus `[-b:a]` has a genuine transverse base locus, not a second affine
selector chart.  Blowing up `(a,b)` resolves this projective coefficient map
and records its first normal direction, in the sense suggested by
[THM-2985](THM-2985-multiparameter-normal-map-and-arc-factor-separation.md).
It does not select one of the two common `y` roots.

Finally, partial differentiation of `(6)` in `x`, followed by reduction at
`c=c_*`, gives `(15)`.  Hence this same transverse coefficient base locus
lies on a genuine double-root component of the residual discriminant.

## 4. Owner-wall separation

The exact owner resultants in `(16)` have the degrees and squarefree factor
profiles in `(17)`.  Division in `K_i[c]` gives the exact linear quotient
`L_i`; direct `S`-jet evaluation identifies it with the live `q_3` factor.
The same calculation gives `(18)` and the midpoint identity.  Exact gcds in
`K_i[c]` give

```text
gcd(D_i,B_iF_i)=1,                                      (23)
```

so the new degree-36 degeneracy is neither finite-clutch failure nor hidden
owner-boundary contact.

## 5. Audit and normalization record

The primary route constructs the full bivariate objects before passing to
`A_i` and computes the characteristic polynomial with a domain matrix.  The
independent route reconstructs `H_0,H_1,H_2` from the four exact specializations
`c=0,1,2,3`, builds multiplication through a separate dense matrix, and proves
that `c_*` is primitive by a nonzero 36-by-36 Krylov determinant.  It also
recomputes the quadratic row, normal determinant, `H_x/H_xx`, and both owner
resultants without importing or executing the primary.

The two scripts use different local text serializers for diagnostic polynomial
digests.  Those strings are not theorem data.  The load-bearing objects are
the monic characteristic polynomial defined by `(21)`, its exact factor
profile, the quotient identities `(20)`, and the norm identities `(11)`.

Both sources replay identically in normal and optimized modes against their
stored transcripts.  Dependency hashes, source ASTs, zero floating literals,
and zero Python assertion nodes are checked in the transcripts listed in the
frontmatter.

## 6. Scope and nonconsequences

This theorem is sharply fixed-slice:

- only the two explicit THM-3212 response pairs are covered;
- `C=c+x` and `E'=1`, hence `d=k=1`, are fixed;
- `B=1`, the THM-3289 localization, and its degree-44 boundary are retained;
- the degree-119 residual part of the full repeated-`H` eliminant is not
  classified or projected to parameter space;
- the blow-up resolves only the projective coefficient pair `[-b:a]`; it does
  not choose between the two common `y` roots;
- no deformation in `d,k,V,A`, no second-coordinate cofactor, no Keller mate,
  and no inverse cover is constructed.

In particular this proves neither `JC(2)` nor `DC(2)`, both of which remain
open.

QED.
