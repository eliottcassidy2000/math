---
id: THM-3891
title: "Split quadratic C8 carriers pay two pencil branch values"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  Every split quadratic binary-cubic leading row
  of factor-degree type (0,1,1) whose discriminant is one eighth power has
  exactly two coordinate classes.  The moving class has the THM-3889 Newton
  boundary, reproduced here self-containedly.  The apparently new constant
  carrier C^2 F0 is governed by the branch-value divisor of a binary-cubic
  pencil; Riemann--Hurwitz forces at least two weighted tangent addresses
  unless the total discriminant is visibly reducible.  Hence an irreducible
  discriminant in this entire split grammar has at least two normalization
  places at its unique projective infinity point.
source: jc_sparse_direct_search / post-THM-3889 split-leading classification, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The exact companion verifies the product
  discriminant and both coordinate normal forms, replays every moving-class
  Newton edge and coefficient seam, proves both weighted-initial-form
  identities in the constant-carrier class, and checks the proportional
  reducible exits.  It also contains explicit gcd-degree 0, 1, and 2 pencil
  controls and a separately typed 6,561-row FINITE-EXACT hostile census.
  Normal and -O replays are byte-identical to the frozen 13157-check
  transcript.  Independent audit must recheck the pencil branch-value lemma,
  weighted blow-up-to-normalization implication, and coordinate
  classification.
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3889-maximally-confluent-quadratic-binary-cubic-two-place-obstruction
script: 04-computation/jc2_split_quadratic_c8_pencil_branch_values_thm3891.py
output: 05-knowledge/results/jc2_split_quadratic_c8_pencil_branch_values_thm3891.out
script_sha256: 08c228b39ac0e3156a519a2ffdd0c0209d43c432b65b101b00828c7a595fed98
output_sha256: e6e4671e4b691b00ebf5589bb037e8b26f2a35323eb6e317a9fd0fec6c0f3466
semantic_sha256: dc60618edd264533fae608c7ab43c0db86a47939b99f473de53534c69a576305
hash_basis: raw LF bytes
---

# THM-3891 -- a split eighth power cannot hide one infinity place

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Let `A,C` be target coordinates and `U,V` binary-form
coordinates.  Consider a binary cubic

```text
Phi=Phi_2+A Phi_A+C Phi_C,                                  (1)
```

where `Phi_A,Phi_C in k[U,V]_3` and the quadratic coefficient row factors as

```text
Phi_2=ell_0(U,V) ell_1(A,C;U,V) ell_2(A,C;U,V).             (2)
```

Here `ell_0` is a nonzero constant binary linear form, while `ell_1,ell_2`
have base degree one and binary degree one.  Suppose

```text
Disc(Phi_2)=kappa C^8,                         kappa!=0.    (3)
```

Put `Delta(A,C)=Disc(Phi)`.  Then:

1. `Delta` has degree eight and degree-eight part `kappa C^8`, so its
   projective closure meets the line at infinity only at `[1:0:0]`, with
   total intersection multiplicity eight;
2. the Delone--Faddeev binary index form is contained in `(A,C)`, so the
   index form represents no unit and the corresponding cubic order is
   globally nonmonogenic; and
3. if `Delta` is irreducible, its normalization has at least two distinct
   points over `[1:0:0]`.

The third assertion is the new content.  It distinguishes a single
projective support point from a single normalization place and closes the
whole split factor-degree `(0,1,1)` eighth-power grammar.  It does not cover
nonsplit quadratic rows or a factor-degree pattern not equivalent to `(2)`.

## 1. Exact coordinate classification

After a constant `GL_2(k)` change in `U,V`, take `ell_0=U` and write

```text
ell_i=L_i U+M_i V,                    L_i,M_i in k[A,C]_1. (4)
```

The product-of-determinants formula gives

```text
Disc(Phi_2)=M_1^2 M_2^2 (L_1M_2-L_2M_1)^2.                (5)
```

Unique factorization in `k[A,C]` and `(3)` force

```text
M_1=m_1 C,          M_2=m_2 C,
L_1M_2-L_2M_1=d C^2,              m_1m_2d!=0.             (6)
```

Absorb `m_1,m_2` into the two factors.  There is then a linear form `L` and
a scalar `e!=0` such that, up to a nonzero overall scalar,

```text
Phi_2=U(LU+CV)((L-eC)U+CV).                                (7)
```

There are exactly two cases.

* If `L,C` are independent, a base change preserving the line `C=0` and a
  diagonal source change put `(7)` into

  ```text
  U(AU+CV)((A+C)U+CV).                                    (8)
  ```

* If `L` is proportional to `C`, then

  ```text
  Phi_2=C^2 F_0(U,V),                                     (9)
  ```

  where `F_0` is a squarefree split constant binary cubic.

Both alternatives occur.  Thus `(9)` is the only split coordinate class not
already visible in the moving normal form `(8)`.

## 2. The moving class has two exact Newton addresses

For completeness, the moving calculation is included rather than imported
from a provisional theorem.  Expanding `(8)` and adding arbitrary linear
terms gives

```text
a=A^2+AC+alpha A+alpha_1 C,
b=2AC+C^2+beta A+beta_1 C,
c=C^2+gamma A+gamma_1 C,
d=delta A+eta C.                                           (10)
```

In the infinity chart `x=C/A,z=Z/A`, set

```text
H(x,z)=z^8 Delta(1/z,x/z).                                 (11)
```

Then `H(x,0)=x^8`.  Direct expansion has three exhaustive cases.

* If `delta!=0`, the compact Newton boundary has vertices

  ```text
  (0,2),(3,1),(8,0),                                      (12)
  ```

  and its two nonzero edge equations are

  ```text
  delta w(4-27delta w),                 1+4delta w.        (13)
  ```

  They give distinct branches of orders three and five.

* If `delta=0,gamma!=0,eta!=gamma`, the first edge under `z=wx^2`
  has equation

  ```text
  w[-4gamma^3w^2
    +(-27eta^2+36eta gamma-8gamma^2)w
    +4(eta-gamma)],                                        (14)
  ```

  whose bracket has a nonzero root, while the separate order-four edge is

  ```text
  1+4(eta-gamma)w.                                         (15)
  ```

* If `delta=0,gamma!=0,eta=gamma`, the order-two edge is

  ```text
  gamma^2w^2(1-4gamma w),                                  (16)
  ```

  while the order-three edge is

  ```text
  1+2(2beta+gamma-2gamma_1)w+gamma^2w^2.                   (17)
  ```

  Both have a nonzero root.

Finally, `delta=gamma=0` makes `C` divide `Delta`, so an irreducible
discriminant never reaches that seam.  This proves assertion 3 in the moving
class without any bounded parameter restriction.

## 3. A pencil cannot have only one branch value

The constant-carrier proof uses the following intrinsic lemma.

**Binary-cubic pencil lemma.**  Let `F` be a squarefree binary cubic and `G`
any binary cubic.  The projective quartic divisor

```text
D_F,G(s,t)=Disc(sF+tG)                         on P^1_[s:t] (18)
```

has at least two distinct support points unless `G` is proportional to `F`.

To prove it, put `E=gcd(F,G)` and `r=3-deg E`.  If `r>=2`, cancel `E` and
consider the degree-`r` rational map

```text
-F/E : G/E.                                                (19)
```

Its branch values are support points of `(18)`.  Riemann--Hurwitz gives total
ramification `2r-2`.  Ramification over one target value is at most `r-1`, so
there are at least two branch values.

If `r=1`, the squarefree quadratic `E` has two distinct roots.  The moving
linear residual `(F+tG)/E` collides with those roots at two parameter values.
If the values were equal, one linear form would vanish at both roots of `E`,
forcing it to vanish identically and hence making `F,G` proportional.  Thus
the two collision values are distinct.  The remaining case `r=0` is exactly
the proportional case.  This proves the lemma, including common-factor and
root-at-infinity boundaries.

## 4. Weighted blow-up of the constant carrier

Write the arbitrary linear perturbation of `(9)` as

```text
Phi=C^2F_0+A F_1+C F_2.                                   (20)
```

The form `F_0` is squarefree.  In the same infinity chart, the regular binary
cubic whose discriminant is `H` is

```text
x^2F_0+zF_1+xzF_2.                                        (21)
```

Give `x,z` weights `1,2`.  The first two summands of `(21)` have weight two
and the last has weight three.  Hence the weighted initial form of `H` is

```text
Disc(x^2F_0+zF_1).                                         (22)
```

On the exceptional weighted projective line, `(22)` is precisely the pencil
divisor `(18)` for `(F_0,F_1)`, with `[s:t]=[x^2:z]`.

If `F_1` is not proportional to `F_0`, the pencil lemma supplies at least two
distinct exceptional points.  The proper transform of every reduced branch
meets the exceptional divisor in one point; conversely every point in the
support of the weighted initial form is reached by a branch.  Thus two
exceptional points give two distinct normalization places.

It remains to treat `F_1=lambda F_0`.

* If `lambda=0`, then `Phi=C(CF_0+F_2)`, so

  ```text
  Delta=C^4 Disc(CF_0+F_2)                                (23)
  ```

  is reducible.

* Suppose `lambda!=0` and use the regular local coordinate

  ```text
  s=x^2+lambda z.                                          (24)
  ```

  Equation `(21)` becomes

  ```text
  sF_0+(xs/lambda)F_2-(x^3/lambda)F_2.                    (25)
  ```

  With weights `wt(x,s)=(1,3)`, its discriminant has initial form

  ```text
  Disc(sF_0-(x^3/lambda)F_2).                              (26)
  ```

  If `F_2` is not proportional to `F_0`, the pencil lemma again gives at
  least two exceptional points and hence two normalization places.  If
  `F_2=mu F_0`, then globally

  ```text
  Phi=(C^2+lambda A+mu C)F_0,
  Delta=Disc(F_0)(C^2+lambda A+mu C)^4,                    (27)
  ```

  so `Delta` is nonreduced.

Every irreducible constant-carrier discriminant therefore lies in one of the
two weighted-pencil cases and has at least two infinity places.  Together
with Section 2, this proves assertion 3.

## 5. Index gate, exact boundary, and finite hostile evidence

All four coefficients of `(1)` lie in `(A,C)`.  For
`T=X omega+Y theta` in the Delone--Faddeev order, the determinant of
`(1,T,T^2)` is `-Phi(X,Y)`, hence also lies in `(A,C)` and cannot be a unit.
Adding a scalar to `T` does not change this determinant.  This proves
assertion 2.

The companion also declares the finite universe obtained from the constant
carrier `F_0=UV(U-V)` by choosing each of the eight coefficients of
`F_1,F_2` from `{-1,0,1}`.  Among all `3^8=6,561` labelled rows there are no
coordinate-factor-free rows whose unique infinity point has a raw single
lattice-primitive Newton edge.  This is **FINITE-EXACT hostile evidence**,
not the universe of the theorem; Sections 3--4 are the all-parameter proof.

The theorem closes only split quadratic rows with factor degrees `(0,1,1)`
and a pure eighth-power leading discriminant.  Nonsplit quadratic rows,
higher-degree coefficient profiles, non-binary cubic orders, the later
affine-atlas conditions, and JC(2) remain **OPEN**.
