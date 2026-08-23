---
id: THM-3838
title: "The nonlinear cubic root ratio has numerator and denominator degree at least five"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In every dominant
  plane atlas of the THM-3811 surface, both
  intrinsic row entries h and k have total degree at least five.  THM-3827's
  dual genus-three floor first gives degree at least four; equality would
  make the primitive generic fibre both a smooth plane quartic and a
  hyperelliptic genus-three curve, which is impossible.  Since (h,k) is
  unimodular, these are the reduced numerator and denominator degrees of
  z=h/k.
source: jc_quartic_c3_construct / dual genus-floor degree boundary, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED by root on 2026-08-23.  The audit checked that
  total degree composes exactly, that the generic projective plane-fibre genus
  defect includes every affine and infinite singularity, that equality in the
  primitive genus floor transfers through the forced linear outer polynomial,
  and that adjunction makes a smooth plane quartic canonically embedded and
  hence nonhyperelliptic.  Unimodularity gives the claimed reduced numerator
  and denominator without cancellation.  The companion verifies the plane-curve genus
  table, the unique degree-four composition boundary, both sidecar degrees
  and squarefree controls, the genus-three values, opposite infinity parity,
  and canonical degrees.  The generative-polynomial and equality statements
  are inherited from independently audited THM-3827.  Normal and optimized
  runs byte-match the frozen transcript.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3835-polynomial-marked-root-ratio-nonentry
  - THM-3836-cubic-factor-cofactor-darboux-packet
script: 04-computation/jc2_root_ratio_degree_five_floor_thm3838.py
output: 05-knowledge/results/jc2_root_ratio_degree_five_floor_thm3838.out
script_sha256: 48ba84f2ebd51c51e3136878a932c7c26584159d645b1e3e82ada064fefc6dd7
output_sha256: 48bef28b33ee6d8138e0a02cae306ccd5c3255263fd98d37a596da13ea798197
semantic_sha256: 86617f87ada99bc0e1b51334b244b406172f646a4dd4c7375369e598ed9bb233
hash_basis: raw LF bytes
---

# THM-3838 -- the reduced root ratio starts in degree five

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.
Let

```text
psi:A2_(x,y) -> U                                                (1)
```

be a dominant etale plane atlas of the THM-3811 nonlinear cubic surface, and
let `h,k in K[x,y]` be the pulled-back intrinsic first row.  Then, for total
polynomial degree,

```text
deg h >= 5,                       deg k >= 5.                    (2)
```

The determinant-one law

```text
Ck-mh=1                                                           (3)
```

makes `gcd(h,k)=1`.  Hence `(2)` is also the numerator/denominator floor for
the reduced marked-root ratio

```text
z=h/k.                                                            (4)
```

No support restriction is used.  Degree five and the existence of the atlas
remain open.

## 1. The inherited dual primitive fibrations

THM-3827 applies the closed-polynomial, or generative-polynomial,
factorization separately to the two nonconstant pullbacks `h,k`.  It gives

```text
h=p(g),                         k=q(ell),                         (5)
```

where `p,q in K[T]` are nonconstant and `g,ell in K[x,y]` are closed.  If
`Gamma_g,Gamma_ell` denote the smooth projective geometric generic fibres,
then

```text
genus(Gamma_g)>=3,              genus(Gamma_ell)>=3.             (6)
```

At equality, THM-3827 proves more: `Gamma_g` is isomorphic to its explicit
degree-eight hyperelliptic sidecar, while `Gamma_ell` is isomorphic to the
dual degree-seven hyperelliptic sidecar.  Both have genus three; the former
has two points at infinity and the latter one.

These are proved inheritance statements, not conclusions of the small exact
companion attached here.

## 2. The coarse degree-four floor

For a nonconstant polynomial `f in K[x,y]` of degree `d`, the projective
closure of the generic plane fibre `f=t` has degree `d` and arithmetic genus

```text
p_a=(d-1)(d-2)/2.                                                (7)
```

Its normalization has geometric genus at most `p_a`.  Since the values for
`d=1,2,3` are `0,0,1`, equation `(6)` gives

```text
deg g>=4,                         deg ell>=4.                    (8)
```

Composition with a nonconstant one-variable polynomial has exact degree

```text
deg p(g)=(deg p)(deg g),          deg q(ell)=(deg q)(deg ell).   (9)
```

Thus `deg h,deg k>=4`.

## 3. Why equality four is impossible

Suppose first that `deg h=4`.  Equations `(8),(9)` force

```text
deg p=1,                          deg g=4.                        (10)
```

The degree-four projective plane closure of the generic fibre has arithmetic
genus three.  By `(6)` its normalization has geometric genus at least three,
so equality holds.  The genus defect is the sum of all delta invariants,
including those at infinity; it is therefore zero.  The projective closure
is a smooth plane quartic.

But a smooth plane quartic is nonhyperelliptic.  Indeed, adjunction identifies
its canonical bundle with `O(1)`, so its canonical map is the given plane
embedding.  A hyperelliptic genus-three curve instead has canonical map of
degree two onto a conic.  This contradicts THM-3827's equality statement,
which makes `Gamma_g` hyperelliptic.

The same argument applied to `k=q(ell)` uses the degree-seven sidecar and
excludes `deg k=4`.  Together with the coarse floor this proves `(2)`.

The result sharpens THM-3835: not only must `z=h/k` retain a genuine
denominator, but both entries of its reduced presentation begin at degree
five.  It does not yet classify the degree-five plane models or solve the
factor/cofactor and Keller equations of THM-3836.  No Jacobian counterexample
is claimed.  **QED.**
