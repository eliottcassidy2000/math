---
id: THM-3835
title: "A polynomial marked-root ratio cannot support the nonlinear cubic plane atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  dominant morphism from the polynomial plane to the THM-3811 surface, the
  marked-root ratio z=h/k is genuinely rational and its reduced denominator
  is k up to a scalar.  If z were polynomial, the intrinsic Bezout law would
  make k a scalar unit, and the THM-3832 chart equation would make z,C
  algebraically dependent, contradicting dominance.  In fact neither z nor
  any PGL2(K) transform of z is polynomial or integral over K[x,y].
  Etaleness and the Keller equation are not needed.  A denominator-free homogeneous
  polynomialization/Keller passport survives for the rational branch.  No
  plane atlas or Jacobian counterexample is constructed.
source: jc_quartic_c3_construct / triangular root-ratio polynomialization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / audit_incoming_3830, 2026-08-23).
  The audit reduced the proof to the polynomial-unit implication from the
  intrinsic Bezout row, checked dominance/injectivity and nonconstancy using
  the exact Laurent arm, and derived the symmetric, normality, and PGL2
  strengthenings.  It separated the non-dominant constant-denominator
  hostile and confirmed that neither etaleness nor characteristic-zero
  valuation arithmetic is load-bearing.  The deterministic companion checks the
  determinant factorization, polynomial-unit boundary, non-dominant
  polynomial-ratio hostile, genuinely rational unimodular-row control,
  nonzero chart dependence relation, exact homogenizations of r and s, and
  the denominator-free Jacobian chain.  Normal and optimized replay agree.
depends_on:
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
script: 04-computation/jc2_polynomial_marked_root_ratio_nonentry_thm3835.py
output: 05-knowledge/results/jc2_polynomial_marked_root_ratio_nonentry_thm3835.out
script_sha256: a98cf391f4bca82e89f9a10d4b404325ebcd9c7d355a3ab96cac15384cd7cf0b
output_sha256: 3aa5626786e0498db37b2a7625c6cf4e5a4b37d0ea365c1e9a5a1f3e8e39286c
semantic_sha256: 436eb2336cafb2ead5cc67ffc4ae28050315804a2651bf8133be532d63884c2b
hash_basis: raw LF bytes
---

# THM-3835 -- the marked-root ratio must retain its denominator

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.  Let

```text
psi:A2_(x,y) -> U
```

be a dominant morphism to the THM-3811 nonlinear cubic surface.  Write again
`h,k,m,C in K[x,y]` for the pullbacks of the intrinsic functions.  Then

```text
z:=h/k in K(x,y)                                                (1)
```

does not belong to `K[x,y]`.  More precisely, `(1)` is already reduced in the
UFD `K[x,y]`, so its denominator is `k` up to a nonzero scalar.  It is not
integral over `K[x,y]`, and no constant projective change of the row makes it
polynomial.

This theorem is stronger in hypotheses than an etale obstruction: it uses
dominance but neither etaleness nor a constant Jacobian.  It does not exclude
the genuinely rational root-ratio chart or construct a plane atlas.

## 1. The determinant fixes the reduced denominator

The intrinsic row law pulls back to

```text
Ck-mh=1.                                                       (2)
```

Thus `(h,k)=K[x,y]`; in particular `gcd(h,k)=1`.  Since polynomial rings over
fields are UFDs, `h/k` is in lowest terms.  Any other reduced presentation of
the same rational function differs by a scalar unit, proving the denominator
claim.

Suppose now that `z=h/k` were polynomial.  Then `h=zk`, and `(2)` would give

```text
k(C-mz)=1.                                                     (3)
```

Hence `k=kappa in K*`.

The non-dominant hostile shows why the last hypothesis matters: the row

```text
h=x,                    k=1,                    m=0, C=1
```

satisfies `(2)` and has polynomial ratio `z=x`.  What fails is not the row
law but dominance over the nonlinear cubic surface.

## 2. The triangular chart makes the dominance contradiction explicit

Use the THM-3832 polynomials

```text
r=3z^3+7z^2+1,
q=7z^2+3,
b=6z^3+7z^2-1,
s=z^2qC-b.                                                     (4)
```

Its exact formula `k=r/(Cs)` and `(3)` give

```text
kappa C(z^2(7z^2+3)C-(6z^3+7z^2-1))
 -(3z^3+7z^2+1)=0.                                            (5)
```

The left side is a nonzero polynomial in two formal variables `z,C`: its
`C^2` coefficient is `kappa z^2(7z^2+3)`.  But THM-3832 proves

```text
K(U)=K(z,C).                                                   (6)
```

Dominance makes `K(U)->K(x,y)` injective, so the images of the transcendence
basis `z,C` cannot satisfy `(5)`.  This contradiction proves `(1)`.

Equivalently, once `(3)` makes `k` constant, dominance already fails because
the nonconstant intrinsic function `k` has constant pullback.  Argument `(5)`
records that failure entirely inside the new root-ratio chart.

## 3. Every projective row direction retains a denominator

Since `K[x,y]` is normal, a fraction in its field that is integral over it is
already polynomial.  The result above therefore also says that `z` is not
integral over the source ring.  Interchanging `h,k` in the unit-ideal argument
shows similarly that `k/h` is not polynomial.

More generally, let

```text
zeta=(a h+b k)/(c h+d k),              ad-bc!=0.                (7)
```

The two linear forms still generate the unit ideal because their constant
coefficient matrix is invertible.  If `zeta` were polynomial, its denominator
`c h+d k` would therefore be a scalar unit.  But THM-3827 proves that `h,k`
are algebraically independent in every dominant plane pullback, so no nonzero
constant linear combination of them is scalar.  Hence every `PGL_2(K)`
transform `(7)` retains a genuine denominator.  This is a projective row
statement, not invariance under a nonlinear target change.

## 4. The denominator-free surviving system

The rational branch should not be studied by cancelling `k`.  Homogenize `(4)`
with `z=h/k`:

```text
R(h,k)=3h^3+7h^2k+k^3,
S(h,k,C)=C h^2(7h^2+3k^2)-k(6h^3+7h^2k-k^3).                 (8)
```

The chart polynomialization law and the THM-3832 weighted-area equation become

```text
CS=R,                                                         (9)
k Jac(h,C)-h Jac(k,C)=lambda R,             lambda in K*.    (10)
```

Together with `(2)`, equations `(8)--(10)` retain the `k=0` divisor erased by
the affine coordinate `z`.  They are necessary conditions only.  Solving
them, extending the remaining intrinsic generators polynomially, and meeting
the cubic two-arm and nonproperness passports all remain open.  No planar
Jacobian counterexample is claimed.  **QED.**
