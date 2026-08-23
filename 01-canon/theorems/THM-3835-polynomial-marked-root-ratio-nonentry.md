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
  Etaleness and the Keller equation are not needed.  An independent
  denominator-free resultant proof gives the same contradiction.  The
  intrinsic divisor k=0 is exactly G_m, so an etale plane atlas would pull it
  back to genuine simple, multi-ended poles.  A denominator-free homogeneous
  polynomialization/Keller passport survives for the rational branch.  No
  plane atlas or Jacobian counterexample is constructed.
source: jc_quartic_c3_construct / triangular root-ratio polynomialization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / audit_incoming_3830, 2026-08-23).
  The audit reduced the proof to the polynomial-unit implication from the
  intrinsic Bezout row, checked dominance/injectivity and nonconstancy using
  the exact Laurent arm, and derived the symmetric, normality, and PGL2
  strengthenings.  It separated the determinant-only constant-denominator
  hostile and confirmed that neither etaleness nor characteristic-zero
  valuation arithmetic is load-bearing.  The deterministic companion checks the
  determinant factorization, polynomial-unit boundary, determinant-only
  polynomial-ratio hostile, genuinely rational unimodular-row control,
  nonzero chart dependence relation, exact homogenizations of r and s, and
  the denominator-free reconstruction/Jacobian chain.  This supersedes the provisional
  square-sidecar valuation proof whose zero-polynomial branch was untyped;
  see MISTAKE-456.  A second independent audit rederived a denominator-free
  syzygy, every exceptional resultant branch, the sharp characteristic-three
  failure, and the quotient B/(k)=K[C,C^-1].  Its 52 gates and the primary
  companion agree under normal and optimized replay after LF normalization;
  stored-file raw hashes agree.
depends_on:
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
related:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3845-nonlinear-cubic-keller-atlas-total-degree-contradiction
  - THM-3836-cubic-factor-cofactor-darboux-packet
script: 04-computation/jc2_polynomial_marked_root_ratio_nonentry_thm3835.py
output: 05-knowledge/results/jc2_polynomial_marked_root_ratio_nonentry_thm3835.out
script_sha256: 884f087731f9ab751c49d236e149002d2f26eefe5fab0f4f6858af4f0fcd6007
output_sha256: a70cf307e867a3f73b029f400b9bbf04a49c910bd19ca80b6bd29abffdbda77c
semantic_sha256: 7d1ad90998513133e4c8f017378a3bf35f978a564c9bd1b090ff6a5ae9b0fb9b
second_script: 04-computation/jc2_polynomial_marked_root_ratio_second_proof_thm3835.py
second_output: 05-knowledge/results/jc2_polynomial_marked_root_ratio_second_proof_thm3835.out
second_script_sha256: 5782d25a7ee3d3e4ab273ff2ea55c5b153109c33193846d9b5ec4a494f513df9
second_output_sha256: 3b2770eedf8673b66af216b9fce187794f7db03bb89271adcdbc9b8c878025ec
second_semantic_sha256: eea5967757bc0b8d7e2d482ec61323149a5a7ec31c1d8be5e859aa5f410b6316
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

The determinant-only hostile shows why the intrinsic row law is insufficient:

```text
h=x,                    k=1,                    m=0, C=1
```

satisfies `(2)` and has polynomial ratio `z=x`.  This row is not asserted to
extend to a morphism into the nonlinear cubic surface; it isolates exactly the
extra target-chart content used below.

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
`c h+d k` would therefore be a scalar `kappa in K`.

This last possibility can be excluded directly in the THM-3832 chart, without
importing an etale hypothesis.  Since `h=zk`, the scalar equation is

```text
k(cz+d)=kappa.                                                (8)
```

If `kappa=0`, then `cz+d=0` in `K(z,C)`; the transcendence of `z` forces
`c=d=0`, contradicting the invertibility in `(7)`.  If `kappa!=0`, substitute
`k=r/(Cs)` into `(8)` to obtain

```text
r(z)(cz+d)=kappa C s(z,C).                                   (9)
```

The left side has `C`-degree zero, whereas the coefficient of `C^2` on the
right is `kappa z^2(7z^2+3)`, which is nonzero in `K[z]`.  Since dominance
injects the transcendence basis `z,C` into `K(x,y)`, `(9)` is impossible.
Hence every `PGL_2(K)` transform `(7)` retains a genuine denominator.  By
normality, none is integral either.  This is a projective row statement, not
invariance under a nonlinear target change.

## 4. The denominator-free surviving system

The rational branch should not be studied by cancelling `k`.  Homogenize `(4)`
with `z=h/k`:

```text
R(h,k)=3h^3+7h^2k+k^3,
S(h,k,C)=C h^2(7h^2+3k^2)-k(6h^3+7h^2k-k^3).                (10)
```

The chart polynomialization law becomes

```text
CS=R.                                                        (11)
```

It also clears all denominators in the THM-3832 reconstruction formulas.  If
`A,omega,D` denote the intrinsic functions, then

```text
D S     = k^2+C h^2,
A S     = h(k^2+C h^2),
omega S = k(k^2+C h^2).                                    (12)
```

Finally the weighted-area equation becomes

```text
k Jac(h,C)-h Jac(k,C)=lambda R,             lambda in K*.    (13)
```

Together with `(2)`, equations `(10)--(13)` retain the `k=0` divisor erased
by the affine coordinate `z`, as well as the divisibility obligations needed
to reconstruct `A,omega,D`.  They are necessary conditions only.  Solving
them was the local construction problem; THM-3845 subsequently proves that
the factor/area subpacket is already inconsistent for a dominant polynomial
Keller atlas.  That later etale conclusion does not weaken this theorem's
dominance-only denominator statement.

## 5. Independent denominator-free proof

There is a shorter proof that uses neither the triangular dependence
calculation `(5)` nor the Keller equation.  Under the contradiction
hypothesis of Section 1, write `k=u in K*` and put

```text
q=7z^2+3,                    r=3z^3+7z^2+1,
b=6z^3+7z^2-1.                                                (14)
```

The determinant and the two intrinsic lift laws have the literal syzygy

```text
u C [z^2 q C-b]-r=0.                                          (15)
```

No division by `z,q,C`, or a spectral factor occurs.  The first lift and
`A=hD` also give

```text
u q A-z(1+2uC)=0.                                              (16)
```

Eliminating `z` from `(15),(16)` yields

```text
7u(1+2uC)^2 R_u(A,C)=0,                                       (17)
deg_A R_u=4,                 deg_C R_u=5.
```

The coefficients of `A` and `C^2` in `R_u` are respectively `3` and `-1`,
so `R_u` is nonzero after every specialization `u in K*`.  Since `K[x,y]`
is a domain, `(17)` either makes `C=-1/(2u)` scalar or makes `A,C`
algebraically dependent.  Dominance fails in both cases.

The abstract row/lift lemma works over every field of characteristic not
three.  This boundary is sharp for the mechanism: in characteristic three,

```text
k=1, h=z=x, D=y, C=1-x^2y, m=-xy, A=xy                        (18)
```

satisfies the determinant, both lift laws and `A=hD`, and the target pair is
birational.  This is a hostile boundary for the proof, not a
characteristic-three realization of the full intrinsic surface.

## 6. The denominator is a genuine pole divisor

In the intrinsic coordinate ring `B` of `U`, setting `k=0` gives exactly

```text
B/(k)=K[C,C^-1],                                               (19)
h=3/(7C^2),       m=-7C^2/3,       D=7C^4/9,
A=C^2/3,          omega=0,          theta=-7C^2/3.              (20)
```

The displayed assignments satisfy every defining relation, and `D` is a
unit on the whole arm.  Thus `(19)` is scheme-theoretic: the divisor is
prime, reduced and smooth, with no extra component.

For an actual etale plane atlas, the theorem makes `k` nonconstant.  Etale
base change pulls `(19)` back to a reduced divisor.  Modulo `k`, the determinant
gives `mh=-1`, so `h` is a unit and `h/k` has a genuine simple pole on each
component.  The pulled-back `C` is a nonconstant unit on every component;
on a projective completion it has a zero and a pole at distinct boundary
points.  Hence each pole component is multi-ended.  This records the exact
geometry of the now-closed THM-3811 candidate lane: THM-3841 and THM-3845
exclude its plane atlases by independent global mechanisms, while `JC(2)`
remains OPEN.  **QED.**

## 7. Exact verification

Run

```text
python -B 04-computation/jc2_polynomial_marked_root_ratio_nonentry_thm3835.py
python -B -O 04-computation/jc2_polynomial_marked_root_ratio_nonentry_thm3835.py
python -B 04-computation/jc2_polynomial_marked_root_ratio_second_proof_thm3835.py
python -B -O 04-computation/jc2_polynomial_marked_root_ratio_second_proof_thm3835.py
```

and compare, after LF normalization, with the two frozen outputs named in
the metadata.  Both companions are assertion-free and use exact symbolic
identities; no finite-field inference proves the characteristic-zero claim.
