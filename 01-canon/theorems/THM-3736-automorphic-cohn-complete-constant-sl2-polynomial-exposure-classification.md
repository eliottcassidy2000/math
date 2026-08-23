---
id: THM-3736
title: "Automorphic Cohn complete constant SL2 polynomial exposure classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For an arbitrary
  constant R in SL2 and an arbitrary polynomial exposed left parameter h,
  every closed row of M0R is classified.  Constant h gives exactly the
  two-sheet Broughton locus of THM-3726.  A nonconstant lower closure forces
  R upper triangular and a nonconstant upper closure forces R lower
  triangular; after an explicit source shear and constant translation of h,
  these are exactly the paired binomial towers of THM-3734.  The exceptional
  nondiagonalizable gauge slope has no closure.  No classified closed row has
  a polynomial Jacobian mate, so no final opposite left shear completes it.
source: root + jc_sparse_direct_search / 2026-08-22
audit: >
  PASS.  An independent audit rederived both raw PDEs and their unique
  highest-degree layers, both valuation descents, both exceptional Jordan
  operators, all four triangular-gauge cancellations and potentials, the
  r=1/r>=2 boundary, and the arbitrary-mate syzygy consequence.  It found no
  lower-degree or zero-entry leak.  Normal, optimized, and frozen output
  agree; script/output/semantic hashes and CHECKS=1436 match.
depends_on:
  - THM-3726-automorphic-cohn-constant-sl2-orbit-broughton-classification
  - THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers
related:
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3723-automorphic-cohn-c3-two-right-resonance
  - THM-3725-automorphic-cohn-opposite-two-right-hyperbolic-resonance-nonentry
script: 04-computation/jc2_automorphic_cohn_complete_constant_sl2_exposure_thm3736.py
output: 05-knowledge/results/jc2_automorphic_cohn_complete_constant_sl2_exposure_thm3736.out
script_sha256: f2d54586b17c0a33e0e44fbb10b03fdcb7171d23f76c68d29c9f45956895ec18
output_sha256: acabffd6173cc1ab8d85fb4f760c16955eac8f233d6510dc06b991ba5b4d421a
semantic_sha256: 2792a981c87a992894cb4422a9d4b2313194ff49cd5ef6a6ed19d0cb99762c34
hash_basis: raw LF bytes
---

# THM-3736 -- the whole constant right orbit is now closed

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3726
classified constant exposed parameters over every constant right-`SL_2`
matrix.  THM-3734 then found genuinely polynomial closures on diagonal
slices.  The apparent remaining problem—nonconstant exposed parameters over
a general constant matrix—collapses by a highest-degree valuation argument:
every survivor is triangular, and every triangular survivor is a source-
sheared member of the diagonal towers.

Work over a characteristic-zero field `k` and pass to an algebraic extension
when a displayed square root is needed.  Put

```text
M_0=[4T^2,2XT-1;1+2XT,X^2],
R=[a,b;c,d] in SL_2(k),
N=M_0R=[alpha;beta],                 z=XT.             (1)
```

Use `curl(U,V)=partial_T U-partial_X V` and
`J(f,g)=f_Xg_T-f_Tg_X`.

## 1. The complete classification

For constant `h`, the answer is exactly THM-3726:

```text
beta+h alpha closed  iff  ch=d-a,       (4a-d)h=b,     (2)
alpha+h beta closed  iff  (a-d)h=-c,    bh=4a-d.       (3)
```

Every row in `(2)--(3)` is a Broughton component `L+L^2S` in
Jacobian-one linear source coordinates and has no polynomial mate.

All **nonconstant** solutions are as follows.

### Lower exposure

The row `beta+h alpha` is closed for a nonconstant polynomial `h` if and only
if

```text
c=0,                      d=a^{-1},
a^2=(r+1)/(2r)             for an integer r>=2,        (4)
D=2a+d,
Y=X+(2b/D)T,
h=b/D + Y^2 H_r^L(YT),                              (5)
H_r^L(z)=[(1+2z/r)^r-1-2z]/(4z^2).
```

The closed row is the gradient of

```text
Q=a Y Phi_r^L(YT),
Phi_r^L(z)=r[(1+2z/r)^(r+1)-1]/[2(r+1)z].             (6)
```

### Upper exposure

The row `alpha+h beta` is closed for a nonconstant polynomial `h` if and only
if

```text
b=0,                      d=a^{-1},
a^2=r/[2(r+1)]             for an integer r>=2,        (7)
D=2a+d,
V=T+(c/D)X,
h=2c/D + V^2 H_r^U(XV),                              (8)
H_r^U(z)=[-(1-2z/r)^r+1-2z]/z^2.
```

The closed row is the gradient of

```text
Q=d V Phi_r^U(XV),
Phi_r^U(z)=-r[1-(1-2z/r)^(r+1)]/[2(r+1)z].            (9)
```

The formal `r=1` instances of `(4)--(9)` have `H_1^L=H_1^U=0`; they are
constant-parameter boundary points already contained in `(2)--(3)`.  Thus
the nonconstant statement correctly starts at `r=2`.

No potential in `(2)--(9)` has a polynomial Jacobian mate.  Consequently no
polynomial final opposite left shear turns the corresponding exposed matrix
into a Jacobian matrix.

## 2. The two raw closure equations

The rows in `(1)` are

```text
alpha=(4aT^2+c(2z-1), 4bT^2+d(2z-1)),
beta =(a(1+2z)+cX^2,  b(1+2z)+dX^2).                 (10)
```

Direct differentiation gives the full lower equation

```text
[4aT^2+c(2z-1)]h_T-[4bT^2+d(2z-1)]h_X
 +2[cX+(4a-d)T]h+2[(a-d)X-bT]=0,                     (11)
```

and the full upper equation

```text
[a(1+2z)+cX^2]h_T-[b(1+2z)+dX^2]h_X
 +2[(a-d)X-bT]h+2[cX+(4a-d)T]=0.                     (12)
```

These identities also retain all zero-entry boundaries of `R`.

## 3. Highest-degree valuation forces triangularity

Suppose `m=deg h>0` and let `h_m` be its highest homogeneous part.  Define
the linear vector field

```text
W=(2aT+cX) partial_T -(2bT+dX) partial_X.             (13)
```

The degree-`m+1` part of `(11)` is

```text
2T W(h_m)+2[cX+(4a-d)T]h_m=0.                        (14)
```

For `F=2T h_m`, equation `(14)` is exactly

```text
W(F)+(2a-d)F=0.                                       (15)
```

Because `F` is divisible by `T`, take its finite `T`-valuation `s>=1`.  If
`c!=0`, the term `cX partial_T F` in `(15)` has `T`-valuation `s-1`, with
nonzero leading coefficient, while every other term has valuation at least
`s`.  This is impossible.  Therefore every nonconstant lower closure forces
`c=0`.

Likewise, the degree-`m+1` part of `(12)`, with `G=2Xh_m`, gives the same
eigen-equation

```text
W(G)+(2a-d)G=0.                                       (16)
```

Now `G` is divisible by `X`.  If `b!=0`, the term
`-2bT partial_X G` has uniquely lower `X`-valuation.  Hence every
nonconstant upper closure forces `b=0`.

This is the decisive reduction: a general four-entry constant matrix cannot
support a nonconstant exposed closure.  The conclusion is caused by the
*off-triangular entry opposite the exposed row*, not by a low-degree search.

## 4. The exceptional Jordan slope is empty

In the lower case `c=0`, the determinant gives `d=a^{-1}`.  The source shear
used below has denominator `D=2a+d`.  If `D=0`, then `d=-2a` and
`a^2=-1/2`.  On homogeneous polynomials of degree `m+1`, `(15)` becomes

```text
[2a(m+3) I -2b T partial_X]F=0.                       (17)
```

The operator `T partial_X` is nilpotent, whereas the scalar `2a(m+3)` is
nonzero in characteristic zero.  Thus the bracket in `(17)` is invertible,
contradicting `F!=0`.

For the upper case `b=0`, the same boundary gives

```text
[2a(m+3) I +c X partial_T]G=0,                        (18)
```

again nonzero scalar plus nilpotent.  Therefore `D=0` hides no exceptional
closure in either orientation.

## 5. Every triangular survivor is a diagonal tower in disguise

Assume first `c=0` and `D!=0`.  Set

```text
Y=X+(2b/D)T,                    kappa=b/D,
h(X,T)=kappa+g(Y,T).                                   (19)
```

Substitution into `(11)` cancels both extra terms because

```text
4a(2b/D)-4b=-2d(2b/D),
(4a-d)kappa-(a-d)(2b/D)-b=0.                          (20)
```

What remains is precisely the diagonal lower closure equation for `g` in
the Jacobian-one source coordinates `(Y,T)`.  THM-3734 proves that all its
polynomial solutions are `(4)--(6)`, with no homogeneous additions.

Assume instead `b=0`.  Put

```text
V=T+(c/D)X,                      kappa=2c/D,
h(X,T)=kappa+g(X,V).                                   (21)
```

The two cancellation identities are

```text
c-d(c/D)=2a(c/D),
(a-d)kappa+c-(4a-d)(c/D)=0.                           (22)
```

Equation `(12)` becomes the diagonal upper equation in `(X,V)`, and
THM-3734 gives exactly `(7)--(9)`.

The source changes have

```text
J(Y,T)=1,                        J(X,V)=1.              (23)
```

Hence the triangular potentials `(6)` and `(9)` are source-automorphic
images of the diagonal divided-power potentials.  Existence of a polynomial
Jacobian mate is invariant under such a source automorphism.  THM-3734
proves that the diagonal potentials have no mates, while THM-3726 proves the
same for the constant Broughton locus.  This completes both the
classification and the nonentry assertion.  **QED.**

## 6. Exact audit surface and remaining counterexample route

Reproduce with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_complete_constant_sl2_exposure_thm3736.py
python3 -B -O 04-computation/jc2_automorphic_cohn_complete_constant_sl2_exposure_thm3736.py
```

The companion derives `(11)--(16)` symbolically, verifies both gauge
cancellations, checks the raw matrix gradients for both towers through depth
eight and four triangular parameters per orientation, and runs the
top-degree valuation systems through degree six over all 116 integer `SL_2`
matrices in the `[-3,3]` box.  It also checks both Jordan orientations through
homogeneous degree eight.  The transcript contains 1,436 exact gates;
normal and optimized runs byte-match.

This theorem completes the **constant right matrix** orbit only.  It does not
classify nonconstant right matrices or longer alternating Cohn words.  A
counterexample in this framework must now acquire its charge coupling from a
genuinely polynomial right factor (or leave the Cohn orbit), rather than from
an arbitrarily complicated polynomial left exposure of a constant right
matrix.
