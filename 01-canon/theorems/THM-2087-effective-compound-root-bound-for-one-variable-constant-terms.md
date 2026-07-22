---
id: THM-2087
title: "An effective compound-root bound for one-variable constant terms"
status: >
  PROVED. Let f be a complex Laurent polynomial with exact extreme exponents
  -M<0<N and nonzero coefficients at both extremes. Put a=min(M,N),
  b=max(M,N), d=a+b, C=binom(d,a), and K=binom(d-1,a-1). Then
  CT(f^m) is nonzero for some 1<=m<=C+K-1. Since K=Ca/d<=C/2,
  this is at most floor(3C/2)-1. The proof is self-contained: the logarithmic small-root
  identity converts L initial zero constant terms into contact of order L+2
  between the small-root product and c*t; the compound polynomial of all
  a-subset root products has t-pole order at most K by weighted degree plus
  complementary-subset duality; and its evaluation at
  c*t is a nonzero polynomial of degree at most C+K by the same transitive
  Galois orbit-product obstruction as THM-2067. For a=1 the bound is the
  sharp d=M+N. This supplies an unconditional effective seed bound for
  THM-2022, but it does not prove the open sharp d bound for general a,b.
source: codex-2026-07-22-gmc-paper-audit
depends_on: []
related:
  - THM-1550-an-exact-criterion-for-the-toral-nullcone
  - THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2
  - THM-1650-newton-polygon-of-the-effective-dvdk-bound
  - THM-2067-galois-orbit-product-closes-one-variable-dvdk
  - THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation
script: 04-computation/tnc_compound_root_effective_bound_codex_20260722.py
output: 05-knowledge/results/tnc_compound_root_effective_bound_codex_20260722.out
script_sha256: 8b585a8ec96360bcf9d597a8d1d1c9d84f43cec7e10de44bcac88c441ebd967f
output_sha256: 01f27494811037e542ea7581886ddba9a1042e50bf7f05c047ec98c7b7f0d048
hash_basis: repository blobs with LF line endings
---

# THM-2087 -- an effective compound-root constant-term bound

THM-2067 proves the bare one-variable seed needed by NC2/GMC(2): a Laurent
polynomial with exponents of both signs has a nonzero constant term in some
positive power.  Its closing paragraph deliberately claims no bound on the
first such power.  The same root object in fact gives a coarse but explicit
bound once one retains not just the chosen small-root product, but the
polynomial of **all subset products of the same size**.

## 1. Statement

Let

```text
f(u)=sum_(q=-M)^N c_q u^q,       c_(-M)c_N != 0,
M,N>=1.
```

Put

```text
a=min(M,N),       b=max(M,N),       d=a+b,
C=binom(d,a),
K=binom(d-1,a-1)=C*a/d,
B=C+K-1.                                                   (1)
```

> **Theorem.** There is an integer `m` with
>
> ```text
> 1<=m<=B
> ```
>
> such that `CT(f^m)!=0`.

Since `a<=b`, one has `K<=C/2` and hence

```text
B<=floor(3C/2)-1.                                         (2)
```

If `a=1`, then `C=d`, `K=1`, and the bound is

```text
B=d=M+N.                                                   (3)
```

Thus the theorem proves the sharp Sturmfels bound on the entire
`min(M,N)=1` boundary, while giving a finite exponential-size bound in every
bidegree.  It does **not** prove the conjectural bound `m<=M+N` when
`min(M,N)>=2`.

Replacing `u` by `u^(-1)` preserves every constant term and exchanges `M`
and `N`.  We may therefore orient the polynomial so that its negative width
is `a` and its positive width is `b`.  Write

```text
f(u)=u^(-a)R(u),
R(X)=r_0+r_1X+...+r_dX^d,       r_0r_d!=0.                (4)
```

If `c_0!=0`, then `CT(f)!=0` and the theorem is immediate.  The argument
below includes that case but is aimed at the cancellation regime.

## 2. Small roots and the exact contact order

Introduce

```text
Phi(X,t)=X^a-tR(X).                                       (5)
```

For small nonzero `t`, exactly `a` roots tend to zero and `b` roots tend to
infinity.  Let `S_0` be the set of small roots and put

```text
Pi(t)=product_(alpha in S_0) alpha.                       (6)
```

The small roots form the roots of the degree-`a` Weierstrass factor of
`Phi` at `(X,t)=(0,0)`.  Consequently `Pi` is an ordinary convergent power
series in `t`, not merely a Puiseux series.  Its leading term is

```text
Pi(t)=c*t+O(t^2),       c=(-1)^(b+d+1)r_0!=0.             (7)
```

For completeness, factor on an annulus separating the two root sets:

```text
1-tf(u)
 =(-tr_d)
   product_(alpha in S_0)(1-alpha/u)
   product_(beta notin S_0)(u-beta).                      (8)
```

The logarithms of the first product have only negative powers of `u`.
After writing `u-beta=-beta(1-u/beta)`, the nonconstant logarithms of the
second product have only positive powers.  Vieta gives

```text
product_(all roots gamma) gamma=(-1)^d r_0/r_d.           (9)
```

Taking the Laurent constant term of the logarithm of (8) therefore yields
the exact germ identity

```text
log(Pi(t)/(c*t))
   =sum_(m>=1) CT(f^m)t^m/m.                              (10)
```

In particular, if

```text
CT(f)=CT(f^2)=...=CT(f^L)=0,                              (11)
```

then exponentiating (10) gives

```text
Pi(t)-c*t=O(t^(L+2)).                                     (12)
```

The extra power of `t` in (12) is important: (10) controls the ratio
`Pi/(ct)`.

## 3. The compound polynomial and its pole order

Let `Omega` be the `d` roots of `Phi` in a splitting field over `C(t)`.  For
each `a`-subset `S` of `Omega`, define

```text
p_S=product_(alpha in S) alpha.
```

There are `C=binom(d,a)` such products.  Their compound polynomial is

```text
H(Y,t)=product_(|S|=a)(Y-p_S)
      =sum_(j=0)^C h_j(t)Y^(C-j).                         (13)
```

Each `h_j` is symmetric in the roots, so `H` belongs to `C(t)[Y]`.  We now
bound its denominators exactly enough for effectivity.

Divide (5) by its leading coefficient `-tr_d`.  The resulting monic root
polynomial is

```text
R(X)/r_d-X^a/(tr_d).                                      (14)
```

If `e_i` denotes the `i`th elementary symmetric function of the `d` roots,
then every `e_i` is constant in `t` except

```text
e_b=A+B/t,       B!=0,                                    (15)
```

because `e_b` corresponds to the coefficient of
`X^(d-b)=X^a` in (14).

The coefficient `h_j` is a symmetric homogeneous polynomial of total root
degree `ja`: it is a sum of products of `j` distinct `a`-subset products.
By the fundamental theorem of symmetric polynomials, it has a unique
weighted-homogeneous expression in `e_1,...,e_d`, with `e_i` assigned weight
`i`.  Any monomial in that expression containing `e_b^s` has

```text
s*b<=j*a.
```

Using (15), this proves

```text
h_j(t) in C[t^(-1)],
pole_order_0(h_j)<=floor(ja/b).                           (16)
```

Every root occurs in exactly

```text
eta=binom(d-1,a-1)=C*a/d                                  (17)
```

of the `a`-subsets, so in particular

```text
h_C(t)=(-1)^C e_d^eta in C*.                              (18)
```

There is a stronger use of complements.  Let `H_dual` be the analogous
compound polynomial of all `b`-subset root products, and write its
coefficients as `h_dual_l`.  Complementation is a bijection between the
`a`-subsets and the `b`-subsets.  If `S_dual=Omega\S`, then

```text
p_S^(-1)=p_(S_dual)/e_d.
```

The elementary identity

```text
E_j(x_1,...,x_C)
  =(product_i x_i) E_(C-j)(x_1^(-1),...,x_C^(-1))
```

therefore gives, up to the displayed harmless sign,

```text
h_j=(-1)^C e_d^(eta-C+j) h_dual_(C-j).                   (19)
```

The power of `e_d` is a nonzero constant in `t`, even when its written
exponent is negative.  The coefficient `h_dual_(C-j)` has total root degree
`(C-j)b`.  Since the only moving elementary symmetric coordinate is still
`e_b`, its pole order is at most `C-j`.  Combining this with (16) yields

```text
pole_order_0(h_j)<=min(floor(ja/b),C-j).                  (20)
```

Now `C-eta=Cb/d` is an integer.  For `j<=C-eta`, the first entry in the
minimum in (20) is at most `eta`; for `j>=C-eta`, the second entry is at most
`eta`.  Hence every coefficient has pole order at most

```text
K=eta=binom(d-1,a-1).                                     (21)
```

It follows that

```text
Htilde(Y,t):=t^K H(Y,t) in C[t,Y].                        (22)
```

## 4. The evaluation at the forbidden line is nonzero

Set

```text
G(t)=Htilde(c*t,t)=t^K H(c*t,t) in C[t].                 (23)
```

First, `G` is not the zero polynomial.  Here is the full orbit argument.

The polynomial `Phi` is irreducible over `C(t)`.  Indeed, in `C[X,t]` it is
linear in `t`; in any factorization one factor is independent of `t` and
must divide both `X^a` and `R(X)`.  Since `R(0)!=0`, their gcd is one.  Gauss
then gives irreducibility over `C(t)`, and characteristic zero gives
separability.  Hence the Galois group acts transitively on `Omega`.

If `G` vanished identically, then one factor of the product (13) would give
an `a`-subset `S` with

```text
p_S=c*t in C(t).                                          (24)
```

Let `O` be the Galois orbit of `S`, of size `r>0`.  Since (24) is fixed by
the Galois group, every subset in `O` has the same product `c*t`.  By
transitivity, every root lies in the same positive number `mu` of members
of `O`.  Multiplication over `O` and (9) give

```text
(c*t)^r=((-1)^d r_0/r_d)^mu.                              (25)
```

Here `mu>0` is the common orbit-incidence number; it need not equal the
full-family incidence `eta` in (17).  The two sides of (25) have respective
`t`-adic valuations `r>0` and zero,
a contradiction.  Therefore

```text
G(t)!=0.                                                  (26)
```

Second, (20)--(22) show that every coefficient
`t^K h_j(t)` has ordinary `t`-degree at most `K`.  Substituting `Y=c*t` in
(13) therefore gives

```text
deg_t G<=K+C.                                             (27)
```

No height estimate, coefficient sign, genericity, or saddle selection enters
(26)--(27).

## 5. Contact cannot exceed polynomial degree

The chosen small-root product is one of the roots in (13), so

```text
Htilde(Pi(t),t)=0.                                        (28)
```

Because `Htilde` is a polynomial in `Y` with analytic `t`-coefficients, the
difference

```text
G(t)
 =Htilde(c*t,t)-Htilde(Pi(t),t)                           (29)
```

is divisible as an analytic germ by `c*t-Pi(t)`.  Under (11), equation (12)
therefore implies

```text
ord_(t=0) G>=L+2.                                         (30)
```

On the other hand, a nonzero polynomial cannot vanish to order greater than
its degree.  Equations (26)--(27) give

```text
L+2<=K+C,
L<=K+C-2.                                                 (31)
```

Consequently the first `B=K+C-1` positive powers cannot all have zero
constant term.  This proves the theorem.

## 6. Boundaries, use, and what remains open

1. **All complex coefficients are allowed.**  The proof works in the exact
   support after zero coefficients are deleted.  Dihedral cancellation such
   as THM-2070 may kill infinitely many return lengths, but cannot postpone
   every nonzero constant term beyond (1).
2. **The one-sided boundary is excluded by hypothesis.**  Both exact extreme
   coefficients must be nonzero and `M,N>=1`.  A nonzero constant coefficient
   is detected already at `m=1`.
3. **The bound is sharp when `min(M,N)=1`.**  The binomial
   `u^(-1)+u^N` has first return `N+1=d`, so (3) cannot be improved on that
   boundary.
4. **The general estimate is deliberately crude.**  At `(M,N)=(2,2)`, for
   example, (1) gives `B=8`, whereas the conjectural sharp bound is `4`.
   Complement duality improves the one-sided weighted-degree estimate, but
   the actual compound coefficients are sparser still: in this bidegree they
   are only linear in the moving coordinate `e_2`.  Understanding that
   sparsity uniformly is a precise algebraic route toward the Sturmfels
   bound.
5. **NC2/GMC(2) paper dependency.**  On a lowest balanced Wick face with
   charge extremes `-M,N`, THM-2022 may choose its seed exponent `m0` with
   `m0<=B`.  Thus THM-2067's project-internal replacement of the external
   DvdK citation is not only logically complete; it can be made effective at
   the seed stage.  The later good prime may still depend on the algebraic
   coefficients, so this is not a coefficient-uniform bound on the final
   Gaussian moment order.

The surviving sharp problem is now concrete: replace the weighted-degree
pole bound (16) by enough compound-coordinate structure to force contact
degree at most `d+1`.  The theorem above proves a finite bound without
claiming that open sharpening.

## 7. Paper-proof dependency audit

This argument was found while auditing whether THM-2022 still imports the
Duistermaat--van der Kallen theorem under another name.  It does not.  The
paper chain can be made entirely project-internal as follows.

1. The lowest balanced face in THM-2022 has no projection collisions: equal
   charge and equal exposed height force equal `(a_i,b_i)`.  Its Laurent face
   polynomial is therefore an exact Laurent polynomial, not a quotient with
   merged coefficients.
2. If that face contains a neutral charge, its coefficient is already the
   nonzero constant term at `m=1`.  Otherwise its balanced convex hull
   contains charges of both strict signs, exactly the hypothesis of this
   theorem (or of THM-2067).
3. Equations (8)--(10) are the only analytic-germ input behind THM-1550.  The
   small product is an ordinary Weierstrass coefficient; the logarithm is
   normalized to have constant value zero; and (10) is an identity of
   convergent germs.  No asymptotic continuation or DvdK critical-value
   theorem is used.
4. Irreducibility in Section 4 uses only `gcd(X^a,R)=1`, supplied by
   `R(0)!=0`.  The selected small-root subset need not be Galois-invariant:
   if its product were the rational function `ct`, its *orbit* supplies the
   uniform-incidence norm contradiction (25).  This is precisely why no
   unproved monodromy classification is hidden in the argument.
5. THM-2022 then uses only the resulting nonzero seed, rational face height,
   ordinary algebraic specialization, Kummer/Lucas congruences, and
   Frobenius.  Those steps neither cite nor imply the stronger DvdK
   limsup/critical-value conclusion.

Thus the external DvdK dependence is genuinely removed from the **paper
proof** of NC2 and GMC(2).  What remains external only if one wants it is the
stronger published analytic theorem or historical attribution.  Lean
formalization of the root/compound argument is a separate engineering task
and is not asserted here.
