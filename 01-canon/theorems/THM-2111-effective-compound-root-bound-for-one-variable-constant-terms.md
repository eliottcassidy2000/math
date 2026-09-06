---
id: THM-2111
title: "An effective compound-root bound for one-variable constant terms"
status: >
  PROVED. Let f be a complex Laurent polynomial with exact extreme exponents
  -M<0<N and nonzero coefficients at both extremes. Put a=min(M,N),
  b=max(M,N), d=a+b, and C=binom(d,a). Then CT(f^m) is nonzero
  for some 1<=m<=C. The proof is self-contained: the logarithmic small-root
  identity converts L initial zero constant terms into contact of order L+2
  between the small-root product and c*t; the compound polynomial of all
  a-subset root products has coefficient pole order at most C-j in its
  Y^(C-j) coefficient by complementary-subset duality. Consequently its
  evaluation at c*t is already a nonzero polynomial of exact degree C;
  subtracting the chosen small-root product loses at most one order because
  the divided difference has an exact simple t-pole. In fact, if m_* is the
  first nonzero return, then ord_0 Q=m_* and its leading coefficient is an
  explicit nonzero multiple of CT(f^m_*). Beyond the small-root
  logarithmic identity the proof is Galois-free. For a=1 the bound is the
  sharp d=M+N. This supplies an unconditional effective seed bound for
  THM-2022. The separate parabolic-critical argument in
  THM-4417-width-two-laurent-first-return-parabolic-critical-bound now proves
  the sharp d bound for a=2. The general sharp d bound for a>=3 remains open;
  this compound proof continues to give the stated binomial bound.
source: codex-2026-07-22-gmc-paper-audit
depends_on: []
related:
  - THM-2101-additive-orbit-residue-dvdk-bypass
  - THM-1550-an-exact-criterion-for-the-toral-nullcone
  - THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2
  - THM-1650-newton-polygon-of-the-effective-dvdk-bound
  - THM-2067-galois-orbit-product-closes-one-variable-dvdk
  - THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation
  - THM-4417-width-two-laurent-first-return-parabolic-critical-bound
script: 04-computation/tnc_compound_root_effective_bound_codex_20260722.py
output: 05-knowledge/results/tnc_compound_root_effective_bound_codex_20260722.out
script_sha256: 41aa720913c8bffeced16384b1dd6556537c0c129686ed16a7dba927bca75f3f
output_sha256: 21ab80cf081412f089ccfa455de312b269c6c31f497440bcbb4b15d81e0ece5e
hash_basis: repository blobs with LF line endings
---

# THM-2111 -- an effective compound-root constant-term bound

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
C=binom(d,a).                                              (1)
```

> **Theorem.** There is an integer `m` with
>
> ```text
> 1<=m<=C
> ```
>
> such that `CT(f^m)!=0`.

If `a=1`, then `C=d`, and the bound is

```text
C=d=M+N.                                                   (2)
```

Thus the theorem proves the sharp Sturmfels bound on the entire
`min(M,N)=1` boundary, while giving a finite exponential-size bound in every
bidegree. This compound proof does **not** reach `m<=M+N` when
`min(M,N)>=2`. The separate argument in
[THM-4417, width-two Laurent first-return parabolic critical bound](THM-4417-width-two-laurent-first-return-parabolic-critical-bound.md)
now proves it for `min(M,N)=2`, including a stronger distinct-critical-point
bound. The remaining general sharp-width problem is `min(M,N)>=3`.

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

The second entry in (20), aligned with the exponent of `Y`, is the decisive
one.  Put `k=C-j`.  Since `h_j` is a polynomial in `t^(-1)` and has pole
order at most `k`,

```text
t^k h_j(t) in C[t],       deg_t(t^k h_j)<=k.             (21)
```

Thus evaluation on a line `Y=constant*t` clears each coefficient by its own
`Y`-power.  A uniform clearing factor is unnecessary.

## 4. The evaluation at the forbidden line is nonzero

Set

```text
Q(t)=H(c*t,t)
    =sum_(j=0)^C c^(C-j) t^(C-j) h_j(t).                 (22)
```

Equation (21) shows directly that

```text
Q(t) in C[t],       deg_t Q=C.                           (23)
```

The degree is exactly `C` because the monic term `Y^C` contributes
`c^C t^C`, while every term with `j>0` has degree at most `C-j<C`.
In particular `Q` is nonzero.  No irreducibility, Galois group, height
estimate, coefficient sign, genericity, or saddle selection enters (23).

## 5. Contact cannot exceed polynomial degree

The chosen small-root product is one of the roots in (13), so, as an identity
of algebraic germs,

```text
H(Pi(t),t)=0.                                             (24)
```

Although the coefficients of `H` may have poles at `t=0`, complement duality
shows that its divided difference along the two order-one germs has at worst
a simple pole.  Indeed, putting `k=C-j` and using (24),

```text
Q(t)=H(c*t,t)-H(Pi(t),t)
    =(c*t-Pi(t)) A(t),                                    (25)

A(t)=sum_(k=1)^C h_(C-k)(t)
       sum_(r=0)^(k-1) (c*t)^(k-1-r) Pi(t)^r.             (26)
```

Both `c*t` and `Pi(t)` vanish to order one.  The inner sum in (26) therefore
vanishes to order at least `k-1`, while (20) gives
`h_(C-k)(t)=O(t^(-k))`.  Hence

```text
A(t)=O(t^(-1)).                                           (27)
```

Under (11), equations (12), (25), and (27) imply

```text
ord_(t=0) Q>=L+1.                                         (28)
```

On the other hand, the nonzero polynomial `Q` has degree `C` by
(23).  If the first `C` constant terms vanished, (28) with `L=C` would give
`ord_0 Q>=C+1`, a contradiction.  Consequently one of the first `C` powers
has nonzero constant term.  This proves the theorem.

## 6. The compound determinant records the exact first return

The simple pole in (27) is always exact. This gives a stronger interpretation
of the same polynomial `Q`.

Let

```text
m_*=min{m>=1:CT(f^m)!=0}.                               (29)
```

Existence follows from the theorem just proved. Give the algebraic closure of
`C((t))` its additive valuation `v(t)=1`. The Newton polygon of (5) has
vertices

```text
(0,1), (a,0), (a+b,1),                                 (30)
```

so the `a` small roots have valuation `1/a` and the `b` large roots have
valuation `-1/b`. If an `a`-subset `S` contains `ell<=a-1` small roots, then

```text
v(p_S)=ell/a-(a-ell)/b
      =1-(a-ell)(1/a+1/b)<1.                           (31)
```

Thus `v(c*t-p_S)=v(p_S)` for every `S!=S_0`. Every root belongs to `eta` of
the `a`-subsets, so with

```text
E=e_d^eta=product_(|S|=a) p_S in C*,                   (32)
```

the sum of all subset-product valuations is zero. Since `v(Pi)=1`, equations
(25)--(26) give the exact factorization

```text
A(t)=product_(S!=S_0) (c*t-p_S),                        (32a)
```

and therefore

```text
v(A)=-1.                                                (33)
```

More precisely, each `ct/p_S` tends to zero for `S!=S_0`, and `Pi/(ct)` tends
to one. Therefore

```text
lim_(t->0) tA(t)=(-1)^(C-1) E/c.                       (34)
```

The logarithmic identity (10) now gives

```text
ct-Pi(t)=-c*CT(f^m_*)/m_* * t^(m_*+1)
          +O(t^(m_*+2)).                               (35)
```

Multiplying (34)--(35) in (25) proves the exact formulas

```text
ord_(t=0) Q=m_*,
[t^m_*]Q=(-1)^C E*CT(f^m_*)/m_*.                       (36)
```

Hence the first constant-term return is the `t`-adic order of the finite
coefficient polynomial `Q`, not merely bounded by its degree. If `T(t)` is a
companion matrix for the monic normalization of `Phi`, then

```text
H(Y,t)=det(Y I - exterior_power^a T(t)),                (37)
```

so (36) is also a compound-determinant certificate. It exposes a sharp no-go:
the pole in the divided-difference argument is structural, not a loose bound;
trying to prove `A=O(1)` cannot improve `C` to `a+b`.

## 7. Boundaries, use, and what remains open

1. **All complex coefficients are allowed.**  The proof works in the exact
   support after zero coefficients are deleted.  Dihedral cancellation such
   as THM-2070 may kill infinitely many return lengths, but cannot postpone
   every nonzero constant term beyond (1).
2. **The one-sided boundary is excluded by hypothesis.**  Both exact extreme
   coefficients must be nonzero and `M,N>=1`.  A nonzero constant coefficient
   is detected already at `m=1`.
3. **The bound is sharp when `min(M,N)=1`.**  The binomial
   `u^(-1)+u^N` has first return `N+1=d`, so (2) cannot be improved on that
   boundary.
4. **The general estimate is deliberately crude.**  At `(M,N)=(2,2)`, for
   example, (1) gives `C=6`, whereas the sharp bound is `4`.
   [THM-4417, width-two Laurent first-return parabolic critical bound](THM-4417-width-two-laurent-first-return-parabolic-critical-bound.md)
   proves the width bound on the entire `min(M,N)=2` boundary by converting
   small-root contact to parabolic petal cycles; it does not change this
   compound-family proof.
   Complement duality removes the earlier uniform pole-clearing loss, but the
   compound family still has `C` members.  Replacing that family-size degree by
   the original root degree `d`, or proving that contact beyond `d+1` forces
   a lower-dimensional degeneration, is the remaining algebraic route toward
   the Sturmfels bound.
5. **NC2/GMC(2) paper dependency.**  On a lowest balanced Wick face with
   charge extremes `-M,N`, THM-2022 may choose its seed exponent `m0` with
   `m0<=C` directly from this theorem.  Beyond the small-root logarithmic
   identity, the effective proof is Galois-free; THM-2067 is a related
   historical route to bare existence, not a dependency.  The later good
   prime may still depend on the algebraic coefficients, so this is not a
   coefficient-uniform bound on the final Gaussian moment order.

For `min(M,N)>=3`, the surviving general sharp problem is concrete: replace
the full compound-family degree `C` by enough structure to force contact
degree at most `d`. The theorem above proves a finite bound without claiming
that remaining open sharpening.

## 8. Paper-proof dependency audit

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
4. Complementary-subset duality is a symmetric-polynomial identity.  It
   clears the forbidden line coefficient by coefficient, and the unique
   leading term `c^C t^C` makes `Q` nonzero.  Thus the effective proof uses
   neither irreducibility nor a Galois-orbit or monodromy classification.
5. THM-2022 then uses only the resulting nonzero seed, rational face height,
   ordinary algebraic specialization, Kummer/Lucas congruences, and
   Frobenius.  Those steps neither cite nor imply the stronger DvdK
   limsup/critical-value conclusion.

Thus the external DvdK dependence is genuinely removed from the **paper
proof** of NC2 and GMC(2).  What remains external only if one wants it is the
stronger published analytic theorem or historical attribution.  Lean
formalization of the root/compound argument is a separate engineering task
and is not asserted here.
