---
id: THM-3056
title: "Product-Gamma reciprocal hypergeometric dual and Hankel reversal"
status: >
  PROVED + VERIFIED-EXACT, with one explicitly OPEN FINITE-EXACT extension.
  For a J-factor product-Gamma moment sequence w_n=c^n product_j(alpha_j)_n,
  its coefficientwise reciprocal has generating function
  _1F_J(1;alpha_1,...,alpha_J;z/c), exact entire order 1/J and type
  J*c^(-1/J).  It is strictly log-concave, so its first Hankel determinant is
  negative and reciprocalization exits the Hamburger/Stieltjes moment cone.
  For J=1, every generalized Hankel minor has strict sign
  (-1)^(r choose 2), with a closed contiguous determinant.  For THM-3047 at
  positive integer t, the width values form an integer support whose harmonic
  mass is an explicit hypergeometric value, whose support Dirichlet abscissa
  is zero, and whose prime-divisor shadow is exactly the primes not dividing
  t.  All-order sign regularity for J>1 is FINITE-EXACT through the frozen
  scout only and remains OPEN.
source: codex-2026-08-01-product-gamma-reciprocal-dual
depends_on:
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3050-rational-product-gamma-radial-nullcone-and-critical-borel-order
related:
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary
  - THM-3053-beta-gamma-prefix-transport-and-multiplicative-holotopy-cone
script: 04-computation/gmc_product_gamma_reciprocal_dual_thm3056.py
output: 05-knowledge/results/gmc_product_gamma_reciprocal_dual_thm3056.out
script_sha256: 6d2c70618c4f180017c3f18f50b771ef84b98bcba5ac1ed88191a1783b3e8be1
output_sha256: 8dd3df2e852a7acaa8dc92b4b264bdc143c5dbdea59c148edbcb8cdbf7b8b734
hash_basis: LF-normalized bytes
---

# THM-3056 -- product-Gamma reciprocal hypergeometric dual and Hankel reversal

**PROVED + VERIFIED-EXACT, with the `J>1` sign-regular extension explicitly
OPEN / FINITE-EXACT.**

THM-3047 puts the formal width flag inside the strict Stieltjes cone, and
THM-3050 identifies its critical Borel order.  Coefficientwise reciprocalization
does something sharply different: it creates a hypergeometric entire function
and an extremely sparse harmonic support, but it reverses the first Hankel
wall.  Thus the reciprocal sequence is analytically excellent and
moment-theoretically hostile at the same time.

## 1. The reciprocal product-Gamma transform

Fix `J>=1`, `c>0`, and positive real shapes

```text
alpha_1,...,alpha_J>0.
```

Put

```text
w_n=c^n product_(j=1)^J (alpha_j)_n,
r_n=1/w_n.                                                           (1)
```

Then the ordinary generating function of the reciprocal coefficients is

```text
R(z)=sum_(n>=0) r_n z^n
    =_1F_J(1;alpha_1,...,alpha_J;z/c).                               (2)
```

Indeed the `n`th hypergeometric coefficient is

```text
(1)_n/(product_j(alpha_j)_n) * (z/c)^n/n!
 =z^n/[c^n product_j(alpha_j)_n].                                   (3)
```

The coefficient ratio

```text
r_(n+1)/r_n=1/[c product_j(n+alpha_j)]                               (4)
```

tends to zero like `1/(c n^J)`, so `(2)` is entire.  Stirling gives

```text
log(1/r_n)=J n log n+n log c-Jn+O(log n).                            (5)
```

The coefficient formulas for order and type therefore yield

```text
rho(R)=limsup n log n/log(1/r_n)=1/J,                                (6)

tau(R)=1/(e rho) limsup n r_n^(rho/n)=J c^(-1/J).                    (7)
```

Both are exact.  THM-3050 proves that `(w_n)` has critical Borel order `J`;
equation `(6)` gives the reciprocal identity

```text
(critical Borel order of w) * (entire order of R)=1.                 (8)
```

This is an order duality between two coefficient transforms, not a claim that
the corresponding functions are compositional or multiplicative inverses.

## 2. Reciprocalization crosses the first Hankel wall

Equation `(4)` is strictly decreasing in `n`.  Hence `(r_n)` is strictly
log-concave:

```text
r_n^2>r_(n-1)r_(n+1), n>=1.                                         (9)
```

In particular,

```text
det [[r_0,r_1],[r_1,r_2]]
 =r_0r_2-r_1^2<0.                                                    (10)
```

So `(r_n)` is not even a Hamburger moment sequence, much less a Stieltjes
moment sequence.  This is the exact boundary with THM-3047: the original
`(w_n)` has every generalized Hankel minor positive, while its coefficientwise
reciprocal fails at order two.  The map preserves positivity, support, and
hypergeometric holonomy; it destroys positive-semidefinite Hankel geometry.

## 3. One Gamma factor is strictly Hankel sign-regular

For `J=1`, write `alpha=alpha_1`.  A stronger statement survives the negative
minor.  For every `h>=1` and all strictly increasing nonnegative index sets

```text
p_1<...<p_h, q_1<...<q_h,
```

one has

```text
sign det [r_(p_i+q_j)]_(i,j=1)^h=(-1)^(h(h-1)/2).                    (11)
```

Thus the reciprocal single-Gamma Hankel kernel is strictly sign-regular with
the alternating reversal signature.

### Falling-factorial proof

The positive scale `c` factors from rows and columns, so it does not affect
the sign.  Factor `1/(alpha)_(p_i)` from row `i` and put

```text
x_i=alpha+p_i.
```

The remaining matrix is `[1/(x_i)_(q_j)]`.  Let `C=q_h`.  Multiplication of
row `i` by the positive number `(x_i)_C` turns column `j` into

```text
(x_i+q_j)_(C-q_j).
```

With `y_i=x_i+C-1` and `d_j=C-q_j`, this is the falling-factorial alternant

```text
det [y_i falling_(d_j)].                                             (12)
```

The degrees satisfy `d_1>...>d_h=0`.  Put

```text
lambda_j=d_j-(h-j),                                                   (13)
```

which is a partition.  Dividing `(12)` by the base alternant
`det[y_i falling_(h-j)]` gives its factorial Schur polynomial.  Its tableau
formula is strictly positive here: `p_i>=i-1` and `y_i=alpha+p_i+C-1`, while
every shift occurring in the tableau is at most `C-1`.  Every tableau factor
is therefore positive.  The base alternant is

```text
(-1)^(h(h-1)/2) product_(i<j)(y_j-y_i),                              (14)
```

which proves `(11)`.

For consecutive rows and columns, `(12)` has `lambda=0` and gives the closed
formula

```text
det [r_(m+i+j)]_(0<=i,j<h)
 =(-1)^(h(h-1)/2) c^(-hm-h(h-1))
   * product_(j=1)^(h-1) j!
   / product_(i=0)^(h-1) (alpha)_(m+i+h-1).                          (15)
```

This is the exact reciprocal-Pochhammer determinant evaluation sought by the
Hankel scout.

## 4. Integer formal-width supports

Now specialize to THM-3047.  For `k>=2`, put

```text
A=k!(H_k-1), B=k!(k+1-2H_k), I=A+B=k!(k-H_k),                        (16)

F_M=product_(s=1)^(M-1)(1+st)^I (1+Mt)^B.                           (17)
```

Take a positive integer `t`.  Then every `F_M` is a positive integer.  With
`a=1/t`, `(2)` becomes

```text
sum_(M>=0) z^M/F_M
 =_1F_I(1; a repeated A times, a+1 repeated B times; z/t^I).         (18)
```

The quotient

```text
F_(M+1)/F_M=(1+Mt)^A(1+(M+1)t)^B                                   (19)
```

is greater than one everywhere except at `M=0,k=2`, where it equals one.
Therefore the value support

```text
H_(k,t)={F_M:M>=0}                                                   (20)
```

has no repetitions when `k>=3`, and has exactly the duplicate
`F_0=F_1=1` when `k=2`.  Its subset-harmonic mass is consequently

```text
sum_(h in H_(k,t))1/h
 =_1F_I(1;a^A,(a+1)^B;t^(-I))-1_(k=2).                              (21)
```

At `k=2,t=1`, `(21)` is `e-1`, exactly THM-2000's factorial-support row.
The subtracted one is the complete collision tax; it must not be discarded
when one moves between indexed and support semantics.

Stirling sharpens the Abel--Dini location.  One has

```text
F_M=C_(k,t) t^(IM)(M!)^I M^(I/t-A)(1+O(1/M)),                       (22)

C_(k,t)=1/[Gamma(1/t)^A Gamma(1/t+1)^B].                            (23)
```

If `N_(k,t)(x)=#{h in H_(k,t):h<=x}`, inversion of `(22)` gives

```text
N_(k,t)(x)~log x/[I log log x].                                     (24)
```

Thus the support lies far beyond the Bertrand boundary.  Its support
Dirichlet series

```text
D_(k,t)(s)=sum_(h in H_(k,t))h^(-s)                                 (25)
```

has exact abscissa of absolute convergence `Re(s)=0`: it converges for
`Re(s)>0`; at `Re(s)=0` its absolute terms equal one; and for `Re(s)<0` they
grow.

Since `(21)` converges, THM-2438 applies without a tail caveat.  The extended
limiting mean weak and strict multiplicative-fibre losses caused by deleting
`H_(k,t)` both exist and equal the hypergeometric mass `(21)`.  This gives the
formal-width sequence an exact divisor-scar interpretation.

## 5. The prime shadow crosses back to the divergent boundary

Although `(20)` has finite reciprocal mass, its prime-divisor shadow is huge:

```text
{p prime: p divides some F_M}={p prime:p does not divide t}.          (26)
```

If `p|t`, every factor `1+st` is `1 mod p`.  If `p` does not divide `t`,
there is a unique `r_p in {1,...,p-1}` with

```text
1+r_p t=0 mod p.                                                      (27)
```

For `B>0`, this factor first occurs in `F_(r_p)`; for `B=0` (equivalently
`k=2`), it first occurs in `F_(r_p+1)`.  This proves `(26)`.

The complete multiplicity sidecar is also explicit.  For each `ell>=1`, let
`r_(p^ell)` be the unique integer in `{1,...,p^ell-1}` satisfying `(27)`
modulo `p^ell`.  Then, for `p` not dividing `t`,

```text
v_p(F_M)
 =I sum_(ell>=1) max(0,1+floor((M-1-r_(p^ell))/p^ell))
  +B v_p(1+Mt).                                                       (28)
```

Only finitely many summands are nonzero.  Formula `(28)` simply counts the
unique root of `1+st` in every complete `p^ell` residue block.

Removing the finitely many primes dividing `t` does not repair Euler's
divergent prime-reciprocal series.  Hence the prime shadow of `(20)` has
divergent reciprocal mass even though `(20)` itself has finite mass.  The map

```text
integer support -> union of its prime divisors                              (29)
```

preserves eventual prime occurrence but destroys magnitude and term
multiplicity; `(28)` is the sidecar that restores the latter.  This is an
exact example of a quotient crossing the Abel--Dini/Bertrand boundary.

## 6. Open multi-factor Hankel sign question

For `J>1`, the sequence `(r_n)` is the coefficientwise product of `J`
single-factor reciprocal-Pochhammer sequences.  Equation `(10)` still proves
the negative order-two sign, but `(11)` does **not** follow from a generic
Hadamard-product closure: the available Stieltjes/Hankel multiplier theorems
apply to positive Hankel kernels, while these factors are sign-regular and
already indefinite.  THM-3051 and THM-3053 make this transport distinction
load-bearing.

The exact companion finds the same signature `(-1)^(h choose 2)` in `4,250`
multi-factor generalized minors, including the seven-factor `k=3` width
families, with two independent determinant evaluations.  This is a
**FINITE-EXACT SCOUT, not an all-order theorem**.  The live obligation is to
prove a special reciprocal-Pochhammer Hadamard closure (for example through a
positive factorial-Schur expansion) or find its first hostile minor.

## 7. Exact companion and scope

The dependency-free rational companion checks:

- `162` hypergeometric coefficient ratios and `144` negative Hankel-two cells;
- `270` closed determinants `(15)` and `26,460` generalized single-factor
  sign-regular minors, with elimination and Leibniz paths agreeing;
- `135` THM-3047 integer-width cells and the exact `k=2` collision boundary;
- `5,040` prime-valuation cells for `(26)--(28)`;
- the separately labelled `4,250`-minor multi-factor scout.

Reproduce with

```text
python 04-computation/gmc_product_gamma_reciprocal_dual_thm3056.py
python -O 04-computation/gmc_product_gamma_reciprocal_dual_thm3056.py
```

Both modes byte-match the stored eight-line transcript.  No claim here turns
the reciprocal coefficients into a moment law, a physical moving-lower
resultant, or a new GMC/NC2 proof.  The only `J>1` all-order sign statement is
explicitly open.
