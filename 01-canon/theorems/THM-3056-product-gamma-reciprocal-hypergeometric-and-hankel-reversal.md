---
id: THM-3056
title: "Product-Gamma reciprocal hypergeometric dual and Hankel reversal"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For a J-factor product-Gamma moment sequence w_n=c^n product_j(alpha_j)_n,
  its coefficientwise reciprocal has generating function
  _1F_J(1;alpha_1,...,alpha_J;z/c), exact entire order 1/J and type
  J*c^(-1/J).  It is strictly log-concave, so its first Hankel determinant is
  negative and reciprocalization exits the Hamburger/Stieltjes moment cone.
  For every J, every generalized Hankel minor of order h has strict sign
  (-1)^(h choose 2); a nested-prefix-polynomial reduction factors its
  column-reversal through a generalized Vandermonde matrix and a totally
  nonnegative weighted-Pascal matrix.  For J=1 there is also a closed
  contiguous determinant, and the sign law extends from integer row indices
  to arbitrary increasing nonnegative real row nodes.  For THM-3047 at
  positive integer t, the width values form an integer support whose harmonic
  mass is an explicit hypergeometric value, whose support Dirichlet abscissa
  is zero, and whose prime-divisor shadow is exactly the primes not dividing
  t.
source: codex-2026-08-01-product-gamma-reciprocal-dual
audit: >
  An independent hostile audit ACCEPTED the row/column reduction, the
  descending-layer prefix ordering, the weighted-Pascal tail-Jacobi
  factorization, and the strictly positive M=T term in the Cauchy--Binet
  expansion proving all-J strict sign regularity.
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
script_sha256: 66f50dc1bb117b053908add0ce2689f8f1533930d4deaf9e40069b7d427c2690
output_sha256: 1c4c8cb29bfce0d6ef028ee5a55464bb30b50259f8cdceeed241b719ac659af1
hash_basis: LF-normalized bytes
---

# THM-3056 -- product-Gamma reciprocal hypergeometric dual and Hankel reversal

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3047 puts the formal width flag inside the strict Stieltjes cone, and
THM-3050 identifies its critical Borel order.  Coefficientwise reciprocalization
does something sharply different: it creates a hypergeometric entire function
and an extremely sparse harmonic support, but it reverses the first Hankel
wall.  Thus the reciprocal sequence is analytically excellent and
moment-theoretically hostile at the same time.  Surprisingly, the entire
reciprocal product-Gamma family still has a rigid alternating sign law at
every generalized Hankel order.

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

## 3. Every finite product is strictly Hankel sign-regular

A stronger statement survives the negative order-two minor.  For every
`J>=1`, every `h>=1`, and all strictly increasing nonnegative index sets

```text
p_1<...<p_h, q_1<...<q_h,
```

one has

```text
sign det [r_(p_i+q_j)]_(i,j=1)^h=(-1)^(h(h-1)/2).                    (11)
```

Thus every reciprocal finite product-Gamma Hankel kernel is strictly
sign-regular with the alternating column-reversal signature.  This is not a
generic closure theorem for coefficientwise products of sign-regular kernels;
it comes from a special nested-prefix factorization.

### Positive prefix-polynomial proof

The positive powers of `c` factor from rows and columns, so they do not affect
the sign.  For each shape and row put

```text
x_(ell,i)=alpha_ell+p_i.
```

Factoring `product_ell 1/(alpha_ell)_(p_i)` from row `i` leaves

```text
product_ell 1/(x_(ell,i))_(q_j).
```

Let `C=q_h`.  Multiplication of row `i` by the positive number
`product_ell(x_(ell,i))_C` turns column `j` into the polynomial

```text
Q_j(p_i)=product_(u=q_j)^(C-1) product_(ell=1)^J
          (p_i+alpha_ell+u).                                        (12)
```

Choose `0<delta<min_ell alpha_ell` and write `z_i=p_i+delta`.  Order all the
positive numbers

```text
beta=alpha_ell+u-delta
```

first by the layer `u=C-1,C-2,...,0`, and arbitrarily but consistently within
each layer.  If

```text
P_t(z)=product_(s=1)^t(z+beta_s),
t_j=J(C-q_j),                                                        (13)
```

then `Q_j(p_i)=P_(t_j)(z_i)`.  The original `q_j` increase, so the `t_j`
strictly decrease.  Reversing the `h` columns makes their prefix lengths
strictly increase and contributes exactly `(-1)^(h(h-1)/2)` to the original
determinant.

It remains to show that the reversed evaluation determinant is positive.
Let

```text
B_(m,t)=[z^m]P_t(z).
```

The infinite upper-unitriangular coefficient matrix `B` is totally
nonnegative.  Indeed

```text
B_(m,t)=beta_t B_(m,t-1)+B_(m-1,t-1),                               (14)
```

with the out-of-range entries zero.  There is a direct finite Jacobi
factorization which fixes the network orientation.  Truncate at any
`0<=m,t<=T`, transpose, and write

```text
M_(n,k)=B_(k,n)=e_(n-k)(beta_1,...,beta_n),
L_i=I+beta_i sum_(r=i)^T E_(r,r-1).                                  (14a)
```

The prefix recurrence gives exactly

```text
M=L_T L_(T-1) ... L_1.                                               (14b)
```

Each `L_i` is a nonnegative lower-bidiagonal tail Jacobi matrix, hence is
totally nonnegative; equivalently, it is the path matrix of one oriented
layer with a stay edge of weight one and a downward edge of weight `beta_i`.
Cauchy--Binet makes their ordered product `M`, and therefore `B=M^T`, totally
nonnegative.  This factorization is also the precise weighted-Pascal planar
network behind `(14)`.

Now let `E_(i,m)=z_i^m`.  Every minor of `E` with increasing rows and degrees
is a positive generalized Vandermonde determinant because
`0<z_1<...<z_h`.  The reversed matrix factors as

```text
[P_(t_j)(z_i)] = E B[:,T].
```

where `T={t_1<...<t_h}` after reversal.  Cauchy--Binet makes its determinant
a sum of nonnegative terms.  It is strictly positive: take the degree set
`M=T`.  The selected coefficient minor `B[T,T]` is upper triangular with
diagonal one, while the generalized Vandermonde minor `E[:,T]` is positive.
Undoing the column reversal proves `(11)`.

### Mixed continuous--discrete extension

Nothing in the proof used integrality of the row nodes.  Define, for `x>=0`
and `n` a nonnegative integer,

```text
K(x,n)=c^(-(x+n)) product_(ell=1)^J
       Gamma(alpha_ell)/Gamma(alpha_ell+x+n).                        (15a)
```

For arbitrary strictly increasing real `x_1<...<x_h` in `[0,infinity)` and
strictly increasing integer `q_1<...<q_h`, positive row and column factors
reduce `[K(x_i,q_j)]` to

```text
[product_ell 1/(alpha_ell+x_i)_(q_j)].
```

The same prefix-polynomial argument therefore gives

```text
sign det[K(x_i,q_j)]=(-1)^(h(h-1)/2).                               (15b)
```

This half-discrete kernel is the natural interpolation of the reciprocal
sequence on one coordinate.  Allowing arbitrary real nodes on **both**
coordinates is not proved here: the integer prefix lengths in `(13)` are the
mechanism, so erasing that coordinate would require a different argument.

### Closed one-factor contiguous determinant

For `J=1`, consecutive rows and columns also admit the closed formula

```text
det [r_(m+i+j)]_(0<=i,j<h)
 =(-1)^(h(h-1)/2) c^(-hm-h(h-1))
   * product_(j=1)^(h-1) j!
   / product_(i=0)^(h-1) (alpha)_(m+i+h-1).                         (15c)
```

After the same positive row clearing, this is the ordinary
falling-factorial alternant; its Vandermonde product gives the displayed
evaluation.  This is the exact reciprocal-Pochhammer determinant sought by
the original Hankel scout.

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

## 6. Why the multi-factor closure is special

For `J>1`, `(r_n)` is the coefficientwise product of `J` single-factor
reciprocal-Pochhammer sequences.  The proof of `(11)` deliberately does not
invoke a generic Hadamard-product closure: the available Stieltjes/Hankel
multiplier theorems apply to positive Hankel kernels, while these factors are
sign-regular and already indefinite.  Such a generic closure would be false
without further hypotheses.

The actual preserved object is narrower and more informative.  After positive
row and column scalings and one global column reversal, every finite minor is
an evaluation matrix of **nested polynomials whose zeros are negative and
whose new linear factors arrive in positive layers**.  The weighted-Pascal
coefficient matrix remembers the nesting sidecar that an undifferentiated
Hadamard-product description loses.  THM-3051 and THM-3053 make the same kind
of distinction elsewhere: transport works because a concrete positive kernel
or prefix structure is retained, not merely because each factor separately
has a sign label.

There is also a precise Maclaurin--Gregory--Newton bridge.  The polynomials
`P_t` are the ordered block-Newton/falling-factorial basis forced by the Gamma
increments, while `B` changes that basis to Maclaurin monomials.  The tail
Jacobi factors in `(14b)` say that this basis change is an ordered path
incidence map.  THM-3053's Beta--Gamma cone uses the same grammar on exponent
inventory: prefix sums are the cut coordinates of an oriented path, and
adjacent Beta edges are its elementary path moves.  Here ordered prefix
**products** force total nonnegativity; there prefix **sums** characterize
feasible multiplicative transport.  The shared mechanism is path incidence,
not a metaphorical resemblance between addition and multiplication.

This suggests a reusable extension criterion.  Any kernel reducible, by
positive diagonal scalings and reversal of its parameter order, to evaluations
`P_t(z_i)` with positive increasing nodes and nested products
`P_t=P_(t-1)(z+beta_t)`, `beta_t>0`, inherits strict total positivity after
the permutation.  The reciprocal product-Gamma family is one exact instance;
arbitrary sign-regular Hadamard products are not.

## 7. Exact companion and scope

The dependency-free rational companion checks:

- `162` hypergeometric coefficient ratios and `144` negative Hankel-two cells;
- `270` closed determinants `(15c)` and `26,460` generalized single-factor
  sign-regular minors, with elimination and Leibniz paths agreeing;
- `2,550` mixed continuous--discrete sign minors at exact rational row nodes,
  with both determinant paths agreeing;
- `8,820` minors of three weighted-Pascal coefficient matrices, all
  nonnegative;
- `135` THM-3047 integer-width cells and the exact `k=2` collision boundary;
- `5,040` prime-valuation cells for `(26)--(28)`;
- `4,250` generalized multi-factor sign minors, including the seven-factor
  `k=3` width families, with two determinant paths agreeing.

Reproduce with

```text
python 04-computation/gmc_product_gamma_reciprocal_dual_thm3056.py
python -O 04-computation/gmc_product_gamma_reciprocal_dual_thm3056.py
```

Both modes byte-match the stored eight-line transcript.  No claim here turns
the reciprocal coefficients into a moment law, a physical moving-lower
resultant, or a new GMC/NC2 proof.  Equation `(11)` is an all-order theorem for
every **finite** positive-shape product; no infinite-product limit or generic
Hadamard sign-regular closure is asserted.
