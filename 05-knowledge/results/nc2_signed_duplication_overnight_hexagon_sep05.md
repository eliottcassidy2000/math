# Signed duplication by an explicit SOS, and real-rooted Laurent return cores

Status: **PROVED ANALYTICALLY / INDEPENDENTLY AUDITED.** Root and the
independent `observer_collision` referee audited the complete written
proof, all consumers and boundaries, and the Hermite-limit corollary.
No external theorem is required for the duplication inequality. The
arithmetic-progression complex-phase consumer additionally uses the proved
repository theorem THM-4436. The Hermite specialization is explicitly
credited to classical Mehler inheritance. The finite controls are
**FINITE-EXACT**, not the proof of the unbounded statements. Canon
THM-4440 is root-owned; its promotion is coordinated separately.

## 1. The exact signed duplication identity

Let `n>=1`, `1<=k<=n` be integers, and let `r_1,...,r_n` be arbitrary
real numbers: signs and zeros are allowed. Write `e_j=e_j(r_1,...,r_n)`,
`D=binom(n,k)`, and `r_S=prod_(i in S) r_i`. Put

```text
C = binom(2n,2k)/binom(n,k)^2,
beta_s = binom(n-s,k-s)/binom(n,k),
c_(k,s) = 2^(2k-2s+1) (2s-2)! (k-s)! k!
          / [(s-1)! (2k)!],                      1<=s<=k.
```

The duplicated list means exactly
`r^dup=(r_1,r_1,r_2,r_2,...,r_n,r_n)`. Then

```text
C e_k(r)^2 - e_(2k)(r^dup)
 = sum_(s=1)^k c_(k,s)
     sum_(|S|=s) [r_S e_(k-s)(r outside S)-beta_s e_k(r)]^2.       (1)
```

Every coefficient `c_(k,s)` is strictly positive; the last is
`c_(k,k)=1/(2k-1)`. In particular,

```text
e_(2k)(r^dup)/binom(2n,2k) <= [e_k(r)/binom(n,k)]^2.             (2)
```

The constant in (2) is sharp, since all `r_i` equal gives equality.
The exact equality criterion in (1) is that all `k`-subset products
`r_I` are equal. For `1<=k<n`, this is equivalent to

```text
fewer than k of the r_i are nonzero, OR r_1=...=r_n.             (3)
```

For `k=n`, equality always holds. Zeros must not be omitted from (3).
The last square level alone gives the quantitative strengthening

```text
C e_k^2-e_(2k)(r^dup)
 >= [sum_(|I|=k) r_I^2-e_k^2/D]/(2k-1).                        (4)
```

### Inheritance, adversaries, and the live board

The closest proved mechanisms are
[THM-4436, complete factorial-row simple negative roots and phase wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
[THM-4432, two first channels with all doubling carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md),
and the [contiguous Euler/interlacing and actual-mass hostiles](nc2_channel_contiguous_overnight_hexagon_sep05.md).
Individual first-channel real-rootedness does not imply cross-mass
interlacing or coprimality. The genuine actual-mass root word `PQQPQPQ`
and the positive individual Fourier-pair products in that note are retained.
The corrected near miss is replacing `CT(q^2)` by the positive Hermitian
quantity `CT(q(u)q(u^(-1)))`; those are different coefficient maps.

The least-used relevant operation is **lifting products to the full space
of k-subsets**, keeping their intersection sizes. The live board is:
complete factorial rows; real versus complex carrier roots; ordinary
coefficient duplication; subset-distance kernels; centering the constant
mode; actual Laurent congruence; coefficient-phase normalization. The new
positive mechanism concerns the real roots of an ordinary carrier in its
variable, not merely the real roots of a different channel polynomial in
the phase parameter.

Targeted searches for the exact duplication statement, square-coefficient
inequality, signed elementary-symmetric forms, subset kernels, and the
arithmetic-progression consumer found no existing repository proof of (1)
or its unbounded consumers. This is not an external-priority claim. For
comparison, [Nikolov, *On certain inequalities for real-root polynomials*
(2022), equation (1) and Theorem 2.1](https://www.math.bas.bg/smb/2022_PK/tom_2022/pdf/077-085.pdf)
concerns alternating derivative sums; the ordinary-square coefficient here
is different. Neither those inequalities nor a general named kernel
theorem is being imported into the proof below.

## 2. Elementary proof: coefficient multiplicity, centering, and squares

Index coordinates by the `D` subsets `I` of `[n]` of size `k`. Set

```text
v_I=r_I,
ell(I,J)=|I\J|=|J\I|,
d_ell=4^ell/binom(2ell,ell),
K_(I,J)=d_(ell(I,J)).                                           (5)
```

For a fixed monomial with `k-ell` squared variables and `2ell` singly
occurring variables, the number of ordered decompositions as `r_I*r_J`
is `binom(2ell,ell)`. Its weight in `v^T K v` is therefore `4^ell`.
That is also its coefficient in `e_(2k)(r^dup)`, since each singleton
variable can be chosen from either of its two copies. Consequently

```text
sum_I v_I=e_k(r),       v^T K v=e_(2k)(r^dup).                  (6)
```

Permutation symmetry makes every row sum of `K` equal. Substituting
`r_i=1` in (6) proves that the common row sum is

```text
R=binom(2n,2k)/D=C*D.                                          (7)
```

Let `w=v-(e_k/D)*1`, so `sum_I w_I=0`. Equation (7) cancels both cross
terms and the constant-mode contribution, giving

```text
C e_k^2-v^T K v=-w^T K w.                                      (8)
```

The exact increment of (5) has a positive integral representation:

```text
d_(ell+1)-d_ell
 = 4^ell (ell!)^2/(2ell+1)!
 = integral_0^1 [4t(1-t)]^ell dt.                               (9)
```

The integral equality follows, for example, by repeated integration by
parts of `t^ell(1-t)^ell`. Put `q=4t(1-t)`. Summing the finite geometric
series gives

```text
d_ell=1+integral_0^1 (1-q^ell)/(1-q) dt.                        (10)
```

At `t=1/2`, the quotient has the continuous value `ell`; there is no
singularity in (10). Define `Q(q)_(I,J)=q^(ell(I,J))`. On a vector of
coordinate sum zero, (8)--(10) yield

```text
C e_k^2-e_(2k)(r^dup)
 = integral_0^1 w^T Q(q) w/(1-q) dt.                           (11)
```

Here is the decisive finite expansion, valid including `q=0,1` with
the polynomial conventions:

```text
q^(k-|I intersect J|)
 = sum_(S subset I intersect J) q^(k-|S|)(1-q)^|S|.
```

Indeed the right side is
`q^(k-|I intersect J|) [q+(1-q)]^|I intersect J|`.
Therefore

```text
w^T Q(q)w
 = sum_(s=0)^k q^(k-s)(1-q)^s
     sum_(|S|=s) [sum_(I containing S) w_I]^2.                 (12)
```

The term `s=0` is zero, since `sum w_I=0`. Dividing the remaining terms
by `1-q` makes (11) a finite integral of nonnegative squares, with no
endpoint problem. Moreover,

```text
sum_(I containing S) w_I
 = r_S e_(k-s)(r outside S)-beta_s e_k.
```

Finally the substitution `u=|2t-1|` gives

```text
integral_0^1 q^(k-s)(1-q)^(s-1) dt
 = integral_0^1 (1-u^2)^(k-s)u^(2s-2) du
 = c_(k,s).                                                   (13)
```

The final evaluation follows by integration by parts (or expanding the
integer power); it is the positive factorial expression stated above.
This proves (1) and (4). The `s=k` squares vanish exactly when all
`v_I=e_k/D`; this condition also makes every earlier square vanish.

To finish the equality audit, suppose `k<n`. If fewer than `k` numbers
are nonzero, all products are zero. If some but not all numbers are zero
and at least `k` are nonzero, one product is nonzero and a product using a
zero is zero, so equality is impossible. If all numbers are nonzero,
compare products using any common `k-1` indices outside a chosen pair
`i,j`; such indices exist because `k<n`. Equality forces `r_i=r_j`.
This proves (3), including every zero-root boundary.

## 3. An interior zero coefficient forces a strictly negative square coefficient

Let `H` be a nonzero real-rooted real polynomial of exact degree `n`,
with `H(0)!=0`. Factor

```text
H(s)=H(0) prod_(i=1)^n (1+r_i*s),
```

where every `r_i` is real and nonzero. For `1<=k<n`, if `[s^k]H=0`,
equations (1) and (4) give

```text
[s^(2k)]H(s)^2
 <= -H(0)^2/(2k-1) * sum_(|I|=k) r_I^2 < 0.                  (14)
```

Repeated real roots are allowed. The full sharp coefficient inequality is

```text
([s^k]H)^2/binom(n,k)^2
 >= [s^(2k)]H^2/binom(2n,2k),                                (15)
```

with equality, for `1<=k<n` and `H(0)!=0`, exactly when all roots of `H`
coincide. Formula (14), not a presumed sign of every summand, is the
needed cross-row information.

The weakest interior condition for the strict zero-coefficient consequence
is `ord_0(H)<k<deg(H)`: factor `H=s^r G`, apply (14) to `G` at `k-r`,
and shift the squared coefficient by `2r`. Thus a nonzero constant term
is convenient, but not necessary. If the coefficient is zero outside
this interior range, the doubled coefficient is zero as well. For example,
`H=s^2,k=1` has both coefficients zero. Dropping real-rootedness is fatal:
`H=1+s^2,k=1` has first coefficient zero and squared coefficient `+2`.

Nor is termwise negativity valid even for real-rooted `H`:

```text
H=(s+1)(s+2)(s-10)(s-4)=s^4-11s^3+92s+80,
[s^2]H=0,
[s^0]H*[s^4]H=80>0,
[s^1]H*[s^3]H=-1012,
[s^4]H^2=2*(80-1012)=-1864.
```

The subset centering retains the exact weights that make the combined
quantity negative; it does not erase the positive individual pair.

## 4. Arbitrary-support real-rooted Laurent cores: a linear width bound

Let integers `d>=1`, `N>=1`, `0<a<d*N` be given. Let `R` be a real
polynomial of exact degree `N`, with `R(0)!=0`, all of whose roots are
real. Its coefficients need not be positive, and internal coefficients
may vanish. Define

```text
f(u)=u^(-a)R(u^d),
e=gcd(a,d),       m0=d/e,       k0=a/e.
```

The endpoint width is exactly `d*N`. At every mass `m`, literal
coefficient extraction gives

```text
CT(f^m)=0                          if d does not divide a*m,
CT(f^m)=[s^(a*m/d)]R(s)^m          otherwise.                 (16)
```

The congruent masses are exactly the positive multiples of `m0`.
At `m0`, take `H=R^m0`: it is real-rooted, has nonzero constant term,
and `0<k0<N*m0`. If its indicated coefficient is zero, (14) says

```text
CT(f^m0)=0  ==>  CT(f^(2m0))<0.                              (17)
```

There is no intermediate congruent mass. Hence the first nonzero moment
is **exactly either `m0` or `2m0`**, not merely at most some later mass.
If `N>=2`,

```text
first nonzero mass <= 2m0=2d/e <= d*N=width.                  (18)
```

If `N=1`, the binomial `R^m0` has no zero interior coefficient, so the
first nonzero mass is `m0<=d`, again within the width. The hypothesis
does not assert that the first congruent support coefficient is nonzero:
its possible vanishing is exactly what (17) handles.

The same implication holds at **every** positive congruent mass, with
`m` and `2m` replacing `m0,2m0`. It follows by applying (14) to `R^m`.

For a complex Laurent polynomial `F`, these conclusions about vanishing
also hold if there exist nonzero complex `lambda,sigma` with
`lambda F(sigma*u)=u^(-a)R(u^d)` satisfying the real hypotheses above.
Indeed `CT[(lambda F(sigma*u))^m]=lambda^m CT(F^m)`. The strictly negative
sign belongs to the declared real gauge; it is not an order statement
about an unnormalized complex number.

## 5. All-mass arithmetic-progression trinomial coprimality

Let `d>=2`, `1<=a<d`, and let `alpha,beta,gamma` be arbitrary nonzero
complex numbers. For

```text
F(u)=alpha*u^(-a)+beta*u^(d-a)+gamma*u^(2d-a),
tau=alpha*gamma/beta^2,
m0=d/gcd(a,d),
```

the support-return masses are exactly the multiples of `m0`. For a
congruent mass `m`, put `k=a*m/d`, so `1<=k<m`. The complete channel row is

```text
P_m(t)=sum_(j=0)^floor(k/2)
       m! t^j / [(m-k+j)!(k-2j)!j!],
CT(F^m)=alpha^(m-k)*beta^k * P_m(tau).                        (19)
```

All roots of each nonconstant `P_m` are strictly negative, by
[THM-4436](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md)
with `A=1,B=2,x=m-k,z=0,h=floor(k/2),r=k mod 2`. The constant term of
`P_m` is positive, so zero is not a root either.

For every real `t<0`, the ordinary quadratic

```text
R_t(s)=1+s+t*s^2
```

has two distinct real roots of opposite signs. Moreover
`P_m(t)=[s^k]R_t(s)^m` and
`P_(2m)(t)=[s^(2k)]R_t(s)^(2m)`. Therefore (14) proves

```text
P_m(t)=0  ==>  P_(2m)(t)<0,                                 (20)
```

at every root of `P_m`. It follows that **`P_m` and `P_(2m)` are coprime
over C for every congruent positive mass `m`**, without a bound on their
degrees or the number of channels. Constant rows cause no exception.
The identification in (19) is essential: the doubled Laurent moment
uses the square of `R_t^m`, not the square of the channel polynomial
`P_m(t)`.

For completeness, the real gauge is exact for arbitrary complex
coefficients. Choose `sigma^d=alpha/beta` and
`lambda=sigma^a/alpha`. Then

```text
lambda F(sigma*u)=u^(-a)[1+u^d+tau*u^(2d)].
```

Thus the first nonzero moment is exactly `m0` unless `P_m0(tau)=0`,
in which case it is exactly `2m0`. In particular the width bound
`first mass<=2d=width` holds for every such complex trinomial.
Both alternatives are attained precisely when `floor((a/gcd(a,d))/2)>=1`;
otherwise the first row is constant and only `m0` occurs. Here "both"
means among nonzero coefficient choices for the fixed support.

### Unbounded multiplicity and width-sharpness

For every integer `h>=1`, take

```text
a=2h,       d=2h+1,
support=(-2h,1,2h+2),       m0=2h+1,
P_m0(t)=sum_(j=0)^h (2h+1)! t^j / [(1+j)!(2h-2j)!j!].        (21)
```

This has exactly `h` distinct strictly negative roots. Each gives actual
nonzero coefficients `(alpha,beta,gamma)=(1,1,t)` whose first detection
is `2m0=4h+2`, exactly the endpoint width. The first-return fibre has
`h+1` channels, unbounded. For example, at `h=5`,

```text
support=(-10,1,12),       first mass=11,       width=22,
P_11(t)=11+495t+4620t^2+11550t^3+6930t^4+462t^5.
```

Each of its five negative roots gives first nonzero moment exactly 22.
For `h>=5`, these examples lie outside the incoming min-width-eight
strip and have more than two first channels. The family is not presented
as the first known width-sharp example; its contribution is the
unbounded-channel analytic classification (20)--(21).

## 6. The exact boundary of the real-rooted-core method for trinomials

The core hypothesis is a genuine restriction, not another description
of the general negative phase wall. The primitive support

```text
f=u^(-3)+u-u^9,
R(s)=1+s-s^3,
CT(f^4)=[s^3]R^4=0,
CT(f^8)=[s^6]R^8=-224
```

has compressed cubic discriminant `-23`, so its core is not real-rooted.
Its noncancellation is already explained by THM-4432; the present
core proof does not cover it. It is a scope hostile, not a counterexample
to the desired general two-rung conclusion.

There is a precise small surviving extra sector. A real-rooted polynomial
with nonzero constant term cannot have two consecutive internal zero
coefficients. To see this, for any nonconstant real-rooted polynomial
`G` and real `x` with `G(x)!=0`,

```text
(G'/G)'(x)=-sum_roots multiplicity/(x-root)^2<0.
```

Thus `G'(x)=G''(x)=0` forces `G(x)=0`. Apply this successively to
derivatives to propagate two consecutive coefficient zeros down to the
constant term, a contradiction. All relevant derivatives have real
roots by Rolle's theorem.

Consequently a primitive compressed trinomial
`R=alpha+beta*s^A+gamma*s^B`, `0<A<B`, `gcd(A,B)=1`, with real nonzero
coefficients and real roots must have `A<=2` and `B-A<=2`.
Therefore `B<=3`; the sole possible `B=4,A=2` is not primitive.
The real-rooted-core route cannot directly cover primitive `B>=4`.

For `B=3`, in the negative phase gauge, the precise real-rooted sector is

```text
-4/27 <= tau < 0.                                             (22)
```

For `A=1`, normalize the core to `1+s+tau*s^3`; its discriminant is
`-tau*(4+27tau)`. For `A=2`, normalize to `1-v*s^2+s^3`, `v>0`, with
`tau=-1/v^3`; its discriminant is `4v^3-27`. A real cubic has all roots
real, with repetitions allowed, exactly when its discriminant is
nonnegative, as also follows directly from its two extrema. These
calculations give (22) in both cases. The negative sign in (14) remains
strict at the repeated-root endpoint. Thus (17) also certifies every
congruent cubic-core moment in this sector. The rest of the general
trinomial negative ray remains OPEN.

## 7. The correctly scaled Hermite limit has an exact negative-square identity

This section recovers a classical specialization, not a new literature
identity. The primary reference is [DLMF18.18.28, Mehler's formula](https://dlmf.nist.gov/18.18#E28),
with [DLMF18.18.23, Hermite linearization](https://dlmf.nist.gov/18.18#E23)
as a related coefficient identity. Both formulas and their normalizations
were checked directly. The new connection here is the self-contained
finite-subset SOS limit and its precise all-root sign consumer.

Define the probabilists' Hermite polynomials by the formal generating
function

```text
exp(x*z-z^2/2)=sum_(j>=0) He_j(x)*z^j/j!.
```

For every integer `k>=1` and real `x`, the full limiting identity is

```text
He_(2k)(sqrt(2)*x)
 = 2^k He_k(x)^2
   - sum_(s=1)^k binom(k,s) 2^(k-s) (2s-3)!! He_(k-s)(x)^2,   (23)
```

where `(-1)!!=1`. In particular,

```text
He_k(x)=0  ==>  He_(2k)(sqrt(2)*x) <= -(2k-3)!! < 0.           (24)
```

To see the classical inheritance explicitly, diagonal Mehler in this
normalization and elementary binomial expansion give respectively

```text
sum_(j>=0) 2^j He_j(x)^2 z^j/j!
 = (1-4z^2)^(-1/2) exp[2x^2*z/(1+2z)],
sum_(j>=0) He_(2j)(sqrt(2)*x) z^j/j!
 = (1+2z)^(-1/2) exp[2x^2*z/(1+2z)].
```

Multiply the first series by
`sqrt(1-2z)=1-sum_(s>=1)(2s-3)!! z^s/s!` and compare coefficients
to recover (23). These are also formal series identities, so no endpoint
convergence is being assumed. The independent SOS derivation follows.

The factor `sqrt(2)` is essential: it is the actual scale induced by
squaring the Gaussian generating function, not an arbitrary comparison
of `He_k` and `He_(2k)` at the same argument. Formula (24) answers the
scaled Hermite hostile suggested by large-offset arithmetic-progression
rows; no unproved asymptotic sign inference is needed.

Here is a self-contained coefficient limit from (1). For a fixed real
`x`, use the real-rooted polynomial approximants

```text
H_L(z)=(1+x*z/L)^L (1-z^2/(2L))^L.
```

Regard these as products of `n=3L` factors `1+r_i*z`, allowing `r_i=0`
when `x=0`. Their reciprocal-root list consists of `L` copies of `x/L`
and `L` copies each of `1/sqrt(2L)` and `-1/sqrt(2L)`. Thus

```text
max_i |r_i| -> 0,       sum_i r_i^2 = 1+x^2/L -> 1.
```

For every fixed coefficient index, the binomial formulas show directly
that `H_L` converges coefficientwise to `exp(x*z-z^2/2)`. Write its
limiting coefficients as `a_j=He_j(x)/j!`. Then the ordinary squared
coefficient converges to
`2^k He_(2k)(sqrt(2)*x)/(2k)!`, and

```text
binom(2n,2k)/binom(n,k)^2 -> 4^k/binom(2k,k).
```

For a fixed `s`, the nonnegative elementary sum `e_s(r_1^2,...,r_n^2)`
tends to `1/s!`. Indeed the difference between
`(sum r_i^2)^s` and `s! e_s(r_i^2)` consists of tuples with a repeated
index, bounded by a constant depending only on `s` times
`max r_i^2 * (sum r_i^2)^(s-1)`, which tends to zero.

Removing any fixed `s` factors from `H_L` changes each fixed coefficient
by `o(1)`, uniformly in the removed indices. This follows by finite
formal division by `prod_(i in S)(1+r_i*z)`: the coefficients of its
reciprocal through any fixed order tend uniformly to those of 1.
Consequently

```text
sum_(|S|=s) r_S^2 e_(k-s)(r outside S)^2 -> a_(k-s)^2/s!.
```

The centering terms in (1) vanish in this limit. More explicitly,
`sum_(|S|=s) r_S e_(k-s)(r outside S)=binom(k,s)e_k`,
`beta_s=O(n^(-s))`, and `beta_s^2 binom(n,s)=O(n^(-s))`.
All relevant fixed coefficients stay bounded. Passing to the limit in
the finite sum (1) therefore gives

```text
[4^k/binom(2k,k)] a_k^2
 - 2^k He_(2k)(sqrt(2)*x)/(2k)!
 = sum_(s=1)^k c_(k,s) a_(k-s)^2/s!.
```

Substituting the factorial expression for `c_(k,s)` proves (23).
Its final square is the nonzero constant `(2k-3)!!`, proving (24).
The weaker estimate (24) also follows immediately by passing to the
limit in the terminal bound (4).

All Hermite roots are real and simple: differentiating the generating
function gives `He_(j+1)=x*He_j-j*He_(j-1)`, with `He_0=1,He_1=x`.
The recurrence first excludes common roots of consecutive terms. At
the roots of `He_j`, the signs of `He_(j+1)=-j*He_(j-1)` alternate;
the two exterior signs and the monic leading coefficient complete the
usual elementary interlacing induction. Thus (24) applies to every
Hermite root, and also proves that `He_k(x)` and
`He_(2k)(sqrt(2)*x)` are coprime as polynomials over C. For `k>=2`
the first inequality in (24) is strict: the lower Hermite squares cannot
all vanish (for `k=2` use `He_1` at roots of `He_2`; for `k>=3` use
`He_1=x` and `He_2=x^2-1`). No best-constant claim is made for that
quantitative bound at `k>=2`.

## 8. Incoming connection and limits of the advance

The fresh incoming
[independent path / full-support Hadamard proof and two-row sign bank](overnight3_20260906_moments_root_geometry.md)
was read in full. It expressly recovers the same all-row simple-negative
conclusion as THM-4436, independently, and leaves the unbounded doubled-row
sign problem open. Its complete auxiliary support is essential: the
valid `(-5,2,9)` row has a non-real-rooted truncated auxiliary prefix
`21+20t+5t^2`, although the full shifted product gives the correctly
real-rooted first row `21+140t+105t^2`.

The present map is different: it takes the actual ordinary carrier
`R^m`, lifts its real reciprocal roots to **all** `k`-subset products,
and preserves the ordinary-square coefficient through exact ordered-pair
multiplicities. It then centers the constant subset mode and retains all
incidence levels `s=1,...,k` in (1). Deleting those weights would lose
the signed target. No independent network coupling or arbitrary-row
coprimality is inferred.

Equation (20) proves the `(A,B)=(1,2)` part of the incoming two-row sign
bank for all heights, all offsets, all admissible masses, and both possible
right-carry states; `A=1` has no left carry. Equation (22) covers an
additional cubic phase sector. It does not prove the other sampled
support classes in that bank. Similarly it covers an infinite strip
outside the [incoming width-eight theorem](overnight_20260906_moments_width8.md),
while that theorem covers many non-real-rooted-core supports outside the
present scope. Neither broader conclusion subsumes the other.

The square-coefficient proof does not use real-rootedness of a channel
polynomial as a substitute for real-rootedness of `R^m`. The general
unbounded first-root/doubled-row negativity target for primitive `B>=3`
outside (22), the general trinomial sharp bound, and the general
Laurent min-width-at-least-three sharp bound remain OPEN.

## 9. Exact checker and audit manifest

The [standalone checker](../../04-computation/nc2_signed_duplication_overnight_hexagon_sep05.py)
uses only standard-library integer and rational arithmetic. Its
[stored output](nc2_signed_duplication_overnight_hexagon_sep05.out)
is generated by:

```text
python3 04-computation/nc2_signed_duplication_overnight_hexagon_sep05.py
python3 -O 04-computation/nc2_signed_duplication_overnight_hexagon_sep05.py
```

The declared finite universes are independent of a favorable sampler:

- All 41 distance increments through `ell=40`, integrated by expanded
  polynomials; all 528 SOS weights for `1<=s<=k<=32`.
- All 4,353 Gram entries by allowed intersection size for `1<=n<=35`,
  every `1<=k<=n`; independent Johnson row counts also check the constant.
- All 3,829 interior `(root multiset,k)` rows for `n=2..7` and roots in
  `{-2,-1,0,1,2}`, including exactly 546 equality cases.
- 140 direct polynomial/kernel/SOS comparisons for `n=2..8`, all interior
  `k`, and five specified signed root lists, not just an eigenvalue check.
- 2,913 exact rational-root carriers tuned to an interior zero coefficient
  at `n=3..16`; the strict sign and terminal-SOS quantitative bound are
  checked independently with integer coefficients.
- 360 arithmetic-progression rows: every `d=2..16`, `1<=a<d`, and the
  first three admissible masses. Literal negative/positive charge counts
  reconstruct both full rows; an independent primitive integer
  pseudoremainder algorithm proves each sampled gcd is constant.
- 558 literal arbitrary-support real-rooted-core cases, degree `1..6`,
  three declared real root lists, `d=1..4`, every `0<a<d*N`; 86 have first
  cancellation. Direct Laurent multiplication and ordinary multiplication
  agree, including all noncongruent intermediate masses.
- All sharp-family rows `h=1..20`, and explicit positive/hostile controls
  for nonreal roots, zero-root equality, mixed pair signs, and nontrivial
  versus constant polynomial gcds.
- All scaled Hermite polynomial identities (23) for `k=1..32`, with
  positive square coefficients and constant gcds. The generating-function
  coefficient formula independently checks the recurrence through degree64.

The producer raises explicit exceptions at every failed gate; optimization
cannot disable checks. Normal and optimized runs are byte-identical and
pass **29,613 active gates**. Root and `observer_collision` independently
audited the full signed inequality, equality cases, actual core consumer,
AP coprimality/sharpness, and cubic method boundary. Both additionally
audited the complete Hermite limit and both classical generating-function
normalizations. Source comparisons were checked directly in the primary
Nikolov and DLMF documents. No large finite search is a proof dependency.

```text
script SHA-256:
654cabc326769d3852f35524ce143444f4c03aad7edf52fffa5b9bb4049a6c96
normal/optimized/stored output SHA-256:
60ee438bc8dd80355b4ad16a17e210d7a64f634b32951015b3bf4be0135e519a
internal exact-control trace SHA-256:
8cf72b05bb0c4c3123a7e1727afdbc0af0e6f32f9ece04be46cba6abcbd53503
```
