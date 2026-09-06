# Every trinomial return fibre has simple negative cancellation roots

**Status: PROVED ANALYTICALLY / INDEPENDENTLY AUDITED, relative to the
CITED finite hyperbolicity-preserver theorem.** Root and observer separately
audited the factorial splitting, Laguerre parameters, simplicity argument,
complete all-mass Laurent consumer, and its nonvacuity; root also checked
the primary Borcea--Branden statement and exact symbol hypothesis. Exact
Sturm checks are **FINITE-EXACT** corroboration. Canonical promotion and
maintained navigation are root-owned.

## 1. Statement with weaker-than-Laurent hypotheses

Let `A,B,h,r,z,x` be integers satisfying

```text
0<A<B,     h>=0,     x>=0,     0<=r<B,     0<=z<A,
delta=B-A,                  m=x+B*h+r+z.
```

Define the complete factorial-row polynomial

```text
P(t)=sum_(j=0)^h m! t^j /
      [(x+delta*j)! (B*h+r-B*j)! (z+A*j)!].          (1)
```

If `h=0`, this is a nonzero constant. If `h>=1`, it has exactly `h`
**distinct, strictly negative real roots**. Neither `gcd(A,B)=1`, positive
Laurent charges, first-return minimality, nor any no-carry assumption is
needed for this factorial statement. All parameters are integers; the
canonical residue inequalities are load-bearing in the factorization.

### Laurent consequence, at every mass

Let `a>=1`, `1<=b<c` be arbitrary integers, with no content normalization
required, and let

```text
f(u)=alpha*u^(-a)+beta*u^b+gamma*u^c,
alpha*beta*gamma!=0,
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g,
tau=alpha^(B-A)*gamma^A/beta^B.                     (2)
```

At **every** positive mass `m` with a nonempty support-return fibre,
`CT(f^m)` is a nonzero coefficient monomial times a polynomial `(1)` in
`tau`. Consequently:

- It can vanish only at finitely many strictly negative real values of
  `tau`, all simple as scalar roots.
- If `tau` is not strictly negative real, every nonempty support moment is
  nonzero. In particular the first nonzero moment is the first support
  return, and is at most `(a+b)/gcd(a,b) < a+c`.
- Each individual moment's zero set on the coefficient torus is reduced
  and smooth. Every scalar root is genuinely attainable there.

The global negative ray is only a **necessary locus** for cancellation.
Different masses can have different finite root sets on it. Simple negative
roots alone do **not** prove that two different moment polynomials are
coprime. The general sharp trinomial return problem on that remaining ray
is not closed here.

### All-charge proportional central-resonance consequence

For every pair of positive integers `p,q` and every `m>=0`, define

```text
A_m^(p,q)(kappa)=sum_(j=0)^floor(m/(p+q))
                 m! kappa^j/[(p*j)!(q*j)!(m-(p+q)*j)!].           (9)
```

When its degree is positive, this polynomial also has only distinct,
strictly negative real roots. Set `A=q`, `B=p+q`, `delta=p`, `x=z=0`,
`h=floor(m/B)` and `r=m mod B` in `(1)`. This proves `(9)` immediately,
without a coprimality requirement on `p,q`.

These are precisely the scalar factors in
[THM-2021, proportional central-resonance factorization, Section A](../../01-canon/theorems/THM-2021-gmc2-legendre-finite-recurrence-closure.md).
Its Section C locates zeros only in the symmetric primitive case by a
Legendre transform; `(9)` extends that zero-location conclusion to all
positive charges and adds simplicity. Thus its entire exceptional scalar
set `union_m {kappa:A_m^(p,q)(kappa)=0}` lies on the strictly negative
real ray. It does **not** show that a fixed negative parameter vanishes
at only finitely many levels, nor exclude common roots of different
levels. Those are different predicates, exactly as THM-2021's correction
lineage warns. Full NC2/GMC(2) is already proved elsewhere and is not
newly claimed here.

## 2. Inheritance, duplication and the recovered operation

The closest Laurent mechanism is
[THM-2023, hyper-Bessel sharp-boundary Laguerre--Polya theorem](../../01-canon/theorems/THM-2023-gmc2-hyperbessel-boundary-laguerre-polya.md):
factorial splitting there gives a positive-parameter entire hypergeometric
function with negative zeros. Its theorem is a limiting GMC(2) boundary
statement, not the present finite trinomial-row statement. Its explicit
opposite-boundary warning is retained: not every root-of-unity filter or
factorial series has the same zero geometry.

The less-used, more exact bridge is in a different problem:
[THM-2760, exact-prefix even Faber flux gcd and smooth-boundary exclusion, Section 4](../../01-canon/theorems/THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion.md).
That proof factors a finite coefficient polynomial by degree-bounded
Schur--Szego composition of Jacobi transforms, then uses simple roots to
separate a polynomial from a derivative response. The source is the
finite composition operation, not an analogy between the original objects.
Here the target factors are Laguerre transforms instead of Jacobi ones.

The coefficient map keeps every factorial weight and the complete affine
index; composition is coefficientwise multiplication **after division by
`binom(h,j)`**. It preserves hyperbolicity. It forgets relations between
different return masses, so a two-mass separation needs an additional
response, interlacing or signed-remainder sidecar.

The canonical hostile `(-13,1,8)` has two first channels but five second
channels. All five remain in this theorem. Another hostile is the first
row for `(-4,1,6)`: `5+30t+10t^2` has simple negative roots, whereas
deleting its middle channel leaves `5+10t^2`, whose roots are nonreal.
Neither positivity alone (`1+t+t^2`) nor a chosen sub-fibre suffices.

Targeted theorem/result/reference searches for trinomial/channel/return
real-rootedness, negative-axis cancellation, finite multipliers, Laguerre
and Schur--Szego found the two mechanisms above, not this all-row theorem.
The [incoming width-eight report](overnight_20260906_moments_width8.md)
proves a different fact: finite symbolic resultants separate enough early
moments throughout an endpoint strip. The
[all-carry two-first-channel theorem](nc2_two_rung_overnight_hexagon_sep05.md)
separates two levels for arbitrary endpoints when the first row has two
channels. Neither separation follows merely from the present zero-location
theorem. No external-priority claim is made.

The live board is: full affine fibres; factorial splitting; finite stable
symbols; open simple-root locus; scalar phase wall; cross-mass response.
The new operation connects the first four. The last remains the missing
ingredient for arbitrary higher-multiplicity two-rung separation.

## 3. Exact factorial splitting

Write rising and falling factorials as

```text
(u)_j=u(u+1)...(u+j-1),
u_down_j=u(u-1)...(u-j+1).
```

For `h>=1` introduce the multisets

```text
eta_s=h+(r-s)/B,         0<=s<B, s!=r,
mu_i=(x+i)/delta,        1<=i<=delta,
nu_i=(z+i)/A,            1<=i<=A, i!=A-z,
kappa=B^B/(delta^delta*A^A)>0.                       (3)
```

Every `eta_s>h-1`, while all `mu_i,nu_i>0`. There are exactly `B-1`
falling parameters and `B-1` inverse-rising parameters. If `p_j` is the
coefficient in `(1)`, direct grouping of the consecutive factorial factors
gives

```text
p_j/p_0=binom(h,j)*kappa^j
         product_(s!=r) (eta_s)_down_j /
         [product_i (mu_i)_j product_(i!=A-z) (nu_i)_j].             (4)
```

Indeed the numerator factorial ratio contributes
`B^(B*j) product_(s=0)^(B-1) (h+(r-s)/B)_down_j`.
The two increasing factorial ratios contribute
`delta^(delta*j) product_i(mu_i)_j` and
`A^(A*j) product_(i=1)^A((z+i)/A)_j`.
The omitted numerator parameter `s=r` is exactly `h`; the omitted
denominator parameter `i=A-z` is exactly `1`. Their quotient is
`h_down_j/(1)_j=binom(h,j)`, proving `(4)` without an asymptotic formula.

## 4. The precise finite-preserver input and its Laguerre symbols

For degree bound `h`, let `T_lambda(t^j)=lambda_j*t^j` and set

```text
J_lambda(t)=sum_(j=0)^h binom(h,j)*lambda_j*t^j.
```

If `J_lambda` has degree `h` and only negative real roots, then
`T_lambda` preserves all real-rooted polynomials of degree at most `h`.
Here is the exact cited sufficient condition: Borcea--Branden,
[*Polya-Schur master theorems for circular domains and their boundaries*, arXiv:math/0607416v6, Theorem 2(b), printed page 5](https://arxiv.org/pdf/math/0607416#page=5),
states preservation when `T[(s+w)^h]` is real stable. In our case this
symbol is `w^h J_lambda(s/w)`, a nonzero constant times
`product_i(s+rho_i*w)`, `rho_i>0`; it cannot vanish when both variables
have positive imaginary parts. The proof of the cited sufficiency is
Lemma 4 and Theorem 4 on printed pages 12--13. Thus the required property
is preservation on the **full** real coefficient space, not just on
positive-coefficient inputs.

For the two kinds of diagonal factors in `(4)` the symbols are explicitly

```text
sum_j binom(h,j)*eta_down_j*t^j
  =h!*t^h*L_h^(eta-h)(-1/t),                       eta>h-1,       (5)

sum_j binom(h,j)*t^j/(mu)_j
  =h!/(mu)_h * L_h^(mu-1)(-t),                    mu>0.          (6)
```

These identities follow coefficientwise from

```text
L_h^alpha(s)=sum_(k=0)^h binom(h+alpha,h-k)*(-s)^k/k!.
```

For completeness the needed Laguerre zero fact has a short proof.
Rodrigues' formula is

```text
L_h^alpha(s)=s^(-alpha)*exp(s)/h! *
             (d/ds)^h[exp(-s)*s^(h+alpha)].
```

For `alpha>-1`, integration by parts `h` times makes this polynomial
orthogonal to every polynomial of degree less than `h` under the positive
weight `s^alpha*exp(-s)` on `(0,infinity)`. Boundary terms vanish: at zero
the differentiated bracket in every boundary term has a factor of order
at least `s^(alpha+1)`, and at infinity the exponential dominates every
power. If there were fewer than `h` sign-changing roots in the interval,
multiplication by the product of those linear factors would produce a
nonzero constant-sign integrand against the same positive weight, contrary
to orthogonality. Degree `h` therefore forces `h` distinct positive roots.
This also excludes other or multiple roots.

The parameters in `(5)` and `(6)` are both greater than `-1`. Reflection
and negative reciprocal therefore give simple negative roots for every
symbol in `(4)`. There is no zero or degree loss, since all finite
diagonal entries are strictly positive.

## 5. Why simplicity also survives

An invertible real linear operator on the degree-at-most-`h` coefficient
space that preserves real-rootedness sends the open set of degree-`h`
polynomials with simple real roots into itself, provided it preserves
degree. Indeed its image of that open set is open and consists of
real-rooted polynomials. A repeated-root polynomial `Q`, with repeated
root `rho`, is not an interior point of that set: the arbitrarily close
polynomial

```text
Q(t)+epsilon*Q(t)/(t-rho)^2,
epsilon>0,
```

has the nonreal factor `(t-rho)^2+epsilon`. This is an explicit boundary
witness, valid for every multiplicity at least two.

Each factor operator in `(4)` is invertible on the full coefficient space
and preserves degree, because all its `h+1` diagonal entries are positive.
Start with `(1+t)^h`, apply one falling-factor operator, and obtain its
simple-negative Laguerre symbol `(5)`. Apply the remaining factors in any
order. The cited preservation theorem and the openness argument retain
simple real roots at every step. Positive coefficients force them to be
negative; nonzero constant terms rule out zero. Finally the positive
rescaling `t -> kappa*t` and nonzero factor `p_0` give exactly `(1)`.
There is at least one falling factor because `B>=2`. The case `h=0`
was already handled. This proves the factorial theorem.

This simplicity mechanism is also the elementary basic-case perturbation
behind [Kostov--Shapiro, *On Schur--Szego composition of polynomials*, arXiv:math/0605377, Theorem 6(ii)](https://arxiv.org/pdf/math/0605377#page=2),
used in THM-2760. No general multiplicity classification from that paper
is required here.

## 6. Completeness and the Laurent consumer

A return channel `(x',y',z')` at mass `m` obeys

```text
x'+y'+z'=m,                 a*x'=b*y'+c*z',
A*y'+B*z'=a*m/g.                                      (7)
```

If the right side of the last equation is not integral, the support
fibre is empty. Otherwise `gcd(A,B)=1` gives a unique `z` with
`0<=z<A` and `B*z=a*m/g (mod A)`. Set

```text
y=(a*m/g-B*z)/A=B*h+r,  0<=r<B,
x=m-y-z.
```

Whenever the fibre is nonempty, every nonnegative solution of the last
equation in `(7)` is uniquely

```text
y'=B*h+r-B*j,             z'=z+A*j,              0<=j<=h.
```

In particular `h>=0`. For these solutions the forced negative count
`x'=m-y'-z'` is strictly positive, since
`a*x'=b*y'+c*z'>0`. The strict inequality follows from `a*m/g>0`, which
prevents both positive counts being zero. Thus neither a hidden left
channel nor a negative-count truncation is omitted. The full fibre is
exactly `(x+delta*j,B*h+r-B*j,z+A*j)`, and `x>0` at its left end.

Its coefficient monomial factors as

```text
alpha^(x+delta*j)*beta^(B*h+r-B*j)*gamma^(z+A*j)
 =alpha^x*beta^(B*h+r)*gamma^z * tau^j.                (8)
```

Equations `(1)` and `(8)` prove the all-mass claim. The scalar `tau` is
unchanged both by multiplying all coefficients by a common nonzero
constant and by changing the Laurent variable `u -> lambda*u`: the
corresponding exponents sum to zero and satisfy
`-a*(B-A)+c*A-b*B=0`. It is the faithful monomial quotient needed here,
not a guessed phase coordinate.

For every nonzero scalar value of `tau`, set `alpha=gamma=1` and take
any `B`th root for `beta` of `1/tau`; hence every negative scalar root
occurs for actual nonzero coefficients. The derivative
`partial tau/partial gamma=A*tau/gamma` never vanishes on the torus.
Since `P'(tau)!=0` at its zeros, `(8)` gives a nonvanishing gradient
there and proves the reduced, smooth zero-set assertion.

Finally the pair of charges `(-a,b)` always supplies a support return at
mass `(a+b)/gcd(a,b)`. Off the negative `tau` ray no earlier nonempty
support moment can cancel, which proves the stated detection bound.

## 7. Exact controls and honest boundary

The independent checker uses only integer/rational arithmetic. It splits
every row coefficient as in `(4)`, but tests root location and simplicity
by a separate Euclidean Sturm chain. It also reconstructs complete weighted
rows from the original negative-count equation before comparing them with
the canonical parameter producer. No floating-point root finder or
repository theorem import is used.

Finite universe:

- `2<=B<=5`, `1<=A<B`, `0<=h<=5`, all canonical residues, `0<=x<=2`,
  with **no** coprimality or charge-positivity filter: 1,530 rows, including
  144 with noncoprime `A,B`.
- All `1<=a<=5`, `1<=b<=4`, `b<c<=7`, all `1<=m<=18`, with no
  primitive-content filter: 1,045 nonempty rows and 575 empty rows.
- 112 pairs of exact Laguerre symbols, degree `1..8`, with
  `eta=h-1+j/7`, `mu=j/7`, `j=1..14`.
- Thirty wider rows, including degree 12, `x=0`, and noncoprime controls.
- All 304 proportional-resonance rows `1<=p,q<=4`, `0<=m<=18`, produced
  directly from `(9)` and compared with the factorial specialization.

Selected full coefficient lists are

```text
(-4,1,6),     m=5:   [5,30,10]
(-13,1,8),    m=7:   [42,210]
(-13,1,8),    m=14:  [14,6006,120120,210210,18018]
(-15,33,49),  m=16:  [4368,7280]
(-15,33,49),  m=32:  [64512240,2356099200,294512400].
```

The shifted second row for `(-13,1,8)` explicitly keeps both floor
carries. Its left channel differs from twice the first left channel; the
new monomial prefactor in `(8)` accounts for that, while `tau` stays fixed.

Hostile controls retain the nonreal zeros after deleting the middle
channel of `5+30t+10t^2`, and the fact that the simple-negative
polynomials `1+t` and `(1+t)(1+2t)` share a root. The latter is not a
Laurent counterexample; it refutes only the proposed logical implication
from individual zero geometry to two-polynomial coprimality.

Reproduce with

```bash
python3 04-computation/nc2_channel_realrooted_overnight_hexagon_sep05.py
python3 -O 04-computation/nc2_channel_realrooted_overnight_hexagon_sep05.py
```

The unbounded conclusion depends on the audited proof and its cited
finite-preserver input, not on the size of this control box. The final
normal/optimized/stored replay manifest is recorded below.

All 27,847 explicit gates pass. Normal and optimized runs are byte-identical
to the stored output. Raw-LF SHA-256:

```text
394c1fe5acadf47fdd1e1b4cbd104efc2fe2e5f3e9e72bcb274d0d0459b3c5bc  04-computation/nc2_channel_realrooted_overnight_hexagon_sep05.py
a51cbae1b91c40718db22284bc221b577501f28c52a5d41b3d4fbb90739c8248  05-knowledge/results/nc2_channel_realrooted_overnight_hexagon_sep05.out
```

Audit separation: root checked the primary finite-preserver theorem and
independently audited the full mathematical proof; observer independently
audited the factorial identities, Laguerre orthogonality, full-space
openness argument, all-mass completeness, quotient invariances and torus
consumer. The additional proportional-resonance specialization is the
literal substitution recorded after `(9)` and is checked on all 304
declared rows. These finite checks assert neither cross-level coprimality
nor a uniform signed remainder.
