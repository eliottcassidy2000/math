# Short intervals and the resolution of the singular series

19 October 2026. Continuation of
[the power-depth rough-pair theorem](power-depth-rough-pairs.md).
All statements below concern rough integers. Almost-prime conclusions
require the separate size comparison in Section 5.

There are two useful findings. The rough-count theorem is uniform in
arbitrarily large targets and every sufficiently long interval. The full
singular series, however, can contain significant information from primes
larger than the sieve cutoff. A sharp tail bound quantifies that loss, and
primorial targets attain its limiting envelope.

## 1. Definitions and the uniform short-interval theorem

For m>=2 and a in {1,...,m-1}, let

  F_m(a)=(6a-1)(6a+1)(6(m-a)-1)(6(m-a)+1).

All four factors lie strictly between 1 and 6m. Let J=(A,A+L] be a real
half-open interval with A>=0, L>0, and A+L<=m-1. Write

  R(m;J,z)=#{a in J intersect Z : no prime p<z divides F_m(a)},
  rho_p(m)=2 if p|m, 3 if p|9m^2-1, and 4 otherwise (p>=5),
  V_m(z)=product_(5<=p<z)(1-rho_p(m)/p).

Primes 2 and 3 never divide F_m(a).

### Theorem 1. PROVED, with the cited classical sieve input.

Fix u>37. For every epsilon>0, there is L_0(u,epsilon) such that for every
real L>=L_0, every integer m, and every admissible interval J of length L,

  | R(m;J,L^(1/u)) / [L V_m(L^(1/u))] - 1 |
       <= exp(37-u)+epsilon.                             (1)

In particular, R(m;J,L^(1/u)) >>_u L/(log L)^4 uniformly in m and J.
There is no upper bound on m in terms of L in this statement.

Proof. Fix K>1. The fixed pre-sieve Q from the preceding note makes the
strong product condition of
[Thorner--Zaman, Lemma 6.2, pp. 1058-1060](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-s.pdf)
hold with dimension four and this K, uniformly in m. Retain all r_0(m)
admissible classes modulo Q. Counting them jointly with the roots modulo
a squarefree d gives

  A_d=L [r_0(m)/Q] rho_d(m)/d+r_d,
  |r_d|<=Q 4^omega(d).

These estimates use only the interval length, not its location. Applying
the same upper and lower weights with D=L/(log L)^10 and z=L^(1/u) gives

  |R(m;J,z)-L V_m(z)|
    <= exp(37-s) K^10 L V_m(z)+Q D(1+log D)^3,
  s=log D/log z = u-10u log log L/log L.                  (2)

Here s>=37 eventually. Also V_m(z)>>1/(log z)^4 uniformly, so the final
error divided by L V_m(z) is O_(K,u)((log L)^(-3)). First let L tend to
infinity, then let the fixed K decrease to 1. This proves all the uniform
quantifiers in (1). QED.

The strict p<z convention is handled by approaching z from below inside
its preceding prime gap, as in the preceding note. Multiplicities in the
indexed sequence remain present; a and m-a are separate indices.

## 2. The omitted singular-series tail has a sharp envelope

Define

  b_p(m)=(1-rho_p(m)/p)/(1-1/p)^4,
  S_4(m)=81 product_(p>=5) b_p(m),
  T_m(z)=product_(p>=z) b_p(m)      (z>=11).              (3)

Thus T_m(z) is exactly the factor omitted when S_4 is truncated below z.
It should not be confused with an error in the short-interval count.

### Theorem 2. PROVED uniform tail envelope.

There is an absolute C>0 such that, uniformly for m>=2 and z tending to
infinity,

  exp(-C/z) <= T_m(z)
     <= [1+o_z(1)] max{1, log log(9m^3)/log z}^2.         (4)

The o_z(1) can be chosen independently of m. The exponent 2 in this
bound is attained: for every fixed c>0, let

  m_y=product_(p<=y) p,   z_y=y^(1/c).

Then, as y tends to infinity,

  T_(m_y)(z_y) -> max{1,c}^2.                            (5)

Proof of the upper and lower bounds. Put

  b_p^0=(1-4/p)/(1-1/p)^4,
  N_m=m(9m^2-1).

The generic factors satisfy 0<b_p^0<1 and log b_p^0=O(p^-2). Their product
for p>=z therefore lies between exp(-C/z) and 1. The exact remaining
correction is

  b_p(m)/b_p^0=(p-rho_p(m))/(p-4).

It equals 1 unless p|N_m, and always lies between 1 and (p-2)/(p-4).
This proves the lower bound in (4).

For the upper bound split the exceptional primes at

  Y=max{z, log(9m^3)}.

The product of their corrections for z<=p<Y is no larger than

  product_(z<=p<Y) (p-2)/(p-4)
     = (log Y/log z)^2 [1+o_z(1)],                       (6)

uniformly for Y>=z. To see uniformity, the ratio of
product_(11<=p<t)(p-2)/(p-4) to (log t)^2 tends to a positive constant:
Mertens gives the logarithmic factor and the remaining Euler product
converges absolutely. Its relative error is uniformly small for every
t>=z. Taking the ratio at Y and z proves (6).

For p>=Y, use

  log[(p-rho_p(m))/(p-4)] <= 4/p,
  sum_(p|N_m, p>=Y) 1/p
       <= log N_m/(Y log Y) <= 1/log Y.                  (7)

The first inequality uses p>=11. The second follows from
1/p<=log p/(Y log Y) and sum_(p|N_m)log p<=log N_m.
Thus the large-prime corrections contribute at most exp(4/log z)=1+o_z(1).
Since log Y/log z is exactly the maximum in (4), the upper bound follows.

Proof of sharpness. For m=m_y, every prime 5<=p<=y divides m, so rho_p=2.
Chebyshev's estimate log m_y is of order y gives
log log(9m_y^3)=log y+O(1). Applying (4) just above y shows that all factors
beyond y have product 1+o(1). If c>1, the remaining product from z_y to y is

  product_(z_y<=p<=y) (1-2/p)/(1-1/p)^4
       = (log y/log z_y)^2 [1+o(1)] -> c^2.              (8)

Again the quotient after removing the squared Mertens product converges
absolutely. If c<=1, (4) directly gives T_(m_y)(z_y)->1. QED.

This proof uses the logarithmic size of the exceptional-prime budget,
not just the number of exceptional primes. Bounding all corrections by
1/z alone would miss the sharp squared envelope.

## 3. When the full singular series is an honest normalization

Mertens and the definition (3) give the exact asymptotic comparison

  V_m(z)
    = [e^(-4 gamma) S_4(m)/(log z)^4]
         [1+o_z(1)]/T_m(z),                             (9)

uniformly in m; gamma denotes Euler's constant. There is no suppressed
assumption that T_m(z) tends to 1.

By (4), the full and truncated normalizations do agree whenever

  log log(9m^3)/log z <= 1+o(1),   z -> infinity.          (10)

For the short-interval cutoff z=L^(1/u), a convenient sufficient condition
is L>=(log m)^u. More generally (10) permits
L>=(log m)^(u-o(1)); this is why an exact threshold among all possible
functions L(m) would be too strong a claim.

Combining (1), (9), and (10) proves the same relative-error bound
exp(37-u)+epsilon with denominator

  L e^(-4 gamma) S_4(m)/(log z)^4,                       (11)

uniformly in every admissible interval of that length. In particular this
holds when L>=m^theta for any fixed theta>0, but power length is far more
than the normalization itself needs.

## 4. A concrete obstruction below the power-of-log boundary

### Corollary 3. PROVED sharpness at the fixed exponent scale.

Fix u>37 and A>0. Along the primorial targets m=m_y, set

  L_y=(log m_y)^A,   z_y=L_y^(1/u).

Then

  V_(m_y)(z_y)
   / [e^(-4 gamma) S_4(m_y)/(log z_y)^4]
      -> min{1,(A/u)^2}.                                (12)

Thus A=u is the sharp exponent boundary for uniform replacement by the
full singular series in the family L=(log m)^A.

Proof. Chebyshev's estimate gives
log z_y=(A/u)log y+O(1). The proof of (5), with this harmless O(1) change
in logarithms, gives T_(m_y)(z_y)->max{1,(u/A)^2}. Use (9). QED.

For example take u=50 and A=25. For every choice of admissible intervals
J_y of length L_y, every limit point of

  R(m_y;J_y,z_y)
       / [L_y e^(-4 gamma)S_4(m_y)/(log z_y)^4]

lies between (1-exp(-13))/4 and (1+exp(-13))/4. The full-series prediction
is approximately four times too large, while the truncated-density
prediction is controlled by Theorem 1. This is a counterexample to a
uniform short-interval use of the full singular series, not a counterexample
to its usual global-scale role.

The mechanism is explicit. Many primes between z_y and y divide the target.
They change the full arithmetic factor, but the sieve at depth z_y does not
remove them. The truncated product records exactly which local tests were
actually imposed.

## 5. Bounded prime-factor counts require a power-sized interval

For an almost-prime consequence fix theta in (0,1], c_0>0, and an integer
r with r theta>37. Consider any admissible interval with

  L>=c_0 m^theta,   z=(6m)^(1/r).

The proof of (2), with D=L/(log L)^10, now has
s>=r theta-o(1). It yields, uniformly over these intervals,

  R(m;J,z) >>_(theta,c_0,r) L/(log m)^4.                 (13)

Indeed the pre-sieve K can be chosen so that
K^10 exp(37-r theta)<1, and the same remainder is negligible. Each of the
four forms is below 6m and has no prime factor below z. Hence each has at
most r-1 prime factors, counted with multiplicity.

For instance every sufficiently long interval of length at least sqrt(m)
inside {1,...,m-1} contains many such representations with at most 74 prime
factors in each form: take theta=1/2 and r=75. These counts are deliberately
unoptimized consequences of the conservative sieve.

The polylogarithmic intervals of Section 3 have no analogous fixed P37
conclusion: their cutoff grows only as a power of log m while the forms
can be of size m. Their proved conclusion is about rough representations
and the validity of the arithmetic normalization.

## 6. Practical research consequence

A comparison of observed short-interval counts against S_4(m) should
record the cutoff z and the interval length, and either verify (10) or
retain T_m(z) explicitly. The local model, its singular-series completion,
and the prime endpoint have different information content. Here the first
two are related by an exact tail and a sharp uniform bound; the prime
endpoint still requires information beyond this sieve.

## 7. Supporting controls and a finite-tail diagnostic

The companion script checks 600 finite tail factorizations and 212
translated-interval discrepancies exactly. It also illustrates the
primorial limit 4 at several cutoffs, using ordinary floating arithmetic
rather than certified numerical intervals.

The residual analytic bound behind that diagnostic is explicit. For an
integer B>=11,

  -13/B <= log product_(p>B) b_p(m)
           <= 4 log N_m/(B log B).                      (14)

Indeed the power series gives

  |log b_p^0| <= 8/[p^2(1-4/p)] <= 13/p^2

for p>=11; summing over integers greater than B bounds the generic tail.
The exceptional tail satisfies (7) with Y=B. Formula (14) is proved; the
script's decimal evaluations are illustrative and are not used in any
proof above.
