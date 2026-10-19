# Every-target rough twin-center sums at power sieve depth

19 October 2026. A classical-sieve application, with the target uniformity
and local factors proved here. This continues the
[structural note](prime-pairs-structural-note.md), Section 5, beyond its
logarithmic-depth CRT construction.

## Result first

Write P^-(n) for the least prime factor of n>1, and Omega(n) for its number
of prime factors counted with multiplicity. For a positive integer m put

  I_m = {a in Z : m/3 < a <= 2m/3},
  F_m(a) = (6a-1)(6a+1)(6(m-a)-1)(6(m-a)+1),
  R(m,z) = #{a in I_m : P^-(F_m(a)) >= z}.

Thus a counted representation m=a+(m-a) has two centers, each with both
neighbors free of primes below z. This is a condition on four integers,
not merely a condition on their centers.

For each prime p>=5 define

  rho_p(m) = 2 if p divides m;
             3 if p divides 9m^2-1;
             4 otherwise,
  V_m(z) = product_(5<=p<z) (1-rho_p(m)/p).                 (1)

The two divisibility cases are disjoint for p>=5.

### Theorem 1. PROVED, using the cited combinatorial sieve.

Fix any real u>37 and set z_m=(6m)^(1/u). Then

  limsup_(m -> infinity)
    | R(m,z_m) / [(m/3) V_m(z_m)] - 1 |
      <= exp(37-u).                                      (2)

The limit is over all positive integers m, with no exceptional target set.
In particular, there are absolute c>0 and m_0 such that, for every m>=m_0,
at least c m/(log m)^4 integers a in I_m have all four numbers

  6a-1, 6a+1, 6(m-a)-1, 6(m-a)+1

with at most 37 prime factors each, counted with multiplicity.

For the last assertion take u=38. Every factor is greater than 1 and
strictly less than 6m for m>=2. If any factor had 38 prime factors, all at
least (6m)^(1/38), it would be at least 6m, a contradiction. This choice of
scale avoids losing a constant in the almost-prime bound.

Two convenient quantitative specializations of (2), valid for all
sufficiently large m, are

  0.6 (m/3)V_m(z_m) <= R(m,z_m) <= 1.4 (m/3)V_m(z_m)
      when u=38;
  |R(m,z_m)/[(m/3)V_m(z_m)]-1| < 0.00001
      when u=50.                                         (3)

The second statement follows because exp(-13)<0.00001. No numerical onset
m_0 is asserted. Neither 37 nor 50 is an optimized sieve parameter, and the
almost-prime conclusion is not claimed as a literature record or a new
sieve theorem.

## 1. The analytic input and its precise scope

We use a conservative consequence of Jesse Thorner and Asif Zaman,
[A unified and improved Chebotarev density theorem, Section 6 and Lemma 6.2, pp. 1058-1060](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-s.pdf).
To match their printed hypothesis (6-2) literally, we impose a slightly
stronger product condition than the usual sieve-dimension condition.

Let 0<=g(p)<1/2 and h(p)=g(p)/(1-g(p)). Assume, for every 2<=w<=z,

  product_(w<=p<z) (1-h(p))^(-1)
       <= K (log z/log w)^kappa,                          (4)

where K>1. Let A_d=Xg(d)+r_d be divisibility sums of nonnegative indexed
weights, with P(z)=product_(p<z) p and g,h extended multiplicatively on
squarefree integers. For s=log D/log z>=beta=9 kappa+1, the cited lower and upper
beta-sieve coefficients satisfy |lambda_d|<=1 and vanish for d>=D.
Writing theta=1*lambda, their Lemma 6.2 bounds sum_b theta(b)h(b) between
1-exp(9 kappa-s)K^10 and 1+exp(9 kappa-s)K^10, with the corresponding
choice of lower or upper coefficients. Set their g'=g,g''=0 for the lower
estimate, and reverse them for the upper estimate. This is a literal
specialization of (6-1), (6-2), and (6-6), not an optimal-beta assertion.

The finite Euler-product identity

  sum_(d|P(z)) lambda_d g(d)
    = product_(p<z)(1-g(p)) sum_(b|P(z)) theta(b)h(b)

follows by inserting theta(b)=sum_(d|b)lambda_d and interchanging the two
finite sums. Applying the lower or upper divisor inequality to each index
and bounding its remainder by sum |r_d| therefore gives the weaker bound

  |S-X product_(p<z)(1-g(p))|
   <= exp(beta-s) K^10 X product_(p<z)(1-g(p))
        + sum_(d<D, d|P(z)) |r_d|.                       (5)

We use kappa=4 and beta=37. The coefficient of the exponential in (5) is
literally 1; there is no hidden O_kappa constant. The published lemma
cites Friedlander and Iwaniec, Opera de Cribro, Lemma 6.8. The exponent
37-s deliberately weakens its 36-s estimate.

For the source's p<=z cutoff convention, take cutoffs t<z approaching z
inside the preceding prime gap. Their prime set is exactly p<z, and the
bound is continuous in log t. This justifies our strict cutoff even when z
is itself prime; no endpoint prime is silently discarded.

The input allows multiplicities: our sequence is indexed by a, so the
identity F_m(a)=F_m(m-a) does not merge two representations.

## 2. Local admissibility, including the target-dependent bad primes

Modulo p>=5 the roots of F_m are

  1/6, -1/6, m+1/6, m-1/6.                              (6)

Each displayed pair has two distinct entries. If m=0 modulo p the pairs
coincide. Otherwise they intersect in one entry exactly when
m=1/3 or m=-1/3 modulo p. Those are the roots of 9m^2-1, giving (1).
There are no further collisions. Primes 2 and 3 divide none of the forms.
In particular rho_p(m)<=4<p for every prime p>=5 and every target m.

For a squarefree d composed of primes at least 5, the Chinese remainder
theorem gives exactly

  rho_d(m) = product_(p|d) rho_p(m)

roots modulo d. Thus rho_d(m)<=4^omega(d), uniformly in m. Primes dividing
m(9m^2-1) reduce the number of roots; they never create a local obstruction.
There is no need to exclude targets with many such prime divisors.

This root formula is inherited from the earlier note and from repository
[THM-2447, twin-gap local prime harmonics](../../../01-canon/theorems/THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md).
The new step here is its uniform lower-sieve application in the growing
interval I_m.

## 3. A fixed pre-sieve makes K as close to 1 as needed

Fix K>1, independently of m. Mertens' product theorem and convergence of

  product_(p>=11) [(p-4)/(p-8)] (1-1/p)^4

show that

  H(t) = product_(11<=p<t) (p-4)/(p-8)
       = C (log t)^4 (1+o(1))

for a positive constant C. Consequently there is a fixed w_0>=11 such that

  H(z)/H(w) <= K (log z/log w)^4
       for every z>=w>=w_0.                             (7)

For completeness, choose the relative error in the asymptotic smaller
than delta with (1+delta)/(1-delta)<K, simultaneously for all t>=w_0.
Taking the ratio proves (7), including short intervals and endpoints.

Let Q=product_(5<=p<w_0) p, a fixed integer, and first retain all residue
classes a modulo Q for which F_m(a) is coprime to Q. Their number is

  r_0(m) = product_(5<=p<w_0) (p-rho_p(m)),
  1 <= r_0(m) <= Q.                                     (8)

For the subsequent sieve set g_m(p)=rho_p(m)/p when p>=w_0, and g_m(p)=0
otherwise. Then g_m(p)<=4/p<1/2 on its support. Moreover,

  (1-g_m(p)/(1-g_m(p)))^(-1)
    = (1-g_m(p))/(1-2g_m(p)) <= (p-4)/(p-8).

Together with (7), this proves exactly (4), with kappa=4 and the fixed K.
If w<w_0<=z, replace w by w_0 on the product side; log z/log w only becomes
larger. If z<w_0, the product is 1. Thus (4) holds for all required w,z,
with constants independent of the target.

The pre-sieve is fixed after K is selected. It is not an m-dependent
primorial whose size is hidden inside an error term.

## 4. The full local density survives the remainder estimate

Put L=m/3, the real length of the half-open interval defining I_m, and

  X_m = L r_0(m)/Q.

For squarefree d with prime factors >=w_0, count a in I_m satisfying both
(F_m(a),Q)=1 and d|F_m(a). There are r_0(m)rho_d(m) admissible residue
classes modulo Qd. Counting each class in a half-open real interval incurs
error at most 1. Hence, exactly in the form required by (5),

  A_d = X_m rho_d(m)/d + r_d,
  |r_d| <= r_0(m)rho_d(m) <= Q 4^omega(d).               (9)

This is uniform in d and m; it does not invoke distribution of primes in
arithmetic progressions. Also,

  X_m product_(w_0<=p<z)(1-g_m(p)) = L V_m(z).          (10)

Keeping all small-prime residue classes, rather than choosing one of
them, is what preserves the full target-dependent product in (10).

Let tau_4(d) count ordered factorizations of d into four positive integers.
On squarefree d, 4^omega(d)=tau_4(d). The elementary bound

  sum_(d<=D) tau_4(d)
   <= D [sum_(n<=D) 1/n]^3
   <= D(1+log D)^3

follows by summing floor(D/(abc)) over a,b,c. Therefore

  sum_(d<D, d|P(z)) |r_d| <= Q D(1+log D)^3.            (11)

Take D=m/(log m)^10 and z=(6m)^(1/u). Then

  s=log D/log z -> u.

For fixed u>37, s>=37 eventually. The generic product also gives

  V_m(z) >= product_(5<=p<z)(1-4/p) >> 1/(log z)^4,     (12)

uniformly in m. Dividing (11) by L V_m(z) gives O_(K,u)((log m)^(-3)),
which tends to zero. Applying (5) and (10) yields

  limsup_m |R(m,z)/[L V_m(z)]-1| <= K^10 exp(37-u).

This holds for every fixed K>1. Let K decrease to 1 after taking the
limsup. The quantity on the left does not depend on the pre-sieve, so (2)
follows. There is no claim that K and Q vary with m in one application.
Finally (12) and u=38 prove the asserted c m/(log m)^4 lower bound. QED.

## 5. The full singular series really controls this rough problem

Define the positive singular series

  S_4(m) = 81 product_(p>=5)
                 (1-rho_p(m)/p)/(1-1/p)^4.              (13)

For fixed m, all but finitely many factors have rho_p=4 and are 1+O(p^-2),
so the product converges. It also has an absolute positive lower bound,
because replacing every rho_p by 4 decreases every factor. In the notation
of the earlier structural note, S_4(m)=(12 C_2)^2 Gamma(m).

### Proposition 2. PROVED uniform truncation.

For fixed u>0 and z=(6m)^(1/u),

  V_m(z) = [e^(-4 gamma) S_4(m)/(log z)^4] (1+o_u(1)),  (14)

uniformly over all positive integer targets as m tends to infinity;
gamma is Euler's constant.

Proof. In the tail of (13), the logarithms of the generic rho_p=4 factors
sum to O(1/z). An exceptional prime divides N_m=m(9m^2-1), and replacing
its generic factor changes the logarithm by at most 2/(p-4). The number
of exceptional primes p>=z is at most log N_m/log z. Their total change
is therefore

  O(log N_m/(z log z)) = O_u(1/z).

The full singular series and its truncation thus have ratio 1+O_u(1/z).
Lastly Mertens gives

  product_(5<=p<z)(1-1/p)^4
       = 81 e^(-4 gamma)/(log z)^4 (1+o(1)).

Combining this with (13) proves (14). QED.

Consequently the denominator in (2) can equivalently be replaced by

  (m/3) e^(-4 gamma) S_4(m)/(log z_m)^4.                (15)

This is an unconditional, all-target comparison between actual rough
representations and the same arithmetic local factor that appears in the
four-prime heuristic. It supplies a controlled adjacent model for studying
the local symmetry identities.

## 6. What remains at the prime boundary

The theorem proves representations with uniformly bounded numbers of
prime factors. It does not prove that both centers belong to
K={k>=1:6k-1 and 6k+1 are prime}.

To force the four positive forms below 6m to be prime merely by removing
small factors, one would sieve above sqrt(6m), corresponding roughly to
u=2. Our elementary distribution estimate has D=m/(log m)^10, so s would
also approach 2. The input used here only supplies the stated bound for
s>=37; its lower estimate gives no information at s=2. The number 37 is a
conservative parameter of this chosen sieve, not the optimal sieving limit
and not a claim that every improvement must fail there.

More importantly, neither the exact local roots, the singular-series
floor, nor the Fourier positivity proved earlier supplies an estimate for
the unsieved prime-pair error at that boundary. The transition from rough
integers to primes still requires additional arithmetic information.

Useful next steps are concrete:

1. Replace the conservative combinatorial weights with a precisely sourced
   dimension-four lower sieve or weighted almost-prime sieve, and quantify
   the improved factor count without losing uniformity in m.
2. Retain the all-target local factor while allowing short subintervals of
   I_m; their length determines the admissible distribution level directly.
3. Compare measured rough counts at increasing depth with (15), separating
   the rigorously controlled range from the conjectural prime endpoint.

The supporting script checks local roots, mixed CRT residue counts, and
interval discrepancies exactly. It does not numerically certify the
asymptotic onset or replace the cited fundamental lemma.
