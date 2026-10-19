# Prime pairs: local symmetry, additive surplus, and the missing correlation

19 October 2026. Status: PROVED elementary propositions, explicitly marked
CONDITIONAL transfers, and supporting FINITE-EXACT tests. The statements are
self-contained; no claim of literature priority is attached.
Primary literature and exact inherited theorem paths are in the
[companion source ledger](literature-and-inheritance.md).

## What this note adds

A twin center is an integer k >= 1 for which 6k-1 and 6k+1 are both prime.
Write K for the set of twin centers in these divided-by-six coordinates.
The question whether every sufficiently large m belongs to K+K is a strong
Goldbach-adjacent problem. It remains open, as does the infinitude of K.

There is nevertheless a precise tractable layer between modular pictures and
those open questions:

1. A centrally symmetric forbidden set makes sum counts and difference counts
   identical in every finite residue model. This extends the repository's
   established local twin-gap factor to the Goldbach direction.
2. The twin-center local factor has an explicit positive Fourier spectrum,
   a universal positive floor, and arithmetic moments. Its mean over integer
   targets is exactly 1 in the limit after the correct normalization.
3. Requiring the target itself to survive the same finite sieve produces an
   exact local representation advantage. The reason is that the forbidden
   pair of residues is sum-free. The full product advantage is about 1.41189.
4. Every sufficiently long interval contains the predicted number of rough
   twin-center representations, with an elementary uniform error. This is an
   unconditional adjacent theorem, although its sieve depth is far short of
   primality.
5. Prime-power corrections are quantitatively negligible at the conjectured
   prime-pair and four-prime main-term scales. The real missing input is a
   joint prime correlation, not the deterministic square correction in a
   one-point prime race.
6. An exact Boolean-poset model separates full one-coordinate information
   from prime-pair information. It explains a genuine information loss without
   pretending to prove a universal impossibility theorem about sieves.

The prior [THM-2447 local root count](../../../01-canon/theorems/THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md)
is inherited, not relabelled as new. Its repository-relative path is also
recorded in the source ledger for relocation with this packet. The
main additions here are the general symmetry and sum-free-defect formulations,
the target-conditioned product, and the rigorously separated transfer gates.

## 1. A finite-group theorem, before specializing to primes

Let G be a finite abelian group of order q, written additively. Let B be a
proper subset of G, let s=|B|, and put A=G\B. Assume B=-B; then also A=-A.
For t in G define

  C_A(t) = #{a in A : a-t is in A},
  S_A(t) = #{a in A : t-a is in A},
  J_B(t) = #{b in B : b-t is in B},
  gamma_A(t) = q C_A(t)/(q-s)^2.

All pairs here are ordered, and elements can repeat.

### Proposition 1. Symmetry, spectrum, and exact variance. PROVED.

For every t,

  S_A(t) = C_A(t) = q-2s+J_B(t).                         (1)

The average of gamma_A over G is 1. With normalized Fourier coefficients

  hat f(chi) = (1/q) sum_x f(x) conjugate(chi(x)),

one has, for every nontrivial character chi,

  hat gamma_A(chi) = |sum_(b in B) chi(b)|^2/(q-s)^2.   (2)

The trivial coefficient is 1. In particular gamma_A is positive definite.
If E(B)=sum_t J_B(t)^2 is the additive energy of B, then

  average_t (gamma_A(t)-1)^2
       = [q E(B)-s^4]/(q-s)^4.                         (3)

If q>2s, the uniform pointwise bound is

  gamma_A(t) >= q(q-2s)/(q-s)^2
             = 1-s^2/(q-s)^2 > 0.                    (4)

Proof. Reflection gives t-a in A iff a-t in A, proving the first equality.
For the second, exclude the two translates B and t+B from G and use
inclusion-exclusion. Summing C_A over t counts all ordered pairs of A and
therefore gives |A|^2. Fourier transformation of autocorrelation gives the
squared modulus in (2); the nontrivial Fourier sum of 1_A is the negative of
that of 1_B. Finally sum_t J_B(t)=s^2. Taking the variance of
q(q-2s+J_B)/(q-s)^2 yields (3). Dropping J_B>=0 gives (4). QED.

A translate of this theorem is equally valid. If A is invariant under
x -> 2c-x, then S_A(t)=C_A(t-2c). The center and the shift must be retained.

This is an intrinsic symmetric relation. A tournament adds no information:
there is no canonical winner between t and -t, whose counts are equal.
The useful object is a positive-definite correlation kernel, not an arbitrary
orientation of tied residues.

### Proposition 2. The sum-free-defect identity. PROVED.

Let

  T(B) = #{(b1,b2) in B^2 : b1+b2 is in B}.

Then

  #{(a1,a2) in A^2 : a1+a2 is in A}
       = q^2-3qs+3s^2-T(B).                            (5)

Consequently the average local pair factor conditional on the target t being
in A is

  average_(t in A) gamma_A(t)
       = 1 + [s^3-q T(B)]/(q-s)^3.                    (6)

Proof. Among the q^2 pairs (x,y), remove the three events x in B, y in B,
and x+y in B. Each has size qs. Every intersection of two has size s^2,
and their triple intersection has size T(B). This proves (5), even without
symmetry. In the symmetric setting, summing C_A(t)=S_A(t) over t in A gives
the left side of (5). Multiply by q/(q-s)^3 and expand
(q-s)^3=q^3-3q^2s+3qs^2-s^3 to obtain (6). QED.

Thus sum-free forbidden sets produce a positive internal-target surplus.
A surplus is not automatic: B={0} gives T(B)=1 and conditional mean
1-1/(q-1)^2. This negative example rules out a generic claim that conditioning
a target to be admissible always helps.

### A genuine graph bridge: cubic bias and quartic fluctuations

Assume additionally 0 is not in B, so the Cayley graph with edges x~y iff
y-x belongs to B is a simple undirected s-regular graph. Let M be its
adjacency matrix. Its eigenvalues are real, indexed by characters, with
lambda_trivial=s. Closed-walk counting gives

  trace(M^3)=q T(B),  trace(M^4)=q E(B).

Equations (3) and (6) therefore become

  conditional mean = 1 - sum_(chi nontrivial) lambda_chi^3/(q-s)^3,
  variance = sum_(chi nontrivial) lambda_chi^4/(q-s)^4.         (6a)

This is a precise graph/number-theory bridge: additive-triple bias is a
signed third spectral moment, while local fluctuation energy is a fourth
spectral moment. The relation is undirected; arbitrary tournament
orientations would destroy the relevant self-adjoint spectrum.

For a simple graph trace(M^3) is six times its triangle count. Thus the
conditional mean exceeds 1 exactly when the graph has fewer than s^3/6
triangles. This is a numerical threshold for this regular Cayley graph,
not a comparison with an unspecified random-graph model.

### Explicit equal-energy, opposite-bias example. PROVED.

In F_13 choose B_1={1,5,8,12}={+/-1,+/-5} and
B_2={1,2,11,12}={+/-1,+/-2}. Both sets have size 4 and energy 36.
Their J_B(t) values for t=0,...,12 are respectively

  B_1: 4,0,1,1,2,0,2,2,0,2,1,1,0;
  B_2: 4,2,1,2,1,0,0,0,0,1,2,1,2.

Squaring and adding either row gives 36. Directly, T(B_1)=0 and T(B_2)=6;
for B_2 the six ordered pairs are (1,1),(-1,-1),(2,-1),(-1,2),
(-2,1),(1,-2). Consequently both local kernels have variance 212/6561,
but their internal-target means are 793/729 and 715/729 respectively.
One is above 1, the other below.

Their Cayley graphs both have 13 vertices, degree 4 and 13 four-cycles,
but have respectively 0 and 13 triangles. Indeed their fourth trace is
13*36=468, and the number of four-cycles is
[468-13*4*(2*4-1)]/8=13. Their third traces are 0 and 78.
Thus even degree, size, and exact four-cycle count fail to determine the
sign of the target-conditioning effect. This example closes a concrete
obstruction question rather than merely proposing a search for one.

## 2. The twin-center wheel and the sum/difference bridge

For every prime p>=5, set u=6^(-1) in F_p and

  B_p={u,-u},  A_p=F_p\B_p.

A residue k belongs to A_p exactly when neither 6k-1 nor 6k+1 vanishes
modulo p. This is a condition on roughness. A prime equal to p fails that
condition, so finite-prime exceptions must not be silently identified with
primality itself.

For a sum target m, the four forms are

  6a-1, 6a+1, 6(m-a)-1, 6(m-a)+1.                     (7)

For a difference target m, they are

  6a-1, 6a+1, 6(a+m)-1, 6(a+m)+1.                     (8)

Primes 2 and 3 divide none of these forms. At p>=5, Proposition 1 gives the
same count in (7) and (8):

  C_p(m)=p-c_p(m),
  c_p(m)=2  if p divides m,
          3  if p divides 9m^2-1,
          4  otherwise.                              (9)

The cases are disjoint. The overlap J_(B_p)(m) is 2 at m=0, is 1 at
m=+/-2u=+/-1/3, and is zero elsewhere. This re-proves the inherited
THM-2447 formula and explains why it also governs sums.

For a finite set P of primes >=5, let W=product_(p in P) p. Define A_W by
requiring membership in every A_p, and let

  delta_W=|A_W|/W=product_p (1-2/p),
  gamma_p(m)=p C_p(m)/(p-2)^2,
  Gamma_W(m)=product_p gamma_p(m).

The Chinese remainder theorem gives the exact equality

  #{a mod W : a and m-a both belong to A_W}
       = W delta_W^2 Gamma_W(m).                      (10)

The same right side counts a and a+m. No independence of actual primes is
being assumed: CRT independence is a theorem for uniform residue classes.

### The source/target map and its loss

Source: residue pairs (a,b) with a,b in the symmetric set A_W.
Map: (a,b) -> (a,-b).
Target: a difference representation a-(-b)=m instead of a+b=m.
Preserved: membership, multiplicity, every local sieve constraint.
Lost after returning to positive integers: the positivity region and its
cutoff. Goldbach uses 1<=a<=m-1; fixed-gap counting uses 1<=a<=X with m fixed.
The map sends the second positive summand to a negative integer. It therefore
cannot turn a fixed-shift prime-tuple theorem into a Goldbach theorem.

This is the decisive sidecar: the arithmetic local factor is the same,
but the region where all four forms are positive is different.

## 3. A positive floor, exact moments, and the most difficult local target

For an integer m>=1 put

  G = product_(p>=5) [1-4/(p-2)^2],
  w(m) = product_(p>=5,p|m) (p-2)/(p-4)
         times product_(p>=5,p|9m^2-1) (p-3)/(p-4),
  Gamma(m)=G w(m).                                    (11)

Only finitely many factors occur in w(m). The product G converges to a
strictly positive number because sum_p 4/(p-2)^2 converges. Equivalently,
Gamma(m)=lim_(z->infinity) Gamma_(product_(5<=p<=z)p)(m).
The target m=0 is excluded: its local product diverges, as it should for a
repeated rather than a genuinely four-form pattern.

### Proposition 3. Sharp universal floor. PROVED.

For every positive integer m,

  Gamma(m)>=G,

and equality holds exactly at m=1.

Proof. Every correction factor in w is >=1. Equality requires that all prime
factors of m(3m-1)(3m+1) lie in {2,3}. Since neither 3m-1 nor 3m+1 is divisible
by 3, both must be powers of 2. They differ by 2 and are at least 2; the only
such pair is 2,4. Hence m=1. Conversely m=1 has no exceptional prime >=5. QED.

This floor is local. For example m=701 has no representation in K+K in the
known finite census, although Gamma(701)>G. The smallest local factor and the
absence of actual representations are very different predicates.

### Proposition 4. Exact arithmetic moments. PROVED.

For every fixed real r, the following mean exists:

  lim_(N->infinity) (1/N) sum_(1<=m<=N) Gamma(m)^r
    = product_(p>=5) (1/p) sum_(t mod p) gamma_p(t)^r.  (12)

In particular the mean of Gamma is 1, and

  lim (1/N) sum_m (Gamma(m)-1)^2
    = product_(p>=5) [1+(6p-16)/(p-2)^4] - 1.          (13)

The variance product converges. At each prime, Proposition 1 has s=2,
E(B_p)=6 and therefore variance (6p-16)/(p-2)^4.

Proof of the passage from CRT averages to integer averages. This step is
important: finite-wheel independence alone is not a proof of (12).
Write g_p=p(p-4)/(p-2)^2. The factor gamma_p(m)^r/g_p^r is

  1 + alpha_p 1_(p|m)
    + beta_p 1_(p|3m-1) + beta_p 1_(p|3m+1),           (14)

where alpha_p=((p-2)/(p-4))^r-1 and
beta_p=((p-3)/(p-4))^r-1. At p>=5 the three events are disjoint.
For fixed r there is a constant C_r such that
|alpha_p|+2|beta_p|<=C_r/p for all p>=5.

Expand the finite product (14). A term chooses disjoint squarefree divisors
of m, 3m-1, and 3m+1. For their product d, each particular assignment is one
residue class modulo d, and its coefficient has absolute sum at most
C_r^omega(d)/d across all assignments. For 1<=m<=N, every relevant d divides
m(3m-1)(3m+1)<9N^3, so the complete expansion is captured by d<=9N^3.
The number of m in a particular residue class is N/d+O(1).

After division by N, the total rounding error is at most

  (1/N) sum_(d<=9N^3, squarefree) C_r^omega(d)/d
       <= (1+log(9N^3))^k/N -> 0,                     (15)

where k is any integer >=max(C_r,1). To justify the bound without an analytic
prime estimate, dominate k^omega(d) on squarefree d by the k-fold divisor
function and sum its harmonic series over a containing k-dimensional box.
The main coefficient series is absolutely convergent, bounded by
sum_d C_r^omega(d)/d^2. Letting N grow therefore gives the Euler product of
the exact local means. Multiplication by G^r proves (12).
The local mean is 1 for r=1; r=2 gives (13). QED.

The same expansion proves the arithmetic moments for negative r; the
coefficients then need not be positive, but the absolute bounds still hold.
These are moments of the explicit arithmetic local weight, not moments of
actual twin-prime representation counts.

### Fourier form

For every nonzero frequency k modulo p,

  hat gamma_p(k)=4 cos^2(2 pi k u/p)/(p-2)^2 >=0.       (16)

CRT tensorizes these nonnegative spectra. This is a complete finite-wheel
spectral description. It does not make the Fourier error of the prime
indicator nonnegative, small, or independent.

A useful probabilistic restatement is available but not needed in the proof.
For independent uniform p-adic residue coordinates, Gamma_W is a mean-one
martingale with uniformly bounded second moment. Its second-moment limit is
the product in (13), and its values have the common lower bound G. This
profinite random model and the arithmetic averages agree at the level of
(12); neither identifies the prime process with independent coordinates.

## 4. Why internal twin-center targets are locally easier

For p>=5 the forbidden pair B_p={u,-u} is sum-free. Its sums are 0 and +/-2u;
none belongs to {u,-u}, since p is neither 2 nor 3. Proposition 2 therefore
gives

  average_(m in A_p) gamma_p(m)=1+8/(p-2)^3.           (17)

Equivalently, the exact number of local additive triples is

  #{(a,b) mod p : a,b,a+b all belong to A_p}
       = p^2-6p+12.                                  (18)

### Proposition 5. Exact finite-wheel target advantage. PROVED.

For every finite wheel W as above,

  (1/|A_W|) sum_(m in A_W) Gamma_W(m)
       = product_(p|W) [1+8/(p-2)^3].                 (19)

Proof. Both the conditioning m in A_W and the multiplicative observable
Gamma_W factor by CRT. Apply (17) in each factor. QED.

The right side tends to a finite number B>1. At p=5 alone, the surviving
classes are {0,2,3}; their seven additive triples yield factor 35/27.
At p=7 the factor is 133/125. This gives a quantitative version of the
repository's observation that internal twin-center targets avoid some of
the weakest mod-5 parent channels.

For truncation at z>4,

  0 <= log(B/B_z)
     <= 8 sum_(n>z) 1/(n-2)^3 <= 4/(z-2)^2.          (20)

The final inequality uses comparison with an integral; the bound sums over
all integers and is therefore valid for the smaller set of prime factors.
Thus the product through one million determines B to a relative tail below
4.00002e-12. A directed-rounding interval is recorded with the supporting
calculations. Explicitly,

  1.41188824329665943969327923841603421673627157
    <= B <=
  1.41188824330230701525676813013037419105105626.

The rounded value is B=1.41189 to five decimal places. The separately labelled
binary floating-point product is not used to certify this interval.

Scope. Equation (19) conditions on survival of a finite residue sieve. It is
not the assertion that averaging over actual m in K multiplies their parent
count by B. That assertion needs a further joint-primality theorem involving
the target forms 6m+/-1 as well as the four parent forms. The six affine forms
include parallel pairs, so one must not casually invoke a finite-complexity
linear-forms theorem for primes. The ordinary twin-prime obstruction is still
inside that family.

## 5. An unconditional Goldbach theorem for rough twin centers

Fix any finite wheel W from primes >=5. Let R_W(m) count ordered pairs of
positive integers a,b with a+b=m and a,b in A_W modulo W. These are pairs for
which all four integers 6a+/-1,6b+/-1 avoid the selected small primes.

### Proposition 6. Uniform periodic representation theorem. PROVED.

For m>=2 write

  r_W(m)=W delta_W^2 Gamma_W(m).

Then

  |R_W(m)-(m-1)delta_W^2 Gamma_W(m)| <= r_W(m),         (21)

and more concretely

  R_W(m) >= floor((m-1)/W) r_W(m).                    (22)

In particular every integer m>W has such a representation, and for m-1>=2W,

  R_W(m) >= (m-1)delta_W^2 G/2.                       (23)

Proof. The indicator of a and m-a belonging to A_W is periodic in a with
period W and has exactly r_W(m) successes in each complete period. Partition
1,...,m-1 into complete periods and a final segment. The final segment
contains between 0 and r_W(m) successes. Every local factor is positive by
(9), so r_W(m)>=1. The lower bound follows using Gamma_W>=G. QED.

There is a useful growing-depth version. Standard prime number and Mertens
estimates, used only for this corollary, give

  log W_z ~ z,   delta_(W_z) asymptotic to c/(log z)^2

for W_z=product_(5<=p<=z)p and an explicit positive c. Taking
z=(1-epsilon) log m for a fixed 0<epsilon<1 gives W_z=o(m). Uniformly for
large m,

  R_(W_z)(m) >= c_epsilon m/(log log m)^4              (24)

for some positive constant c_epsilon. Thus every sufficiently large target
has many twin-center pairs whose four entries have no prime factor <=
(1-epsilon)log m.

The primality threshold is completely different. To force integers up to
6m+1 to be prime by trial division requires removing factors up to their
square roots, far beyond logarithmic depth. The period then dwarfs m and
(22) becomes empty. This identifies exactly where this elementary theorem
stops. It is not a shortcut around the sieve parity problem.

## 6. The correct conditional transfer and its sharp threshold

Let

  R_K(m)=#{(a,b) in K^2 : a+b=m},
  C_2=product_(p>2) [1-1/(p-1)^2],
  S_2=12 C_2,
  S_4(m)=S_2^2 Gamma(m).

The factor 12 is essential because K is parameterized by k while the two
prime entries are about 6k. Direct local multiplication gives

  S_4(m)=81 product_(p>=5) (1-c_p(m)/p)/(1-1/p)^4.     (25)

The factors at 2 and 3 are 16 and 81/16; their product is 81.
For m>2 define the positive singular integral

  I_4(m)= integral from 1 to m-1 of
    dt/[log(6t-1) log(6t+1)
        log(6(m-t)-1) log(6(m-t)+1)].                 (26)

As m grows, I_4(m)~m/log^4(6m). This follows by separating fixed-relative-
size interior intervals from the endpoints and then allowing those relative
sizes to shrink; small power-size endpoint regions are negligible.

The Hardy-Littlewood-type expected main term is

  M(m)=S_4(m) I_4(m).

Its use as the actual main term is CONDITIONAL. Define the actual discrepancy
E(m)=R_K(m)-M(m); this definition itself makes no conjecture.

### Proposition 7. Pointwise gate. PROVED implication, OPEN input.

If for every sufficiently large m,

  E(m) >= -S_2^2 G I_4(m),                            (27)

then every sufficiently large m lies in K+K. A uniform estimate
E(m)=o(m/log^4 m) is sufficient.

Proof. Proposition 3 gives Gamma(m)>G for every m>1, since the unique
equality target is m=1. Therefore

  R_K(m)>=S_2^2 [Gamma(m)-G] I_4(m)>0

throughout the m>2 domain. QED.

The non-strict error bound in (27) is sufficient because the local floor
itself is strict here. For a quantitative margin, the stronger premise
E(m)>=-(1-eta)S_2^2 G I_4(m), with eta>0, yields
R_K(m)>=eta S_2^2 G I_4(m). A generic floor argument that only proves R>=0
would not imply existence; the unique equality case is load-bearing in
allowing eta=0 in this particular theorem.

For the weakest possible assumption at each target, it is enough that
E(m)>-M(m). Formula (27) is the useful uniform gate obtained using only the
universal local floor. A fixed-shift prime-quadruplet asymptotic for each
fixed m in (8) supplies neither (27) nor its necessary uniformity in (7).

### Mean-square gate and its limitation

On a dyadic interval N<=m<=2N, the floor gives M(m)>=cN/log^4 N for large N.
If Z_N is the number of targets with R_K(m)=0, then

  Z_N <= (log^8 N)/(c^2 N^2)
          sum_(N<=m<=2N) |E(m)|^2.                   (28)

Hence a total mean-square error o(N^3/log^8 N) yields only Z_N=o(N).
To force Z_N=0 from this inequality one needs the much stronger total bound
o(N^2/log^8 N), or some other mechanism controlling individual exceptions.
This factor-N distinction is exactly why an average result does not establish
an every-target Goldbach assertion.

One Fourier version makes the missing estimate concrete. Zero-pad two real
sequences f,g supported on 1,...,N into Z/MZ with M>2N. Their unnormalized
Fourier transforms satisfy, for every target t,

  |(f*f)(t)-(g*g)(t)|
    <= (1/M) sum_chi |hat f(chi)-hat g(chi)|
                           (|hat f(chi)|+|hat g(chi)|).       (29)

This follows by factoring hat f^2-hat g^2 and using Fourier inversion.
Use f=1_K and a calibrated rough-number g. Proposition 6 supplies a lower
bound for g*g; a bound in (29) below that margin would transfer existence.
The nonnegative spectrum of the uncalibrated wheel alone supplies no such
bound on f-g. No estimate in this note closes that gap.

## 7. Prime powers: the square bias is real, but too small to be the main obstacle

Define f(n)=1 when n is prime, and zero otherwise. Put

  j(n)=1/k if n=p^k for a prime p and k>=1,
  j(n)=0 otherwise,
  e(n)=j(n)-f(n),   H(X)=sum_(n<=X) e(n).

Thus 0<=f<=j<=1 and e is supported on proper prime powers. Uniqueness of
prime factorization makes j well defined. Exactly,

  H(X)=sum_(k>=2) pi(X^(1/k))/k
      = pi(sqrt X)/2 + O(X^(1/3) log X).               (30)

Using the prime number theorem, H(X)~sqrt X/log X. An elementary bound
H(X)=O(sqrt X) also suffices for the qualitative conclusions below.

### Proposition 8. Exact pair correction and four-form error. PROVED.

For ordinary Goldbach convolutions,

  (j*j)(N)-(f*f)(N)=2(f*e)(N)+(e*e)(N),
  0 <= (j*j)(N)-(f*f)(N) <= 2H(N).                    (31)

For fixed h>=1 and 1<=n<=X,

  0 <= sum_n j(n)j(n+h)-sum_n f(n)f(n+h)
        <= 2H(X+h).                                 (32)

For the four forms L_i(a) in (7), with 1<=a<=m-1,

  0 <= sum_a product_(i=1..4) j(L_i(a)) - R_K(m)
        <= 4H(6m+1).                                (33)

Proof. Expand j=f+e for the first equality. For numbers 0<=x_i<=y_i<=1,

  0<=product_i y_i-product_i x_i<=sum_i(y_i-x_i),

by telescoping. Apply this with x_i=f(L_i),y_i=j(L_i). Each form is injective
on the summation interval and its positive values are at most 6m+1, giving
(33). The same argument gives (31)-(32). QED.

Consequences. At the conjectured X/log^2 X prime-pair scale or m/log^4 m
four-prime scale, the right sides are lower order. An asymptotic formula at
those scales is equivalent with f replaced by j. A prime-power weighted
four-form lower bound exceeding 4H(6m+1) would already give a genuine prime
representation. Obtaining that lower bound remains a hard correlation
problem; weighted prime powers do not supply it automatically.

This explains how a deterministic square term can dominate the mean bias
of a one-point residue race at scale sqrt X/log X while being negligible at
the two-point main-term scale. The race and the pair problem ask for different
statistics. Removing the one-point square bias does not remove the sieve
parity obstruction.

### Deconvolution must precede addition

The inherited exact prime-power transform can be extended tensorially. For
n>=2,

  f(n)=sum_(k>=1, a>=2, n=a^k) mu(k) j(a)/k.           (34)

Thus

  (f*f)(N)=sum_(k,l>=1) mu(k)mu(l)/(kl)
                 sum_(a,b>=2; a^k+b^l=N) j(a)j(b).   (35)

All variables in (34)-(35) are integers; the bases satisfy a,b>=2 and
the exponents are positive. Thus all sums are finite. To prove (34), if n=p^r, its coefficient is
(1/r)sum_(k|r)mu(k), which is 1 at r=1 and zero otherwise; non-prime-powers
contribute zero. Substitute into the convolution to prove (35).

This is Möbius inversion on the product of two exponent-divisibility posets.
It is not a scalar inversion after compressing (a,b) to a+b. For N=10,
(j*j)(10)=11/3 while (f*f)(10)=3. Applying a scalar root correction at
floor(sqrt 10)=3 and floor(cuberoot 10)=2 gives zero corrections and leaves
11/3. Taking a root of the sum does not recover roots of its summands.
The missing coordinate is the separate exponent of each summand.

## 8. A precise poset obstruction to marginal-only arguments

This finite model is deliberately not asserted to be an additive prime
sequence. Its purpose is to say exactly what certain data fail to determine.

Choose r>=2 distinct primes q_1,...,q_r. Associate a subset S of {1,...,r}
with n_S=product_(i in S)q_i, taking n_empty=1. Prime numbers correspond
exactly to singleton subsets. Put epsilon(S)=(-1)^|S|. On ordered pairs of
subsets define two probability measures

  P_sigma(S,T)=4^(-r)[1+sigma epsilon(S)epsilon(T)],
  sigma in {+1,-1}.                                  (36)

### Proposition 9. Exact marginal blindness. PROVED.

Both measures are nonnegative probability measures and have the same full
uniform marginal distribution on each coordinate. For any U,V subseteq
{1,...,r}, the two measures give the same probability to

  U subseteq S and V subseteq T                       (37)

unless U and V are both the full set. Yet under P_-, the probability that
both n_S,n_T are prime is zero; under P_+ it is 2r^2/4^r.

Proof. The parity sum sum_S epsilon(S) is zero. Also
sum_(S containing U) epsilon(S)=0 unless U is the full set: pair subsets
using one element outside U. Factor the difference between the probabilities
of (37) into the corresponding two parity sums. For singleton S,T their
parity product is +1, so their P_- weight vanishes and their P_+ weight is
2/4^r. There are r^2 such pairs. QED.

Interpretation. The poset order is genuine divisibility, and even complete
one-variable information does not determine a two-variable prime event.
The two systems also agree on a large family of mixed divisibility tests.
This is a concrete warning against deriving pair existence from residue
margins or one-point bias corrections alone.

Boundary. The model does not enforce n_T=n_S+h or n_S+n_T=N. It therefore
does not prove that all conceivable sieves or all uses of those constraints
are powerless. Analytic results such as [Ford-Maynard's prime-producing sieve](https://arxiv.org/abs/2407.14368)
framework specify their admissible information much more precisely. This
example only refutes the particular implication from the displayed data to
prime-pair positivity.

## 9. Research roadmap with actual stopping criteria

### A. Internal-target surplus as a diagnostic

Compute parent counts R_K(m), not merely consecutive gap counts. Compare
ordinary targets and actual targets m in K after normalizing by the full
Gamma(m) and the correct logarithmic integral. The exact finite-sieve factor
B is a baseline to be tested, not fitted or assumed. Use disjoint ranges and
explicitly separate the exceptional small primes equal to one of the forms.

Outcome worth obtaining: an observed residual pattern after local effects
are removed, with a preregistered null. Stop before claiming an asymptotic
from persistent finite-range signs. The earlier repository warning about
correlated octaves applies directly.

### B. A roughness-depth experiment with a theorem at its base

At fixed m, increase z in W_z and count surviving representations. At
z<=roughly log m the periodic theorem predicts a strong uniform margin.
Record where the period-based certificate stops and where actual survivors
thin faster than this certificate can track. Include a hostile target such
as 701 and a target with many small prime factors.

Outcome: an exact loss ledger between modular abundance and final primality.
This experiment does not itself resolve the parity barrier. It can locate
which extra statistic a proposed proof must control.

### C. Replace bare support by a quantified bilinear estimate

Any serious prime claim should name its candidate estimate before running a
larger census. Options are a pointwise error like (27), the Fourier norm in
(29), or an explicit Type I/Type II package adequate for a cited sieve
framework. Verify which of these estimates is actually proved, rather than
using the word equidistribution as a substitute.

Outcome: a theorem whose hypotheses are independently checkable, or an
explicit obstruction/countermodel. If only a mean-square error at the scale
in (28) is accessible, state the resulting exceptional-set theorem and stop
short of an every-target conclusion.

### D. General forbidden-pattern classification

Propositions 1 and 2 work for any finite symmetric forbidden pattern.
Classify families where T(B)=0, or more generally qT(B)<|B|^3, so that
internal targets have a provable surplus. The threshold qT(B)=|B|^3 is exact.
Study how additive energy controls fluctuations while T(B) controls the
conditional mean. These are different invariants and can vary independently.

This is a self-contained additive-combinatorial question that can make real
progress without proving anything new about primes. The F_13 example above already separates equal energy from opposite
conditional-surplus sign. A useful next deliverable would classify the
smallest pattern sizes admitting this separation and determine which
restrictions on B make the cubic sign follow from other graph data.

## 10. Claims deliberately not made

No infinitude of twin primes, no binary Goldbach theorem, no all-target
four-prime lower bound, no independence law for primes, and no universal
sieve impossibility theorem is proved here. The exact finite-group and
arithmetic-moment propositions are actual theorems with proofs. The rough
representation theorem is unconditional. The prime transfers identify the
missing analytic hypotheses explicitly.
