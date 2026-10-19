# Scaled integers, cusp windows, and arithmetic correlations

19 October 2026. Status: PROVED elementary transfer results and obstructions;
the cusp application imports the separately proved threshold-width theorem;
the rough-number application imports the repository's classical-sieve result.
No twin-prime, binary Goldbach, normality-of-pi, or tournament-balance claim is
made. No historical priority is claimed.

## Results first

There are three mathematically different ways to put large integers in a
small real interval. They have different arithmetic consequences.

1. Exact lattice scaling k -> k/Q is lossless when Q is retained. Its
   microscopic resolution is 1/Q. An unmarked continuum limit is not lossless.
2. Reading the prefixes floor(b^n alpha) samples a sparse integer sequence.
   Its finite-sieve behavior has an exact digit-driven Markov model when
   alpha is sampled uniformly. A deterministic constant needs a separate
   arithmetic hypothesis. The model has an explicit, sharp mixing error.
3. Varying a quadratic cusp parameter over a window that traverses consecutive
   escape counts transfers ordinary short-interval arithmetic with vanishing
   relative distortion. This gives a proved rough twin-center Goldbach
   statement in actual cusp-parameter windows.

A fourth connection concerns the repository's poset correlation work.
Independently normalizing several integer counts into (0,1) can reverse a
correlation inequality. Exactly which normalization gauges preserve all
diamond ratios can be classified, and the required exponent correction is
explicit.

## 1. Exact scaling and its sharp resolution boundary

Let Q be a positive integer and encode 0<k<Q as x=k/Q. Then Qx=k exactly,
so every residue, factorization, primality predicate and additive relation is
preserved. For example,

    a+b=m  iff  a/Q+b/Q=m/Q,
    a-b=h  iff  a/Q-b/Q=h/Q.

The scale Q is part of the data. In particular, a fixed integer shift h
becomes the microscopic real shift h/Q, not a fixed nonzero real distance.

**Proposition 1, sharp worst-case decoding threshold.** If an observation y
is known to satisfy |y-k/Q|<1/(2Q), rounding Qy to the nearest integer
recovers k. No decoder can uniformly recover even k modulo an integer W>=2
from observations with error allowance >=1/(2Q), provided the permitted
integer range contains two consecutive integers.

Proof. The first claim is the disjoint-nearest-lattice-point property. For
the second, (k+1/2)/Q is compatible with both k and k+1 at that error
allowance; these have different residues modulo every W>=2. QED.

The lower bound is a statement about uniform guarantees under an error
allowance. A particular observation with a wider allowance can still happen
to identify an integer, or a residue, from extra information.

For Q=b^n the worst-case absolute precision is therefore of order b^(-n).
Placing the number in (0,1) moves the information to fractional digits; it
does not remove the required n base-b digits. If the encoding is instead
x=floor(alpha Q)/Q, the exact certification test is that the entire interval
Q times the certified enclosure for alpha lies in one half-open integer
cell. Proximity to a cell boundary matters, not merely a relative error.

The same issue occurs in the inherited odd-order reconstruction theorem.
For B=-B in a finite abelian group of odd order q, the repository proved

    1_B(x) = J_B(2x) mod 2,
    J_B(t)=#{(b,c) in B^2 : b-c=t}.

If only the normalized values J_B(t)/q are stored, an absolute error below
1/(2q) is sufficient to recover every exact integer count and hence B.
An unspecified o(1) error in the normalized kernel does not provide this
integer-rounding guarantee. This is a resolution requirement for that
reconstruction procedure, not a claimed lower bound for every conceivable
method using all structural constraints.

### What an ordinary continuum limit erases

Let Q tend to infinity through multiples of four and put

    A_Q={0<=k<Q : k is even},
    B_Q={0<=k<Q : k mod 4 is 0 or 1}.

Both measures Q^(-1) sum_(k in A_Q) delta_(k/Q) and the analogous measure for
B_Q converge weakly to one half of Lebesgue measure on [0,1]. This follows
by grouping a continuous test function into consecutive blocks of four and
using uniform continuity. Nevertheless their cyclic shift-one correlations
are exactly 0 and 1/4:

    Q^(-1) #{k mod Q : k,k+1 in A_Q}=0,
    Q^(-1) #{k mod Q : k,k+1 in B_Q}=1/4.

Thus the same macroscopic density has two different microscopic pair
statistics. A useful repaired state retains (k/Q, k mod W), or equivalently
the relevant labelled residue/correlation data, before passing to a limit.

## 2. Prefix scaling has an exact finite-state arithmetic model

Fix integers b>=2 and W>=2. Let U be uniform on [0,1), set

    M_n=floor(b^n U),  R_n=M_n mod W.

Write L=b^h. The next h digits form D, uniform on {0,...,L-1} and independent
of M_n. Therefore

    M_(n+h)=L M_n+D,
    P^h(r,s)=L^(-1) #{0<=d<L : Lr+d=s mod W}.             (1)

This is a theorem about Lebesgue measure on real inputs. It does not assume
anything about the digits of a named constant. A fixed integer part can be
added to U without changing any of the bounds below.

**Theorem 2, exact mixing error.** The uniform distribution on Z/WZ is
stationary. Write L=qW+t with 0<=t<W. Every row of P^h has total-variation
distance exactly

    d_W(L)=t(W-t)/(WL) <= W/(4L)                         (2)

from uniform. In particular the row is exactly uniform iff W divides L.

Proof. The bW pairs (r,d), with 0<=r<W and 0<=d<b, enumerate the integers
br+d from 0 to bW-1 bijectively; every residue occurs b times, proving
stationarity. A translated interval of L integers contains q+1 entries in
t residue classes and q in the other W-t. Its half l1 discrepancy from
uniform is

    (1/2)[t((q+1)/L-1/W)+(W-t)(1/W-q/L)]
      =t(W-t)/(WL).

Translation changes which classes receive the extra entry, not this value.
The final bound is t(W-t)<=W^2/4. QED.

For 0<=f,g<=1 on Z/WZ, put bar f=W^(-1)sum f and similarly for g. The actual
prefix sequence starts at R_0=0, and Theorem 2 gives

    |E[f(R_n)g(R_(n+h))]-bar f bar g|
       <= d_W(b^n)+bar f d_W(b^h).                     (3)

Indeed replacement of the initial R_n law by stationarity costs at most
d_W(b^n); conditional expectation of g differs from bar g by at most
d_W(b^h), and averaging f costs only its mean bar f. This also proves the
stationary bound bar f d_W(b^h) without the initial term.

For a rare finite-sieve event of density delta, relative pair decorrelation
requires controlling [d_W(b^n)+delta d_W(b^h)]/delta^2, not just showing
b^h grows. This explicitly exposes the scale cost of increasing the wheel.

There is also an exact spectral explanation. If gcd(b,W)=1, put
Df(r)=f(r)-f(r+1) and Uf(r)=f(br). Telescoping gives bPD=DU.
On mean-zero functions D is invertible and U permutes the Fourier modes
j->bj, so P=b^(-1)DUD^(-1). Hence every nonconstant eigenvalue has modulus
1/b. The change of coordinates becomes poorly conditioned as W grows;
this spectral statement does not eliminate the explicit W/b^h mixing cost.

For arbitrary gcd(b,W), replace U by U_0, its mean-zero projection. Then
bPD=DU_0 still holds. On nonzero Fourier modes U_0 sends j to bj mod W,
with an arrow into zero killed. Thus cyclic components contribute 1/b times
roots of unity to the spectrum of P, while transient components contribute
zero eigenvalues and nilpotent parts. If every prime divisor of W divides b,
all nonconstant modes are nilpotent. This matches exact uniformity after
any h for which W divides b^h, directly verified by (2).

### The exact twin-center factor before averaging the digits

Let W be squarefree with prime divisors >=5, and define

    A_W={r mod W : gcd((6r-1)(6r+1),W)=1}.

For a prescribed digit block d and multiplier L, define

    C_W(L,d)=#{r mod W : r in A_W and Lr+d in A_W}.

CRT gives C_W(L,d)=product_(p|W) C_p(L,d). Put u=6^(-1) mod p.

**Proposition 3, affine local-factor classification.**

- If L=0 mod p, then C_p(L,d)=(p-2) 1_(d not in {u,-u}).
- If L=1 or -1 mod p, then C_p=p-2 at d=0; C_p=p-3 at
  d=2u or -2u; and C_p=p-4 at all other d.
- If L is nonzero and not 1 or -1 mod p, then C_p=p-3 when
  d belongs to {+(L-1)u,-(L-1)u,+(L+1)u,-(L+1)u}; otherwise C_p=p-4.

In the last case all four displayed residues are distinct.

Proof. For invertible L the forbidden r residues are {u,-u} together with
{L^(-1)(u-d),L^(-1)(-u-d)}. Their intersection counts the solutions of
Lx+d=y with x,y in {u,-u}. The four candidate values d=y-Lx are
the displayed residues. They coincide twice at zero precisely for L=+/-1;
otherwise they are distinct. When p divides L the second condition depends
only on d and the first has p-2 choices. QED.

The L=-1 case is exactly the inherited four-form Goldbach local count:

    C_p(-1,m)=p-rho_p(m),
    rho_p(m)=2 if p|m, 3 if p|9m^2-1, and 4 otherwise.

The L=+1 case is its difference analogue. Prefix scaling uses L=b^h, so
the actual affine factors must be retained rather than silently substituting
the sum/difference factor. Primes dividing b can even give a zero local
factor for a particular digit block; averaging all blocks restores the
uniform-residue model in Theorem 2. This is a concrete place where the carry
data matter.

Formula (3) also applies directly to the Goldbach wheel event

    f_m(r)=1_(r in A_W) 1_(m-r in A_W),
    bar f_m=product_(p|W)(1-rho_p(m)/p).

It therefore supplies a quantitative, unconditional cross-scale correlation
model for the exact local object already used in the repository.

This is a full-wheel estimate. It must be kept separate from the
dimension-four sieve used below: that sieve controls sums over truncated
divisors and does not require a sample longer than the complete primorial
period. A fixed-W mixing limit cannot be substituted for growing-W
equidistribution, nor can its exponential-in-cutoff period be silently paid
by a power-depth sieve theorem.

### What changes for one fixed real number

For deterministic alpha, M_n mod b is its nth fractional digit; M_n mod b^j
is its final j-digit block, once n>=j. Consequently

    M_n is uniformly distributed mod b^j for every j>=1
      iff alpha is normal to base b.                    (4)

This is the overlapping-block definition of normality, with at most j-1
initial terms omitted. Irrationality alone is insufficient: a base-b real
whose nonzero digits occur only at factorial positions is irrational, but
M_n mod b is zero with density one.

Thus a theorem asserting the stationary arithmetic distribution for the
particular lacunary sample floor(pi b^n) must not use the uniform-input
calculation as its proof. Even the subfamily of moduli b^j contains the
normality question. Equation (4) does not say that a fixed wheel coprime to
b requires full normality, nor that normality by itself establishes prime
correlations.

As a useful contrast, for irrational alpha the dense Beatty sequence
floor(alpha n+beta) is uniformly distributed modulo every fixed W. A short
proof is available without any normality assumption: for each nonzero integer
k, the normalized exponential sum of k(alpha n+beta)/W is a geometric sum
of size O_(k,alpha,W)(1/N), so fractional parts of (alpha n+beta)/W are
uniform; membership in residue r corresponds to the interval
[r/W,(r+1)/W). This dense sampling result does not apply to n=b^j.

The primary literature goes much further for single primes in finite-type
Beatty sequences: Banks--Shparlinski, “Prime numbers with Beatty sequences,”
https://arxiv.org/abs/0708.1015. Additive restrictions are studied by
Banks--Guloglu--Nevans, https://arxiv.org/abs/math/0701285. These citations
identify a promising dense-parameter route; no prime theorem from either
paper is needed in the proofs above. The distinction between almost every
input and a fixed familiar constant is also explicit in Khoshnevisan,
https://www.math.utah.edu/~davar/ps-pdf-files/NN2.pdf, Sections 2--3.

## 3. A real cusp window transfers integer arithmetic

Consider x_0=0, x_(j+1)=x_j^2+1/4+epsilon and fixed escape radius R>1/2.
Let N_R(epsilon)=min{j:x_j>R}. For each sufficiently large n there is a
unique positive threshold epsilon_n with x_n(epsilon_n)=R, and

    N_R(epsilon)=n iff epsilon_n<epsilon<=epsilon_(n-1).

The [sibling cusp analysis, Section 4](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md)
proves the stronger threshold statement

    n-pi/sqrt(epsilon_n) -> C_R,
    w_n:=epsilon_(n-1)-epsilon_n ~ 2pi^2 n^(-3).         (5)

The phase limit, not merely N_R~pi/sqrt(epsilon), is load-bearing in the
second assertion. Indeed put t_n=pi/sqrt(epsilon_n)=n-C_R+o(1).
Then t_n-t_(n-1)=1+o(1), and factor the difference
pi^2[t_(n-1)^(-2)-t_n^(-2)] to obtain (5).

The next result is an elementary general transfer theorem once these level
widths have been established.

**Theorem 4, relative transfer of arbitrary nonnegative statistics.** Let
M,L be positive integers with M tending to infinity and L=o(M), and let

    J={M+1,...,M+L},
    I=(epsilon_(M+L),epsilon_M].

For every nonnegative function F on J, uniformly over its choice,

    [integral_I F(N_R(epsilon)) d epsilon]/|I|
      = [1+o(1)] L^(-1) sum_(a in J) F(a),              (6)

where zero on the right means exact zero on the left. Also

    |I| ~ 2pi^2 L/M^3.                                (7)

Proof with an explicit distortion bound. Let

    eta_M=sup_(n>M) |w_n/(2pi^2 n^(-3))-1| ->0,
    K_(M,L)=((M+L)/(M+1))^3 (1+eta_M)/(1-eta_M).

For large M, eta_M<1, and max_(n in J) w_n/min_(n in J) w_n is at most
K_(M,L). Thus the ratio of the two sides in (6), whenever nonzero, lies
between 1/K_(M,L) and K_(M,L). This tends to one, proving a relative bound
even when the arithmetic event has very small density. The same pointwise
weight estimates prove (7). QED.

It is important that this transfers ALL nonnegative statistics, including
products of indicators, without replacing a joint statistic by products of
its marginals. No independence or arithmetic equidistribution was inserted.

### An unconditional rough twin-center Goldbach corollary

Fix u>37 and 0<theta<1. For each sufficiently large integer target m, put

    M=floor(m/3),  L=floor(m^theta),  z=L^(1/u),
    F_m(a)=(6a-1)(6a+1)(6(m-a)-1)(6(m-a)+1),
    V_m(z)=product_(5<=p<z)(1-rho_p(m)/p).

Then J lies in {1,...,m-1} eventually. Let E_m consist of epsilon in I for
which no prime p<z divides F_m(N_R(epsilon)). The repository's uniform
short-interval rough-count theorem and (6) imply

    limsup_(m->infinity)
      | |E_m|/[|I| V_m(z)]-1 | <= exp(37-u).            (8)

In particular |E_m| is positive for every sufficiently large m. The factor
V_m(z) is bounded below by a constant times (log z)^(-4), while an upper
bound of that same order is NOT uniform in m: exceptional targets can
enhance the local factor. Accordingly the universally valid conclusion is
the lower bound

    |E_m| >>_(u,theta,R) m^(theta-3)/(log m)^4,          (9)

and the more precise two-sided comparison is |E_m| asymp_(u,theta,R)
m^(theta-3) V_m(z). Keeping that local factor prevents a spurious uniform
two-sided estimate of order m^(theta-3)/(log m)^4.

Every epsilon counted in E_m supplies a decomposition m=a+(m-a) with both
centers having rough neighbors. Moreover, for any integer r>u/theta, all
four neighbors have at most r-1 prime factors, counted with multiplicity,
for sufficiently large m: they are less than 6m whereas z^r>6m eventually.

For a direct numerical scale, u=50 gives the same exp(-13) asymptotic relative
error as the inherited integer theorem. The constants and onset are not
optimized. This is a genuine arithmetic theorem about cusp windows, but its
arithmetic strength is precisely the inherited rough-number theorem.

The singular-series tail warning in the fourth earlier note survives this
change of coordinates: V_m(z) must not be replaced by the full S_4(m) without
its tail criterion. Scaling does not erase primes beyond the actual sieve
cutoff.

### Why the width is the relevant precision

At N of order M, one level cell has width asymptotic to 2pi^2 M^(-3), and
the parameter itself is of order pi^2 M^(-2). Resolving adjacent typical
levels therefore needs absolute parameter resolution of order M^(-3), or
relative resolution of order M^(-1). There is no positive uniform error
tolerance at every input: an input can be arbitrarily close to a threshold.

For epsilon=b^(-2n), this scale is b^(-3n). The point sample retains only one
level at each n. The window theorem concerns many neighboring levels and
does not assert that those geometric point samples hit the arithmetic set.

There is a precise deterministic certification condition. If
|pi/sqrt(epsilon_k)-(k-C_R)|<=eta for all relevant thresholds and
dist(pi/sqrt(epsilon)+C_R,Z)>eta, the threshold intervals give

    N_R(epsilon)=ceil(pi/sqrt(epsilon)+C_R).

Near the excluded boundary, the integer must be determined from the actual
orbit or a sharper threshold enclosure. A leading asymptotic alone cannot
determine even parity: replacing any integer estimate K by 2 floor(K/2)
changes it by at most one while forcing an even output.

## 4. Which normalizations preserve the poset correlation inequalities?

The repository's fixed-interface attachment packet uses deletion counts
H(S)>0, including diamond ratios

    D_H(S;a,b)=H(S)H(S+a+b)/[H(S+a)H(S+b)].

For actual available-maxima diamonds, the imported Fishburn and Chan--Pak
bounds put this ratio between n/(n-1) and 2. The present statement is about
normalization of such data, not a new correlation inequality.

**Theorem 5, exact scaling-gauge classification on a Boolean lattice.** Let
c(S)>0 for every S subset of a finite set T, and set H'(S)=c(S)H(S).
Every labelled diamond ratio is unchanged, for any positive H, iff

    c(S)=c(empty) product_(i in S) t_i                 (10)

for positive numbers t_i.

Proof. The transformed ratio is the old ratio multiplied by
c(S)c(S+a+b)/[c(S+a)c(S+b)]. Thus invariance is equivalent to all these
multipliers being one. Put v(S)=log c(S). The resulting diamond identities
say the increment v(S+a)-v(S) is unaffected by adjoining any other element.
Deleting the elements of S one at a time shows that this increment is
v({a})-v(empty). Summing increments proves (10). The converse is immediate
by cancellation. QED.

Common scaling is allowed; coordinatewise arbitrary mantissa scaling is
not. In base b, if H(S)=b^(e_S) h(S), then the exact repair is

    D_H(S;a,b)=b^(e_S+e_(S+a+b)-e_(S+a)-e_(S+b)) D_h(S;a,b).  (11)

**Natural exact hostile.** For a four-element antichain and two specified
maxima, the counts are (H00,H10,H01,H11)=(24,6,6,2), with ratio 4/3.
Independently placing each integer in (0,1) by its decimal digit count gives

    (h00,h10,h01,h11)=(0.24,0.6,0.6,0.2),
    D_h=2/15.

The normalized array appears to violate even the lower bound one. The
exponent vector (2,1,1,1) gives correction 10^(2+1-1-1)=10, restoring 4/3.
This is a literal loss of the required scale coordinate in an actual
correlation array. Retaining common scale or the exponent correction makes
the operation lossless.

This complements the earlier packet's different hostile: arbitrary mixtures
of normalized arrays can destroy nonlinear realizability. Exact rescaling,
independent normalization, and mixing are three different operations and
should not be substituted for each other.

## 5. Inheritance, sources, and next usable hypotheses

Inherited repository inputs (read in full at their mathematical sections):

- [prime_pairs_local_symmetry_20261019/prime-pairs-structural-note.md](../prime_pairs_local_symmetry_20261019/prime-pairs-structural-note.md):
  symmetric sum/difference correlation, wheel root counts, rough-number
  theorem, and the separate actual-prime correlation gate.
- [prime_pairs_local_symmetry_20261019/forbidden-pattern-classification.md](../prime_pairs_local_symmetry_20261019/forbidden-pattern-classification.md):
  odd-order parity reconstruction and its even-order failure.
- [prime_pairs_local_symmetry_20261019/power-depth-rough-pairs.md](../prime_pairs_local_symmetry_20261019/power-depth-rough-pairs.md):
  classical sieve at u>37 and explicit relative-error envelope.
- [prime_pairs_local_symmetry_20261019/short-interval-local-resolution.md](../prime_pairs_local_symmetry_20261019/short-interval-local-resolution.md):
  location-uniform short-interval theorem and sharp full-series tail warning.
- [boundary_attachment_correlation_20261019/PASSIVE_ATTACHMENT_HULL_AND_FKG.md](../boundary_attachment_correlation_20261019/PASSIVE_ATTACHMENT_HULL_AND_FKG.md):
  fixed-interface deletion arrays, nonlinear realizability, and their
  normalization conventions.

These paths are relative to math/05-knowledge/results. No arbitrary
tournament orientation has been added: the genuine objects are labelled
correlation kernels and deletion diamonds.

Three concrete next questions, with the missing input named:

1. For a particular constant such as pi, bound the discrepancy of
   floor(pi b^n) modulo the chosen finite wheel, and its avoidance of the
   shifted cusp thresholds. This is a precise deterministic substitute for
   the random-digit model, not a consequence of the escape asymptotic.
2. Quantify the phase-tail eta_M in (5), to make the cusp-window distortion
   and exact-integer certification effective at a declared scale.
3. If the intended goal is actual primes, prove the already explicit
   four-form correlation discrepancy estimate. All the maps here preserve
   that obligation; none manufactures it from one-point density.

The companion exact-control script verifies the finite formulas independently
by enumeration and rational arithmetic. Its inverse-square weight controls
test the transfer algebra, not the analytic proof of the actual cusp phase.
Both normal and optimized Python runs pass 148,202 active guards and produce
the identical JSON receipt. Run `python verify_scaling_correlations.py` and
`python -O verify_scaling_correlations.py` from this directory.
