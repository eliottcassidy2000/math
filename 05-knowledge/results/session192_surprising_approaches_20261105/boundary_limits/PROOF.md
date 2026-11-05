# Fixed-feedback growth has two boundary laws, not concentration

**PROVED by the elementary arguments below; independent audit pending.**
Session 192, 2026-11-05. No literature-priority claim. The exact `{7,21}`
Hamilton-path spectrum remains open.

## 0. What is inherited, and what changes

For a fixed tournament Q on r labelled apex vertices, attach an ordered
transitive chain x_1,...,x_m. A letter c_i in {0,1}^r records which apexes
point to x_i. Write T_Q(w) for this physical tournament and H_Q(w) for its
number of directed Hamilton vertex words. The parameter m is the **chain
length**, not the total order m+r. Throughout r>=1 is fixed.

The exact run-assignment representation, the bound H<=r!(r+1)^m, the mean
E H ~ r!(r+1)^m/4^r, and the asymptotic strong probability (1-2^-r)^2 are
inherited from Session26's fixed-apex and random-attachment packets. Session26
also already proves the one-apex strong product-of-uniforms law. We do not
claim these again as new results. Session191's 36-state two-apex transfer and
its exact minimal rational dimension 17 are the nearest current mechanism.
Its factorial-growth strong law concerns **uniform random tournaments**,
not the fixed-feedback model here.

New results in this packet:

1. A uniform deterministic rank-one boundary factorization with an explicit
   error, valid for every word and every fixed core.
2. A general fixed-r product-boundary limit, a three-coordinate affine
   system for r=2, and an exact finite moment solver.
3. Explicit new two-apex limiting moments, including under strongness.
4. Almost-sure exponential growth but persistent normalized oscillation for
   one nested random attachment sequence.
5. Exact two-apex lower-tail and negative-moment thresholds, separating
   neutral padding from intrinsic strong boundary behavior.
6. A sharp core-edge perturbation bound and a hostile showing precisely
   where this asymptotic compression loses exact spectrum information.

Concept board: run slots; coupon/transport coupling; random affine systems;
exact versus approximate Hankel rank; rare boundary words. The unrelated
repository prompts are the information-loss/merge law in THM-4481
(`entropy-merge-law-sign-strategies`), factorial-Hankel faithfulness in
THM-2810 (`factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go`),
and the finite-versus-global warning in THM-4487
(`dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence`). Their
particular Collatz or Gaussian theorems are not asserted for tournaments.
The actual transported mechanism is made explicit below.

## 1. The boundary statistic

Put s=r+1 and alpha=r/s. Independently assign every chain name to one of
s run labels 0,...,r, uniformly. Label 0 is free. For each j>=1, require the
first name assigned to label j to have c_{i,j}=1. Define:

* F_r(u): probability that every constrained label j>=1 appears in the
  finite word u and its first name has the required bit 1;
* G_r(u): the same probability, additionally requiring free label 0 to appear.

These are rational numbers with denominator s^|u|, obtained by a 2^s-state
occupancy process. For an infinite attachment word omega, define

    U_r(omega) = lim_k F_r(omega[1..k])
              = lim_k G_r(omega[1..k]).                       (1)

Both finite probabilities are nondecreasing as k grows. Once a successful
allocation has seen all required labels, every continuation remains
successful. The probability that a named allocation label is still unseen
is alpha^k. A union bound therefore gives, **for every infinite word**,

    0 <= U_r-F_r(prefix_k) <= r alpha^k,
    0 <= U_r-G_r(prefix_k) <= s alpha^k.                       (2)

This is a coupling bound, not a stationary-distribution assertion: an
unseen label is the only location at which the tail can change acceptance.
In particular U_r is a continuous function of the infinite letter sequence.
Let u* be reversal of u with every bit complemented.

## 2. Uniform rank-one boundary theorem

Take any factorization w=u z v, with |u|=p, |v|=q and p+q<=m. Then

    H_Q(w)/s^m = r! G_r(u) G_r(v*) + E_Q(w;u,v),
    0 <= E_Q(w;u,v) <= r! s (alpha^p+alpha^q).               (3)

**Proof.** Fix an order a_1,...,a_r of the apex names. Every Hamilton word
has a unique representation

    L_0, a_1, L_1, ..., a_r, L_r,

where the chain names assigned to each L_j are in increasing order. Suppose
every slot meets both u and v. The first-name boundary requirements are
entirely in u, and the last-name requirements entirely in v. No two apexes
are adjacent, so no core arc is consulted. Middle names can be assigned
arbitrarily. The initial acceptance probability equals G_r(u), independent
of the apex order, by permuting the constrained slot labels. Reversing and
complementing v makes its terminal acceptance probability G_r(v*), again
independent of the apex order. Thus these good assignments give the displayed
rank-one term, for each of r! apex orders.

Every remaining assignment omits some slot from u or from v. For each named
slot these probabilities are alpha^p and alpha^q. There are s slots and r!
orders. Some remaining assignments fail to be paths; this only decreases the
error. The union bound proves (3). QED.

Three exact consequences:

* Changing an arbitrary middle block while retaining boundary words of length
  k changes H by at most 2 r! s r^k s^(m-k). The two values lie above the same
  leading integer and below the same upper endpoint, so no extra factor two
  is needed.
* For prefixes all of length p and suffixes all of length q, the matrix
  H_Q(uv)/s^(p+q) is within delta=r!s(alpha^p+alpha^q) entrywise of one
  nonnegative rank-one matrix. Every 2x2 minor has absolute value at most
  4r! delta, by expanding one factor difference at a time. This is an
  approximate-rank statement; the exact two-apex Hankel dimension remains 17.
* The leading matrix is independent of Q. This is deterministic core
  universality, stronger in quantifiers than independence of a leading mean.

## 3. Random letters: the product law

Now let omega_i be iid uniform letters in {0,1}^r. Write U for U_r(omega),
and let V be an independent copy. For every fixed core Q,

    H_Q(omega_1...omega_m)/s^m  ==>  r! U V.                (4)

All positive integer moments converge too. Indeed all variables lie in
[0,r!]. Split the word into lengths p=floor(m/2), q=m-p. The initial word
and reverse-complement terminal word are independent uniform words. Extend
each independently to an infinite word, giving independent U,V. Combining
(2)-(3) supplies an explicit coupling

    |H_Q/s^m-r!UV| <= 2 r! s (alpha^p+alpha^q).             (5)

Consequently q-th moment error is at most q(r!)^(q-1) times the right side
of (5). Constants are deliberately conservative. This is a genuine bounded
coupling, hence does not need an unknown density or an anti-concentration
hypothesis for U.

For any fixed allocation, the first appearances of the r constrained labels
occur at distinct indices, and their required bits are independent fair
bits. Therefore

    E U = 2^-r,       E(r!UV)=r!/4^r.                     (6)

U is positive almost surely: an infinite random word contains r letters
with all bits 1, and a finite successful assignment can use them for the r
first constrained appearances. Also U is not constant. Conditional on the
first letter being all zeros, U=U_tail/s; this event has positive probability
and U_tail has the same positive mean. A positive constant could not satisfy
that identity. Hence the limiting law in (4) is nondegenerate for every r.

### Strong tournaments

Let B_m be the event that the first letter is nonzero and the last letter is
not all ones. This is necessary for strongness, and has probability
(1-2^-r)^2 for m>=2. The inherited Session26 elementary argument gives

    P(B_m but not strong) <= r(r-1)(3/4)^m:                (7)

with this failure probability all ordered apex pairs have a two-step path
through a chain name. Then the apexes share an SCC, and the two endpoint
conditions make that SCC the whole tournament.

Thus (4), conditioned on strongness, converges to r! U_+ V_+, where U_+,V_+
are independent copies of U conditioned on their first letter being nonzero.
There is no assertion of endpoint independence conditional on strongness at
finite m. Equation (7) controls the replacement of that conditioning.

If M_q=E U^q, the zero-first-letter identity gives exactly

    E U_+^q = M_q (1-2^-r s^-q)/(1-2^-r).                 (8)

The strong limit is therefore explicit at the level of every moment.

## 4. Bernoulli convolution and a three-coordinate two-apex system

For r=2 define A_j as the probability that the first occurrence of one
specified constrained slot has bit j equal to 1. Then

    A_j(omega) = (1/3) sum_(i>=1) (2/3)^(i-1) omega_{i,j}.

A_1 and A_2 are independent Bernoulli-convolution coordinates with contraction
2/3. We make no absolute-continuity or smoothness claim about these measures.
If c=(c_1,c_2) is the first letter, primes denote the tail, and U is the
joint two-slot boundary probability, then

    A_1=(2 A'_1+c_1)/3,
    A_2=(2 A'_2+c_2)/3,
    U  =(U'+c_1 A'_2+c_2 A'_1)/3.                         (9)

The first allocation chooses the free slot (leaving both requirements), or
one of the two constrained slots (checking its bit and leaving the other
requirement). This proves (9), a concrete affine iterated system, rather than
an analogy between spectra. The U coordinate records joint success rather than multiplying the two
single-slot probabilities. For example, the deterministic infinite word
3 followed by zeros has A_1=A_2=1/3 but U=0: one chain name cannot be the
first occurrence of two distinct allocation slots. No claim that U is or
is not a function of A_1,A_2 alone is required.

For general r, retain one coordinate U_S for each subset S of still-required
apex bits. Put U_empty=1. The exact recursion is

    U_S(c omega) = [(s-|S|)U_S(omega)
                    +sum_(j in S)c_j U_(S minus j)(omega)]/s.             (10)

This is a 2^r-dimensional affine boundary system (including the constant).
It is not an exact 2^r-state replacement for H. The full count still needs
its lower-scale continuation coordinates.

## 5. Every limiting moment is a finite rational computation

For the q-th moment, take q independent allocation streams reading the same
random letters. A state is (S_1,...,S_q), with S_a the still-unseen constrained
labels in stream a. At a step choose independently j_a in {0,...,r}.
If j_a is an unseen constrained label, remove it and request that bit of the
current random letter to equal 1. If D is the set of distinct requested bits,
that transition survives with probability

    s^-q 2^(-|D|).                                        (11)

Let P(S) be the probability of eventually satisfying all remaining requests.
P(empty,...,empty)=1. Every non-self transition strictly decreases the total
number of missing labels. The self-loop probability is

    d(S)=product_a (s-|S_a|)/s^q < 1

away from the absorbing state. Therefore

    P(S) = sum_(S' proper descendant) p(S,S')P(S')/[1-d(S)].              (12)

Equations (11)-(12) are a finite triangular rational algorithm with at most
2^(rq) states, and M_q=P(full,...,full). Coincident requests share one bit;
replacing their probability by independent 1/2 factors would be wrong.

For two apexes the exact values are

    M_1=1/4,
    M_2=29/320,
    M_3=1177/29440,
    M_4=237365069/11787776000.

Hence the unconditioned path limit Z=2UV has

    E Z=1/8,
    E Z^2=841/25600,
    Var Z=441/25600 > 0.                                (13)

For strong tournaments the corresponding limit has

    E Z_strong=121/648,
    E Z_strong^2=41209/746496,
    Var Z_strong=136625/6718464 > 0.                      (14)

The positive variances disprove any transfer of Session191's relative
concentration to this ensemble. The appearance of 21^2 in one rational
variance numerator has no established relation to the missing Hamilton count
21 and supplies no spectrum obstruction.

Controls: for r=1 the algorithm gives M_q=1/(q+1), recovering the inherited
uniform boundary law. For r=3 it gives M_1=1/8, M_2=79/3360 and
M_3=951449/167336960. The two-apex second-moment constant is independently
checked by the tensor square of the exact 36-state path transfer.

### 5.1 Exact second-moment spectrum and a full-vector certificate

Let V_m=E H_Q(w)^2 for a two-apex core and uniform length-m letters. Reversing
the sole core arc merely relabels the iid attachment ensemble, so either core
has the same V_m. In addition to the limiting calculation above, for m>=3:

    V_m = (841/25600)9^m
          -[5(10m-109)/3072]6^m
          +[(14255m-71716)/345600]4^m
          -[(3m-89)/144]3^m
          -[(4m^3-24m^2+68m-687)/3072]2^m
          +(1536m^2-5656m+9365)/27648.                  (23)

Thus V_m/9^m approaches 841/25600 with error O(m(2/3)^m), much sharper than
the general coupling bound. The polynomial factors record genuine repeated
spectral channels, not a single-run Hamilton-count recurrence.

The proof is finite exact linear algebra with an all-order certificate. Let
M_c be the 36-state integer matrices defined in verify.py, let v be the tensor
square of its initial vector, and let K=sum_c M_c tensor M_c. The supplied
verify_spectral.py checks the **entire 1296-coordinate vector identity**

    K^3(K-36I)(K-24I)^2(K-16I)^2(K-12I)^2
       (K-8I)^4(K-4I)^3 v = 0.                         (24)

Since 4^m V_m is the accepting row applied to K^m v, (24) gives its scalar
recurrence at every order. Dividing the spectral roots by 4 gives the channels
in (23). The root zero has multiplicity three, so only the first three terms
need a separate head. Fourteen consecutive exact values at m=3,...,16
identify the coefficients of the fourteen-dimensional positive-root solution;
the script checks them and further terms through m=100. Therefore (23) is not
extrapolated from a fitted scalar recurrence. All tests use integers and
fractions, with no numerical eigenvalue or rank threshold.

## 6. Nested random growth: an exponent, but no normalized strong law

Fix one iid infinite attachment sequence, and enlarge the tournament by
successive chain prefixes. Almost surely,

    log H_Q(omega[1..m])/m -> log s,                      (15)
    liminf_m H_Q(omega[1..m])/s^m = 0,
    limsup_m H_Q(omega[1..m])/s^m = r! U_r(omega).          (16)

These statements hold simultaneously for every core on the r apex names,
since there are finitely many such cores. They do not condition the whole
nested sequence to be strong.

**Growth proof.** Almost surely a finite initial block has G_r>0. For a
terminal window of length k_m=ceil(C log m), the chance that it contains fewer
than r all-zero letters is at most a constant times
k_m^(r-1)(1-2^-r)^(k_m-r+1). Choose C sufficiently large to make this summable.
Borel-Cantelli implies that all sufficiently late terminal windows contain r
such letters. In the reverse-complement window, assign those r all-one
letters to the constrained slots, and some other name to the free slot;
assign all remaining names to the free slot. This is a valid first-boundary
allocation, so G_r>=s^-k_m. Apply the nonnegative main term of (3) with the
fixed successful prefix, this terminal window, and the intervening middle.
Thus H_Q/s^m>=c s^-k_m. The upper bound H_Q<=r!s^m gives (15).

**Oscillation proof.** Fix a prefix length p and let the suffix grow.
Equation (3), with G_r<=1, gives limsup H_Q/s^m<=r!G_r(prefix_p)+r!s alpha^p.
Let p tend to infinity. This gives the upper bound r!U.

Every finite all-zero terminal block occurs infinitely often almost surely,
by independent disjoint block trials. Let its length tend to infinity along
a subsequence. Its reverse-complement is all ones, so G_r tends to 1 by the
coupon union bound. The lower half of (3), first with p fixed and then
p tending to infinity, gives limsup>=r!U.

Similarly arbitrarily long all-one terminal blocks occur. In any accepted
Hamilton word, at least one of the r slots immediately preceding an apex
must be absent from such a block, since a last boundary there would point in
the wrong direction. Thus H_Q/s^m<=r! r alpha^k for a terminal all-one block
of length k. This tends to zero. QED.

Since U>0 almost surely, (16) is a strict failure of convergence, not merely
a failure of our proof method. In contrast, the distributional law (4) still
holds: the right boundary is renewed as the chain grows, while the left
boundary converges to a fixed random value.

## 7. Exact core edits and the hostile lower exponential scale

If Q and Q' differ in t named core arcs, then, for every attachment word w,

    |H_Q(w)-H_Q'(w)| <= t (r-1)! r^m.                     (17)

It suffices to reverse one arc a->b. Hamilton words not using that adjacency
are unchanged. Words lost must have a,b adjacent in that order; there are
(r-1)! possible apex orders with that adjacency, and at most r^m assignments
of chain names to the remaining run slots. Words gained obey the same bound
with b,a. The absolute difference of the gained and lost counts is bounded
by their common upper bound. Telescope over the t edits.

This is sharp for r=2 and m>=1: take w=1^m, meaning a->chain->b. With a->b the whole
tournament is transitive and H=1. Reverse the core arc: the cyclic blow-up
C3(1,m,1) has H=2^m+1, by assigning an arbitrary subset of the chain before
one apex and its complement after the other, plus the one remaining path.
The difference is exactly 2^m.

Both normalized counts tend to zero on the 3^m scale. The leading boundary
factor is zero, yet the subleading difference is unbounded. This is the
canonical hostile for turning boundary universality into exact-H equivalence.
Likewise the bound (3) cannot resolve unit or +2 gaps: with nonoverlapping
boundary words, its guaranteed absolute error is exponentially large.

## 8. A Mellin/renewal bridge: two different lower-tail exponents

There is a further exact distinction between neutral padding and intrinsic
strong boundary difficulty. Specialize to r=2 and put

    beta=log(4)/log(3),      gamma=log(2)/log(3/2).

Then beta<gamma, since (3/2)^2<3. For the two-apex boundary variables,

    P(U<=epsilon)=Theta(epsilon^beta),
    P(U_+<=epsilon)=Theta(epsilon^gamma)                  (18)

as epsilon tends to zero. Equivalently, for positive t,

    E U^(-t)<infinity  iff t<beta,
    E U_+^(-t)<infinity iff t<gamma.                      (19)

The critical moments diverge. This is not a statement about exact attained
integer counts.

**Proof of the strong boundary exponent.** Let A be one of the one-slot
Bernoulli coordinates in (9), and alpha=2/3. If its first 1-bit occurs at
position L, then

    alpha^(L-1)/3 <= A <= alpha^(L-1),
    P(L=l)=2^-l.

Therefore P(A<=epsilon)=Theta(epsilon^gamma). Conditional on any nonzero
first letter, (9) gives U_+>=A_tail/3 for one appropriate coordinate with
the law of A. Conversely, conditional on first letter 1 (probability 1/3
among the nonzero letters), the joint success U is at most the single-slot
success A_2=(2/3)A'_2. Hence

    (1/3)P(A<=3epsilon/2) <= P(U_+<=epsilon)
                           <= P(A<=3epsilon).

This proves the second tail estimate and its moment threshold by integrating
P(U_+<=epsilon) epsilon^(-t-1).

**Neutral-prefix renewal.** Let K be the number of initial zero letters.
Then P(K=k)=(3/4)4^-k, and deleting that prefix leaves an independent word
conditioned on a nonzero first letter. The first-letter equation gives the
exact distributional identity

    U = 3^(-K) U_+,                                     (20)

with independence on the right. Consequently, wherever finite,

    E U^(-t) = [(3/4)/(1-3^t/4)] E U_+^(-t).             (21)

Since beta<gamma, this proves (19), including divergence at beta. For the
first tail in (18), split the geometric sum in (20) at
k=floor(log_3(1/epsilon)). Its large-k tail is Theta(epsilon^beta), while
the small-k part is bounded by a geometric sum using the exponent gamma>
beta. A lower bound follows from U_+<=1. This proves (18).

The real Mellin singularity in (21) is therefore produced by source padding,
not by a mysterious property of the Hamilton spectrum. The exponent beta
is exactly the ratio of four-letter input entropy to ternary bulk growth.
Deleting neutral padding changes the tail exponent to gamma. This provides
a precise entropy/renewal connection with a missing-coordinate warning.

For independent copies, an elementary product-tail integration also gives

    P(2UV<=epsilon)=Theta(epsilon^beta log(1/epsilon)),
    P(2U_+V_+<=epsilon)=Theta(epsilon^gamma log(1/epsilon)).               (22)

For completeness, if c x^a<=P(X<=x)<=C x^a near zero, condition on the
second factor Y. On Y>=epsilon the integrand is comparable to
(epsilon/Y)^a; integration by parts gives
E[Y^(-a);Y>=epsilon]=Theta(log(1/epsilon)). The omitted part Y<epsilon
is O(epsilon^a). This proves (22).

A concrete surviving route is now visible: because the strong boundary has
a finite negative beta-moment, a lattice renewal count across chain lengths
may control bounded-H strong two-apex presentations at scale X^beta. The
uniform additive approximation (3) alone does not justify that passage at
the rare-event boundary, so no such counting asymptotic is asserted here.

## 9. Connection contract and surviving research directions

* **Transport/coupling -> path growth:** couple allocation histories with the
  same first-seen labels. The preserved object is boundary acceptance and
  its leading path mass; missing slots carry all error. Exact integer H and
  rare low-growth words are lost. The necessary sidecar is the full
  continuation response, 17-dimensional in the two-apex language.
* **Information-loss laws -> grammar design:** old input bits become invisible
  to leading normalized H once all relevant first/last slots are occupied.
  This is a proved exponential influence bound, not an imported Collatz
  entropy inequality. An exact spectrum grammar must exploit the small
  channels instead of mistaking asymptotic compression for exact state loss.
* **Affine probability systems -> a reusable model:** (9)-(12) compute the
  genuine boundary distribution and all rational moments. They preserve
  large-scale distributional growth; they do not prove an interval of
  attained integers or the absolute continuity of U.
* **Hankel faithfulness -> a numerical warning:** uniformly near-rank-one
  continuation matrices can still have exact rank 17. A spectral compression
  is a fast coarse estimator, not an exact low-state ternary constructor.

Concrete next questions: classify the small-channel state needed after
subtracting the rank-one term; derive a useful anti-concentration bound for
U_2 without assuming a Bernoulli-convolution density; and use the sharp core
edit (17) to search for cancelling 2^m-scale changes between equal leading
boundaries. None is answered by the product law alone.

## 10. Reproduce and source scope

Run `python3 verify.py`, `python3 -O verify.py`, `python3 verify_spectral.py`
and `python3 -O verify_spectral.py`; deterministic receipts
must match. The verifier is self-contained and uses only the standard library.
It checks physical Hamilton counts by subset DP, literal run allocation,
all split boundary inequalities in its declared small universe, the affine
recursion, the moment solver against a separate finite-time absorption
iteration, tensor-square path moments, and the sharp core-flip hostile.

Relevant inherited sources:

* [Session 191 two-apex transfer](../../session191_tournament_growth_20261104/two_apex/TWO_APEX.md).
* [Session 191 uniform-tournament limit](../../session191_tournament_growth_20261104/random_growth/PROOF.md), a contrasting ensemble not used to prove the fixed-apex results here.
* [Session 26 fixed-apex and slot work](../../session26_random_bridges_20261022/session26-fixed-apex/proof.md), [random-attachment transform](../../session26_random_bridges_20261022/session26-attachments/proof.md), and [one-apex product law](../../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md).

The statements above were independently reviewed. No claim of global spectrum completeness follows.
