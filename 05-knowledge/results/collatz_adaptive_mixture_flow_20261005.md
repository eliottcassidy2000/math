# Adaptive Collatz weights from a positive difference array

**Status: PROVED scoped constructions and bounds; FINITE-EXACT checks;
positivity at every integer and universal Collatz remain OPEN.**

The main result is a summable incoming flow with root weight 1, total mass
at most **16/3**, and no chosen fixed exchange rate between its two costs.
Its weights have a factorial formula, exact infinite-fibre sums, and at most
a polynomial loss against the best fixed exchange rate for each certified
source. It is positive exactly on the rooted component. This improves the
weights; it does not establish that the rooted component contains every input.

## 1. Inheritance, portfolio, and the retained obstruction

The closest mechanism is **C8** in
[three bits and critical flow](collatz_three_bits_critical_flow_20261005.md):
the canonical sibling-base flow is unconditionally summable, and its positive
support is exactly the rooted component. The corrected near miss is the
distinction between this support statement and a proof of positive support
everywhere. The hostile is the same note's **C7**: 107->161->121->91 with an
extra positive predecessor 429 defeats the specified coarse golden observer.
The least-used relevant sidecars here are the removed root loop, its first
nontrivial child 3, and the full source/word accompanying a compressed bill.

The concept board is:

| Lane | Object | Decisive question |
|---|---|---|
| Anchor | Critical incoming flow | Can its positive support be proved universal? |
| Niche | Mixtures of legal price systems | Can one avoid exponential loss from a poorly chosen price? |
| Wildcard | Positive finite differences / Pascal array | Can the entire incoming fibre be paid by a telescoping identity? |
| Bridge | Checked common-future controller | Does its literal certificate transport the new weight? |
| Hostile | Binary density and discarded order | Does a proposed transfer retain the original integer and legal word? |

The binary-cylinder idea was recovered, rather than counted as a new result:
[THM-4501, recursive motif families and frequency](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
already proves that the positive basin of any target coprime to 3 is dense
in the 2-adics. Its proof extends a prescribed shortcut parity prefix by
solving an exponent congruence modulo a power of 3. It constructs another,
possibly enormous source in the cylinder; it does not certify the given
source. This is the relevant failure boundary for a continuity shortcut.

## 2. Objects and the fixed-price flow

Work on positive odd integers with U(n)=oddpart(3n+1), S(n)=4n+1.
Every n uniquely has n=S^j(b), where b is a base with v2(3b+1) in {1,2}.
For b!=1 write U(b)=S^k(G(b)). Kill the actual edges entering 1; the
incoming operator on V={3,5,...} is

    (K v)(m) = sum_(n in V : U(n)=m) v(n).

The root value below is bookkeeping, not an inequality on its self-loop.
For 0<r<1, put d_k=(1-r)r^k. On a base whose first-hit G-chain reaches 1,
let g_r be the product of its edge factors d_k, with g_r(1)=1. On every
other base put g_r=0. Extend f_r(S^j b)=r^j g_r(b).

For a rooted n define L(n) to be its number of base edges, and K(n) to be
the total sibling depth, including the initial j. Then

    f_r(n) = (1-r)^L(n) r^K(n).                              (1)

The counters can also be read from a strict actual ROOT word
(a_1,...,a_tau): for n>1,

    L=tau-1,       K=sum_i floor((a_i-1)/2).                 (2)

At n=1 both are zero. Every rooted n>1 has K>=1: its last nonroot odd
predecessor of 1 has valuation at least 4. The counters require a verified
completed certificate; they are not independently known functions of the
input size. On an unrooted input (1) is replaced by zero, with no invented
finite counters.

Inherited fibre elimination proves K f_r<=f_r. In fact equality holds at
every nonroot target not divisible by 3; targets divisible by 3 have no
incoming odd edge. This also holds on unrooted components, where all these
weights vanish.

## 3. Two root generations sharply improve the mass bound

**P1 — PROVED.** For every 0<r<1, without assuming universal convergence,

    sum_b g_r(b) <= 2+2r-r^2,
    sum_(positive odd n) f_r(n) <= (2+2r-r^2)/(1-r).          (3)

Proof. Put D=1+r+r^2. At a base c the inverse child of depth k exists
exactly when 3 does not divide S^k(c). Since S^k(c)=c+k modulo 3, one
class k=a mod3 is forbidden. The sum of its child factors is

    1-r^a/D <= kappa := 1-r^2/D < 1.

At the root omit k=0, which would be the self-loop. Its first generation
has total mass

    A = r(1+r^2)/D.

Its k=1 child is exactly 3, with weight e=(1-r)r. At 3 the forbidden
depth class is 0, so its outgoing factor sum is B=(r+r^2)/D.
All other first-generation vertices still have factor sum at most kappa.
Thus the second generation has mass at most

    A_2 <= e B + (A-e) kappa.

Every subsequent generation loses a factor at most kappa. Summing gives

    sum_b g_r(b) <= 1+A+A_2/(1-kappa) = 2+2r-r^2.

Finally each base's sibling fibre contributes g_r(b)/(1-r). The proof
counts only the inverse tree rooted at 1 and makes no assertion that it
contains every base. Its two exceptional starting generations are why the
bound stays bounded as r tends to zero.

For comparison, the earlier valid bound was D/[r^2(1-r)]. At r=1/16 the
new full-flow bound is 181/80; at r=1/4 it is 13/4. This sharpens that
bound rather than retracting it.

## 4. A universal positive array and its arithmetic realization

Mix the existing valid flows using the probability density 2(1-r) dr:

    W(n) = integral_0^1 2(1-r) f_r(n) dr.                   (4)

**P2 — PROVED.** W(1)=1, K W<=W, and sum_n W(n)<=16/3. At a rooted
source the exact value is

    W(n)=w(L,K),       w(L,K)=2 K! (L+1)!/(L+K+2)!.          (5)

It is zero at every unrooted source. Its positive support is exactly the
rooted component.

Proof. The integral in (4) is 2 B(K+1,L+2), yielding (5). Nonnegative
integration preserves each incoming inequality. Tonelli and (3) give

    sum_n W(n) <= integral_0^1 2(2+2r-r^2) dr = 16/3.

No interchange needs an unproved orbit bound. A rooted monomial is positive
throughout (0,1); an unrooted source has zero integrand throughout.

The formal array w is strictly positive on **every** pair (L,K) of
nonnegative integers. It has the exact two-way split

    w(L,K)=w(L+1,K)+w(L,K+1),                              (6)
    w(L+1,K)/w(L,K)=(L+2)/(L+K+3),
    w(L,K+1)/w(L,K)=(K+1)/(L+K+3).

This is an adaptive price rule: the two next fractions sum to 1 and depend
on accumulated counter values. Their product recovers the same weight in
any formal order. Order independence pays a bill; it does not preserve a
Collatz source guard.

**P3 — PROVED exact fibre payment.** Summing (6) along K gives

    sum_(j>=0) w(L+1,K+j)=w(L,K),                          (7)
    sum_(j>=J) w(L+1,K+j)=w(L,K+J).

For an actual nonroot target m coprime to 3 with counters (L,K), its
minimal inverse base has counters (L+1,K), and its successive siblings
have (L+1,K+j). Thus (7) is exactly its entire incoming Collatz row,
including an explicit remainder after any finite truncation.

There is also a literal difference-family reading. If Delta acts in K by
Delta h(K)=h(K+1)-h(K), then

    w(L,K)=(-Delta)^L [2/((K+1)(K+2))].                    (8)

Every iterated signed difference is positive because it is the integral
of r^K(1-r)^L against 2(1-r) dr. The initial sequence is completely
monotone in precisely this elementary sense. This supplies a rigorous
version of a difference family with two compatible operations: shifting
K and taking its positive difference. Its Pascal identity is not a claim
about the fixed-seed Fourier Pascal tower elsewhere in the repo.

In particular, on the actual infinite root ray S^j(1),

    W(S^j(1))=2/((j+1)(j+2)),      sum_(j>=0) W(S^j(1))=2.

These are polynomial tails in sibling depth, in place of the exponential
r^j. The abstract array also has w(L,0)=2/(L+2), but (L>0,K=0) is not a
completed positive Collatz certificate. That boundary must not be turned
into a claim that such actual rooted corridors exist.

## 5. Polynomial loss against the best fixed price

**P4 — PROVED.** Write t=L+K and

    M(L,K)=sup_(0<r<1) (1-r)^L r^K.

Then, including zero counters,

    [2(L+1)/((t+1)(t+2))] M(L,K) <= w(L,K) <= M(L,K).        (9)

For t>0 the supremum is (L/t)^L(K/t)^K, with 0^0=1. Its product with
binom(t,L) is a binomial probability and is at most 1. Substitute

    w(L,K)=2(L+1)/[(t+1)(t+2) binom(t,L)]

to obtain the lower bound. The upper bound follows because (4) averages
the monomial over a probability density. At t=0 both sides equal 1.

Thus one valid, summable flow loses only a quadratic factor against the
best fixed price chosen separately for each certified source. A badly
matched fixed price can lose exponentially. This is the useful coding
connection: averaging compatible history likelihoods gives a quantified
adaptation guarantee, with no appeal to randomness of an integer's orbit.
It removes a pricing obstruction, not the obligation to obtain its ROOT
certificate.

More generally a beta(alpha,beta) prior with alpha>0 and beta>1 gives

    W_(alpha,beta)(n)=B(K+alpha,L+beta)/B(alpha,beta),
    sum_n W_(alpha,beta)(n)
      <= 3(alpha+beta-1)/(beta-1)-beta/(alpha+beta).          (10)

This follows by writing (3)'s full bound as 3/(1-r)-(1-r) and integrating.
The condition at r=1 is real: a uniform prior gives weight 1/(j+1) on
S^j(1), so even that single fibre has infinite mass. A price mixture must
be checked at its parameter boundaries.

## 6. The exact injection measure lives at multiples of three

**P5 — PROVED.** On V, the defect is

    W-KW = lambda,       lambda(n)=W(n) if 3|n, else 0,
    sum_n lambda(n)=1.                                    (11)

Equality away from multiples of 3 is (7). At multiples of 3 there are
no incoming odd predecessors. Sum over V, using finite mass: the total
defect equals the mass of the sources whose next step is 1. Those are
S^j(1), j>=1, and their mass is 2-W(1)=1.

Consequently lambda is an unconditional probability measure supported on
rooted multiples of 3. W|_V is its occupation measure before hitting 1:

    W|_V=sum_(t>=0) K^t lambda,
    E_lambda[odd steps to 1]=sum_(n in V) W(n)<=13/3.        (12)

To justify the series, iterate (11); its remainder K^T W has norm tending
to zero by dominated convergence, since W has finite mass and is supported
on finite rooted paths. The measure was built from that rooted tree, so
its small mean is not a conclusion about the fixed source priors mu or nu.

The inherited three-divisible inverse section makes this a focused
positivity target. For a unit n>1 modulo 3, choose a from its residue:

| n mod9 | 1 | 2 | 4 | 5 | 7 | 8 |
|---|---|---|---|---|---|---|
| a | 6 | 5 | 4 | 1 | 2 | 3 |

Then z=(2^a n-1)/3 is positive odd, divisible by 3, less than 22n, and
U(z)=n. This section is proved in
[fusion helpers, section 3](collatz_fusion_helpers_20261005.md).
For a certified n of counters (L,K), put j=floor((a-1)/2)<=2. The new
weight transport is exact:

    W(z)/W(n) = (L+2)(K+1)_j / (L+K+3)_(j+1),             (13)

where (x)_j is the rising factorial and (x)_0=1. Its counters are
(L+1,K+j). Positivity of lambda at **every** positive odd multiple of 3
would imply positivity at every unit, and hence universal Collatz.
Equation (11) establishes its normalization, not that missing full support.

## 7. Checked route switching transports these weights

The actual six-letter word v=(1,2,1,1,1,2), on n=155+2048t, ends at
S(h), where h=111+1458t<n. This is the proved guarded controller in
[recursive dependency kernel](collatz_recursive_dependency_kernel_20261004.md).
All six letters add zero sibling depth, and S adds one. Therefore a
certificate for h with counters (L,K) gives

    counters(n)=(L+6,K+1),
    W(n)/W(h)=w(L+6,K+1)/w(L,K).                           (14)

No root search for n is needed once the child's verified certificate is
supplied. The existing substitution proves the literal first-hit word;
the new array evaluates its paid weight. For h=111 the counters are
(23,7), and n=155 has (29,8), giving weight ratio 5800/131461.
At m legal repetitions the counter increment is (6m,m). Reaching an
uncertified exit still leaves an uncertified source.

| Map | Preserved predicate | Lost data | Required sidecar / hostile |
|---|---|---|---|
| ROOT word -> (L,K) | Exact weight polynomial | Order and source guard | 53 and113 both have (1,3), first valuations5 and2 |
| Fixed flows -> positive mixture | Incoming inequalities and rooted support | Chosen price | Parameter-integrability test; uniform prior fails |
| Target -> inverse sibling fibre | Actual common future | Source magnitude changes | Exact base and depth; telescoping tail (7) |
| Child certificate -> guarded controller source | ROOT and computable weight | Ordinary source changes | Guard and literal splice (14) |
| Rooted tree -> lambda | Unit injection mass | Whether an arbitrary atom is present | Positivity at each multiple of3 stays OPEN |

Neither counter can simply be dropped. Ignoring K assigns the same weight
to infinitely many siblings of a base, giving an infinite incoming sum.
Using only the root profile in K gives the same weight to 3 and5, while
13 is an additional positive predecessor of5, violating the incoming row.
Keeping both counters still does not replace the word: 53 and113 have
equal counters and different first operations.

## 8. What is solved and the next admissible targets

The enriched object is (n, checked ROOT word, L, K, moment curve
r^K(1-r)^L). It retains the integer and guard, compresses its price family,
and supports exact common-future composition. Its positive difference
array covers all formal counter pairs; arithmetic realization for every
integer is precisely the unresolved part.

W(n) is a uniformly computable real without assuming convergence. Follow
a finite G-prefix, retaining its accumulated counters. If it has rooted,
(5) is exact. Otherwise zero is a valid lower bound and w(L,K) is a valid
upper bound, at most 2/(L+2). As the prefix length grows this interval
shrinks effectively. Approximating a value that might be zero is not
certifying it positive.

The immediate productive target is a positive lower bound on W(n), or
lambda(3m), that does **not** use an already completed ROOT word. Bounds
in terms of the unknown counters do not accomplish that. Three helper
tests are now precise:

1. Can the retained source/ternary carry bound the counter bill on the
   unresolved exit of a paid controller? The necessary output is a finite
   source-dependent lower bound, not a formal array entry.
2. Can an independently grounded family of three-divisible sources cover
   the section in (13)? Its factor22 gives a fixed ordinary-size cost;
   residue density by itself is insufficient.
3. Can one improve mixtures using a checked price reset at a macro boundary?
   The reset must preserve the full fibre payment, not just a single edge.
   Equations (6)-(7) are exact benchmarks for such a proposal.

No known negative cycle was identified with a positive ROOT. No value from
finite verification is treated as an unbounded convergence theorem.

## 9. Reproduction

[Exact program](../../04-computation/experiments/collatz_adaptive_mixture_flow_20261005.py)
and [JSON output](collatz_adaptive_mixture_flow_20261005.json):

```sh
python3 04-computation/experiments/collatz_adaptive_mixture_flow_20261005.py
python3 -O 04-computation/experiments/collatz_adaptive_mixture_flow_20261005.py
```

**672,721 explicit checks**, identical outputs normally and with -O.
The universe is all positive odds below 2^15; an independent reader follows
actual U-words while the other follows sibling-base edges. The full array
rectangle 0<=L,K<=60 checks splitting, tails and the regret inequality;
0<=L,K<=18 also integrates by expanding the polynomial rather than using
factorials. Inverse target rows below4096 retain six actual siblings and
their exact infinite remainder. The controller tests t=0..127; the
three-divisible section tests unit odd targets below8192. Fixed-price root
generations and beta-prior formulas have separate exact checks.

The finite source head has W-mass approximately2.599557 and lambda-mass
approximately0.572152; JSON retains exact fractions. These are finite lower
bounds, not estimates that the unknown remainder is small. The proved
global quantities are W-mass<=16/3 and lambda-mass=1. Controls include the
divergent uniform-prior root fibre and equal-cost, different-word sources.
All load-bearing calculations use integers and rational fractions; decimals
are display only. This session includes self-audit and separate code paths,
not an independent external audit or a formal proof assistant build.
