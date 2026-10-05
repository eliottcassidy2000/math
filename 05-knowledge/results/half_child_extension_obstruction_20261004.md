# The exact half-child decoder has an all-length obstruction

2026-10-04. **PROVED:** the rank-charge criterion, the finite supplied-child
decision bound, the all-length obstruction at7, and the repaired families.
**FINITE-EXACT:** the explicitly bounded certificate census below.
Universal positive Collatz convergence and universal grounded selection remain
**OPEN**. This is an obstruction to a specified representation of a join,
not an obstruction to reaching the root.

Artifacts: [script](../../04-computation/experiments/half_child_extension_obstruction_20261004.py)
and [saved output](half_child_extension_obstruction_20261004.out).

## 1. Inheritance and the precise extension being tested

The closest proved mechanism is the ordered-carry decoder in
[the eight-row debt search](debt_word_join_search_20261004.md), extended to
sixteen selected families by
[the incoming binary/ternary fusion](collatz_binary_ternary_guard_fusion_20261004.md).
Those identities preserve the child map `(n-1)/2`; comparing affine slopes
forces equal numbers of odd steps and a halving-cost difference of one.
The incoming [general family lift, section2](collatz_join_shields_and_lifts_20261004.md)
already allows different word lengths and changes the child parameter map.
That broader mechanism supplies the repair below; it is inherited.

Canonical hostile: source7 and its half-child3. Corrected near miss: the
failure of one short word shape did not exclude longer joins, but a genuine
all-length invariant can. Least-used sidecar: the first-hit odd rank and
halving cost separately, including what happens if root loops are appended.

Anchor: determine whether the same half-child decoder can become universal.
Niche: replace an unbounded word search by a bound supplied by the existing
child certificate. Wildcard: classify the attainable affine slopes after
root padding. The concept board is **source / fixed child / ordered word /
odd rank / halving cost / root boundary**. No literature-priority claim or
new general convergence theorem is made. The targeted inherited statements
do not contain the rank-charge classification below.

Let `U(n)=oddpart(3n+1)` on positive odd integers. For a word w of positive
valuations, write

    length(w)=r, cost(w)=A, F_w(X)=(3^r X+B_w)/2^A.

The identity under examination is an equality of affine functions:

    F_w(X) = F_v((X-1)/2).                             (1)

Its actual seed must satisfy `n>=3, n=3 mod4`, so `m=(n-1)/2` is positive
odd, and both words must be the actual valuation prefixes of those sources.
Pointwise equality only at X=n is a weaker condition and is not (1).
Comparing slopes in (1), and using unique factorization, gives

    length(w)=length(v),  cost(w)=cost(v)+1.            (2)

When (2) holds, a common endpoint at the one seed already forces the
constant terms in (1) to agree. This supplies both directions of the test.

## 2. A complete criterion for supplied home certificates

For a source x with a supplied first-hit route to1, let `r(x)` count odd
steps and let `A(x)` count halvings. Define the integer

    I(x)=A(x)-2r(x),       I(1)=0.                     (3)

Equivalently, I is ordinary first-hit rank minus three times odd first-hit
rank. The formal root self-edge has valuation2, adds one odd step and two
halvings, and leaves I unchanged. Formal padding is arithmetic bookkeeping;
a stored first-hit home certificate stops before that self-edge.

**PROVED strict criterion.** Suppose both n and m have supplied home
certificates. There is a join of the form (1) using actual prefixes that
do not continue past either first root visit if and only if

    r(n)=r(m),  A(n)-A(m)=1.                           (4)

Necessity: by (2) the two prefixes have equal lengths and a cost difference
of one. Their common endpoint has the same remaining first-hit suffix,
so appending it preserves both differences and gives (4). Sufficiency:
use the two complete first-hit words. Their endpoints are1; (4) and the
seed equality imply (1). Empty prefixes are allowed, but do not furnish
a half-child identity on distinct positive sources.

**PROVED padded criterion.** If actual `U(1)=1` loops of valuation2 may be
appended to the two routes, there is a join of the form (1) at some finite
length if and only if

    I(n)-I(m)=1.                                      (5)

For necessity, take the equal-length join from (2). If its endpoint is
not1, append the same actual suffix to root on both sides. If it is1,
no suffix is needed. At the resulting equal length L, the two costs are

    A(n)+2(L-r(n)),  A(m)+2(L-r(m)).

Their difference equals both one and `I(n)-I(m)`. For sufficiency, pad the
two supplied home routes to length `L=max(r(n),r(m))`. Equation (5) gives
the required cost difference, and both endpoints are1. Therefore the
slopes and the constant terms match. This proof covers every possible
length; it does not extrapolate from a finite padding search.

A useful local form does not require a home certificate. If n and its
half-child actually coalesce after the same number L of odd steps and
their accumulated halving-cost difference is d, that difference remains
d at every later equal-time endpoint. If d is not one, no earlier or
later half-child affine join can exist: an earlier join with difference
one would retain difference one through time L. This is a finite witness
of an all-length obstruction even when the shared future is not known.

## 3. Source7 is blocked at every decoder length

The exact first-hit routes are

    7 --1--> 11 --1--> 17 --2--> 13 --3--> 5 --4--> 1,
    3 --1--> 5 --4--> 1.

They have `(r,A)=(5,11)` and `(2,5)`. Both charges are1, so (5) fails.
At equal time five, the child has added three root loops, making its cost
`5+3*2=11`, exactly the source cost. Every further simultaneous root loop
leaves that difference zero. No longer search through positive valuation
words can produce (1) for source7.

This is the least obstruction in the stated source domain. At source3,
the half-child is1; the source word `(1,4)` and padded child word `(2,2)`
do satisfy (1), with cost difference one. They are not two first-hit
prefixes, because the child starts at root. Source3 can of course be
handled by its direct first-hit route. Root padding must not create a
false certificate rank.

The incoming hostile sources27 and703 also fail (5), with charge
differences -15 and -13 respectively. These are checked properties of
their supplied home routes, not evidence of nonconvergence.

## 4. The existing child certificate gives a finite search bound

**PROVED source-aware decision.** Suppose only the half-child's first-hit
certificate is supplied, with odd rank R and cost D. Every strict join
of type (1) must occur within the first R source steps. Indeed, joining
at time j<=R leaves a shared suffix of R-j steps, so n itself must first
reach1 at exactly time R with cost D+1.

Consequently, replay n for at most R odd steps, stopping immediately at1.
Accept the strict affine ansatz exactly when that first root occurs at
time R and its accumulated cost is D+1. The converse is (4), proved
directly by the two complete words. If the source has not reached1 by R,
the strict ansatz is ruled out, without a conjecture about its later
trajectory. A root reached earlier also rules out this strict ansatz.

The executable `decide_strict_with_supplied_child` implements this bound.
It validates the child source and first-hit word, retains the original n,
observed exact prefix and actual frontier, and returns the scoped verdict.
It does not manufacture a child certificate. A rejected half-child rule
does not reject another child, a different child map, or a join with
unequal word lengths. It also does not by itself decide the padded ansatz
if the source remains unfinished at the bound.

For7, the supplied child3 has R=2. Two actual source edges leave frontier17;
the strict half-child decoder is already exhausted. The explicit five-step
route in section3 supplies the stronger obstruction even to root padding.
The charge I of an unfinished arbitrary source is not assumed computable:
the all-length classification consumes home certificates when it uses I.

## 5. What changing the child map recovers

The seed7 really does share a future with3:

    7 --(1,1,2,3)--> 5 <--(1)-- 3.

The word lengths are four and one, so this pointwise join is not (1).
The inherited least-period lift turns it into the exact family

    7+256t --(1,1,2,3)--> 5+162t <--(1)-- 3+108t,
    t an integer >=0.                                (6)

The source period256 preserves the four-letter exact cylinder; the child
period108 preserves its one-letter exact cylinder. Direct substitution
gives the shared endpoint. At every height the child is positive and
strictly smaller. Its new affine map is

    m=(27n+3)/64,

whose slope27/64 replaces1/2. A supplied certificate for this actual
child can be transported through the shared endpoint. The family does
not independently supply every child certificate, and does not prove
that all sources enter its cylinder.

Changing the intercept instead is another valid repair, distinct from (1):

    7+4096t --(1,1,2,3,4)--> 1+486t <--(2,2,2,2,2)-- 1+2048t.

Here the child map is `(n-5)/2`. At t=0 the child loops are formal padding
and must be trimmed from a first-hit certificate; the original source
still has its genuine five-step route. Both repairs retain the original
seed, and explicitly record which child coordinate has changed.

There is also an exact classification of this freedom. For any two
supplied home sources n,m, the positive slopes of affine child maps
`Y=lambda*X+beta` admitting an actual common-future diagram, with root
padding allowed, are precisely

    lambda = 2^(I(m)-I(n))*(3/4)^k,   k in Z,          (7)
    beta = m-lambda*n.

To prove necessity, append the same suffix from the common endpoint to
root. Let the resulting lengths be L and M. Their costs are `2L+I(n)`
and `2M+I(m)`. Slope comparison gives (7) with k=L-M, and the seed fixes
beta. Conversely every integer k is obtainable by adding sufficiently
many nonnegative root loops to the two existing certificates; the same
calculation gives the slope, and equality at the seed fixes the constant.
This classifies representation choices for already certified sources.
It is not a construction of a home certificate for an arbitrary input.

The half-child slope1/2 requires k=0 and `I(n)-I(m)=1` by unique
factorization. For7 and3, the attainable slopes are `(3/4)^k` instead;
the successful slope27/64 in (6) is k=3. Thus the obstruction and the
repair are consequences of the same exact invariant.

## 6. Exact finite controls and limits

The declared universe is all5000 sources `n=3 mod4`, `3<=n<20000`, with
their actual half-children. Both home routes are explicitly obtained under
a5000-step cap and independently read using ordinary halving loops.
This bounded run succeeds for all declared inputs; it makes no claim
about arbitrary larger inputs.

The results are:

| Exact property in this finite universe | Count |
|---|---:|
| Strict affine half-child join possible | 2870 |
| Join possible with formal root padding | 2871 |
| Impossible at every length, including root padding | 2129 |
| Padding-only instance | 1, namely source3 |

The script checks195199 literal synchronous prefixes against independently
composed affine coefficients,5000 bounded decisions supplied with the
actual child certificate,909 general-slope instances, and260 repaired
family instances. The latter use t from0 through127 and the additional
exact heights10^20 and10^100. Six malformed sources, a word belonging to
the wrong source, and a root-padded alleged first-hit certificate are
rejected. Checks do not disappear under optimization.

Neither these counts nor the single padding-only instance is promoted to
a density theorem or an all-source classification without its stated
certificate hypotheses. In particular, no simplification from equal odd
ranks alone to equal halving costs is assumed.

Reproduce from the repository root:

    python -X utf8 -B 04-computation/experiments/half_child_extension_obstruction_20261004.py
    python -O -X utf8 -B 04-computation/experiments/half_child_extension_obstruction_20261004.py

Both outputs agree after LF normalization; SHA256:
`4b517302d32a49582569ef41e3a7a83d728ba4ccddfb5ddf66798c6759d5957b`.

Connection contract: source = two exact supplied routes and their seed
identities; target = an affine common-future identity; map = composition
of guarded words; preserved = actual endpoint and source identity; lost
by keeping only a pointwise join = slope, carry and two separate odd
lengths; required sidecars = first-hit boundary, child parameter map and
existing rooted suffix. The cheapest decisive hostile is7/3. The strongest
survivor is the exact unequal-length family (6), together with the finite
decision bound that prevents indefinite search inside the failed ansatz.
