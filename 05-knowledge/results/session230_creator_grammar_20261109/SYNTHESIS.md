# First birth is an optimization; a generator is operation-relative

**Status:** PROVED elementary reductions and independently checked finite
certificates, with cited inherited classification clearly separated. The
five-operation grammar through order eleven remains a finite result whose
large census completeness is inherited. Chan–Pak completeness remains OPEN.

## 1. A constructive answer to the question

Keep the **whole birth fiber** for every attained count h: all isomorphism
classes of the least order b(h) realizing h. Session224's H=9 example already
shows that one class or one converse pair is insufficient. No definition of
b(h) for all odd h is assumed; 7 and 21 are excluded and global realization
of every other positive odd number is not established here.

There are three different meanings of a fundamental tournament:

- For **ordinal sum**, the indecomposable objects are strongly connected
  tournaments. The SCC word is unique and H is its multiplicative norm.
- For **substitution**, the indecomposable objects are prime tournaments.
  H alone is not the substitution invariant; retain module/path-cover data.
- For **first appearance**, an object minimizes vertex cost at a fixed H.
  This is an optimization property, not either preceding irreducibility.

These meanings are provably different. The H=5 birth object is already
C3[TT2,1,1], so the full birth library is not an irredundant basis. Conversely,
the prime W7 has H=31, already realized at order six, so even **all** first
creators cannot generate it by substitution and ordinal sum alone.

The independently certified first failure occurs at order seven. There are
exactly eight noncreator primes, with H-multiset
27,31,31,41,41,41,43,43; at order eight there are 214. The canonical-orbit
mass audit exhausts all labelled tournaments, not just a sampled class list.
See [proof and constructions](primes/REPORT.md) and
[independent audit](primes/audit/REPORT.md).

## 2. The 63 exception is a change of operation, not a new prime

The two first creators of 63 have no first-creator one-vertex deletion.
Nevertheless each has three disjoint two-vertex modules. Contracting them
produces W5, the prime first creator of 9. With W5's chain-first/hub-last
convention, their module sizes are

    (1,2,2,2,1) and (2,2,2,1,1).

Every nontrivial child is TT2, itself two singletons in ordinal sum. Thus the
famous creator-ancestry orphan is built directly from smaller birth objects
once substitution is permitted. Its deletion-count multiset is
{1,11,25,25,33,33,37,37}; this does not contain a count new at order seven.

This generalizes beyond 63. The previously known Fibonacci fans F_k with
H(F_k)=2^k-1 have exactly two fixed prime kernels, obtained by compressing
constant-orientation runs of their transitive chain. For every r>=2,

    F_(2r)   = W5[TT_(r-1),TT2,TT_(r-1),1,1],
    F_(2r+1) = W7[TT_(r-1),1,1,1,TT_(r-1),1,1].

The Mersenne realization and fan definition are inherited from THM-4748;
the present result identifies their module kernels. The even-exponent
branch uses a first creator, while the odd-exponent branch uses the
noncreator prime W7. This is an exact illustration of why an arithmetic
family can alternate between two structural primitives.

## 3. A nontrivial paired growth operation

Define W_(2k+1) as a transitive chain of length 2k with an extra hub beating
exactly the even-indexed chain vertices (indices start at zero). Its primality and criticality have
elementary proofs: an interval of length at least two is distinguished by
the alternating hub; for k>=2 every vertex deletion creates a proper nontrivial module. W3 is the seed exception.

The exact all-order count and genuine two-vertex lift are

    H(W_(2k+1)) = (4^(k+1)+6k+5)/9,
    H(W_(n+2))  = 4 H(W_n)-n.

The first values are 3,9,31,117,459,1825. This is not unconstrained vertex
insertion: every result remains prime, while every intermediate one-vertex
deletion is nonprime. No finite seed set generates all tournaments by
substitution alone, since these prime orders are unbounded. No claim is
made about the first-birth status of every later W_n.

The cited Schmerl–Trotter prime-deletion mechanism provides a structural
all-order generation route with prime-preserving one- or two-vertex growth
and proper substitution. Its exact scope and relationship to creator
stopping leaves are recorded in the prime report. This route is distinct
from proving the much tighter finite five-operation grammar at every order.

## 4. The 7/21 proof object literally computes the growth operation

Let H(empty)=1. For distinct original vertices a,b define K_T(a,b) as the
sum, over directed a-to-b paths with an even number of vertices, of
H(T minus the path vertices). Then a new vertex v beating precisely A has

    H(T+v) = H(T) + 2 sum_(a in A,b not in A) K_T(a,b).

This positive directed-cut kernel comes directly from the odd-cycle packing
formula: the unique odd cycle containing v is cut open into the even path.
The other disjoint cycles contribute the residual H. Thus the polymers in
the omission proofs are not only analogous to generation coefficients;
they are those coefficients. Existing directed-cut/submodular calculus in
THM-4097/4114/4118 is credited, not rediscovered without attribution.

Consequences include:

- Neutral cuts are exactly terminal SCC unions.
- Every pair curvature is at most -2H(T-{a,b})<0.
- On a strong parent of order at least three, the cheapest nonneutral
  extension is H(T)+2 mincut(K_T).
- If this exceeds max S(n), every nonneutral extension is a new first
  creator. This is a sufficient certificate, not a necessary condition and
  not an efficient algorithm for first constructing K.
- Equal order and H do not determine the extension spectrum. Explicit
  five-vertex H=3 examples have different possible child counts.

See [positive-kernel proof](cut/PROOF.md) and
[independent audit](cut/audit/AUDIT.md).

For a transitive parent, K(i,j)=0 backwards, 1 at distance one, and
2^(j-i-2) at larger forward distance. A cut word can discard its initial
zeros and final ones. To have cut weight 3 or 10, the residual word must
start 1, end 0, and have length at most six: a longer word's endpoint term
alone is at least16. All31 remaining words avoid both weights.

Therefore H=7 and H=21 are excluded for **every apex-transitive tournament**
by a finite certificate with a proved all-order reduction. This is a direct,
noncircular subclass proof, not a replacement for the known global omission
proofs. The common language is feedback-vertex transversals, odd-cycle
packings and realizable integer cut capacities, rather than numerology.

## 5. Recovered connectivity gives a sharp structural boundary

The cut argument is inherited from THM-4811(B), recovered during the concurrent S17 audit. We restate it with explicit sharp-boundary witnesses, not a novelty claim. For every integer k>=1 and tournament of order 2<=n<=4k+1,

    lambda(T)>=k  iff  minimum indegree and outdegree are at least k.

The proof sums degrees on a cut. A cut with fewer than k reversed arcs
would need at least2k+1 vertices on **both** sides. At n=4k+2 this boundary
is sharp: two regular (2k+1)-blocks joined in one direction with one cross
arc reversed have minimum semidegree k but lambda=1, for every k>=2.
For k=1, the nonstrong ordinal sum of two C3 blocks is the boundary witness.

This recovers S17's n<=13 test at k=3 and identifies its sharp failure at14.
Separately, at every n, if the connectivity-k H-floor exceeds max S(n-1),
all noncreators are within k-1 cross-arc reversals of a proper ordinal sum.
That implication has no degree-threshold restriction; the threshold is only
needed to replace connectivity by a score test.

## 6. What remains open

The finite grammar through eleven uses every birth class, restricted
extensions of birth parents, proper substitution, ordinal sum, and a
coupled repair of at most two cross arcs. Our audit confirms its logical
noncircularity and all97 supplied exceptional witnesses, but does not rerun
the full large extension census. No all-order validity follows from it.

The next strong target is an all-order characterization of prime noncreator
residuals, or a bounded-defect operation that reaches them. During this session, the incoming S17 audit established class-level creator genealogy at order eleven: zero orphan classes among 506,912,096 candidate extensions from noncreator parents. That is an inherited independently audited finite result, not a fresh Session230 replay. The only creator orphans through eleven remain the two63 classes. Its source and limits are in the [updated audit status](audit/REVIEW.md#current-status-update-after-the-baseline-review). Count realization, class generation, and basis minimality
remain three separate problems.
