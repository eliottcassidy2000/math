# Random fair attachments to a transitive base: exact path-cover transform

Status: PROVED, with exact computational controls and independent source-bound review. This result is recorded as a research packet, without a numbered canon reservation or a claim of literature priority.

## Scope and inherited coordinates

Fix a labelled tournament Q on k>=1 core vertices. Add m>=0 labelled vertices carrying the transitive tournament in their fixed total order. Orient each of the km cross edges independently and fairly. All internal core edges remain fixed. Let E_Q(m) be the expected number of directed Hamiltonian paths, counted as vertex words. No strong-connectivity conditioning or isomorphism-class weighting is imposed.

Let p_r(Q) count ordered spanning covers of Q by r nonempty directed paths. Thus p_r is F_Q(r), NOT the unordered p_Q(r), in THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md, section 1. Put P_Q(t)=sum_{r=1}^k p_r t^r. To avoid conflicting with that source's path-colour Q_T, write R_q(m)=q! S(m,q), with S(0,0)=1 and S(m,0)=0 for m>0. Every expression j^0 uses 0^0=1, as required by the exponential series.

Closest inherited mechanism: ordered path covers / permutation cuts, THM-3181, equations (5),(12). Canonical hostile: K1▷C3 and C3▷K1 have identical profiles but are nonisomorphic (section 8). Corrected near miss: scalar Hamiltonian data does not in general recover a path-cover profile (THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md); quotient-walk sidecars are needed for cyclic substitution. Least-used relevant sidecar here: THM-3181's full backward-adjacency distribution and reversal symmetry.

## The exact transform

For every m>=0,

E_Q(m) = sum_{r=1}^k p_r/4^r [R_{r+1}(m)+4R_r(m)+4R_{r-1}(m)].                 (1)

Proof. In any Hamiltonian word, take its maximal consecutive runs of core vertices. These are an ordered cover by r directed paths. There are r-1 compulsory nonempty base runs between core runs and u optional nonempty base runs at the two ends, with u=0,1,2. Choose the occupied end positions in binom(2,u) ways. The q=r-1+u occupied base slots receive the m distinct base labels by a surjection, counted by q!S(m,q). Within each slot the increasing order is forced, so each assignment specifies exactly one candidate word. The candidate uses 2r-2+u cross edges. They are distinct edges, even when a run is a singleton, because a Hamiltonian word visits every vertex once. Independence makes its success probability 2^{-(2r-2+u)}. Sum over candidates and u by linearity of expectation. No independence between DIFFERENT candidate words is needed. This gives (1). When m=0 only r=1,u=0 survives, yielding E_Q(0)=p_1=H(Q), with no singular exception.

## Exponential generating function, finite spectrum, rational OGF

Let G_Q(x)=sum_{m>=0} E_Q(m)x^m/m!, and y=e^x-1. Then

G_Q(x) = ((y+2)^2/y) P_Q(y/4).                                             (2)

The apparent singularity at y=0 is removable because P_Q(0)=0. Equivalently, define the ordinary polynomial

A_Q(z) = (z+1)^2 sum_{r=1}^k p_r (z-1)^{r-1}/4^r = sum_{j=0}^{k+1} a_j z^j. (3)

Then G_Q(x)=A_Q(e^x), and for EVERY m>=0,

E_Q(m)=sum_{j=0}^{k+1} a_j j^m.                                            (4)

Consequently sum_{m>=0} E_Q(m)w^m = sum_{j=0}^{k+1} a_j/(1-jw). The j=0 term is the constant a_0, not a positive-index tail. A tail recurrence has characteristic polynomial dividing product_{j=1}^{k+1}(X-j); cancellation may reduce it. An all-index annihilator can retain X as well. No minimal recurrence claim is made. Always A_Q(-1)=A'_Q(-1)=0; these are two signed linear identities on the channels. Equation (4) follows alternatively from R_q(m)=sum_{j=0}^q (-1)^{q-j}binom(q,j)j^m.

## Exact tomography and information loss

The transform is injective on ordered path-cover profiles. Its formal inverse is

P_Q(t) = t/(1+2t)^2 G_Q(log(1+4t)).                                        (5)

More strongly, only E_Q(0),...,E_Q(h-1), h=ceil(k/2), are sufficient to reconstruct the entire profile when k is known.

First, E_Q(m) depends only on p_1,...,p_{m+1}, and its coefficient of p_{m+1} is m!/4^m, which is nonzero. Thus the first h expectations recover p_1,...,p_h recursively. For the backward-adjacency distribution B_b, the cut identity is

p_r=sum_{b=0}^{r-1} binom(k-1-b,r-1-b) B_b.                                 (6)

This is triangular with diagonal one, hence the recovered p_1,...,p_h give B_0,...,B_{h-1}. Reversing a permutation gives B_b=B_{k-1-b}, recovering all B_b. Equation (6) then gives all p_r. This is a sufficiency statement, not a claim of optimal sampling minimality on the finite set of realized tournaments.

Two fixed-order cores have equal expectation sequences iff they have equal profiles. This is not tournament classification: K1▷C3 and C3▷K1 have sorted outdegree lists (1,1,1,3) and (0,2,2,2), but both have ordered profile (3,18,36,24). The whole expectation curve collides. Source/target/map: core tournament -> fair-attachment expectation curve, via (1); preserved predicate: equality of ordered path-cover profiles; lost data: SCC order and other profile collisions; restoring sidecar for classification is not supplied here. This also says that any existing fixed-order same-H/different-profile pair is eventually distinguished by at least one attachment experiment, despite agreement at m=0.

## Universal leading channels and first triangle sensitivity

Let c_3(Q) be the number of cyclic triples. The top profile coordinates are

p_k=k!,
p_{k-1}=binom(k,2)(k-1)!,
p_{k-2}=(k-2)![binom(k,3)+2c_3(Q)+3binom(k,4)]  (k>=3).                     (7)

The first formula orders the singletons. The second selects one directed two-vertex path and orders the resulting k-1 blocks. For the third, either use one three-vertex path or two disjoint two-vertex paths. A transitive triple has one directed Hamiltonian word and a cyclic triple has three, giving binom(k,3)+2c_3. On any four chosen vertices there are three pair partitions, each with its uniquely directed pair words. Order the k-2 blocks to finish.

In (3), the coefficient a_{k+1}=k!/4^k. The coefficient a_k is

p_k(3-k)/4^k + p_{k-1}/4^{k-1} = (k+1)k!/4^k.

Hence for fixed k>=2, as m tends to infinity,

E_Q(m)=k!/4^k [(k+1)^m+(k+1)k^m]+O_k((k-1)^m).                             (8)

The error can be taken uniformly over k-vertex cores since the universe is finite. Only p_k,p_{k-1},p_{k-2} contribute to a_{k-1}; the first two are universal. For two k-vertex cores Q,Q', k>=3,

E_Q(m)-E_Q'(m) = (2(k-2)!/4^{k-2})(c_3(Q)-c_3(Q'))(k-1)^m
                 + O_k((k-2)^m).                                         (9)

More generally, if s is the largest index with p_s(Q) != p_s(Q'), their difference has leading exponential channel (s+1)^m with coefficient (p_s(Q)-p_s(Q'))/4^s. Thus the largest differing cover coordinate determines eventual sign. No claim about each realization or a strong-tournament census follows.

For k=3, p(C3)=(3,6,6), p(T3)=(1,6,6), and exactly

E_C3(m)-E_T3(m)=(2^m+2+0^m)/2.                                            (10)

For m>=1 this is 2^{m-1}+1, but at m=0 it is 2, not 3/2. The finite head must not be discarded. E_C3(1)=9/2 is an expectation; individual attachment outcomes here have H=3 or H=5.

## Reproducible finite audit

Run python3 audit.py and python3 -O audit.py in this artifact directory; outputs must byte-match. The code uses only the Python standard library and exact fractions, not sampled means or floating point.

Independent mechanisms: core profiles are obtained from every core permutation and its forward/backward cuts; each attachment mean is computed by exhaustive orientation-mask enumeration and a direct Hamiltonian subset DP on the actual enlarged tournament. The DP does not use the slot formula.

Universe: all 1,099 labelled cores of orders 1 through 5 for top-profile identities, exponential channels and half-prefix tomography. Exhaustive attachment audits cover every core with k<=3 and m=0..3, and every k=4 core with m=0..2: 236 core/base pairs, 22,337 attachment masks. Positive controls: singleton and transitive cores, zero-base head, directed triangle. Hostiles: nonisomorphic profile collision; m=0 in the C3-T3 difference; actual-H span versus a nonintegral mean. No inherited strong-connectivity, irreducibility, or class-census filter is present.

## Search and limits

Targeted local reads: THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md and THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md. A targeted local text search for random attachment, expected Hamiltonian, and Stirling/attachment in canon and result files found no competing exact theorem, only unrelated finite-result hits. This is NOT a novelty search and implies no priority claim. No external literature claim is made. Parent owns startup/synchronization, final canon audit, theorem-ID reservation, Git edits and publication.

## Final corollary: spectral reflection and the forced zero

The cut dictionary also gives

P_Q(t)/t = sum_{b=0}^{k-1} B_b t^b(1+t)^{k-1-b}.

Replacing t by -1-t contributes the sign (-1)^{k-1} and exchanges the exponents b and k-1-b. Because B_b=B_{k-1-b},

P_Q(-1-t)/(-1-t) = (-1)^{k-1} P_Q(t)/t.

Substitute t=(z-1)/4 in (3). Under z -> -2-z, t becomes -1-t and (z+1)^2 remains unchanged. Therefore

A_Q(-2-z)=(-1)^{k+1} A_Q(z).

This is reflection of a polynomial about z=-1, not a claim about a dynamical system or reflection of the positive spectrum about the origin. For odd k, A_Q is even about -1 and has a zero of multiplicity at least two there. For even k, A_Q is odd about -1 and is already divisible by (z+1)^2, so its zero has multiplicity at least three. No assertion of exact multiplicity is required.

## Final hostile: positive expectations do not imply positive spectral atoms

For the transitive four-vertex core, the ordered profile is (1,14,36,24). Direct substitution gives

A_T4(z) = [3(z+1)^5-8(z+1)^3]/32,
(a_0,a_1,a_2,a_3,a_4,a_5)=(-5,-9,6,22,15,3)/32.

The negative a_1=-9/32 remains in every positive-index tail; this is not merely the removable zero-base finite-head channel a_0.

Let H be the 5-by-5 shifted Hankel matrix H_{ij}=E_T4(i+j+1), 0<=i,j<=4. Since the zero channel vanishes at positive indices, H=V diag(a_j j) V^T, where V_{ij}=j^i, 1<=j<=5. Thus

det(H) = product_{j=1}^5(a_j j) product_{1<=i<j<=5}(j-i)^2
       = -16238475/1024 < 0.

Exact rational Gaussian elimination independently obtains the same determinant. In particular this positive expectation sequence has no nonnegative-measure moment representation of the form E(m)=integral x^m dmu(x) for all m>=1 with mu supported on nonnegative x: that representation would make H positive semidefinite, contradicting its negative determinant. Here positivity of the scalar observable does not imply positivity of its spectral expansion.

The legitimate cross-topic connection to signed runner-curvature calculations is this diagnostic boundary: a positive observable may be assembled from signed spectral coordinates, and dropping the signs can invalidate curvature or moment arguments. This is not an identification of the underlying runner and tournament measures or an automatic theorem transfer.

Final finite audit adds both the polynomial reflection and forced root multiplicities for all 1,099 cores, the exact T4 coefficient vector, and the independently evaluated shifted Hankel determinant. Normal/-O replay outputs again byte-match. Frozen artifact SHA256: audit.py dd95ac36eba1ea12e8c786516585b0367ba56a4eab87b411a02b50361efcf047; audit.out e0c7a19fc451457f21dd1c7774045e979f5b4cc79d750ae29ceb7083972f03aa.
