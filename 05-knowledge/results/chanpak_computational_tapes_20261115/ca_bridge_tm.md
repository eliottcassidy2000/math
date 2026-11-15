# A tournament tape has an exact observation horizon

**PROVED:** a sharp continuation coordinate theorem and an unbounded linear register obstruction for indivisible one-port modules. **Not proved:** all-target synthesis, or that bounded-state nonlinear encodings are impossible.

## Historical inputs

- [Session 135 synthesis](https://github.com/eliottcassidy2000/math/blob/7f6987db0cc31f79d603cfa58488a23d3e5c4b75/05-knowledge/results/session135_fractran_tm_structure_20261101/SYNTHESIS.md) and [literal machine translation](https://github.com/eliottcassidy2000/math/blob/7f6987db0cc31f79d603cfa58488a23d3e5c4b75/05-knowledge/results/session135_fractran_tm_structure_20261101/machine_translation/PROOF.md): packed one-tape and literal multi-tape implementations preserve an exact **boundary decode** even though intermediate rows can break a macro-level invariant. Failed FRACTRAN tests restore the entire partial write before the next branch; one-hot control can be erased on a certified image; the named-vertex square-sum verifier rejects a repeated-vertex word that passes several aggregate counters. Translate that discipline here: do not nominate a register by what the present H observes; characterize what all future readers can still observe.
- [Observable Hankel clocks](https://github.com/eliottcassidy2000/math/blob/951dfeafe94865ebd360e40e97136a863d07bc86/05-knowledge/results/observable_hankel_clocks_20261020/README.md) separates a reduced recurrence, its Hankel image, and a precision-dependent filtration of observable states. The analogue here is a reader-budget-dependent filtration, not the same recurrence.
- [June tournament automata](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/07-reflections/tournament-three-state-automaton-left-middle-right-s582.md) distinguishes a hidden middle state and its final binary shadow in an LRC application. Here the hidden coordinates are directional ordered path covers and the shadow is H. [Session 197's certified binary–ternary bridge](https://github.com/eliottcassidy2000/math/blob/b49eb615b9238151032033db17bc0265c104a89b/session197/binary-ternary-bridge/PROOF.md) supplies a contrasting positive case: a fixed small protected interface works for a particular family of controlled moves. The theorem below concerns **all** homogeneous-module continuations in a budget, not that restricted family.

- Direct technical ancestors: [THM-1975 on the refined compositional invariant](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md), [THM-3121 on the factorial walk-content kernel](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md), and [THM-3134 on the endpoint jet](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md) already give the complete profile and general substitution expansion. The contribution here is not that composition formula: it is the sharp external-contact costs, the physical unitriangular family in the return-gap-weighted coordinates, the modular equivalence in every characteristic, and the unrestricted reader rank obstruction.

## Machine model

Write P_i(D) for the number of partitions of all vertices of nonempty D into i ordered nonempty lists, each list a directed path; the order of the lists matters. P_1=H. A reader is any tournament R with marked vertex u, into which we substitute D: every exterior vertex v sees all of D uniformly in the direction prescribed by its arc with u. Denote the output by R[u<-D]. A sequence of such one-active-module substitutions flattens to one reader because every added vertex still sees the original D uniformly. Count the total number of exterior vertices rather than the number of macroinstructions. No assertion here covers a write or contact that distinguishes individual vertices of D.

**Theorem.** For any two sources D,E and any budget q>=0, the following are equivalent:

1. Every reader with at most q exterior vertices gives the same H for D and E.
2. P_i(D)=P_i(E) for 1<=i<=floor(q/2)+1.

The same equivalence holds for strong sources and strong readers/outputs only. A one-vertex reader is admitted and reads H. For general sources the separating readers with positive exterior budget still produce strong outputs, even if the input D or E is not strong.

**Upper direction.** Compress the maximal contiguous D intervals in a Hamiltonian path to a repeated marker. Between two such intervals there are at least two exterior vertices: an intervening singleton would have to receive from D and send back to D, which uniform adjacency forbids. Thus i intervals cost at least 2(i-1) exterior vertices. Given the compressed path (ordered exterior vertices and i plug markers), the above conditions completely determine its admissibility; each valid compressed path admits exactly P_i(D) substitutions. Summing partitions genuine paths rather than walks, with no omission or duplication, yields H as a linear combination of the claimed coordinates, with nonnegative integer coefficients dependent only on R.

**Lower direction and physical query design.** Let R_m[D]=C3[D,TT_m,TT_m] for m>=1, with block directions D->A->B->D and two transitive exterior blocks of size m. Take R_0[D]=D. To visit D m+1 times one must spend every exterior vertex in m return gaps, each necessarily of form a_j,b_j, with no exterior prefix or suffix. The relative order of the A vertices across those gaps is unrestricted, and the same is true for B; there are (m!)² literal assignments. (Transitivity of A,B does not order their separated visits.) Hence

    H(R_m[D]) = (m!)² P_{m+1}(D) + sum_{i=1}^m c_{m,i} P_i(D).

The coefficient therefore never vanishes, and m=0,1,...,floor(q/2) is a physical triangular interrogation basis. In particular, if the first different coordinate is P_j, then the smallest external budget at which the sources are distinguishable is exactly 2(j-1); one explicit separator is R_{j-1}. The diagonal (m!)² is not invertible modulo every integer. The sparse-reader result below gives the correct integral coordinates for modular work.

## Arbitrarily large residual rank, even for strong tapes

Use an infinite matrix with one row per source D and one column per literal R_m, entry H(R_m[D]). Over Q its column rank is infinite. For transitive sources TT_s, P_i(TT_s) counts surjections of an s-set onto i ordered boxes. The submatrix (s=1,...,d, i=1,...,d) is triangular with diagonal i!; changing columns to m=0,...,d-1 is the invertible triangular change above. Thus any linear state x(D) in a fixed-dimensional Q-vector space with every such future readout linear in x is impossible; allowing affine readout adds only one dimension.

Unbounded rank persists for the strong sources S_s=C3[TT_s,1,1], s>=1, and already for Q_j(S_s), where Q_j counts ordered j slots of paths with empties allowed. The change of basis Q_j=sum_{i=1}^j binom(j,i) P_i is unit triangular. The two exterior singletons allow at most j+1 maximal TT_s intervals in j slots, because each interior re-entry gap consumes both singletons. There are exactly j compressed skeletons at j+1 intervals: both singletons in one slot making the return, and one TT_s interval in each of the other j-1 slots. Consequently Q_j(S_s)=j P_{j+1}(TT_s)+lower terms in P_i(TT_s), i<=j. Since P_i(TT_s)~i^s as s tends to infinity, Q_j(S_s)~j(j+1)^s; any finite linear combination with top index j has a unique dominant exponential, and cannot vanish on every s. Strong-source and strong-output restriction therefore does not make a fixed-dimensional linear representation possible.

The obstruction is not an information-theoretic prohibition on integer encodings: a single integer can encode an arbitrary polynomial nonlinearly. It applies to finite-rank linear or affine state readout. A proposed bounded initial-segment test for arbitrary-depth interchangeability would still need a separate theorem that its segment determines the entire profile on the physical source family; this rank proof does not assert collisions at every fixed segment depth. A bounded-active-interface architecture can evade it by proving a uniform visit cap, selecting a restricted reader language of finite rank, or using a genuine growing tape. Computing H exactly on a finite set, or producing a Turing-machine translation of a count formula, does not in either case establish that a prescribed target physically realizes.

## How this can change a target search

A continuation budget q has a precise cost of visibility: q=0,1 sees only H, q=2,3 sees (H,P2), q=4,5 sees (H,P2,P3). A search over these homogeneous contexts can quotient by exactly those coordinates without ever coalescing two distinguishable states. To compare two signed-repair candidates, identify the earliest different P_j: no reader spending fewer than 2(j-1) exterior vertices can exploit the difference, while a canonical strong reader at that budget gives the first conclusive test. This is a constructive alternative to counting all possible marked hosts, but not an inverse formula choosing one at every N. Potential next move: design a restricted recontact language with a **constant uniform visit cap** and independent cycle updates onto the surviving coordinate differences, or deliberately use unbounded grades with a fast online compiler. Between those two extremes, an unfounded finite register claim is now falsifiable by literal readers.


## The integral observation state and the finite modular tape

**Return-gap action.** A compressed frame with i occurrences of the module has i−1 interior exterior gaps; distinguish its exterior prefix and suffix. Every gap starts at a vertex receiving from the module and ends at one sending into it. Arbitrarily permuting the interior gaps preserves legality, since no exterior arc connects two distinct gaps directly. Their vertex sets are disjoint and nonempty, so this action of S_(i−1) is free. Thus the coefficient of P_i in *every* reader is divisible by (i−1)!. Define the physical weighted grades W_i=(i−1)! P_i.

At the exact frontier q=2m, the coefficient of P_(m+1) is zero unless the exterior splits into m outbound vertices A and m inbound vertices B. If so, it is

    m! per(M),   M_(a,b) = [a -> b],

where per is the permanent: a maximal-run frame determines a perfect matching of the m two-vertex return gaps and any of its m! orders is allowed. The matching-only orientation (a_i->b_i, b_j->a_i for j≠i) realizes the smallest possible nonzero top coefficient m!, not (m!)². Keep the arcs inside A,B transitive, for definiteness. This reader is still strong after substitution, as every a_i and b_i lies on a directed triangle with the module and every module vertex lies on such a triangle.

**Exact ring version.** The matching-only reader with q=2m therefore gives W_(m+1) plus an integer linear combination of lower weighted grades. The reader basis in W coordinates is lower unitriangular and can be inverted over any commutative ring. In particular, for every modulus M≥2 and budget q, two sources have the same H modulo M in every admissible reader if and only if all W_i agree modulo M through i=floor(q/2)+1. The same equivalence holds with no finite bound on q and equality of all modular W_i. Integer linear combinations of the physical queries reconstruct W_i; this does not mean a single tournament realizes each signed combination.

An unordered i-path cover has i! orderings, so W_i is always a multiple of (i−1)!i!. Consequently all modular grades at or beyond the first i with M dividing (i−1)!i! vanish. For a prime p, the observation coordinates equivalently are the raw P_i modulo p through

    min(floor(q/2)+1, p−1),

since the lower factorials are units. This gives a finite modular tape and an infinite-rank exact tape **for the same physical reader model**. Over F_p the finite-budget matrix rank is exactly that minimum for unrestricted sources, by the same transitive triangular minor (whose weighted diagonal is (i−1)!i!) and the reader unitriangular change. In particular mod 2 the entire machine only sees H and mod 3 it sees H,P2. For M=p^k the sufficient cutoff is obtained from Legendre's formula:

    ν_p((i−1)!i!) = (2i−1−s_p(i)−s_p(i−1))/(p−1),

where s_p is base-p digit sum. Thus the first forced-zero grade is (p−1)k/2+O_p(log(k+2)); exact classification uses W_i whether or not additional universal relations remove other states. This is a precise mechanism for congruence-complete observations disguising an infinite exact continuation profile, not a claim that all residues or all modular states are physically reachable.

## A literal monic tape head, not only a readout

The sparse matching reader M_m with 2m exterior vertices gives a **single strong tournament operation** with all-grade formula

    W_j(M_m[D]) = W_(j+m)(D) + integral linear combination of W_i(D), i<j+m,

for every j≥1 and every nonempty D. More generally, an arbitrary reader R of budget q updates grade j using only weighted grades i≤j+floor(q/2), with integer coefficients.

For the general claim, recover W_j of the output via the integral unitriangular readout with queries of cost at most 2(j−1). Compose each query with R: the resulting homogeneous context has cost at most q+2(j−1), so the previous theorem gives the claimed cut and integrality. Alternatively, decompose into ordered slots; the number of within-slot return gaps needed to exceed j module intervals gives the same upper bound.

For the leading coefficient in the 2m case, count ordered covers of the output by j nonempty paths. With j+m module intervals, every slot must contain a module interval; all exterior vertices must lie in m internal return gaps, each consisting of one A→B pair. There are per(M) matchings, m! linear orders of their pairs, and binom(m+j−1,j−1) ways to distribute the ordered gaps into the j slots. The intervals themselves can be filled in P_(j+m)(D) ways. Multiplication by (j−1)! turns the coefficient into (m+j−1)! per(M), so it is per(M) on the W scale and exactly 1 for the sparse matching reader.

In particular the two-vertex move L(D)=C3[D,1,1] satisfies a simple band formula, with W_i=0 for i≤0:

    W_j(L(D)) = W_(j+1)(D)
              + j(j+2) W_j(D)
              + j(j−1)(2j−1) W_(j−1)(D)
              + j(j−1)^2(j−2) W_(j−2)(D).

For verification before weighting, count compressed frames with k module intervals. For k=j+1 there are j possible locations of the two-vertex return. For k=j, either the two external vertices go independently to the eligible ends of slots (j²), or both form one exterior prefix or suffix (2j). For k=j−1 there is one exterior-only slot: put both external vertices there (j), or put one there and attach the other to one of the j−1 module slots (2j(j−1)). For k=j−2 the singletons make distinct slots (j(j−1)). These cover all cases, with impossible grades at k≤0 contributing zero because D is nonempty. Taking the factorial ratios gives the displayed band.

This makes the machine analogy constructive: two more vertices physically bring the next hidden grade into the observed register with a **unit leading coefficient**. Iterating the operation gives an actual coherent tape evolution. A prescribed H can be solved for a load in the *formal* integer linear recurrence, but the load must be the W-profile of an actual tournament. The exact reachability problem has not become a register-machine universality problem. At fixed modulus M the truncated tape does update without seeing hidden grades, because every omitted W_i vanishes at that resolution; formal solvability at every M still does not produce a physical integer preimage.

## A one-letter strong clock is already superexponential

Start with one center and iterate L(D)=C3[D,1,1] exactly k times, making the literal strong tournament T_k on 2k+1 vertices for k≥1. Index the two vertices of shell i by a_i,b_i, inner to outer. The a_i are transitive in increasing index and the b_i in decreasing index. Across types, a_i→b_i on a matching and b_i→a_j for i≠j. The center beats every a_i and loses to every b_i. Then, with natural logarithms,

    log H(T_k) = 2k log k − 2k log log k + O(k).

**Lower bound.** Let r=floor(k/log k), for sufficiently large k. Choose the largest r indices as r bridges, one fixed bridge a_i→b_i per slot. Reserve the smallest r indices: put their a-vertices all in the first slot, and make their b-vertices one fixed terminal tail per slot. Independently assign each of the remaining k−2r a-vertices, and each of the remaining k−2r b-vertices, to an arbitrary slot. In every slot, order its a-vertices in increasing index (ending at a_i) and its b-vertices in decreasing index (starting at b_i and ending at its reserved tail). Concatenate these words in fixed slot order and prepend the center. A reserved tail beats the first a-vertex of the next slot because no a-vertex of a reserved index is in a slot after the first. Every distinct assignment yields a distinct Hamiltonian word, because the fixed bridges recover each slot. This gives at least r^(2(k−2r)) paths.

**Upper bound.** Remove the center as a separator, keeping the two fragments before and after it. If there are r A→B transitions inside the fragments, each is a distinct matched bridge, and each type has at most r+2 transitive runs across both fragments. The two type patterns have only O((k+1)^2) possibilities (each alternates runs and is determined by starting type and run count, permitting either fragment empty). Order the bridge indices in at most k^r ways. Assign each of the other k−r vertices of each type to its run, at most (r+2)^(2k−2r) choices; a transitive run can have only one internal order. Summing over r adds only a polynomial factor to the maximum. Put B=log k−log log k. If r≤k/log k, the logarithm of this bound is at most 2kB+O(k), since r log k≤k and log(r+2)≤B+O(1). For r≥k/log k put a=(r+2)log k/k≥1. The logarithm is

    2kB + 2k log a − r(log k−2log log k+2log a)
      ≤ 2kB + k(2 log a−a/2) + O(log k)
      ≤ 2kB + O(k).

The lower construction has the same error scale. Thus H(T_k) exceeds any exponential in k and cannot satisfy an eventual constant-coefficient linear recurrence over Q; finite-dimensional linear evolution on this one-letter physical language is impossible even without independent probes. This is an explicit growth mechanism (rare matched transitions, many unconstrained run assignments), not a proof of dense absolute coverage by the O(log X) or fewer values from this particular clock.

## Replay

Run `python ca_bridge_tm_check.py` and `python -O ca_bridge_tm_check.py`. Both emit the same JSON receipt: 1,099 reader contexts, 1,033 permanent-coefficient checks, 31 literal substitutions, 640 modular comparisons, eight frontier query checks, and 435 direct ordered-cover tests of the two-vertex band. The proof of the one-letter asymptotic is the two explicit combinatorial bounds above, not a fit to these finite checks.
