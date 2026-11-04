# Prime coverage requires growing constructor information

Session 176. **Proved within the specified transitive-inflation constructor model.** This is not a Boolean circuit lower bound, a global Hamilton-spectrum obstruction, or a new claim of literature priority.

## 1. Exact prime-relevant library

Let a finite library contain tournaments Q. Its permitted constructions are positive transitive inflation of every quotient vertex, reversal, and arbitrary ordered sums. The library may depend on a target cutoff N. All statements below apply separately at each N; there is no uniformity assumption on how the library is selected.

Decompose each Q into its ordered strongly connected components (SCCs). Call a component nontrivial if it has more than one vertex. A nontrivial tournament SCC has at least three vertices. Retain a quotient precisely when it has exactly one nontrivial SCC, and extract that SCC. Deduplicate the extracted cores up to isomorphism and reversal. Let the resulting cores be C_1,...,C_K with orders q_1,...,q_K. Write d=max q_j when K>0.

**Theorem 1 (exact prime-core reduction).** The prime Hamilton-path counts produced by the entire constructor class are exactly the union of prime counts produced by positive transitive inflations of the extracted cores C_j.

**Proof.** Hamilton-path counts multiply over ordered SCCs. Transitive inflation of a singleton SCC is transitive and contributes one. Positive transitive inflation of a nontrivial strong SCC remains strong: paths through other modules allow travel between any two internal vertices. Every such inflation has at least three Hamilton paths: each quotient vertex belongs to a directed triangle, and selecting one vertex from its three modules yields an induced triangle; insertion monotonicity gives the bound. Thus a quotient with two or more nontrivial SCCs always contributes a composite count, and a quotient with none contributes one. A quotient with exactly one contributes precisely the count of its inflated core. Conversely any inflation of an extracted core is realized with singleton sizes one at every padding vertex of its parent quotient. Reversal preserves counts. Ordered sums cannot create a new prime count because their counts are products of positive integers. Isomorphic cores and reversed cores have identical inflation spectra after permuting sizes. These observations prove both inclusions. ∎

This improves the inherited K_original*d overcount to the actual prime-relevant core inventory; it does not merely extract every component of a composite template and incorrectly allow it independently.

## 2. Finite coverage budget

For N>=3 set L=floor(log_2(N-1)). Let P_library(N) be the number of distinct prime counts at most N supplied by the library.

**Theorem 2.**

    P_library(N) <= sum_(j=1..K) L^(q_j) <= K L^d.

The empty sum is zero. More precisely, the first sum may be replaced by the number of orbits of size tuples in {1,...,L}^{q_j} under the automorphism group of each core, followed by summation; no orbit formula is needed for the stated bound.

**Proof.** The Session175 triangle lemma gives H(C_j[T_(r_1),...,T_(r_q)]) >= 2^(r_i)+1 for every i. Hence H<=N forces 1<=r_i<=L. There are at most L^q size tuples. Sum over the exact eligible core library from Theorem1. Size tuples related by a core automorphism yield isomorphic inflated tournaments and hence the same H, justifying the optional orbit refinement. ∎

For a completely explicit, elementary necessary condition for covering every odd prime other than7 through N, put m=floor(N/2), and let t be the least nonnegative integer such that

    (2m)^t >= binom(2m,m).

Then define B(N)=max(0,t-2). Coverage requires

    sum_j L^(q_j) >= B(N).

Indeed each prime-power contribution to binom(2m,m) is at most 2m, by the floor formula for its p-adic valuation, so binom(2m,m)<=(2m)^pi(2m). Thus t<=pi(2m)<=pi(N). Removing primes2 and7 removes at most two values. This lower bound is deliberately elementary and weaker than the exact prime count, but uses integer arithmetic only.

## 3. Library size versus interface order

**Corollary 3.** Suppose libraries indexed by N cover all odd primes other than7 up to N. If

    log K(N) <= (alpha+o(1)) log N,  with 0<=alpha<1,

then

    liminf_(N->infinity) d(N) log log N / log N >= 1-alpha.

In particular any subpolynomial-size library, including any bounded-size library, requires d(N)>=(1-o(1))log N/log log N.

**Proof.** The elementary estimate binom(2m,m)>=4^m/(2m+1) in the preceding argument gives B(N)>=(log2-o(1))N/log N. Theorem2 therefore yields

    log K+d log L >= log N-log log N+O(1).

Since log L=log log N-log log2+o(1), substitution proves the corollary. If K=0 coverage already fails for N>=3; logs are taken only for positive K. ∎

More generally, covering at least M(N) distinct prime values implies log K+d log L>=log M. If a library has K<=N^alpha and d<=(beta+o(1))log N/log log N with alpha+beta<1, its fraction of all primes through N tends to zero. This follows from K L^d<=N^(alpha+beta+o(1)) and the same elementary prime lower bound. Neither conclusion excludes infinitely many prime outputs.

## 4. Four creative bridges, with exact remaining obligations

1. **Constructor information and circuit resources.** Theorems1–3 turn a fixed-library obstruction into a quantitative tradeoff between library diversity and distinct interface roles. This is a resource bound for one tournament constructor language. To transfer it to B2 circuits would require a sound compiler from arbitrary circuits into this language, controlling K and d in terms of actual gate count. No such compiler is supplied or implied.
2. **Protected gate repair and contextual memory.** S123 saves |U|-nu(U) in a quadratic forest cone whose unprotected products form U. S171 says a compressed state is compositional only if its blind subspace is stable under every allowed continuation. A next concrete target is a deletion-state interface retaining protected-output incidence and future consumers, then an amortized inventory showing protections cannot be reused indefinitely to block distinct repairs. Neither matching counts alone nor a cycle-mean computation proves this global claim.
3. **Ramsey geometry and coordinate-price sidecars.** S117 proves one affine/Walsh orbit spans exact B2 costs n-1 through2n-2. A viable Ramsey-to-circuit bridge must retain the original-coordinate realization cost or a certified cheap change of basis. The constructive target is a boundary-preserving quadratic chart with an explicit XOR implementation ledger; geometric orbit invariants alone provably lose the needed price.
4. **Prime construction and endpoint automata.** S175 supplies exact endpoint-sensitive edge changes and orphan-path corrections; S171 supplies continuation-closed linear response quotients. Grow interface roles at least at the rate required here while minimizing the exact response quotient, then prove realized strong response intervals overlap. Current scalar H, isolated response gcd2, finite cubes, or fixed twin inflation do not establish interval coverage. This is a constructive route rather than a claim of completeness.

## 5. Sources and strongest audited antecedents

- [Session175 clone obstruction](https://github.com/eliottcassidy2000/math/blob/f06078dbe8329477c0f93adfca9e459e2dd8060e/05-knowledge/results/session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md): insertion monotonicity, triangle lower bound, K_original*d*L^d bound, elementary prime counting. The present exact core filter and variable-resource consequences refine that theorem.
- [Session123 protected boundary and forest repair](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session123_realizers_and_binary_capacity_20261031/circuits/RESOURCE_COST_AND_BOUNDARY_KERNEL.md): savings |E\\P|-nu(E\\P), exact simultaneous protected-coordinate kernel, and why semantic compression does not pay physical gates.
- [Session118](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session118_shared_gate_deletions_20261030/README.md): sharp four-deletion AP2, r incomparable constant seeds force2r distinct deletions, exact row frontier b>=max(c,2c-2,4c-8), c<=4. The coefficient-four circuit telescope remains conditional on full state coverage and termination.
- [Session117](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session117_ramsey_circuit_interfaces_20261030/README.md): exact point-defect circuit classification, affine/Walsh counterbenchmark, and skew-double transitive-subset enumerator. No new unrestricted circuit record.
- [Session171 response framework](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session171_knots_information_20261103/framework/NOTE.md): largest continuation-invariant blind subspace, horizon at most dimension minus observer rank, and exact square-zero visited-set Hamilton interfaces. No polynomial-size response theorem.
- [Session175 endpoint interface](https://github.com/eliottcassidy2000/math/blob/50a00fc7f473cde47cd2a360d5d4be4fe03aac92/05-knowledge/results/session175_tournament_spectrum_20261104/interface/PROOF.md): edge response, orphan paths, and exactly bounded seven-vertex response cubes.

The finite checker tests the new SCC prime filter on all labelled quotients through order4 with module sizes1 or2, checks genuine two-strong-component hostile templates, deduplication invariance, and the exact elementary prime bound through N=1000. These are controls on the all-order proofs, not a global spectrum census. No external software is installed.
