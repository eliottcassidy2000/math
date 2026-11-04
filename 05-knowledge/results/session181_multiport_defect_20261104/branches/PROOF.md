# Infinite transitive-module families still cannot close both ternary branches

2026-11-04. **PROVED structural obstruction; FINITE-EXACT controls.** This is not a construction of the missing reusable grammar. It enlarges Session180's finite fixed-suffix-bank obstruction to certain genuinely infinite, unbounded-order suffix families, and identifies a necessary growth coordinate.

## 1. Inheritance and the newly tested gap

The inherited evaluator is H(Z_s(B,w))=h+2^s F, where B is strong, w is a physical marked vertex, s>=1, h=H(B), and F=E_w(B,2). Session180 proves a finite bank of fixed B cannot support both H->3H+2 and H->3H+4 indefinitely. Merely requiring an infinite family of B does not remove that obstruction: allowing arbitrary lengths for finitely many transitive modules is still too small.

The relevant hostile boundary is a nonstrong quotient: TT2[TT_a,TT_b] is transitive and has H=1 regardless of a,b. Strongness is essential in the size argument below. Names and marks are retained in constructions; counting numerical possibilities here deliberately forgets names, which cannot increase the number of different counts.

## 2. Two elementary mechanisms

**Induced monotonicity.** If U is an induced subtournament of T, H(T)>=H(U). Fix an order for the extra vertices. Insert each extra vertex x into each current directed Hamilton path using a deterministic rule: prepend if x beats the first vertex; otherwise insert immediately before the first vertex beaten by x, or append if none exists. Each insertion is valid. Restricting the resulting path to U recovers the original path, so the map is injective. This proof concerns ordinary tournament paths, not native complement guards.

**Triangle blowup.** Let C3[TT_a,1,1] have blocks A->b->c->A. Its Hamilton-path count is 2^a+1. Three paths have A contiguous. Every other path is A_1,b,c,A_2, where A_1,A_2 are nonempty complementary subsets of A, each in its inherited transitive order. There are 2^a-2 such choices.

Every vertex i of a nontrivial strong tournament Q lies in a directed triangle. Indeed, if no edge goes from an outneighbor of i to an inneighbor of i, then the outneighbors cannot reach i, contradicting strongness.

Therefore for every transitive blowup B=Q[TT_{a_1},...,TT_{a_k}] of a strong quotient Q of order k>=3,

    H(B) >= 2^{a_i}+1 for every i.                  (1)

Choose a directed triangle through i, retain its entire i-block and one vertex from each of its other two blocks, and apply induced monotonicity. All a_i are positive integers.

## 3. Fixed strong templates, unlimited module lengths

Fix a finite bank of m nontrivial strong quotient tournaments, each of order at most k (hence each has order at least3). Allow every positive transitive module length, every physical marking w anywhere in B, every transitive prefix clock s>=1, arbitrary switching between templates, and arbitrary finite or infinite observer labels which do not change the physical tournament.

For any represented full count N=H(Z_s(B,w))<=X, h<=X and (1) imply every module length <=L=floor(log_2 X). Also 2^s<=N implies s<=L, since F>=1. A quotient of order j contributes at most L^j length vectors; for each vector there are at most jL marking ranks; and there are at most L clocks. Thus the entire family supplies at most

    m k L^{k+2}                                    (2)

different numerical counts below X. Distinct names, isomorphic presentations, repeated marks with equal responses, and extra observer labels only overcount these possibilities. Singleton strong suffixes contribute at most L additional counts and do not alter the conclusion; a nontrivial transitive suffix is not strong.

Starting at N0, a depth-n word of digits d_i in {2,4} produces

    3^n N0 + sum_i d_i 3^{n-i}.

The 2^n words give distinct numbers, by uniqueness of base-3 expansions after subtracting the common all-2 word, and all lie below X_n=3^n(N0+2)-2. Because L_n=O(n), (2) grows polynomially in n, whereas 2^n grows exponentially. Consequently **no such family closes both ternary branches indefinitely**, even though its suffix order is unbounded and its family is infinite.

The same argument without marks and clocks says that a finite bank of strong quotient templates with arbitrary transitive module lengths cannot itself realize an indefinitely closed two-branch family of Hamilton-path counts.

## 4. Quantitative necessary quotient growth

Suppose a construction uses arbitrary strong quotient templates, not a fixed bank, but at depth n every output suffix can be partitioned into at most K transitive modules with uniform intermodule orientations. There are at most 2^{j(j-1)/2} labelled quotient tournaments of order j. Ignoring strongness only enlarges this enumeration. Applying the previous count bound to every j<=K gives, for K>=3 and L_n>=1,

    2^n <= K^2 * 2^{K(K-1)/2} * L_n^{K+2} + L_n.  (3)

The final L_n covers singleton suffixes. In particular a successful construction cannot have K=o(sqrt(n)). More sharply the necessary maximum K at depth n is at least (sqrt(2)-o(1))*sqrt(n): if K is already larger this is immediate; otherwise K=O(sqrt(n)), and taking base-2 logs in (3) gives n<=K^2/2+o(n).

For a fixed or polynomial-sized menu of templates at each depth, its quotient-order requirement is stronger: the same calculation with its actual menu size M_n gives 2^n<=M_n K L_n^{K+2}+L_n, hence K=Omega(n/log n) when M_n is polynomial in n. These are necessary complexity bounds, not sufficient conditions for realization and not lower bounds on the minimal order of unrestricted tournament witnesses.

The concrete missing coordinate is therefore **growth in the minimum number of transitive uniform modules, or growth of non-transitive internal module structure**, not merely suffix size, added physical names, longer transitive runs, or observer phase. A recursive family with non-transitive modules of unbounded internal depth is outside this theorem. So are arbitrary nonuniform attachments and operations which leave these module families.

## 5. Exact controls and scope

Run `python verify.py` or `python -O verify.py`. All checks are explicit exceptions, not removable assertions. Standard-library subset endpoint DP verifies:

- H(C3[TT_a,1,1])=2^a+1 for a=1,...,11;
- every vertex belongs to a triangle in all 26 labelled strong quotients of orders3 and4;
- strongness and (1) in 1,878 positive module-length assignments, lengths1..3 and total order<=10;
- the nonstrong hostile quotient with lengths7,7 has H=1;
- exact branch uniqueness and maximum through depth12.

These controls audit the stated proofs; they are not all-order evidence in place of proof. The raw `search.cpp` and its C2/C4 files are exploratory unfiltered fixed-port two-birth data and are NOT strong-suffix certificates. In particular the two C4 +2 hits contain a sink. They must not be published as strong continuations. The separate construction worker has the correctly strong-filtered search.
