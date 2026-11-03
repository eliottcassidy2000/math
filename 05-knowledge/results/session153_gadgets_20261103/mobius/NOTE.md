# Hamilton paths via an apex, α-permanents and partition inversion

For any finite directed graph T on nonempty vertex set V, let A be its adjacency matrix, and adjoin a new vertex 0 with both arcs to every vertex, producing matrix B. Denote by H(T[S]) the number of directed Hamilton paths on S, by per the ordinary permanent, and by per(empty)=1. Then

    H(T) = ∑_{π ∈ Π({0} ∪ V)} (−1)^{|π|−1} (|π|−1)! ∏_{C ∈ π} per(B[C]).

This is valid for arbitrary loopless directed graphs, including tournaments; single-vertex paths count as one. The formula applies to the generally nontournament apex extension.

## Proof

Define per_α(B)=∑_{σ} α^{c(σ)} ∏_i B_{i,σ(i)}, where c(σ) is the number of cycles in σ. A one-cycle permutation of B is a directed Hamilton cycle, and cutting at vertex 0 gives a bijection with a directed Hamilton path of T. Thus H(T)=[α]per_α(B).

For a nonnegative integer q, assign each cycle of σ one of q colors. Group vertices by their nonempty color classes, a set partition π, and choose distinct colors for these classes in (q)_{|π|}=q(q−1)...(q−|π|+1) ways. The weighted permutations on each class sum to its principal permanent. Therefore per_q(B)=∑_π (q)_{|π|} ∏_C per(B[C]). Both sides are polynomials, so the identity holds for general α. Since [α](α)_k=(−1)^{k−1}(k−1)!, coefficient extraction proves the formula. These coefficients are the Möbius function μ(π,1̂) of the set-partition lattice.

There is also the direct apex-cycle convolution

    per(B[{0}∪V]) = ∑_{∅≠S⊆V} H(T[S]) per(A[V\S]).

Indeed, the cycle through 0 encodes the path on S, and the other cycles form an arbitrary cycle cover of the complement. Isolating S=V gives a recursion. Both identities use only integer arithmetic.

## Checks and implications

`check.py` independently counts the permanent by direct permutation enumeration, sums over set partitions, and compares with literal Hamilton-path permutation counts. Exhaustive tests cover all 1,099 labeled tournaments on 1–5 vertices. They also cover the 7-vertex cyclic tournament with edges i→j if j−i mod 7 is 1, 2, or 3: the contributions by 1, 2, 3, ... partition blocks are 371, −224, 28, 0, ... and give H=175. Run with both `python check.py` and `python -O check.py`; both transcripts are included. Explicit error raising keeps checks active under optimization.

This formula is an exact alternative and a possible route to bounds on cycle-cover corrections, not a proof that 7 or 21 are absent: signed partition terms cancel and no inequality excluding either value has been derived. The check is not a classification for order 6 or higher. The identity is not claimed as a new literature result.

## Context and sources

The α-permanent and partition decompositions are established in H. Crane, “Some algebraic identities for the α-permanent,” https://arxiv.org/abs/1304.1772 . Original α-permanent context: D. Vere-Jones, “Determinant, permanent, and MacMahon’s master theorem,” Linear Algebra and its Applications 258 (1997), 171–183, https://doi.org/10.1016/S0024-3795(95)00774-1 . The apex specialization and derivation above are elementary applications, not a claim of priority.

Repository audit: `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md` already uses partition-lattice cluster expansions of signed linear forests, a different construction; searches for α-permanent and partition lattice in canon and hypotheses found that related item. This is a targeted audit, not an exhaustive novelty search.
