# Census-free aggregate enumeration of realizable cyclic-triple hypergraphs

**Status:** PROVED consequence of cited structural theorems plus the recovered ordinary/twisted tournament indices; FINITE-EXACT full cycle index through 15; independently verified on actual hypergraph image sets through 6. No new C3 structure theorem or literature-novelty claim.

## Exact object, inherited corrections, and cited inputs

A structure in C on a finite labeled set V is a **distinct 3-uniform hypergraph** H=C3(T) for at least one tournament T on V. Hyperedges are precisely the vertex triples inducing directed triangles. Count H once regardless of the number of realizing tournaments. Include the empty labeled set. This is a quotient by all C3-preserving realization moves, not merely global converse. Let D consist of the singleton and those H whose realizations are strongly connected.

Use the precise module notion of [Boussaïri–Chergui–Ille–Zaidi, Definition 1](https://arxiv.org/html/1805.04607v1): a crossing hyperedge must meet a module M in exactly one vertex, and replacing that vertex by any element of M must again give an edge. A strong module overlaps no module without containment. This is **not** the alternative hypergraph-module convention discussed in their Remark 26.

The needed cited inputs from that paper are: Theorem 10 identifies the strong modules of H with those of every realizing tournament; Theorem 11 equates primeness; Theorem 24 states a prime realization is unique up to converse; Proposition 51 says allowable realization choices at canonical decomposition nodes are independent. Gallai's tournament decomposition is Theorem 16. These structural results are imported, not proved by this packet's finite computation.

The Session67 correction remains mandatory: equal C3 hypergraphs do **not** imply equal tournament automorphism orders. Its order-12 example C3[A,A,A] versus C3[A*,A,A], with A=singleton→directed triangle, has ordinary automorphism orders 81 and 27. The present equations concern hypergraph symmetries and never replace them by those ordinary tournament groups.

## 1. Prime root quotient and its actual action

Let R be the species of prime realizable cyclic-triple hypergraphs, on at least three vertices. The cited prime uniqueness result gives a natural bijection

    R = P / <τ>.

This is a quotient on labeled structures; τ acts freely there because a tournament on at least two fixed labels differs from its converse. It does not mean that every prime tournament is self-converse up to isomorphism.

Let P_e,P_τ be the ordinary/twisted indices recovered in [the converse recovery](CONVERSE_RECOVERY.md). For each relabeling σ, a σ-fixed pair {Q,Q*} contributes exactly two in the combined count of Q fixed by σ and Q fixed by τσ. Averaging yields the full root index

    Z_R = (P_e+P_τ)/2.                                     (1)

Equivalently, for an individual Q the hypergraph automorphism group is Aut±(Q): a permutation preserving C3(Q) must carry Q to Q or Q*, by the cited prime uniqueness. For a non-self-converse Q there are no anti maps on that representative and the two distinct tournament isotypes give one hypergraph isotype. For a self-converse Q its anti coset doubles the ordinary group. Formula (1) handles both cases without enumerating either kind.

## 2. Why the outer construction is SET, not SEQ or global converse

For any realization T, different strong components have no cyclic triple crossing between them. Thus their C3 structures form a disjoint family on the SCC vertex partition. Every order of those components realizes the same H. The partition is nevertheless recoverable from H: at a nontrivial empty root quotient, the common maximal strong modules are precisely the SCCs of every realization, by the cited strong-module identification and the canonical tournament decomposition. Their **order** is forgotten, not their vertex blocks.

Each SCC has singleton or prime root. If H has prime root, the common maximal strong-module blocks, its prime cyclic-triple quotient, and its induced child hypergraphs recover its decomposition uniquely. Conversely, choose any prime root H_Q and arbitrary nonempty realizable child hypergraphs. Realize each child and choose either prime tournament realization Q; tournament inflation constructs the required H. The inflated blocks are canonical maximal proper strong modules because the quotient tournament is prime, and remain so under the cited strong-module identification. Thus the construction neither loses nor adds decompositions.

This gives natural species identities

    C = SET(D),              D = X + R ∘ (C−1).             (2)

No generic hypergraph connectedness assertion is needed: D is defined by singleton/prime canonical root, equivalently by strong tournament realizability. In particular there is no arbitrary nesting or parenthesizing of empty quotient nodes. The outer SET flattens exactly the maximal SCC blocks.

Combining (1)–(2), with ψ_j substituting p_i→p_{ij}, gives

    Z_D = p1 + ((P_e+P_τ)/2) ∘ (Z_C−1),
    Z_C = exp(Σ_{j≥1} ψ_j Z_D/j).                           (3)

All substitutions in (3) are now ordinary species plethysm: the C2 action has already been quotiented at each prime root. This differs from the parity-sensitive tournament converse substitution used to recover P_τ.

Because R has minimum order three, degree n of Z_D depends only on lower degrees of Z_C. Degree n of Z_C then follows from the exponential. This is a unique all-order formal recurrence with explicit partition-based input, not a finite supplied prime census.

At the scalar isotype level, if c(x),d(x) specialize p_i=x^i, the SET identity is the Euler transform

    c(x)=∏_{m≥1}(1−x^m)^(−d_m).

The root term still needs the full Z_R slot action; replacing it by the scalar count of prime root types is invalid.

## 3. Computation, reproduction, and independent image-set census

`recover_c3.py 15` obtains the complete C and D cycle indices using recovered P_e and P_τ, without supplied prime data. To evaluate the exponential exactly, put L=Σ_j ψ_j Z_D/j and use the weighted Euler derivation:

    n (Z_C)_n = [ (Σ_{k≥1} k L_k) Z_C ]_n.

| n | all realizable C3 isotypes | singleton/prime-root C3 isotypes |
|---|---:|---:|
| 1 | 1 | 1 |
| 2 | 1 | 0 |
| 3 | 2 | 1 |
| 4 | 3 | 1 |
| 5 | 9 | 6 |
| 6 | 32 | 22 |
| 7 | 247 | 214 |
| 8 | 3328 | 3074 |
| 9 | 93612 | 90255 |
| 10 | 4774417 | 4680547 |
| 11 | 447114461 | 442336617 |
| 12 | 76601960166 | 76154750809 |
| 13 | 24193986713786 | 24117379959402 |
| 14 | 14176434888156918 | 14152240453790335 |
| 15 | 15496299263240448609 | 15482122751724382239 |

The values through ten agree with the prior Session67 supplied-bank evaluation. Above ten they are outputs of (3), not an independent large hypergraph census.

`check_equivariant.py` performs an independent finite check with no modular decomposition or tournament-automorphism assumption: exhaust all labeled tournaments through six; compute each actual cyclic-triple edge bitset; deduplicate identical bitsets; enumerate vertex relabelings directly on the hypergraphs; compute their relabeling orbits and **hypergraph** automorphism cycle types. It also computes tournament strong connectivity directly and checks that every realization of a given bitset has the same flag. The accumulated full C and D indices agree coefficientwise with (3), not merely after scalar specialization.

The distinct labeled H counts at n=0,…,6 are 1,1,1,2,11,208,8563. At order six the 32768 tournaments give 8563 distinct labeled hypergraphs, 32 hypergraph isotypes, and 22 singleton/prime-root isotypes. This is a census of the realizable image, not all 2^{binom(n,3)} hypergraphs.

    python3 recover_c3.py 15
    python3 check_equivariant.py
    python3 -O recover_c3.py 15
    python3 -O check_equivariant.py

The JSON receipts include full indices. No Python assert is used; normal/optimized outputs are compared. This method counts the C3 quotient and its aggregate symmetries, but supplies no individual canonical representatives, tournament metagraph edge incidence, or ordinary tournament automorphism invariance.


## Independent audit

A second full-index implementation used literal truncated exponential powers
rather than the producer's Euler-derivative recurrence and reproduced every C
and D coefficient through 15. Its separate actual hypergraph-image census
checked all coefficients through five. The public combined checker retains
the producer's independent image-set census through six. The cited structural
theorems were checked in the primary paper with its stated module convention.
