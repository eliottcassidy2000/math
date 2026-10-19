# Twin-signature automorphisms refute hereditary symmetry as a characterization of good trees

19 October 2026. **PROVED, scoped, author checked:** the all-host counterexample, independent-twin divisibility, and rooted orientation-capacity recurrence below. **FINITE-EXACT:** the explicitly delimited accompanying census. Independently [reviewed](review/INDEPENDENT_REVIEW.md); no historical-priority claim is made.

## 1. The inherited open question and its first failure

For an oriented graph S on n vertices and a tournament T on n vertices, let emb(S,T) count bijections carrying every arc of S to an arc of T. Write phi_S(T)=emb(S,T) mod 2. The graph S is parity-rigid when phi_S is constant. An undirected graph U is good when every orientation of U is parity-rigid. A tree U is hereditarily symmetric when every subtree with at least two vertices has a nontrivial automorphism.

THM-4533, [THM-4533](../../../01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md), proves hereditary symmetry implies goodness. HYP-9173(b), [HYP-9173](../../hypotheses/HYP-9173-redei-graphs-growth-and-trees.md), asks for the converse, with finite support through eight vertices. The source is Sections 2–4 of [the original Rédei note](../procgen_shave4_20261001_redei_graphs.md). The later [inheritance bullet](../trees5_tournaments4_erdos_sos_20261005.md) had dropped that open qualifier; this packet corrects it along with the hypothesis and theorem routing.

**Theorem 1 (nine-vertex counterexample).** Let U have vertices 0,...,8 and edges

    01, 02, 03, 04, 05, 16, 27, 68.

Then U is good, every one of its orientations has identically EVEN embedding counts in every tournament, and U is not hereditarily symmetric.

**Proof.** Vertices 3,4,5 are leaves adjacent to 0. In every orientation, at least two of the three incident edges point the same way relative to 0. Swapping those two leaves preserves the oriented graph and has order two. Precomposition by this swap acts freely on all injective embeddings into any host, so their number is even.

The induced subtree on {0,1,2,3,6,7,8} is the spider with leg lengths 1,2,3 from center 0. Its unique degree-three vertex is fixed by every automorphism. Its three branches have different lengths and hence are fixed separately; each rooted path is rigid. This subtree has trivial automorphism group. Thus U is not hereditarily symmetric. ∎

There are only 2^8=256 labelled orientations to inspect, but no orientation census or tournament bound supplies the universal proof: the transposition exists for EVERY orientation and acts on embeddings into EVERY host. The source's through-eight census makes nine the first possible counterexample order. The companion also supplies a self-contained finite verification of that minimality, rather than requiring the old census to be rerun.

This tree has undirected automorphism group S3, acting on the three leaves. No one nonidentity transposition works for every orientation: its choice depends on the direction signature. This is the lost quantifier in attempting to upgrade an undirected symmetry obstruction to a hereditary characterization.

**Corollary 2 (arbitrary hereditary complexity inside good trees).** Every finite tree embeds induced in a good tree for which every orientation has identically even embedding counts. Attach three new leaves to any vertex of the given tree and use the same proof. Consequently goodness cannot be characterized by excluding any collection of subtrees or induced subtrees: goodness is not hereditary. This does not invalidate the proved sufficient condition of hereditary symmetry.

## 2. The general signature mechanism and its exact divisibility budget

Call a set W of m vertices an independent twin class if W is independent and every vertex of W has the same neighbor set N outside W. Put d=|N|. For any orientation S of the entire graph, each w in W has its d-bit signature of directions toward N. There are k=2^d possible signatures. Let their multiplicities be n_sigma, with sum n_sigma=m.

**Theorem 3 (signature subgroup).** There is an orientation-preserving subgroup

    product_sigma Sym(n_sigma) <= Aut(S).

Therefore product_sigma(n_sigma!) divides emb(S,H) for every finite directed host H, where embeddings are injective and need not be spanning. In particular, if m>2^d, every orientation has even embedding count in every host.

**Proof.** Arbitrary permutations within one signature class preserve all incident arcs, because the vertices have no edges among themselves, share exactly the neighbors N, and have matching directions to each neighbor. They fix everything outside W. The resulting direct product acts freely by precomposition on injective embeddings. ∎

The independence hypothesis matters. A tournament on W itself can break every transposition. The m>2^d threshold is exact for forcing duplicate signatures: if m<=2^d, distinct signatures can be assigned. This only says this subgroup obstruction disappears; it does not in general guarantee odd embeddings. For d=1,m=2 the directed three-vertex path is an actual odd-embedding positive control, so the three-leaf threshold itself is sharp.

The subgroup gives more than parity. Define

    V_k(m) = min_{n_1+...+n_k=m, n_i>=0} sum_i v_2(n_i!).

**Proposition 4 (exact orientation-independent subgroup budget).**

    V_k(m) = sum_{j>=1} max(0, ceil((m+k)/2^j)-k).

Thus 2^{V_{2^d}(m)} divides every oriented embedding count whenever the source has such an independent twin class. This is the sharp guarantee obtained from the signature subgroup; other automorphisms or cancellations may increase divisibility.

**Proof.** V_k is nondecreasing: removing one unit from any nonzero part of a partition of m+1 cannot increase its factorial valuation. Write n_i=2b_i+epsilon_i, epsilon_i in {0,1}. If q=sum b_i then

    sum_i v_2(n_i!) = q + sum_i v_2(b_i!),
    q >= q_0 = max(0,ceil((m-k)/2)).

Monotonicity gives the lower bound q_0+V_k(q_0). Conversely choose a minimizing partition of q_0 and distribute m-2q_0 epsilon bits among its k entries; that number lies between zero and k. This realizes the bound. Hence V_k(m)=q_0+V_k(q_0). Iterating yields q_j=max(0,ceil((m-k(2^j-1))/2^j)), and summing q_j proves the formula. ∎

All multiplicity distributions are realizable by orienting the twin-neighbor edges independently, proving sharpness at the level explicitly claimed. The elementary free-action argument is the same mechanism as THM-4533 B3(a); the new point is to retain the changing direction signatures BEFORE reduction modulo two. The quotient loses the size and selected transposition of the stabilizer, so signature multiplicities are the necessary sidecar.

## 3. Exact recursive detection for every tree

The three-leaf test is a first instance of a complete recursive test for whether every orientation has a nontrivial automorphism.

For an un-oriented rooted tree R, let a(R) be the number of isomorphism classes of rooted orientations having trivial rooted automorphism group. Let its pairwise nonisomorphic rooted child-subtree types be T, of multiplicities m_T. A single vertex has a=1.

**Theorem 5 (rooted orientation capacity).**

    a(R) = product_T binom(2a(T), m_T),

with binom(b,m)=0 when m>b.

**Proof.** A rooted automorphism can independently act inside child branches and permute isomorphic oriented child branches. Every child must therefore have trivial rooted automorphism, giving a(T) choices. The parent-child edge has two directions, so there are 2a(T) available branch signatures. Repetition of a signature gives a branch-swapping involution. Conversely, if all selected signatures are distinct and each chosen child is asymmetric, every rooted automorphism is trivial. Since branches of a fixed underlying type are unordered, selecting m_T distinct signatures contributes exactly the displayed binomial coefficient. Different underlying types are independent. ∎

Let A(U) count isomorphism classes of orientations of an unrooted tree U with trivial automorphism group. Every automorphism preserves its center.

- If the center is a vertex v, then A(U)=a(U rooted at v).
- If the center is an edge with nonisomorphic rooted halves R1,R2, then A(U)=2a(R1)a(R2).
- If the rooted halves are isomorphic to R, then A(U)=a(R)^2.

For the last line, orient the central edge from its source half to its target half. An oriented isomorphism must preserve source and target, giving an ordered pair of asymmetric rooted half types. Reversing the edge does not double the count because the un-oriented halves are interchangeable.

**Corollary 6 (complete automorphism-forced even test).** A(U)=0 if and only if every orientation of U has a nontrivial automorphism; this is also equivalent to every orientation having an automorphism of order two. In that case U is good with every parity value zero.

Indeed, recursively the automorphism groups of rooted oriented trees are products of wreath products of symmetric groups. A nontrivial one either has a nontrivial child automorphism, hence recursively an involution, or repeats an oriented branch and has a transposition. For an oriented tree the central edge cannot be swapped by any automorphism. This proves the equivalence without assuming every arbitrary nontrivial finite group has even order.

A(U)>0 is NOT claimed to imply that any orientation is Rédei or nonrigid. It merely removes this automorphism-forced-zero explanation. The property of being good remains an embedding-parity question.

## 4. Corrected frontier

The valid implications are:

    hereditary symmetry => good;
    A(U)=0 => every orientation has identically zero parity => good;
    good does NOT imply hereditary symmetry.

The nine-vertex example is in the second class and not the first. For orientation-specific questions, an even automorphism group forces zero; an asymmetric undirected underlying graph forces nonconstancy by the highest-degree ANF argument in THM-4533 B4. These are different statements and are both retained.

The original tree characterization is therefore REFUTED, not open. Its failure is not a missing tournament sample: it is a universal orbit pairing that turns parity blindness into goodness. The largest Rédei/shaving questions and HYP-9173(a),(c) are untouched, since our counterexample has even rather than odd embedding counts.

The natural next question is a classification after keeping the automorphism-forced-zero branch explicit. No unrestricted replacement equivalence is asserted here.

## 5. Why restricting the host class is a real information loss

For an arbitrary finite oriented graph S,

    emb(S,H) is even for every finite directed host H
    iff |Aut(S)| is even.

The forward implication follows by taking H=S: an injective self-map on the same finite vertex set is a permutation, and preserving all arcs sends its finite arc set bijectively to itself. Thus emb(S,S)=|Aut(S)|. The reverse implication is the free automorphism action. This elementary statement concerns injective embeddings, not noninjective graph homomorphisms.

**Theorem 7 (asymmetric P5 hostile for tournament-only detection).** The oriented path

    1 -> 0, 1 -> 2, 2 -> 3, 3 -> 4

has trivial automorphism group but has even embedding count in EVERY five-vertex tournament.

**Proof.** Its undirected path admits only identity and reversal. Reversal does not preserve the displayed arcs, so its oriented automorphism group is trivial. Deletion–reversal on the first edge says its embedding parity is the sum of the parity of a directed P5 and that of an isolated vertex disjoint from a directed P4. Rédei's odd-path theorem makes the first parity 1. For the second, each of the five choices of isolated-vertex image leaves a four-tournament with an odd number of directed Hamiltonian paths, giving parity 5=1 modulo two. The sum is zero. ∎

The same proof gives an infinite sharp family: reversing only the first edge of a directed P_n has constant tournament parity 1+n modulo two. For n>=4 this orientation is asymmetric, so every odd n>=5 gives an asymmetric universal-even tournament hostile. At n=3 the first-edge-reversed path has a leaf-swapping automorphism; this small exception is retained.

The required Rédei input and the deletion–reversal identity are inherited from THM-4533 Sections B1/B3, with a complete finite control on all 1,024 labelled five-tournaments. Taking this oriented path as its own non-tournament host instead gives exactly one embedding. Accordingly, tournament-only evenness cannot be inverted to an even-automorphism claim, even for trees. No modular homomorphism-count theorem is imported into this injective, restricted-host setting.

## 6. Exact scope and reproduction

Run from this packet directory:

    python3 verify_twin_signature.py
    python3 -O verify_twin_signature.py

Both runs produce the same `twin_signature_receipt.json` data. The file `minimality_witnesses.json` supplies eight explicit orientation/host pairs, one for each nonhereditarily symmetric tree on seven or eight vertices. Every pair has different embedding parity in its listed tournament and the transitive tournament. The verifier checks that these pairs cover the complete generated universe, checks each count both by tree backtracking and by independent permutation enumeration, and verifies both hosts are tournaments. Every tree through six vertices is hereditarily symmetric, completing the finite proof that the nine-vertex counterexample is minimal.

The separate capacity verification enumerates 15,943 orientation masks across all unlabelled trees through nine vertices. It compares the resulting asymmetric orientation-isomorphism classes with the recursive formula. The explicit nine-vertex counterexample also has all 256 labelled orientations checked for even automorphism order. The factorial budget is checked against exact min-plus optimization for 520 pairs 1<=k<=8, 0<=m<=64. The P5 hostile is checked on all 1,024 labelled five-tournaments and on itself as a directed host.

The unlabelled-tree census is generated recursively by adjoining a leaf at every vertex and deduplicating center-rooted canonical forms; every nontrivial tree has a leaf, so this generation is exhaustive. Here H means hereditarily symmetric and Z means A(U)=0, the automorphism-forced-even class:

    n   all trees   H    Z    Z outside H
    1       1      1    0       0
    2       1      1    0       0
    3       1      1    0       0
    4       2      2    1       0
    5       3      3    1       0
    6       6      6    2       0
    7      11     10    4       0
    8      23     16    9       0
    9      47     24   19       3
    10    106     35   46      18
    11    235     49  106      64

This is a complete census of H and Z through eleven, NOT a complete census of all good trees above eight. Failure of both tests does not imply nonrigidity. No random tournament probe contributes a proof in this packet.

The source files named in Section 1 are sufficient dependencies; the new results have direct proofs. The connection to finite groups is the free action on injective embeddings, with direction signatures and rooted branch capacities retained before taking the parity quotient. It is not a connection to planar Jacobian monodromy or quadratic period-three dynamics, and no such transfer is claimed.
