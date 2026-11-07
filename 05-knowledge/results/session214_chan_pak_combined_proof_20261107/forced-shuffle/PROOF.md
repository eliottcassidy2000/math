# A forced shuffle proves the ordered-path half-bound

**Status: PROVED; independent adversarial proof review PASS. Finite checks are recorded separately.** The unrestricted Chan–Pak conjecture remains **OPEN**. No external literature-priority claim is made.

## 1. Inheritance and the missing map

Write H(X) for directed Hamiltonian vertex words in a finite tournament X. Singletons have H=1; the empty tournament has H=1 only as a product convention. Every nonempty tournament has a positive odd H (Rédei).

At the freshly checked main commit `4c63fd58b566f3e006eb87eaedcf48824741cc84`, [HYP-9469](https://github.com/eliottcassidy2000/math/blob/4c63fd58b566f3e006eb87eaedcf48824741cc84/05-knowledge/hypotheses/HYP-9469-minimal-deletion-increment-grows-linearly-small-increment-structures-one-centre-value-set-and-what-remains-for-chan-pak.md), statements (1), (2) and (6), still designate the general positive-increment bound, finiteness of increment levels, and ordered-path lemma as numerical conjectures. THM-4876 proves the endpoint special case and the insertion-slot reduction. Session213 classifies equality in that endpoint special case, but does not give the general cut injection.

Merging two directed words while preserving their internal orders is inherited: see THM-4745, `01-canon/theorems/THM-4745-exact-gluings-of-interval-chains-and-the-reduction-of-busch-equality-to-a-one-way-lemma.md`, statement 1, crediting Busch's Lemma 2. The contribution here is to force a **nonordered** merge whenever any cross arc points in the forbidden direction. Restriction to the two colours retains the complete input pair, repairing the collision issue noted at the end of HYP-9469(6).

Closest mechanism: word merging plus restriction recovery. Hostile: an ordered word whose last A-vertex has no exit into B and whose first B-vertex has no entrance from A, defeating the earlier endpoint-only moves. Corrected near miss: two individually injective endpoint moves can collide. Least-used sidecar: both complete monochromatic subsequences. Concept board: forced shuffle / ordered block / product of induced path counts / insertion slots / strong-component factorization. The anchor is path-spectrum construction; the niche result closes the general extension-growth conjecture; the wildcard is exact low-increment classification.

## 2. Deterministic ordinary merging

Let a=(a_1,...,a_s) and b=(b_1,...,b_r) be disjoint directed words. Compare their current first unused vertices and append whichever beats the other, then advance its pointer. When one word is exhausted, append the other remainder.

This is a directed word preserving both input orders. Within a colour, consecutive output vertices are consecutive in that input. At a colour switch, the previous vertex beat the other current head when it was selected. Empty input words are permitted. This is the standard merge, included only to specify a deterministic map.

## 3. Forced nonordered merging

Let A,B be a nonempty partition of V(X), and suppose there exists an arc from A to B. For **every** pair of Hamilton words

    a=(a_1,...,a_s) of X[A],     b=(b_1,...,b_r) of X[B],

choose

- i greatest such that a_i beats at least one vertex of B;
- j least such that a_i beats b_j.

Let M be the deterministic ordinary merge of (a_1,...,a_(i-1)) and (b_1,...,b_(j-1)). Define

    Phi(a,b) = M, a_i, b_j,...,b_r, a_(i+1),...,a_s.       (1)

**Lemma 1.** Phi(a,b) is a Hamilton word of X, it preserves both input subsequences, and it does not list all B before all A. The map Phi is injective on the set of input pairs.

**Proof.** M is directed. If M is nonempty, its last vertex is either a_(i-1), which beats a_i by the input word, or b_(j-1), which beats a_i by minimality of j. The edge a_i -> b_j is part of the choice. The b suffix is directed. If the final a suffix is nonempty, every vertex of B beats a_(i+1), by maximality of i; in particular b_r -> a_(i+1). The a suffix is directed. These checks also cover i=1, i=s, j=1, j=r and an empty M.

Every vertex appears once, and both colour subsequences have their original order. The occurrence of a_i immediately before b_j makes the word nonordered. Restricting the output word to A and B recovers a and b exactly. Thus different pairs have different outputs. QED.

Only the internal consecutive arcs of the input words and the complete cross relation are used. No claim that arbitrary permutations of A or B remain directed is needed.

## 4. Product surplus and the ordered-path half-bound

Put

    h=H(X),   p=H(X[A])H(X[B]),
    q=#{Hamilton words of X listing every B vertex before every A vertex}.

**Theorem 2.** Whenever there is an A-to-B arc,

    h >= p+q,       q <= min(p,h-p),       2q+1 <= h.      (2)

**Proof.** Lemma 1 supplies p distinct nonordered Hamilton words. They are disjoint from the q ordered words. Also q<=p, because restriction to the two blocks injects every ordered Hamilton word into an input pair. Therefore h>=p+q>=2q. Since h is odd, h>=2q+1. QED.

The last inequality proves HYP-9469(6) in all orders. The more informative bound h>=p+q retains a product sidecar forgotten by the scalar half-bound. It holds even when q=0. If no A-to-B arc exists, X is the dominance join B=>A and h=p=q, so (2) is not asserted; this hypothesis is essential.

The half-bound is sharp: for the reversed-extremal-arc tournament R_n, n>=3, take A to be its penultimate vertex and B the rest. Then h=2^(n-2)+1 and q=2^(n-3). These are the inherited endpoint equalities; they are not a new classification of general-cut equality.

## 5. Universal positive extension growth

Add v with v->A and B->v; call the new tournament Y. For a Hamilton word P of X, insert v in each admissible slot: before the word if its first vertex is in A, after it if its last vertex is in B, or between consecutive B,A vertices. If t(P) is the number of consecutive A,B transitions in P, the number of slots is exactly t(P)+1. This follows by counting the transitions of the binary colour word. It is one precisely for an ordered B-then-A word, and at least two otherwise. The slot map is injective because deleting v recovers P and its slot. Consequently the inherited bound is

    H(Y) >= 2h-q.                                       (3)

**Theorem 3.** Exactly one of the following occurs:

- There is no A-to-B arc. Then Y=B=>{v}=>A and H(Y)=h.
- There is an A-to-B arc. Then

      H(Y) >= max(h+p,2h-p) >= (3h+1)/2,                (4)
      (H(Y)-h)/2 >= ceil((h+1)/4).                      (5)

**Proof.** In the first case, concatenate the two induced Hamilton words; empty A or B is allowed. In the second case, (2) and (3) give (4). The raw difference is even because both counts are odd, giving (5). An A-to-B arc makes a directed triangle v,a,b,v, and (4) gives a strictly positive increment. Thus the two alternatives are exhaustive, not merely sufficient. QED.

This proves HYP-9469(1) with the strengthened bound (5), for all h including h=1 and h=3. It proves the finite-parent part of HYP-9469(2): if a positive deletion increment is exactly 2k, then

    h <= 4k-1,                 H(Y)<=6k-1.              (6)

The number of possible **counts** is finite; this is not a claim that there are finitely many labelled tournaments or orders, because arbitrary singleton components can be padded.

### Sharp asymptotic growth constant

For n>=4, take X=R_n and add v whose sole outneighbor is the penultimate vertex a=n-2. Then

    H(X)=2^(n-2)+1,           H(Y)=3*2^(n-3)+3.          (7)

For clarity, an elementary count of Y is as follows. Write b=n-1 for the last vertex and 0 for the first, with only b->0 reversed in R_n. Paths with v last number H(X)-e_a=2^(n-3)+1. Paths with v first number one. If v is internal, its successor is a. If a is last, the count is H(R_n-a)=2^(n-3)+1. Otherwise the next vertex is b. If b is last there is one path. If b is not last, the next vertex is 0; a nonempty subset of {1,...,n-3} goes before v, and the rest follows 0, each in increasing order, giving 2^(n-3)-1 paths. These disjoint classes sum to (7). Thus H(Y)/H(X) tends to 3/2; no larger universal multiplicative constant is possible. The n=3 boundary is separately H(X)=3,H(Y)=5.

### Arbitrarily many induced additions

If X is an induced subtournament of Z, insert the other vertices in any order. Each step has nondecreasing H by the slot injection. If all steps have zero increment, H(Z)=H(X). Otherwise the first positive step starts at h=H(X), and (4) yields

    H(Z)=h  or  H(Z)>=(3h+1)/2.                         (8)

Thus an induced extension of any size cannot land strictly between h and this floor. Along any such construction ending at H_final, the number of positive-count births is at most log_(3/2)(H_final/H_initial). Count-preserving births are not bounded.

This is a universal induced-growth obstruction, unlike the family-specific density barrier in Session213. It does not forbid noninduced surgeries, rebuilding from a different witness, or realizing the excluded intermediate integers by unrelated tournaments. It gives no ambient forbidden count and no eventual-completeness conclusion.

## 6. Exact increment levels 2,4,6,8

The bound (6), strong-component factorization and a small complete census classify four levels in all orders:

    H(Y)-H(Y-v)=2:     H(Y-v) in {1,3};
    H(Y)-H(Y-v)=4:     H(Y-v) in {1,5};
    H(Y)-H(Y-v)=6:     H(Y-v) in {3,5,9};
    H(Y)-H(Y-v)=8:     H(Y-v) in {1,5,9,11}.             (9)

The first level is inherited as an all-order result. The other three were finite observations in THM-4876/HYP-9469; the point here is the **all-order exhaustion** argument.

Let S be the strong component of Y containing v. If S is a singleton, deletion preserves H, contrary to the positive increment. All other strong components have a product g of their positive odd path counts. Since component order forces Hamilton-word concatenation,

    H(Y)=g H(S),      H(Y-v)=g H(S-v),
    k=g d,           2d=H(S)-H(S-v)>0.                 (10)

For k<=4, the positive odd factor g is either 1 or 3. Both are realizable, by an empty external product or one directed triangle. By (6), H(S)<=H(Y)<=23. Busch's published sharp strong-order bound implies every strong tournament of order at least seven has H>=25. Hence S has order at most six.

The verifier exhausts every labelled tournament on at most six vertices, retains every strong one and every deletion, and then applies exactly the factors g=1,3 compatible with k=gd. Its resulting sets are (9), and it emits a literal marked tournament witness for each value. Completeness is all-order because of (6), (10), and the cited strong-order bound, not because the first six orders happened to stabilize.

Dependency: Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1/Corollary 1, DOI [10.37236/1141](https://doi.org/10.37236/1141). The published bound and source are already pinned in `session197/unrestricted-target/PROOF.md`. Only its consequence H>=25 at strong order>=7 is used here. No new proof or priority claim about that bound is made.

## 7. Reproduction and limits

Run from this directory:

    python3 verify.py > /tmp/session214-shuffle-normal.json
    cmp RESULT.json /tmp/session214-shuffle-normal.json
    python3 -O verify.py > /tmp/session214-shuffle-optimized.json
    cmp RESULT.json /tmp/session214-shuffle-optimized.json

The verifier has explicit exception checks, not assertions. It tests the forced map for every nonzero cross-orientation matrix with 1<=|A|,|B|<=4, with both internal words fixed. This exactly isolates the arcs used by Lemma 1, rather than sampling unrelated tournament orientations. It separately checks every labelled tournament of orders one through six and every cut using endpoint subset dynamic programming; literal permutations check H through five, and explicit pair-image collision checks run through four. A fresh direct subset count of every extension through child order six checks the growth inequalities. The strong-core deletion census and external factors certify (9), with literal witnesses in the receipt.

Finite checks validate implementation and hostiles. Lemma 1 and Theorems 2–3 are all-order combinatorial proofs independent of those computations. The universal spectrum problem remains open: a successful all-odd construction must avoid the induced gap or arrange parent choices that respect it.
