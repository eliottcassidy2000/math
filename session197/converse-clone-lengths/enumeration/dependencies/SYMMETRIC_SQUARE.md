# Marked self-converse rank-two sources: square top coefficients and a new amplitude gap

## Status and scope

This is a separate Session 197 packet. It does not alter the frozen bounded-slack, mirror-star, or symmetric-amplitude-78 packets. Its all-order structural theorem is:

**Theorem 1.** Let a tournament T have an anti-involution tau that fixes the marked vertex r. Suppose the bipartite graph of arcs N+(r) -> N-(r) has matching number two. If a_3 counts root-run words with exactly three occurrences of r, then

    a_3 = 2 q^2  or  a_3 = 4 q^2

for a positive integer q. The proof below gives q explicitly as a path-cover count, not as a square root inferred from examples.

**Theorem 2.** If transitive clones of such a marked self-converse source have H_m = A 3^m - 1, for every positive m or merely for all sufficiently large m, then A >= 242. No attainment at 242 is claimed.

The marked hypothesis is essential. A self-converse tournament with a root not fixed by its anti-involution need not have a square top coefficient. Section 7 gives an odd-order counterexample with an anti-involution-fixed exterior vertex.

Named dependencies are the frozen `MINIMUM_TERNARY.md` (root-run transfer, maximum run degree, root-SCC reduction), `SYMMETRIC_TERNARY.md` (the marked symmetric lower bound 78 and complete order-at-most-nine exclusion), and the peer `session197-orbit-wildcard/EVEN_CONVERSE.md` (even-order self-converse H = 1 mod 4). Busch's bound is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1, DOI 10.37236/1141, official PDF https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf. Only its safe consequence 5^(n-1) <= H^3 is needed here. The square theorem itself needs none of these external results.

## 1. Root-run words and anti-involution normalization

A root-run word uses every vertex other than r exactly once and uses r exactly k times. Consecutive symbols must be tournament arcs, and consecutive copies of r are forbidden. Let a_k be the number of these words. A word with three roots consists of an optional outer prefix, r, a nonempty internal fragment, r, a nonempty internal fragment, r, and an optional outer suffix.

Write O = N+(r), I = N-(r). Each internal fragment starts in O and ends in I. The outer prefix, when nonempty, ends in I; the outer suffix, when nonempty, starts in O.

An anti-involution reverses every arc. It cannot fix two distinct vertices: their connecting arc would have to reverse while its endpoints stayed fixed. Thus, when tau fixes r, **there is no fixed exterior vertex**. It pairs O with I, so |O| = |I| = t and |T| = 2t+1. Write

    O = {u_1,...,u_t},  I = {v_1,...,v_t},  v_i = tau(u_i).

The tournament induced on I is the converse of that induced on O. Define the zero-one matrix M by M_ij = 1 exactly when u_i -> v_j. Anti-invariance gives M_ij = M_ji, including the diagonal entries. Diagonal entries represent arcs between paired vertices; they are allowed and will be treated explicitly.

There is no rank-two case at t = 0 or 1. In particular the one-vertex boundary is vacuous, not an exception to the theorem.

## 2. Complete symmetric matching-support classification

Suppose the bipartite matching number of M is two. View its off-diagonal support as an undirected simple graph on the indices, retaining the diagonal support as loops.

Two vertex-disjoint off-diagonal edges give a bipartite matching of size four, by using each in both directions. An off-diagonal triangle gives a matching of size three. Consequently any two off-diagonal edges intersect and there is no triangle. Such an edge family is a star: after two edges {c,a},{c,b}, any edge not through c must be {a,b}, which would form a triangle.

If the star has at least two leaves, a loop at a leaf a creates the matching (a,a),(c,b),(b,c), where b is a different leaf. A loop at an index outside the star combines with an edge and its reversal to give a matching of size three. Hence only a center loop is permitted.

If there is exactly one off-diagonal edge {a,b}, loops outside {a,b} are forbidden by the same argument. No loops or one endpoint loop give a star with an optional center loop. Two endpoint loops give a full two-by-two block of ones.

If there are no off-diagonal edges, exactly two loops are necessary and sufficient for matching number two.

Thus the list is exhaustive:

1. A nonempty star, with or without its center loop.
2. Two isolated loops.
3. A full two-by-two block of ones.

All other matrix entries are zero. Every listed pattern has matching number two. Self-paired crossing arcs are therefore already included: a star center loop, either isolated loop, or the two diagonal arcs of the full block.

For t labelled indices the number of matrices is

    5 binom(t,2) + 2t (2^(t-1) - t).

The first term accounts for the four single-edge diagonal choices and the diagonal-only pair. The second term accounts for a star with at least two leaves and an optional center loop. For t = 2,3,4,5 the totals are 5,21,62,160.

## 3. The star case: a_3 = 2 q^2

Let c be the center, let L be its nonempty set of leaf indices, and put u = u_c, v = v_c. Delete u and v, leaving X = O - {u} and Y = I - {v}. Every arc between these residual sets points from Y to X.

Each internal root-to-root fragment must cross from O to I. Every available crossing uses u or v. Since the two fragments are disjoint, one contains u but not v and the other contains v but not u. The optional arc u -> v cannot occur in either fragment, because that would leave no crossing available for the other fragment. Thus a self-paired center crossing is not being silently discarded; it is unusable in a three-root word.

The fragment containing u has the exact form

    A, u, D,

where A is a possibly empty X-path ending at an in-neighbor of u, and D is a nonempty Y-path whose first vertex is v_j for a leaf j in L. The fragment containing v has the exact form

    B, v, E,

where B is a nonempty X-path ending at u_j for a leaf j in L, and E is a possibly empty Y-path whose first vertex is an out-neighbor of v in Y.

There are no additional switches between X and Y inside these pieces. Once a path enters Y it cannot return to Y after passing to X without another unavailable crossing through u or v. Before the prescribed crossing it cannot enter Y at all. These observations force precisely the displayed forms.

The outer prefix uses only Y: it contains neither u nor v, and an X vertex could not be followed eventually by an I endpoint. The outer suffix similarly uses only X.

Define q to be the number of ordered triples (A,B,C) of vertex-disjoint directed paths partitioning X, where:

- A may be empty; if nonempty its last vertex points to u;
- B is nonempty and ends at a leaf u_j, j in L;
- C is arbitrary and may be empty.

Here C is the outer suffix. On the Y side, choose an outer prefix P and the two tails D,E partitioning Y. The map that applies tau and reverses each path gives a bijection

    (A,B,C) <-> (E,D,P).

Indeed A ending at x with x -> u corresponds to E beginning at tau(x) with v -> tau(x). B ending at u_j corresponds to D beginning at v_j with u -> v_j. C has no endpoint restriction. This is an actual bijection of complete half-covers, including empty pieces.

The two half-cover choices are independent: all joining requirements have just been listed, and the root separates the two internal fragments. There are q choices on each side and two possible orders of the distinct fragments. Therefore a_3 = 2 q^2. Positivity is explicit: take A empty, B a single leaf, and C any Hamilton path on the remaining X vertices.

## 4. Two loops and the full block

Let the active indices be a and b. Delete u_a,u_b,v_a,v_b, leaving residual sets X and Y = tau(X). Again Y dominates X.

For two isolated loops, the two fragments must use respectively u_a -> v_a and u_b -> v_b. For a full block, they must use a perfect matching between {u_a,u_b} and {v_a,v_b}; there are exactly two choices. Every active vertex must be used in the matching fragments. Otherwise two disjoint O-to-I crossings are impossible.

In either case each fragment has the exact form

    X-prefix, u_i, v_j, Y-tail.

The displayed crossing vertices are adjacent. Inserting an X vertex after u_i would make the required crossing unavailable; inserting a Y vertex before v_j would already require another active crossing vertex. The residual Y -> X direction also rules out additional switches. The outer prefix is entirely Y and the outer suffix entirely X.

Define q to be the number of triples (A_a,A_b,C) partitioning X into directed paths. Either A_i may be empty; if nonempty its last vertex points to u_i. C is arbitrary, possibly empty. Applying tau and reversing paths gives a bijection to the Y-side choices: a tail beginning after v_i must start at tau(x) with x -> u_i, and the outer prefix is the image of C.

For each chosen matching, the X and Y half-covers can be chosen independently. The matching only assigns which Y tail is attached after a particular u_i; it does not change either half-cover condition. This remains true for the self-paired matching i = j, so the fixed crossing support is explicitly covered.

There are two fragment orders. Hence two isolated loops give a_3 = 2 q^2; the full block, with its additional two choices of matching, gives a_3 = 4 q^2. Taking both A_i empty and C any Hamilton path proves q > 0. This finishes Theorem 1.

## 5. From the square theorem to an amplitude gap

For matching number two, the maximum root-run number is three, and the exact transitive-cloning formula is

    H_m = a_1 + a_2(2^m-2) + a_3(3^m-3*2^m+3).

Thus an eventual identity H_m = A 3^m - 1 forces, by independence of the three exponential bases,

    a_3 = A,  a_2 = 3A,  a_1 = 3A-1.

Conversely those equations imply the identity for every positive m. A source with a larger positive top run coefficient cannot satisfy an eventual ternary identity, so this includes sources not initially assumed to have matching number two.

Cloning a tau-fixed root an even number of times preserves self-converseness: extend tau by reversing the transitive clone chain. The clone graph has even order, hence H_2 = 1 mod 4. Since H_2 = 9A-1, this gives A = 2 mod 4. Theorem 1 now forces

    A = 2 q^2 with q odd.

In particular the full-block alternative is impossible for these negative-unit clocks.

The frozen symmetric packet already proves A >= 78. The only values 2 q^2 with q odd in [78,242) are 98 and 162. Remove all SCCs outside the root SCC. Their Hamilton-path product divides both 3A-1 and 9A-1, whose odd gcd is one, so these SCCs are singletons and contribute one. The root SCC is invariant under tau and retains its fixed root, hence has odd order. Its H is at most 3*162-1 = 485. Busch's bound and

    5^12 > 485^3

exclude order at least 13. Thus the strong root core has odd order at most 11. Orders at most nine are excluded by the frozen complete symmetric search. Section 6 excludes order eleven for amplitudes 98 and 162. Therefore A >= 242. This proves Theorem 2, including eventual identities and arbitrary singleton padding.

## 6. Complete order-eleven source check

Normalize the fixed-root source as in Section 1, with five O vertices and their five tau-images. There are exactly 2^10 = 1024 tournaments on the labelled O set. The I tournament is forced. The exhaustive support classification gives 160 possible symmetric rank-two crossing matrices. Hence the physical bank has 163840 sources and covers every order-eleven marked anti-involutive rank-two source, without relying on isomorphism generation or a Hamilton-path normal form.

The half-cover q is evaluated by subset path counts on O. For a star, partition X into A,B,C and sum

    F(A, in_O(u), empty allowed) * F(B, leaves, empty forbidden) * H(C).

For the other cases replace the two endpoint sets by the in-neighbors of u_a and u_b and allow both prefixes empty. F counts paths with endpoint in the indicated set, and H(empty)=1. These are precisely the bijections proved above.

Amplitudes 98 and 162 require q=7 and q=9, respectively, and cannot arise from full blocks. There are exactly 2640 q=7 sources and 5600 q=9 sources. An exact full-graph Hamilton-path dynamic program counts all 8240, without numerical approximation or a Hamilton-count cap. None has H_1 = 6q^2-1. Thus none even meets the first necessary clock equation.

Finite controls are deliberately distinct from this completeness argument. A separate Python augmenting-matching census enumerates every symmetric binary matrix through t=5 and checks the classified support list. The C++ checker tests every physical rank-two source through order nine (4146 sources), comparing the half-cover square formula against a direct word DP with repeated root symbols. It also counts its one-, two-, and three-clone graphs exactly. The direct word DP does not infer its a_3 from the square formula or from clone differences.

The checker emits a deterministic JSON receipt. Normal and Python-optimized runs must agree. These finite checks support, but do not replace, the all-order combinatorial proof.

## 7. Necessary hypothesis and a different positive-offset family

Here is an odd-order hostile to dropping the marked-root condition:

    rows = [2,4,9,19,39,79,31],  root = 2.

A row bit j means an arc to vertex j. The unique anti-involution is reversal [6,5,4,3,2,1,0]. It fixes exterior vertex 3 and sends root 2 to 4. Literal root-run word enumeration gives

    (a_1,a_2,a_3) = (31,30,6).

The first three clone counts are 31,91,247. Since 6 is neither twice nor four times an integer square, unmarked self-converseness is insufficient. This example is not a counterexample to Theorem 1: its unique anti-involution does not fix the root.

A second scope control shows that the earlier ordered mirror-star exclusion of both unit signs must not be extended to all marked self-converse sources:

    rows = [414,260,8,50,294,327,415,318,12],  root = 4.

Reversal fixes the root and reverses every arc. Its run profile is (217,216,72), with first clone counts 217,649,1945,5833. It gives H_m = 72*3^m+1. Its top coefficient is 72 = 2*6^2, entirely consistent with Theorem 1; the odd-q restriction applied only to the negative-unit sign. Sharp positive-unit minimality and arithmetic claims belong to the peer's separate positive-clock packet, not to this proof.
