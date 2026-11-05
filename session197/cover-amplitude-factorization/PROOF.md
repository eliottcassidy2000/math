# Exact minimum-cover factorization of a clone clock's leading amplitude

Session 197. All-order combinatorial theorem, with complete small-source and larger-rank controls. This packet concerns actual rooted tournaments. It is not a realization theorem for arbitrary numerical factors, and it does not prove a unit-offset congruence or the {7,21} spectrum conjecture.

## 1. Root runs and the crossing graph

Let Q be a nonempty tournament with distinguished vertex r. Put O=N^+(r), I=N^-(r), and form the bipartite graph G on O and I whose edges are the directed arcs O->I. Let nu be its maximum matching size.

A root-run word uses each vertex of Q-r exactly once and uses the symbol r exactly k times, with no adjacent occurrences of r. Every consecutive pair must be an arc of Q. Let a_k count these words. The two exterior end fragments may be empty; every fragment between two occurrences of r is nonempty, starts in O, and ends in I.

Replacing r by a transitive tournament of positive order m gives an actual tournament Q_m. Each Hamilton word of Q_m collapses uniquely to a root-run word. For a fixed k-run word, the m distinct clone labels are assigned to k ordered nonempty sets, and within each set appear in transitive order. Hence

    H(Q_m) = sum_(k>=1) a_k k! S(m,k),    m>=1,          (1)

where S(m,k) is the Stirling number of the second kind. In particular, the largest exponential base is the largest k with a_k>0, and its coefficient is a_k, since

    k! S(m,k) = sum_(j=0)^k (-1)^(k-j) binom(k,j) j^m.

There is no assertion here about substituting m=0 into a positive-length exponential formula. Root deletion is a separate boundary value.

Every internal root-to-root fragment contains an O-to-I arc. Selecting one from each fragment produces a matching, since distinct fragments are vertex-disjoint. Thus k-1<=nu and a_k=0 for k>nu+1.

## 2. Two physical half-cover responses

Choose any minimum vertex cover C of G, of size nu. Such a cover exists by the bipartite matching/cover theorem: from a maximum matching, follow alternating paths from unmatched O vertices and take unreached O together with reached I. The usual alternating-path count shows that this set covers every edge and has size nu.

Write C_O=C intersect O, C_I=C intersect I, O'=O minus C and I'=I minus C. Since C covers G, every arc between I' and O' points from I' to O'.

Define L_C as the number of the following collections of directed vertex words, together covering O' exactly once:

- One path P_c for each labelled c in C.
- One additional final path F.
- P_c must be nonempty when c is in C_I, and may be empty when c is in C_O.
- A nonempty P_c must end at a vertex pointing to c.
- F is an arbitrary directed path, possibly empty.

Define R_C dually as the number of collections covering I' exactly once:

- One path S_c for each labelled c in C.
- One additional initial path E.
- S_c must be nonempty when c is in C_O, and may be empty when c is in C_I.
- A nonempty S_c must start at a vertex to which c points.
- E is an arbitrary directed path, possibly empty.

All paths within a collection are vertex-disjoint. Path order matters. Each empty path has one choice. The cover vertices are labels for slots and do not themselves belong to the half-cover paths.

**Theorem.** For every choice of minimum cover,

    a_(nu+1) = nu! L_C R_C > 0.                        (2)

Consequently the exact leading amplitude of H(Q_m), at base nu+1, has this factorization. In particular L_C R_C is independent of the chosen minimum cover, even when its two factors change.

## 3. Exhaustion of the cover and the bijection

Consider a word with nu+1 occurrences of r. There are nu disjoint internal O-to-I fragments. Each fragment contains an O-to-I arc, and hence at least one vertex of C. Since |C|=nu, every internal fragment contains exactly one cover vertex, and every cover vertex is used. In particular no cover vertex lies in an outer fragment.

The internal fragment containing c has the form

    P_c, c, S_c.

Indeed, before c it cannot contain any I' vertex: starting in O, its first O-to-I transition would then meet a cover vertex other than c. Likewise after c it cannot contain an O' vertex: reaching its final I vertex would require another O-to-I transition, again meeting a different cover vertex. Thus the prefix lies entirely in O' and the suffix entirely in I'. The fragment begins in O and ends in I, giving exactly the nonemptiness and incidence conditions in Section 2.

The initial outer fragment lies entirely in I'. If it contained O', its path to its last I vertex would need an uncovered O-to-I edge. Similarly the final outer fragment lies entirely in O'. These assertions also permit empty outer fragments.

It follows that each maximal-run word determines one left collection, one right collection, and an order of the nu distinct cover-labelled internal fragments. Conversely any such two collections and any of the nu! orders produce the valid word

    E, r, P_(c1), c1, S_(c1), r, ...,
          P_(cnu), cnu, S_(cnu), r, F.

No additional compatibility is needed between different internal fragments: they are separated by r. This construction and extraction are inverse, proving (2) without quotienting by graph automorphisms or identifying different path words.

Positivity can also be seen directly. A maximum matching has nu disjoint edges and meets the size-nu cover. Every matching edge must meet exactly one cover vertex; otherwise nu disjoint edges would require more than nu cover vertices. Hence it matches every vertex of C_I to a distinct O' vertex, and every vertex of C_O to a distinct I' vertex. Use these as the mandatory singleton prefix/suffix paths, leave optional paths empty, and choose Hamilton paths on the unused O' and I' vertices for F and E. Every tournament has a Hamilton path, including the empty-path convention. Thus both factors are positive.

## 4. Boundaries, symmetry and scope

For nu=0, C is empty and Q[I]->r->Q[O]. The two factors are L_C=H(Q[O]) and R_C=H(Q[I]), with H(empty)=1. Formula (2) gives a_1=H(Q[I])H(Q[O]), and (1) is constant in m. In particular a singleton Q has a_1=1. Empty sides and singleton components require no separate limiting argument.

At nu=1 the leading amplitude is a product of the two half-cover responses. At nu=2 it is twice their product. If a root-fixed anti-involution preserves a chosen minimum cover setwise, reversal and the anti-involution identify the two collections, so L_C=R_C and a_(nu+1)=nu! L_C^2. The anti-involution swaps O and I, so such a setwise invariant cover has even size. Not every minimum cover is preserved in this way; this observation alone does not replace the separate full rank-two symmetric-support classification.

Formula (2) implies the previously known factorial divisibility, but retains exact positive physical factors. It does not imply that every pair of positive integers occurs as (L_C,R_C), and says nothing by itself about lower exponential coefficients or the constant term. In particular it does not prove any proposed unit-offset mod-8 rule.

A separately frozen and independently audited polynomial strengthening is

    p_nu(x)=x(x+1)...(x+nu-1) L_C(x) R_C(x),

where each unused half is weighted by its labelled-bin polynomial instead of a single Hamilton path. See `session197-odd-cycle-path-polynomials/MINIMUM_COVER_POLYNOMIAL.md`, SHA256 `b03c91bdade816d3b62d5464364bfb6d6fa0913e53a199273eb3dbc3f75bb8a8`, and independent review `session197-bridges/minimum-cover-audit/REPORT.md`, SHA256 `5267ec64e06b78cc5798219601ac46c23b2fe442ad7451d8ef3f7301ee74812d`. The elementary scalar proof above does not require that operator formalism or its code.

## 5. Standalone verification

`verify.cpp` has no source-code or data-bank dependency. It performs two different counts for every tested root and every minimum cover:

1. An exterior-permutation subset DP inserts permitted root gaps and obtains the entire root-run coefficient vector.
2. Separate endpoint path DPs on O' and reversed I', followed by disjoint-subset convolution over the labelled cover slots, count L_C and R_C directly.

It checks maximum-run positivity and cutoff, (2), existence of a minimum cover, and independence of the product from the cover. All minimum covers are enumerated, not just an alternating-path choice.

The tested universes are:

- Every unlabelled tournament through order 7, regenerated by all one-vertex extensions and exact degree-class permutation canonicalization: counts 1,1,2,4,12,56,456; 532 graphs, 3,613 rooted tests, 6,919 minimum covers.
- Every labelled tournament through order 5: 1,099 graphs, 5,405 rooted tests, 8,747 minimum covers.
- Twelve seeded graphs at each order 8 through 13: 72 graphs, 756 rooted tests, 1,336 minimum covers, including matching ranks 0 through 6.

Total: 9,774 rooted tests and 17,002 cover factorizations. The overlap between the two exhaustive universes is intentional. These finite checks supplement the all-order bijection; they are not its proof.

Canonicalization is exact because every isomorphism preserves outdegrees, and all permutations within degree classes are examined. Every tournament has a vertex deletion in the preceding complete bank, so the extension recursion is complete. The labelled universe is an independent enumeration of every orientation of every pair. Seeded larger controls are explicitly not a complete larger-order census.

At the tested orders, all counts fit unsigned 64-bit integers: root-gap words are bounded by (n-1)! 2^n, and intermediate half-cover counts are bounded by q! binom(q+nu,nu), with q+nu<=12. Explicit checks throw exceptions rather than use assertions. Run `python3 verify.py` and `python3 -O verify.py`, or compile the C++ file directly with C++17. Optimized and debug C++ outputs agree exactly.
