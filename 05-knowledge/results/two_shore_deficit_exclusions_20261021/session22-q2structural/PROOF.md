# Two-shore insertion obstruction: Delta >= 2m-2

Status: PROVED; independent analytic and finite-control review passed. Root supplied the insertion construction; this packet independently checks and expands its proof. Earlier degree47 and order10 computational packets remain unchanged.

## Statement

Let T be a strong self-converse tournament on2m vertices, where m>=5 is odd. Suppose a marked Hamilton cycle has exact cut index q=2, so its step-two rotation is an automorphism. Then

    Delta(T)=h(T)-b(T) >= 2m-2.

In particular no deficit-six center has q=2 with m>=5. Combined with q<=8 and odd n/q, this proves that both order10 and order14 have no deficit-six center. The earlier128-frame order10 computation is now an independent finite control, not the main reason for its zero coefficient.

## 1. Intrinsic shores and punctured strongness

Let A,B be the two rotation shores, each of size m. Their induced tournaments are regular and strong. If each A vertex beats k vertices of B, every B vertex beats m-k vertices of A. The full scores on the two shores are

    s_A=(m-1)/2+k,   s_B=(m-1)/2+m-k.

They differ since m is odd, and their sum is2m-1. Hence automorphisms preserve A and B individually, while any anti-isomorphism interchanges them. Interchange the names A,B if necessary so k>m/2. The marked alternating Hamilton cycle supplies cross arcs in both directions. Transitivity within each shore therefore implies

    (m+1)/2 <= k <= m-1.

Every vertex of B has at least one outneighbor and at least one inneighbor in A, and every A vertex has at least k>=3 outneighbors in B.

Deleting one vertex from a regular tournament of odd order m>=5 leaves a strong tournament. Here is a direct proof: if the remainder admits a nonempty split X->Y, regularity and average outdegree comparison force |X|<=1 and |Y|<=1. Then m=|X|+|Y|+1<=3, a contradiction. This is also the special case kappa>=ceil(m/3)>=2 of the inherited regular connectivity lemma.

We use strong-tournament Hamiltonicity, inherited from the previously primary-source-verified Moon theorem. Cutting a Hamilton cycle at any chosen vertex supplies a Hamilton word starting there.

## 2. Source-bound insertion words

Fix v in B and a directed Hamilton cycle C of T[A]. Because v has neighbors in both directions on A, traversing C finds a transition

    a -> v -> a',  with a -> a' consecutive on C.

For each integer i with1<=i<=m-2, cut C so that a is its ith vertex. Insert v between a and a'. This gives a directed word through A union{v}, whose score-shore pattern is

    A^i B A^(m-i).

Let z be its last A vertex. Since z has k>=3 outneighbors in B, at least k-1>=2 of them lie in B\{v}. Choose one, w. The punctured tournament T[B\{v}] is strong, so it has a Hamilton word starting at w. Append that word. The resulting full Hamilton word P_i has pattern

    P_i: A^i B A^(m-i) B^(m-1).                   (1)

Also choose an A vertex u beaten by v. Cut C immediately before u; prepend v. Append a Hamilton word of B\{v} beginning at an outneighbor of the resulting last A vertex, as above. This gives P_0 with pattern

    P_0: B A^m B^(m-1).                           (2)

No common terminal or suffix is assumed. Every appended path is selected after inspecting its literal final A vertex. These are actual directed Hamilton words, not merely feasible binary patterns.

## 3. Orbit separation and mate counting

Every automorphism preserves shores, so different binary shore words certify different Hamilton-word orbits. They do not claim that equal binary patterns imply equivalence.

Fix any anti-isomorphism tau of T. The converse/reversal involution on Hamilton-word orbits sends a word to tau applied to its reversal. On binary shore words this is reverse followed by swapping A,B. Its fixed orbits are precisely the blue orbits. The mates of (1),(2) have patterns

    J(P_i): A^(m-1) B^(m-i) A B^i,  1<=i<=m-2;
    J(P_0): A^(m-1) B^m A.

All m-1 original patterns are distinct: P_0 starts in B; P_i has initial A-run of exactly i. Every mate starts with an A-run of exactly m-1, unlike every original. The mates are pairwise distinct because reversal-and-swap is a bijection on words. Thus all2(m-1) orbits are distinct and none is fixed by J. All are black, proving Delta>=2m-2.

The omitted index i=m-1 is a useful hostile: its pattern A^(m-1) B A B^(m-1) is fixed by reverse-and-swap, so it cannot be credited automatically as another black orbit. The m=3 boundary is also essential: puncturing a regular three-shore leaves a nonstrong two-vertex tournament. The known order6 deficit-six center is not excluded by this argument.

## 4. Consequences and scope

At orders10 and14, q<=Delta+2=8 and n/q odd force q=2, with m=5 and7 respectively. (The other odd-quotient possibility q=n exceeds8.) Therefore both coefficients vanish analytically:

    [z^10]R6 = [z^14]R6 = 0.

This proof is a structural family exclusion, not a full R6 classification. It does not change the degree47 cutoff, settle other residual orders, or reconstruct metagraph adjacency. Literal source/word data remain in the new checker receipts; the earlier order10 deficit histogram is independent supporting computation rather than an ingredient.

## 5. Fresh literal construction controls

`check.py` generates all128 marked q=2 frames with m=5 and all2048 with m=7, independently of whether the frame is self-converse. It explicitly constructs each original Hamilton word, checks every arc and vertex, verifies its stated binary pattern, and checks disjointness from the reverse-complement patterns. The latter become actual mate words only in the self-converse case required by the theorem; no anti-isomorphism is fabricated for non-self-converse controls.

The checker retains twelve literal sample receipts (one for each first-shore orientation at each order), and a SHA256 digest over the canonical JSON of all2176 construction receipts. Each source contains its adjacency rows, high-cross shore, other shore, chosen terminal, original words and pattern/mate-pattern list. The omitted i=m-1 pattern is verified fixed under reverse-complement as a hostile against overcounting.

Run `python3 -B check.py` and `python3 -B -O check.py` in this directory. Both fresh runs perform21,504 always-active checks and have equal non-timing receipt content. Actual UTC start/end and elapsed clocks are retained. This is a complete finite marked-frame construction test at m=5,7, not a class census or an all-m computational proof.

Inherited sources: `effective_black_defect_order_bound_20261020/PROOF.txt` for the cut-index bound and regular connectivity; `metagraph_classification_20261021a/README.md` for the intrinsic two-shore score-word method. Strong Hamiltonicity is the Moon theorem already independently verified during the degree47 proof audit. The all-m insertion and mate-separation argument above, rather than finite computation, proves the new family exclusion.
