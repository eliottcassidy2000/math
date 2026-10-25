# Session 56: an intrinsic tournament LCA query interface

Status: PROVED elementary interface and cut identity, CITED canonical decomposition, FINITE-EXACT certificates through the stated bounds. No novelty claim. No LRC, Hamilton-count complexity, or A000568 speedup claim.

## Inheritance and scope

Source commit: `4f1760fd1eff346767c032524c932d5e317401d2`, recovered read-only from the current publication checkout. Exact inherited files are in `inputs/`, with SHA-256 in `source-hashes.txt`.

- THM-1960, `tournaments-compose-from-regular-seeds-the-spectral-substitution-law`: modular prime seeds differ from SCC/order-join atoms. This session does NOT rediscover modular decomposition.
- THM-1862, `order-join-reduction-principle`: SCC factors compose Hamilton-path counts multiplicatively.
- THM-2183, `order-join-is-an-exact-tournament-metric-product`: fixed matching factor sizes give an exact reversal-distance product on order-joins.
- THM-2195, `transitive-quotients-exactly-control-universal-substitution-products`: that universal product law fails for nontransitive quotients. Do not transfer an SCC metric assertion to a prime-node hierarchy.
- THM-3121, `path-cover-walk-content-substitution-kernel`: a general quotient requires the complete factor path-cover profiles, not just H. This is the relevant compositional sidecar.

Concept board: canonical strong modules; LCA/RMQ navigation; quotient ports; ordered SCC species; Hamilton/path-cover profile; signed cut flux. Anchor: meaningful intrinsic tournament hierarchy. Niche: exact species versus query representation. Wildcard: hostile fixed-Hamilton-order cut signatures.

## 1. Canonical object, not a height decoration

For a tournament T, a nonempty set M is a module if every outside vertex has the same orientation to every element of M. It is strong if it overlaps no module, where overlap means nonempty intersection and neither inclusion. Nonempty strong modules are laminar, contain all singletons and V(T), and hence form an inclusion tree D(T). An internal node's children partition it and are modules. Its quotient Q_M has one vertex per child and inherits the uniform orientation between children.

CITED structural input: Q_M is transitive or prime; a prime here has at least three vertices and no nontrivial module. The tree has at most 2n-1 nodes. Primary pins:

1. McConnell and de Montgolfier, *Linear-time modular decomposition of directed graphs*, §§2.3–2.4 and §4, pp. 6–7, 10–12 in the author PDF: https://www.cs.colostate.edu/pubserv/pubs/McConnell-rmm-digraphDecomp.pdf . Quotient representation and tournament algorithm; O(n+m), thus O(n²) on explicit tournaments.
2. Boussaïri, Chergui, Ille and Zaidi, *3-uniform hypergraphs: modular decomposition and realization by tournaments*, https://arxiv.org/pdf/1805.04607 , accessed 2026-10-25. In the retrieved 33-page version, Theorem 16, Theorem 19, Definition 20, Proposition 21 on printed p.6 (PDF index5) give the tournament decomposition and restriction property. Published record: https://doi.org/10.55016/ojs/cdm.v15i1.67935 . Numbering differs across versions: use this PDF's theorem numbers, not a search-snippet number.

No external algorithm is reimplemented here: `check.py` is deliberately exponential brute force for certificates.

## 2. PROVED exact arc-query interface

Let u≠v be labelled leaves, M=LCA_D(u,v), and C_u,C_v the distinct children of M containing them. Then

    u →_T v  iff  C_u →_(Q_M) C_v.

Proof: M contains both leaves; minimality of M forces distinct children. Uniform adjacency of distinct module children identifies their cross arc with any representative pair, including (u,v). This proves both directions. Reconstructing every pair recovers T exactly.

A sufficient lossless package is:

- the rooted tree, with original-vertex-to-leaf map;
- each child's port identity in its parent's quotient;
- for transitive nodes, the child total order;
- for prime nodes, the child-indexed quotient adjacency relation.

The tree need not be ordered at prime nodes, but an arbitrary indexing must be kept consistently with the quotient. These are named local ports, not decorative node labels.

Euler traversal gives an array of depths and first occurrences; an RMQ over the interval between first occurrences returns the LCA. Ties cause no problem: every minimum-depth occurrence in that interval represents the same LCA node. RMQ does not identify C_u,C_v by itself: retain a child-under-ancestor operation, e.g. binary lifting to depth(M)+1. With a sparse table for Euler RMQ and binary lifting, additional preprocessing is O(n log n) words/time and each arc query O(log n), after decomposition and quotient retention. Stronger known constant-time tree primitives could improve navigation; they are not implemented or claimed as a new algorithm here.

The quotient payload is at most binom(n,2) bits if all quotient pairs are stored: assign each child pair a representative leaf pair, whose LCA is that node, so these assigned pairs never coincide across nodes. If prime arity is at most w, prime payload is at most w(n-1)/2 bits because sum_internal(deg(M)-1)=n-1. Transitive child orders need O(n) machine words. No bounded-w assumption is justified for arbitrary tournaments.

What tree alone loses: a prime tournament has a star strong-module tree with all vertices as leaves. The three explicit n=5 examples in the certificate have that identical tree and H=9,13,15. They are nonisomorphic because H differs. Even retaining a 'prime' tag and child names cannot replace the quotient adjacency. The tree is a navigation skeleton, not a general compression of the tournament's hard part.

## 3. Ordered joins, SCCs, and exact generating functions

Every tournament has a unique ordered list of strongly connected components. The condensation is acyclic and a tournament, hence transitive. Every cross edge goes forward in this order. A directed Hamilton path visits the components consecutively, so H(T)=product H(C_i), recovering THM-1862 rather than claiming a new identity.

Let T be the species including the empty tournament; T_+=T-1; S the nonempty strongly connected tournaments (including the singleton); X the singleton; P the prime tournaments of size≥3. Canonical decomposition gives

    T = SEQ(S),       S = X + P ∘ T_+.

Proof of the second identity: a nontrivial strongly connected tournament cannot have a transitive top quotient, and therefore has a unique prime quotient on its maximal strong modules. Conversely substituting nonempty tournaments into a prime quotient of size≥3 yields a strongly connected tournament, with those blocks the canonical maximal strong modules. For the converse maximality, let a module N touch at least two substituted blocks. Its touched index set is a module of Q: every untouched block has a uniform orientation to N and therefore to those indices. Primeness forces N to touch every block. If N omitted a vertex of block i, uniformity of that vertex toward N would force i to beat all other quotient vertices or lose to all of them, impossible in a prime quotient of size≥3. Thus N is the whole tournament. Every proper module therefore lies in one block, making the blocks exactly the maximal proper modules and strong. This proves uniqueness of the inverse construction.

For the ordinary unlabeled series A(x)=sum a_n x^n and S(x)=sum s_n x^n, the SCC identity gives

    A(x)=1/(1-S(x)),      S(x)=1-1/A(x),
    s_n=a_n-sum_(1≤k<n) s_k a_(n-k).

The A000568 inputs through n=7 recover s_1..s_7=1,0,1,1,6,35,353. This is an arithmetic check of known input counts, not an independently generated n=7 census. For labelled counts the same SEQ identity holds between EGFs, yielding binomial factors in coefficient recurrence.

The prime species identity does NOT become ordinary substitution of unlabeled series. Its correct unlabeled form is

    S(x)=x+Z_P(A(x)-1,A(x²)-1,A(x³)-1,...).

For example the C3 contribution is ((A(x)-1)^3+2(A(x³)-1))/3, because its automorphism group is cyclic of order3. Thus local child names used by queries must be quotiented by prime automorphisms in unlabeled enumeration. Plain prime counts alone do not provide the necessary cycle-index data. This is the exact query/enumeration interface distinction, not a new enumeration shortcut.

For Hamilton counts at prime nodes, use THM-3121's complete path-cover-profile contraction with the quotient walk series, not a product of child H values. The LCA representation supplies the exact Q required by that existing kernel, but it does not make the contraction cheap. An arbitrary prime can remain as a root quotient of order n.

## 4. PROVED signed-cut identity and FINITE-EXACT hostile

Fix any order v_0,...,v_(n-1), Hamiltonian or not. Set

    F(k)=sum_(i<k≤j) (1[v_i→v_j]-1[v_j→v_i]),  0≤k≤n.

Internal pairs cancel in the sum of net outdegrees over the first k vertices, giving

    F(k)=sum_(i<k) (2d⁺(v_i)-(n-1)),
    F(k+1)-F(k)=2d⁺(v_k)-(n-1).

Therefore all cut fluxes with endpoint values F(0)=F(n)=0 carry exactly the ordered score vector, no more. A Cartesian/RMQ tree of these numbers is still weaker, or equal only if the numeric array is retained too.

Minimal H-separating witness within fixed Hamilton order 0→1→2→3→4: both tournaments have scores (3,2,1,2,2), cut vector (2,2,0,0), but H=11 and13. The full matrices, bit convention, leftmost-min Cartesian tree, and independent permutation check are in `certificate.json`. Equivalently start with all forward arcs and reverse:

- H=11 example: (0,3),(1,4),(2,4).
- H=13 example: (0,4),(1,3),(2,4).

All adjacent forward arcs remain. The two examples differ by reversing the directed 4-cycle 0→4→1→3→0, so each affected vertex preserves its degree. This explains the collision instead of merely displaying it.

The n≤4 fixed-Hamilton-order universes were exhausted (1,1,2,8 tournaments respectively), with no H collision on full cut flux. Thus n=5 is minimal for this precisely stated H-separation question. It is not a minimality assertion for arbitrary tree-shape collisions, arbitrary signatures, or an unnamed flux convention.

## 5. Reproducibility and boundaries

Run `python3 check.py` in this directory. It writes `certificate.json`; redirect stdout to `run.out` if desired. No third-party package needed.

- Exhaustive 1099 labelled tournaments at n=1..5: brute-force modules/strong modules, quotient transitive-or-prime classification, all leaf-pair LCA via Euler-depth RMQ against least-containing-set oracle, and every representative-quotient arc reconstruction.
- Explicit prime-star H examples and cut-collision H values: subset DP checked independently by permutation enumeration.
- Fixed Hamilton order cut search: exhausted n≤4, stopped on first n=5 H witness among64 possible inputs.
- A000568/SCC inverse: arithmetic consistency only.

No arbitrary height is assigned to a tournament to force an analogy. The hierarchy is canonical; RMQ accelerates a genuine query; quotient ports and path-cover profiles are the missing sidecars. The species formulas explain the cost of forgetting labels. Neither this packet nor its tests closes unrestricted A000568 enumeration, LRC(14), or general Hamilton-count computation.
