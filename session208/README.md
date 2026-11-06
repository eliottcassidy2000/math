# Session 208: tree transport, native metagraph geometry and count-efficient atoms

This session connects two independently verified research directions to mechanisms in Mohammad Hajiaghayi's [*Logarithmic Embeddings into Probabilistic Spanning Trees for Distances and Capacities*, arXiv:2611.05189v1](https://arxiv.org/abs/2611.05189). The results below distinguish proofs and explicit certificates from conditional applications and unresolved optima. No literature-priority claim is made.

## What comes from the paper

The preprint's Theorem 1.1 states O(log N) expected distance stretch and O(log N) expected fundamental-cut congestion for actual spanning trees of a connected undirected N-vertex graph. It allows DIFFERENT distributions for these two objectives, with support at most the number of original edges. Its main geometric theorem is **cited, not independently certified here**.

The usable elementary mechanisms are exact edge/path incidence, canonical fundamental-cut capacities, weighted tree/edge prices and rational support elimination. Session 208 proves its concrete results independently. On a labelled order-n tournament metagraph, N can be as large as 2^binom(n,2), so a generic O(log N) consequence is O(n²), not O(log n). Work is componentwise unless connectivity is proved.

## Native metagraph: exact structure and certified bounds

The native move preserves a directed Hamilton path and reverses every other tournament arc. It is NOT a single-arc reversal. For n>=3, its labelled degree is the ordinary Hamilton-path count H; class-line degrees involve automorphism factors and must be kept separate.

- [Native tree transport](native-tree-physics/PROOF.md): fundamental capacities dominate cuts and scalar total variation; this does not imply preservation of electrical energy or original adjacency degrees. Bridge and transitive-leaf channels are exact even when their opposite endpoint has a very large H.
- [Bridge factorization](bridge-factorization/PROOF.md): optimal stretch and congestion separately reduce exactly to the bridge-deleted components. Local distributions with support sizes s_i couple using at most 1+sum_i(s_i−1) trees. Logarithmic corollaries remain conditional on the cited paper; the factorization is unconditional.
- [Order-four geometry](core-geometry/PROOF.md): the actual connected 64-vertex labelled graph has 24 transitive leaves; its 40-vertex, 72-edge core consists of two subdivided cubes joined through six squares. The detailed certificates retain physical edges and their multiplicities.
- [Sharp parity invariant](native-capacity-parity/PROOF.md): canonical weighted degree is H plus twice internal transit. Odd-load tree edges form the unique V-join of the tree. For the order-four graph, the EXACT minimum number of odd-load edges is 40, and its EXACT maximum matching size is 24.
- [Strongest certified tree bounds](square-bounds/PROOF.md), for original unit lengths and capacities:

      109/40 <= optimal maximum expected stretch <= 211/66,
      11/4   <= optimal maximum expected congestion <= 269/84.

  These are bounds, NOT exact optima. The lower bounds use weighted adversaries and square boundary data; the upper bounds use explicit rational symmetry mixtures. Distance and capacity certificates use separate distributions. [Earlier refinements](core-improved/PROOF.md) remain valid but are superseded by these intervals.

## Count spectrum: atom mixtures and stronger constructions

Write m(h) for the smallest tournament order realizing exactly h directed Hamilton vertex words. Dominance joins multiply H and add orders. The inherited SCC/min-plus reduction makes prime-exponent mixtures a natural counterpart to the paper's rational support elimination.

- [Rational atom mixtures](rational-atom-mixtures/PROOF.md): for a fixed finite strong-atom dictionary, an optimal rational mixture uses at most the prime-exponent rank in atom types. For powers h^k, its integer optimum differs from k times its LP rate by a bounded additive amount and agrees on denominator multiples. The complete strong-order-at-most-eight dictionary has EXACT rate 40/11 for powers of 15, certified by rational dual prices and the construction 375³·405²=15¹¹ on 40 vertices. This is dictionary optimality, not global optimality. The global asymptotic rate is the infimum over all such finite dictionaries; attainment is not proved.
- [Higher atoms and a lower-bound boundary](higher-atom-constructions/PROOF.md): explicit larger atoms improve the construction rate to 8/3 and refute a global extension of the order-eight dual certificate. The general bound H(T)<=n!/2^floor(n/2) gives a computable lower bound for m(15^k), in particular (log15−o(1))k/log k. This does NOT prove a positive linear asymptotic rate.
- [Exact arc-flip oracle and order 13 atom](balanced-atom-oracle/PROOF.md): complementary-subset endpoint convolution computes every single-arc count change in O(n²2^n) integer arithmetic operations. A verified strong 13-vertex witness has H=15⁵, giving rate 13/5. The discovery graph is distinct from the native metagraph.
- [Strong order 14 atom](balanced-order14/PROOF.md): an explicit independently checked tournament has H=11390625=15⁶. For k=6q+s, 0<=s<6,

      m(15^k) <= 14q+5s,
      lim_k m(15^k)/k <= 7/3.

  This is the strongest construction upper bound in these Session 208 packets. It supersedes 40/11, 8/3 and 13/5 as a global upper bound, without invalidating their witnesses or finite-dictionary results. Neither m(15⁶)=14 nor global optimality of 7/3 is claimed.

## What remains open

No eventual threshold for realizing all sufficiently large admissible odd counts was established. No complete all-order Hamilton-count spectrum, global minimum-order formula, positive linear efficiency lower bound, or exact order-four tree optimum was established. Failed capped searches are not exclusion proofs. The packets include explicit witnesses, standalone verifiers and independent audits; see each README for reproduction.
