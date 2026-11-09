# Inheritance, concepts and source boundaries

## Portfolio and inherited mechanisms

- **Anchor:** test the proposed diagonal/absent-pair implication for ordinary
  Hamiltonian-path totals seven and twenty-one.
- **Niche:** retain reciprocal pairs and endpoint response instead of silently
  replacing a semicomplete digraph by a tournament.
- **Wildcard:** use the unique ordered SCC decomposition as a formal language;
  this produced the general finite-core and rational-generating-function theorem.

Closest proved mechanism: the tournament gaps and strong-order floors in
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Canonical hostile: a triangle with one doubled edge, which remains strong
and keeps exactly three paths because the extra arc is unused.
Corrected near miss: identifying the sole three-path factor with a triangle
omits this second core and undercounts all mixed twenty-one families.
Least-used sidecar: endpoint-resolved path counts, which give the exact
nonnegative arc-addition response and reveal interactions lost by scalar H.

The initial triangle-only inference was corrected before publication; both
producers and the independent auditor verified the repair. It is retained
as a hostile control rather than as a live theorem. No global numerical
theorem IDs were reserved or repurposed.

## Concept board and connection contracts

| source and operation | preserved target | information lost | retained sidecar / decisive check |
|---|---|---|---|
|binary matrix -> erase diagonal|H exactly|loop decorations and some automorphisms|actual loop-subset orbits|
|semicomplete D -> tournament resolution|each selected path remains valid; H is bounded above by H(D)|paths needing the opposite resolution|number of reciprocal pairs used by each path|
|D -> ordered strong components|H becomes the product of core H-values|vertex labels and within-core structure if only counts are stored|ordered core isomorphism types and automorphism groups|
|strong D -> spanning directed cycle|at least |V| distinct paths|all other paths|rotation injection gives a uniform core-size bound|
|fixed-H core words -> rational series|exact isomorphism and labelled counts|individual adjacency matrices|standalone core generator reconstructs them|
|arbitrary D -> union of path edges|the complete Hamiltonian path set|unused arcs|forbidden-hyperedge avoidance recovers exactly the legal additions|

Every useful bridge has a concrete map, not just a shared number. In
particular the strong-subgraph comparison transfers a lower bound, not
an equality or a tournament parity law to general semicomplete graphs.

## Sources and priority

1. **Classical, cited:** Moon, *On Subtournaments of a Tournament*, Canadian
   Mathematical Bulletin 9(3) (1966), 297–301, Theorem 2,
   [DOI](https://doi.org/10.4153/CMB-1966-038-7). The strong-subtournament
   input is inherited from the precisely named THM-1370 file; its finite
   path floors are independently reproduced here.
2. **Published spectrum conjecture:** Chan–Pak, *Computational complexity
   of counting coincidences*, Theoretical Computer Science 1015 (2024),
   114776, Section 6.5, Conjecture 6.2,
   [arXiv:2308.10214v2](https://arxiv.org/abs/2308.10214v2). It credits the
   [2016 MathOverflow discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament).
   Its range assertion remains a conjecture. Source interpretation is
   inherited from the repository's corrected spectrum note and was
   checked against the primary text during final review; no new proof of
   that conjecture is claimed.
3. **Parity context:** the Rédei–Berge theory is routed by
   [CORE-PAPERS](../../reference/CORE-PAPERS.md).
   Grinberg–Stanley is [arXiv:2307.05569](https://arxiv.org/abs/2307.05569);
   Irving–Omar is [arXiv:2412.10572](https://arxiv.org/abs/2412.10572).
   These author/identifier pairs are not interchangeable. No symmetric
   function formula is needed for the finite-core proof.

The all-order gap argument is inherited and credited, not a novelty claim.
The semicomplete classification is a reproducible result of this packet;
no priority claim about the broader graph-theory literature is made.
The [final review](final_review/REPORT.md) records direct source checks of
Moon's original Cambridge PDF and Chan–Pak's arXiv text.

## Reusable research move and boundary

Applied existing meta-patterns: **correct the object before sharpening the
technique**, **search the statement before the method**, and **separate
local support from bounded-height coverage**. The key trigger was a
universal claim whose negation quietly allowed a fourth pair relation.
The small hostile test exposed reciprocal arcs; retaining SCC words then
turned a finite census into a complete classification.

The next generalization is useful but scoped: strong semicomplete cores
have H at least their order, giving finite grammars for every fixed H.
Do not transfer the same core cutoff to arbitrary digraphs merely because
their adjacency matrices are also binary. The missing-pair sector uses
the separate exact path-support generator.

No new global meta-pattern card is promoted from this single research
thread. The concurrent tight-colouring automaton correction was preserved
when synchronizing main; no theorem here depends on that unrelated work.
