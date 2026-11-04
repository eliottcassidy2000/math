# Source lineage for Session 191 tournament growth

This note supplements the first checkpoint without changing its frozen
proofs or receipts. It distinguishes inherited ideas from this packet's
newly checked deliverables; it is not an external priority survey.

## Inherited calculus and operations

* Grinberg--Stanley's odd-cycle collection formula is the all-order input,
  credited in THM-002 and the CORE-PAPERS tournament section. The explicit
  cited consolidation is Irving--Omar, arXiv:2412.10572v3, Corollary 20;
  arXiv:2307.05569 belongs to Grinberg--Stanley.
* THM-070-claim-a-from-ocf.md already gives
  Delta(T,v)=2 sum_(odd C containing v) H(T-V(C)).
* THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md
  already proves the stronger presence-cube mixed-difference positivity
  behind marked-increment heredity. The Session 191 notes restate its simple
  marked specialization rather than claiming a new OCF identity.
* THM-026-deletion-ratio-formula.md already identifies the sum of vertex
  increments with weighted cycle-covered support. The bounds obtained from
  support sizes 3,...,n are immediate corollaries, not a separate new moment
  theorem.
* Moon 1966 Theorem 2 supplies at least two strong vertex-deletions. Its
  exact source and scope are recorded in
  05-knowledge/reference/CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md.
  The familiar odd-cycle lower bound also follows from Moon's vertex-
  pancyclicity; the present proof reconstructs it with rooted descent and
  Camion, without claiming novelty for that cycle statement.
* Moon's interval-chain product and the chain recognition/rigidity argument
  are inherited from THM-4744. The one-arc join and its rooted two-bin
  responses are inherited from Session29, ADDENDUM.md sections 2 and 5, as
  explicitly credited again in the Session190 spectrum source.

## Earlier increment investigations recovered during this session

At default-branch source snapshot
6c1f4697d7200b8d8788d8370df68f1307293129, the GitHub source tool returned

    04-computation/increment_gap_analysis.py
    Git blob 10c4efe180d3f7bf89b7b7913f26ca8bd483c993.

Its header credits kind-pasteur-2026-03-20-S2 and already frames the {7,21}
question in terms of unavailable vertex-addition increments. Its main
experiment groups increments by parent H at orders 5, 6, and 7. It also revisits
the marked OCF formula. An intermediate comment omitting 5-cycles is repaired
later in that same script; the current proof uses all odd cycles throughout.
Neither that script nor a search-result absence is used as proof of a new
all-order statement or of external priority.

THM-331-king-h-increment-bound.md proves a score-dependent lower bound for
a maximum-outdegree vertex from triangle coverage. It is distinct from
the complete rooted capped-increment atlas and is not a dependency here.

## Incoming Session190-era checkpoint

THM-4872 was read directly at pinned commit
3192fab75fd27f262157a67c9e2aa0c17c47a36b, file
01-canon/theorems/THM-4872-unit-start-plus-minus-two-port-pairs-and-the-z3-trivial-symmetry-split-of-three-cycle-tournaments.md.
The source has pre-audit status at that checkpoint. Its Fibonacci numbers
count certain unit-start or few-cycle classes; its +/-2 patches change old
arcs. The present increment is induced vertex arrival/deletion. No identity
between those carriers is assumed, and their published numerical findings
are not claimed as Session 191 discoveries.

## Precisely what this packet adds

The deliverables are the explicitly proved rooted bounded-increment
construction with its independently checked finite atlases; exact order-
by-order minima and the stated all-order cutoff consequences; the exact
all-order marked interval-chain minimum and equality families; and the
labelled neutral-extension count as an elementary SCC corollary. These do
not resolve global spectrum completeness or classify all tournaments.
