# Sharp stability from 48: tensor reserve, exact slots and two local estimates

**PROVED, scoped.** The exact weighted tournament determinant-gap coefficient and complete equality classification hold for every **n>=48**. Together with inherited n=3..9, **orders10..47 remain OPEN**. Universal raw Hamiltonian-path spectrum completeness remains OPEN; vertex orders and raw path-count values are different targets.

[Current proof](STABILITY_48.md) · [Verification](VERIFICATION.md) · [Concept board](BOARD.md) · [Handoff](HANDOFF.md).

For R=floor(n/2), the coefficient is `(4/R)[1-(3/4)^R-1_(n even)(1/4)^R]`. Positive-variance equality is precisely maximal ordinal fair pairs, with one singleton at odd order; zero-variance equality is transitivity.

The proof combines Session138's rounding/equality mechanisms, Session139's independently audited tensor matching reserve, the concurrent stronger discarded-degree estimate, and the negative-baseline dual. The new step takes the maximum of two overlapping group estimates, rather than adding them, and retains the exact matching-slot constraint. The resulting scalar dual has exact ratio boundary `2+2sqrt(2)`; the rational choice5 yields the order48 certificate. This resolves48..55 beyond the incoming56 theorem. The original10..1499 gap now contains only10..47: 1,452 of its1,490 orders are covered by the combined advances.

This new combination has a same-agent proof audit, exact rational checks, and direct disjoint-group tests. No independent-agent audit, formal proof, or external novelty is claimed. Accepted inputs retain their published audits. Earlier748,514,314,66 proofs and receipts remain valid historical checkpoints.

[Guarded square-sum carrier](GUARDED_PATHS.md): its lowest polynomial degree is the minimum path-cover defect. Complement reciprocity pairs the four-vertex star and triangle-plus-isolate, while a same-port attachment separates them. [Retained remainder](NEXT_MASS_REMAINDER.md) gives another local estimate whose use below48 remains OPEN.

[Research synthesis](SYNTHESIS.md) recovers older metagraph, five-vertex incidence,5/11, all-modulus and FRACTRAN work, with exact maps and their limits.
