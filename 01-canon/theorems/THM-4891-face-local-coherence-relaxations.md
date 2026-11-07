---
id: THM-4891
title: "Face-local relaxations of parallel coherence on the permutohedron. Parallel coherence (tournaments, Session 210) is coherence on every square plus every hexagon 2-face. Square coherence alone gives set-context tournaments: the arc a-b may depend on the set of elements before the pair. These have at most n!/2^floor(n/2) sinks (Session 188's bound), with equality at n = 2, 3, 4, 6 and maximum 25 at n = 5. Hexagon coherence alone lets the arc depend on the relative order of the other elements. Exact sink spectra: n = 4, squares {0..6}, hexagons {0..9, 12} (7 realised by 48 acyclic orientations); n = 5, squares {0..25}, hexagons {0..54, 56, 60}; n = 6, squares {0..90}. Acyclic witnesses give 7 and 21 under each relaxation. So neither face type alone excludes 7 or 21. Both together are zone-wide coherence (tournaments), which does. Whether an intermediate relaxation still excludes them is OPEN"
status: >
  PROVED:
  (i) the face-generation equivalence (the classical zonotope fact that zones are connected through 2-faces, written
  out for the permutohedron);
  (ii) and (iii) the dictionaries;
  the recursion f(n) <= C(n,2) f(n-2).
  FINITE-EXACT:
  (iv) at n = 3, 4 by exhaustion of all 2^3, 2^6, 2^12 and 2^24 class orientations, with acyclicity by topological sort;
  at n = 5 by SAT (CaDiCaL 1.5.3 via PySAT, totalizer cardinality encoding, one call per sink count 0..60);
  at n = 6 by SAT for squares (every k = 0..90, plus a witness with at least 90 sinks);
  acyclic witnesses with 7 and 21 sinks by cycle-elimination SAT, each re-verified by an independent rebuild of the
  oriented Cayley graph.
  The independent audit re-derived (i)-(iii) by union-find on geometric opposite-edge relations for n <= 7. It
  recomputed the n = 3, 4 spectra exhaustively in C, and the n = 5 spectra with Glucose4 / sequential counter and
  Minisat22 / kmtotalizer. It found its own acyclic witnesses and checked an acyclic f(6) = 90 witness.
  INHERITED:
  the parallel-coherent carrier and the boundary "dropping coherence admits every count to n!/2" are Session 210
  (session210/coherent-permutohedron/PROOF.md, Theorems 1 and 4);
  the bound n!/2^floor(n/2) and its orbit proof are Session 188's tournament bound
  (05-knowledge/results/session188_structural_bridges_20261104/equal_nine_ambient/PROOF.md (1)) - the orbit argument
  applies verbatim to set-context tournaments;
  the gaps 7 and 21 are THM-1370.
  New: equality in the bound at n = 4 and 6 (6 and 90 sinks, where tournaments reach only 5 and 45), the exact
  relaxed spectra, and the 7/21 witnesses.
  NOT CLAIMED: any new exclusion proof; that zone-wide coherence is necessary for excluding 7 and 21; Chan-Pak (OPEN).
  AUDIT: PASS WITH CORRECTIONS (MAJOR-2 on scope, applied),
  05-knowledge/results/one_apex_collision_balance_20261106/audit/AUDIT.md.
session: opus-2026-11-06-S84
source: 05-knowledge/results/one_apex_collision_balance_20261106/README.md
related:
  - THM-1370
  - THM-207
  - THM-4890
  - session210/coherent-permutohedron/PROOF.md
  - session211/score-permutohedron/PROOF.md
  - 05-knowledge/results/session188_structural_bridges_20261104/equal_nine_ambient/PROOF.md
verification:
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax2.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax3.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_witness.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/audit/aud_perm.py
---

# THM-4891 - face-local relaxations of parallel coherence

## Objects

- **The permutohedron graph.** Its vertices are the arrangements of [n]. Its edges swap adjacent positions: the edge
  that swaps positions i, i+1 of pi exchanges the values {pi_i, pi_(i+1)}.
- **Orientations and sinks.** An orientation directs each edge into one of its two endpoints. A sink is an arrangement
  whose n - 1 incident edges all point into it.
- **Parallel coherence** (Session 210). All edges swapping the same pair {a, b} (one zone) point the same way, so that
  one of a, b always comes first at the sink end. These orientations are exactly the tournaments, and their sinks are
  exactly the Hamiltonian paths.
- **2-faces.** There are two kinds:
  - squares, from commuting swaps s_i, s_j with |i - j| >= 2;
  - hexagons, from the six arrangements of three values {a, b, c} in positions i, i+1, i+2. Opposite edges of the
    6-cycle swap the same pair, with the third value on opposite sides.
- **Partial coherence.** An orientation is square-coherent (resp. hexagon-coherent) if in every square (resp. hexagon)
  opposite edges have equal direction vectors.

## Statement

**(i) Face generation.** An orientation is parallel-coherent iff it is both square-coherent and hexagon-coherent.

**(ii) Squares.** Square-coherent orientations correspond bijectively to *set-context tournaments*. This name differs
from the repository's earlier "context tournament", which means a tournament on exterior vertices.
- A set-context tournament assigns an arc kappa(a, b; S) to each pair {a, b} and each set S of the other elements.
- Sinks are the arrangements pi whose adjacent pairs (pi_i, pi_(i+1)) all follow kappa(.; {pi_1, ..., pi_(i-1)}).
- Equivalently: in every square 2-face of the Boolean lattice B_n, choose one of its two sides. The sinks are the
  maximal chains that use the chosen side at every pair of consecutive steps.
- The maximum number of sinks f(n) satisfies f(1) = f(2) = 1 and f(n) <= C(n,2) f(n-2), so f(n) <= n!/2^floor(n/2).
- Its values are f(3) = 3, f(4) = 6, f(5) = 25 < 30 and f(6) = 90.

**(iii) Hexagons.** Hexagon-coherent orientations correspond bijectively to tables kappa(a, b; sigma), where sigma is
the relative order of the other n - 2 elements. Sinks are the arrangements whose adjacent pairs follow kappa(.; sigma).
At n = 3 these are the tournaments.

**(iv) Exact sink spectra.** Acyclic orientations are in brackets.

| n | square-coherent | hexagon-coherent | parallel-coherent (tournaments) |
|---|---|---|---|
| 3 | {0,1,2,3} [{1,2,3}] | {1,3} | {1,3} |
| 4 | {0..6} [{1..6}] | {0..9, 12} [{1..9, 12}] | {1,3,5} |
| 5 | {0..25} | {0..54, 56, 60} | {1,3,5,9,11,13,15} |
| 6 | {0..90} | not computed | {1,3,5,9,...,33,37,41,43,45} (19 values; 7, 21, 35, 39 missing) |

- At n = 4, hexagon coherence gives 7 sinks for exactly 48 of the 4,096 orientations, all acyclic.
- At n = 5, explicit acyclic orientations with 7 and with 21 sinks exist under each of the two partial coherences
  (`permutohedron/witnesses.out`). The independent check confirms the sink count, acyclicity, and coherence on the
  stated face type and not on the other.
- Square coherence caps the count at f(4) = 6, so 7 first appears there at n = 5.
- 21 exceeds every relaxed maximum at n = 4, so n = 5 is the least order for 21 under each relaxation.
- The square-coherent spectra are full intervals 0..f(n) for n = 3..6.

**(v) Consequence.** Neither face type of local coherence excludes 7 or 21.
- 7 first appears at order 4 under hexagon coherence, and at order 5 under square coherence.
- 21 first appears at order 5 under both, acyclically.
- By (i), the two local relations together are zone-wide coherence, i.e. tournaments, which exclude 7 and 21
  (THM-1370).
- So no coherence condition on 2-faces of a single type carries the exclusion. It is not shown that full zone-wide
  coherence is necessary: intermediate relaxations, such as all hexagons plus some squares, are OPEN.

## Proof

**(i).**
- **Coherence implies both face conditions.** Opposite edges of a square or a hexagon swap the same pair, so they lie
  in the same zone.
- **Both face conditions imply coherence.** Write an edge of the zone {a, b} as (P, ab, Q), where P and Q are the
  sequences before and after the pair.
  - A hexagon on positions (|P|, |P|+1, |P|+2) moves the first element of Q across the pair, and conversely moves the
    last element of P. It keeps the two opposite {a, b}-edges coherent.
  - Squares permute P and Q internally.
  - These moves connect any two edges of the zone: first use hexagons until the prefix set is right, then arrange
    with squares.

**(ii).**
- **Dictionary.** Square moves are exactly the adjacent swaps that avoid the pair. So two {a, b}-edges are square-
  connected iff they have the same set of preceding elements, which gives the dictionary.
- **Boolean-lattice reading.** The prefix sets of pi form a maximal chain of B_n. Steps i and i+1 traverse the square
  of B_n with bottom S_(i-1) and top S_(i-1) + {pi_i, pi_(i+1)}.
- **Bound** (Session 188's orbit argument).
  - A sink's first two entries are the kappa(.; empty)-order of their pair: at most C(n,2) choices.
  - Given (pi_1, pi_2), the tail is a sink of the set-context tournament kappa'(a, b; S') = kappa(a, b; S' + {pi_1, pi_2})
    on the remaining n - 2 elements. The condition at position 2 is dropped.

**(iii).** Hexagon moves carry one element across the pair without changing the relative order of the others. By the
proof of (i), they connect all slots of the pair within one relative order.

**(iv), (v).** These follow from the computations listed in the status. All were re-run in `permutohedron/`, and the
audit repeated them independently. QED.

## Remarks

- **The landscape.** The face-local relaxations sit strictly between Session 210's two extremes:
  - tournaments: the sinks are the Hamiltonian paths, with gaps 7 and 21;
  - arbitrary acyclic orientations: every count 1..n!/2.
  Each relaxation loses 7 and 21 at the first order where the counts fit, so neither face type alone carries the gaps.
- **Analogy with THM-4890.** On the one-apex side, collisions of path counts are decided by a global prefix/suffix
  balance, not by local moves. This is an analogy only.
- **Parity.** Rédei parity is lost too: 0 and even counts occur from order 3 under squares and from order 4 under
  hexagons.
- **Top gaps.** The hexagon spectrum has gaps only near its maximum n!/2: 10 and 11 at n = 4; 55, 57, 58 and 59 at
  n = 5. This is a ceiling effect, like the top band of the tournament spectrum, and not a 7/21-type gap.
- **The n = 4 witness.** All 48 seven-sink hexagon-coherent orientations at n = 4 have exactly four pairs, out of six,
  whose arc flips with the order of the other two (`audit/aud_perm_flips.out`).
