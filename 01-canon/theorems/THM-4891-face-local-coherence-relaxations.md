---
id: THM-4891
title: "Face-local relaxations of parallel coherence on the permutohedron: neither squares nor hexagons alone keep 7 and 21 out. Parallel coherence (tournaments, Session 210) is coherence on every square together with every hexagon. Square coherence alone gives context tournaments (the arc a-b may depend on the set of elements before the pair), with at most n!/2^floor(n/2) sinks (equality at n = 2, 3, 4, 6; 25 at n = 5). Hexagon coherence alone lets the arc depend on the relative order of the other elements. Exact sink spectra: n = 4 squares {0..6}, hexagons {0..9, 12} with 7 realised by 48 acyclic orientations; n = 5 squares {0..25}, hexagons {0..54, 56, 60}; acyclic witnesses for 7 and 21 under each. The exclusion of 7 and 21 needs both face types, i.e. zone-wide coherence"
status: >
  PROVED: (i) the face-generation equivalence (the classical zonotope fact that zones are connected through 2-faces,
  written out for the permutohedron), (ii) and (iii) the dictionaries, and the recursive square bound
  f(n) <= C(n,2) f(n-2).
  FINITE-EXACT: (iv) n = 3, 4 by exhaustion of all 2^3, 2^6, 2^12 and 2^24 class orientations, with acyclicity by
  topological sort. n = 5 by SAT (CaDiCaL 1.5.3 via PySAT, totalizer cardinality encoding, one call per sink count
  0..60). Acyclic witnesses by cycle-elimination SAT (CEGAR), each re-verified by an independent rebuild of the
  oriented Cayley graph. The n = 6 equality f(6) = 90 is a SAT witness for at least 90 sinks.
  INHERITED: the parallel-coherent carrier and the "dropping coherence admits every count to n!/2" boundary are
  Session 210 (session210/coherent-permutohedron/PROOF.md, Theorems 1 and 4). The gaps 7 and 21 are THM-1370.
  NOT CLAIMED: any new exclusion proof; Chan-Pak OPEN.
session: opus-2026-11-06-S84
source: 05-knowledge/results/one_apex_collision_balance_20261106/README.md
related:
  - THM-1370
  - THM-207
  - THM-4890
  - session210/coherent-permutohedron/PROOF.md
  - session211/score-permutohedron/PROOF.md
verification:
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax2.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_relax3.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/permutohedron/perm_witness.py
---

# THM-4891 - face-local relaxations of parallel coherence

## Objects

- **The graph.** The permutohedron graph on S_n has the arrangements of [n] as vertices. Its edges swap adjacent
  positions; the edge swapping positions i, i+1 of pi swaps the values {pi_i, pi_(i+1)}.
- **Orientations and sinks.** An orientation directs each edge into one of its two arrangements. A sink is an
  arrangement whose n - 1 incident edges all point into it.
- **Parallel coherence.** Session 210 calls an orientation parallel-coherent when all edges swapping the same pair
  {a, b} (one zone) point the same way, i.e. one of a, b always comes first at the sink end. Those orientations are
  exactly the tournaments, and their sinks are the Hamiltonian paths.
- **2-faces.** There are two types:
  - squares, from commuting swaps s_i, s_j with |i - j| >= 2;
  - hexagons, from the six arrangements of three values {a, b, c} in positions i, i+1, i+2.
  An orientation is square-coherent (resp. hexagon-coherent) if, in every square (resp. hexagon), opposite edges point
  the same way.

## Statement

**(i) Face generation.** An orientation is parallel-coherent iff it is square-coherent and hexagon-coherent.

**(ii) Squares.** Square-coherent orientations correspond bijectively to context tournaments.
- **Context tournaments.** A context tournament assigns an arc kappa(a, b; S) to each pair {a, b} and each set
  S of the other elements.
- **Sinks.** An arrangement pi is a sink iff each adjacent pair (pi_i, pi_(i+1)) follows kappa(.; {pi_1, ..., pi_(i-1)}).
- **Boolean-lattice form.** Equivalently: in every square 2-face of the Boolean lattice B_n, choose one of its two
  sides. The sinks are the maximal chains that use the chosen side at every pair of consecutive steps.
- **Bound.** The maximum number of sinks f(n) satisfies f(1) = f(2) = 1 and f(n) <= C(n,2) f(n-2), so
  f(n) <= n!/2^floor(n/2).
- **Values.** f(3) = 3, f(4) = 6, f(5) = 25 < 30, f(6) = 90.

**(iii) Hexagons.** Hexagon-coherent orientations correspond bijectively to tables kappa(a, b; sigma), where sigma is
the relative order of the other n - 2 elements. The sinks are the arrangements whose adjacent pairs follow
kappa(.; their sigma). At n = 3 this is exactly the tournaments.

**(iv) Exact sink spectra.** All orientations, with the acyclic ones in brackets:

| n | square-coherent | hexagon-coherent | parallel-coherent (tournaments) |
|---|---|---|---|
| 3 | {0,1,2,3} [{1,2,3}] | {1,3} | {1,3} |
| 4 | {0..6} [{1..6}] | {0..9, 12} [{1..9, 12}] | {1,3,5} |
| 5 | {0..25} | {0..54, 56, 60} | {1,3,5,9,11,13,15} |

- **Counts at n = 4.** 7 sinks occur under hexagon coherence for exactly 48 of the 4,096 orientations, all of them
  acyclic.
- **Witnesses at n = 5.** Explicit acyclic orientations with 7 and with 21 sinks exist under each of the two partial
  coherences (`permutohedron/witnesses.out`). The independent check confirms the sink count and acyclicity, and that
  each witness is coherent on the stated face type and not on the other.
- **Square coherence at n = 4.** It caps the count at f(4) = 6, so 7 first appears at n = 5.
- **Order 5 for 21.** 21 exceeds every relaxed maximum at n = 4, so order 5 is the least order for 21 under each
  relaxation.

**(v) Consequence.**
- Neither face type of local coherence keeps 7 or 21 out. 7 first appears at order 4 under hexagon coherence and at
  order 5 under square coherence. 21 first appears at order 5 under both, acyclically.
- By (i), the conjunction of the two local relations is zone-wide coherence, which gives tournaments.
- So the exclusion of 7 and 21 (THM-1370) is a property of the global zone structure. A coherence condition on 2-faces
  of only one type cannot carry it.

## Proof

**(i).**
- **Coherence implies both face conditions.** Opposite edges of a square or hexagon swap the same pair, so they lie in
  the same zone.
- **Both face conditions imply coherence.** An edge of zone {a, b} is a triple (P, ab, Q), where P and Q are the
  sequences of elements before and after the pair.
  - A hexagon on positions (|P|, |P|+1, |P|+2) moves the first element of Q across the pair, and conversely the last
    element of P. It keeps the two opposite {a, b}-edges coherent.
  - Squares permute P and Q internally.
  - These moves connect any two edges of the zone: first move elements across with hexagons until the prefix set is
    right, then arrange with squares.

**(ii).**
- **Dictionary.** Square moves are exactly the adjacent swaps that avoid the pair. So two {a, b}-edges are joined by
  squares iff they have the same set of preceding elements. That set is S, and the context-tournament dictionary
  follows.
- **Boolean-lattice reading.** The chain of prefix sets of pi is a maximal chain of B_n. Steps i and i+1 traverse the
  square of B_n with bottom S_(i-1) and top S_(i-1) + {pi_i, pi_(i+1)}.
- **Bound.**
  - A sink's first two entries are the kappa(.; empty)-order of their pair: at most C(n,2) choices.
  - Given (pi_1, pi_2), the tail is a sink of the context tournament kappa'(a, b; S') = kappa(a, b; S' + {pi_1, pi_2})
    on the remaining n - 2 elements. The extra condition at position 2 is dropped.
  - Iterating gives n!/2^floor(n/2).

**(iii).** Hexagon moves move one element across the pair without changing the relative order of the others, and
they connect all slots of the pair within a fixed relative order (proof of (i)).

**(iv), (v).** By the computations listed in the status, all re-run in `permutohedron/`. QED.

## Remarks

- **The landscape.** The face-local relaxations sit strictly between Session 210's two extremes:
  - tournaments: sinks are Hamiltonian paths, with the gaps 7 and 21;
  - arbitrary acyclic orientations: every count 1..n!/2.
- **Character of the gaps.** Both relaxations lose the gaps at the first order where the count fits, so the gaps are
  not local. This matches THM-4890 on the one-apex side: collisions of path counts are decided by global prefix/suffix
  balance, not by local moves.
- **Parity.** The relaxed spectra also lose Rédei parity: 0 and even counts occur, from order 3 under squares and from
  order 4 under hexagons. Oddness needs both face types too.
- **Top gaps.** The hexagon spectrum has gaps only near its maximum n!/2: 10 and 11 at n = 4; 55, 57, 58, 59 at n = 5.
  This is a ceiling effect, like the top band of the tournament spectrum, not a 7/21-type gap.
