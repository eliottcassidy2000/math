# Hypohamiltonian criticality and the lonely runner

**Source:** klein-2026-07-10-S238 (owner directive: think of hypohamiltonian
graphs and how they relate to our problem and past thinking)

A graph is **hypohamiltonian** when it is not Hamiltonian but every
vertex-deleted subgraph is — criticality under single deletion. The Petersen
graph is the unique smallest example; the field's method is to let
criticality do the work: every deletion restoring the property forces local
structure (3-connectivity, degree floors, girth constraints), and counting
those forced structures kills candidate orders wholesale.

**LRC(14)'s minimal counterexamples are hypo-objects by construction.** If S
were a 13-runner family with M(S) < 1/14 (or, for the wall, μ(S) = 0), then
LRC(≤13) — the project's axiom — gives every deleted family S∖{v} loneliness
WITH MARGIN: M ≥ 1/13 > 1/14, and Lemma A fattens the margin to a safe set
of measure ≥ 1/(91·max). So a counterexample is precisely "hypo-lonely":
thirteen simultaneous critical-covering constraints, one per runner — each
runner's danger set (v arcs of width 1/(7v), total mass exactly 1/7) must
individually cover the other twelve's fattened safe set. Criticality is not
a hypothesis to add; the induction structure forces it, exactly as in the
graphs.

**Zero-slack covers are rigid — the k = 7 dead-zone theorem is the
paradigm.** Seven arcs of length 1/7 have total mass exactly 1: they cover
the circle only as a perfect net (pairwise disjoint, centers equally
spaced), which confines the configuration to finitely many times — measure
zero, THM-689(A). This is the covering-side shadow of hypohamiltonicity:
a critical cover (remove any arc and a hole opens) with zero slack has no
freedom left. The project has met this rigidity before, every time at the
same place: the tight AP {1..13} sits at ratio exactly 13, where the first
safe window W_P = [1/(14·pmin), 13/(14·pmax)] degenerates to a point
(THM-689(C)), where μ hits zero, and where kps's three faces coincide
([μ = 0] = [ratio 13] = [dilated interval]). The tight AP is our Petersen
graph: the unique minimal critical object, and everything off it has slack
that the machinery converts to witnesses.

**Rédei's inversion, and the parity echo.** In tournaments — this project's
original home — existence is unconditional: every tournament has a
Hamiltonian path (Rédei), and the theorem's content is the REFINEMENT,
parity (the count is odd). In the lonely runner, existence (of a safe time)
is itself the question. But the parity refinement reappeared on the
certificate side: the parity pairing law (S222–S229, LRCParityPairing) —
live multipliers come in exact ± pairs, LM(q) is even for covering
families. A counterexample would have LM(q) = 0 at every modulus: even,
consistent, and silent; every non-counterexample announces itself with a
PAIR of witnesses. Rédei's odd count guarantees a path exists; our even
count means witnesses can only vanish together — the parity laws are the
same species of theorem, one guaranteeing existence, one structuring it.

**What the analogy suggests trying** (logged, not claimed): (i) the
13-fold critical-covering constraints as bulk arithmetic — each runner's
1/7 of danger mass must cover twelve other safe sets simultaneously;
summing the forced incidences against the exhaustive m_P census
(THM-689(B): min m_P ≥ 0.38 already at |P| = 5) may kill criticality by
counting alone, the way degree-sum arguments kill small hypohamiltonian
orders; (ii) hypohamiltonian graphs are built by gluing at critical
vertices — the rational-ray covers (THM-688) are our gluings, and the
adversarial families are glued critical pieces; (iii) girth constraints ↔
the relation lattice: short relations (Schur triples) are the short cycles,
and the extremal families are relation-saturated exactly as the critical
graphs are girth-constrained.

The meta-pattern: the problem keeps rewarding the same move — treat the
extremal object as critical, measure the slack, and let zero slack force
rigidity. The triangle gave the constants; the relation lattice gave the
deviations; criticality gives the endgame its shape.

**Related:** THM-689, THM-687/688, LRCParityPairing (the Rédei transplant),
kps's describing-the-wall reflection, `everything-is-the-triangle.md`.
