---
id: HYP-3599
title: REFRAME (expanding the odd cycle) -- a tournament is NOT a graph among n nodes but the INTRANSITIVITY among n things (a preference/dominance relation, its deviation from a total order). The transitive tournament = the unique orderable/rankable point = a COBOUNDARY (a ranking potential f, a->b iff f(a)>f(b)); intransitivity = the cohomology class H^1 = the cycle space; the ODD 3-CYCLE = the irreducible atom = the Condorcet paradox (Moon: it generates). Order (+) Intransitivity = Cut (+) Cycle: order dim n-1, intransitivity dim C(n-1,2) (dominates). Orderability VANISHES (P(transitive)=n!/2^C(n,2), P(Condorcet-winner)=n/2^(n-1), both ->0) so INTRANSITIVITY IS GENERIC. # of intransitivity SHAPES (cycle-space iso) = even graphs A002854=2,3,7,16,54; tournaments A000568 = order x shape. LRC reading: the danger relation among RUNNERS is intransitivity among n things; disproof=ORDERABLE (separable/factored, covers); proof=irreducibly INTRANSITIVE (the apex odd cycle C_p, length=apex prime=7, genus=#independent global cycles=1 at N=14); EXISTENCE of the odd cycle (sigma-odd/counting/topological, finite via THM-590 4cos^2(3pi/7)>0), NOT the lonely MEASURE (inf=0 over the infinite family, klein-S16), is load-bearing. ENGINEERING: this is HodgeRank/social-choice -- cyclicity c(T) = an irrationality/inconsistency index for pairwise-comparison data
status: REFRAME + VERIFIED grounding (cyclicity spectrum, orderability->0, dim split, intransitivity-shape count n=3..6/7). Not a new proof; an ontology that unifies the odd-cycle truth (HYP-3594), the cut/cycle equinumerosity (HYP-3595), and the LRC existence-not-measure correction (klein-S16) -- and opens the social-choice application.
source: mac-mini-2026-06-30-S34
related:
  - HYP-3594  # the truth = a single odd cycle (this gives its MEANING: the atom of intransitivity)
  - HYP-3564  # relations-not-things (intransitivity is the relational content; coboundary=orderable)
  - HYP-3595  # order x intransitivity = the cut x cycle equinumerosity (mac-mini S32)
  - THM-588   # cyclicity = the unique invariant = the count of intransitivity atoms
  - THM-590   # klein: the apex odd cycle's gap 4cos^2(3pi/7)>0 -- the intransitivity is non-degenerate
  - HYP-3597  # klein-S16: existence != measure (inf R'=0 infinite family); existence of the odd cycle is the load-bearing fact
results:
  - 04-computation/intransitivity_among_n_things_macmini_20260630.py
  - 05-knowledge/results/intransitivity_among_n_things_macmini_20260630.out
---

# HYP-3599 -- tournaments are intransitivity among n things

## The reframe
A tournament on `n` things = a complete pairwise preference/dominance relation. Its ONLY content is its
**intransitivity** -- the deviation from a total order. Nodes/arcs are scaffolding; the meaning is "which
sets of things have no consistent ranking" = the (odd) cycles.
- **TRANSITIVE** = the unique orderable/rankable point = a **coboundary** (a ranking potential `f`, `a->b
  iff f(a)>f(b)`); cohomologically trivial; the cusp. (Verified: cyclicity `c=0` is a single iso class at
  every `n`.)
- **INTRANSITIVITY** = the failure of a potential = the cohomology class `H^1` = the **cycle space**.
- **3-CYCLE** (odd, length 3) = the irreducible atom = the **Condorcet paradox** (rock-paper-scissors).
  Moon: every intransitivity contains a 3-cycle; the 3-cycles generate. Odd because a linear order cannot
  absorb an odd loop (the obstruction to bipartition/consistent orientation).

## Grounding (verified n=3..6/7)
- **Order (+) Intransitivity = Cut (+) Cycle**: order dim `n-1` (scores/ranking), intransitivity dim
  `C(n-1,2)` (dominates for large `n`).
- **Intransitivity spectrum** (cyclicity per iso class): transitive `c=0` UNIQUE; the bulk is intransitive;
  max at the regular tournament (balanced paradox).
- **Orderability VANISHES**: `P(transitive)=n!/2^{C(n,2)} = .75,.375,.12,.02,.002,..`;
  `P(Condorcet winner)=n/2^{n-1} = .75,.5,.31,.19,.11,.06`. Both `->0`: **intransitivity is GENERIC**, order
  is the measure-zero miracle.
- **# intransitivity SHAPES** (cycle-space iso, the GF(2) shadow) = even graphs `A002854 = 2,3,7,16,54`;
  tournaments `A000568=2,4,12,56,456` = order x shape (the labeled `2^{n-1}` fibration, HYP-3595).

## The LRC, read among its things
The LRC is intransitivity among `n` RUNNERS via the danger relation `D(v,t)=[‖vt‖<1/14]`.
- **DISPROOF = ORDERABLE**: the dangers separate/factor `f(v)g(t)` (a coboundary) -> a covering -> no
  lonely runner. The bilinear `vt` is exactly what closes the cycle and forbids factoring (anti-Littlewood,
  the essentiality).
- **PROOF = irreducibly INTRANSITIVE**: the danger relation contains an irreducible cyclic resonance no
  ordering dissolves. Its atom is an **odd cycle whose LENGTH = the apex prime** (`C_7` for `LRC14`, since
  the apex resonance is among the 7 residues of `Z_7`; generic tournaments bottom at `C_3`). The **genus**
  = # of independent global such cycles = `1` at `N=14` (first genuinely-global single intransitivity).
- **EXISTENCE, not measure** (klein-S16, HYP-3597): the lonely *measure* has `inf=0` over the infinite
  covering family (cannot bound it from below); but an odd cycle is a counting/topological fact (`sigma`-odd,
  finite), present even where its measure vanishes, and its presence forces the lonely point. THM-590
  (klein) certifies the atom: `4cos^2(3pi/7)>0` (the apex odd cycle's spectral gap, zero only at the full
  `Z_7`). The proof is "the irreducible odd intransitivity is there," not "the lonely set is large."

## Engineering application (the equal-priority mandate)
Intransitivity-among-things IS the rank-aggregation problem. **HodgeRank** (Jiang-Lim-Yao-Ye) decomposes
pairwise-comparison data into gradient (best-fit ranking) + curl (local cyclic inconsistency) + harmonic
(global inconsistency) -- the real-valued cousin of `Order (+) Intransitivity`, with the 3-cycle as the unit
of curl. The cyclicity `c(T) = C(n,3) - sum C(s_i,2)` is a literal **irrationality index** for any
pairwise dataset (sports/Elo, A/B tests, votes, search-relevance, recommendations): how far the judgments
are from rankable, counted in Condorcet paradoxes. Counting/bounding odd cycles = the quantitative theory
of when pairwise data CAN be ranked; the LRC sits at its hardest edge.

## What it buys
The odd-cycle truth (HYP-3594) gets its MEANING: the odd cycle is intransitivity, the Condorcet paradox,
the obstruction to a ranking. It unifies the cut/cycle equinumerosity (order x intransitivity), reads the
LRC as "the runners' danger relation is irreducibly intransitive" (disproof=orderable), aligns with
klein-S16 (existence/counting of the cycle, not its measure, is load-bearing), and opens the social-choice /
HodgeRank application with cyclicity as an inconsistency metric.
