---
id: HYP-2021
status: OPEN
source: codex-2026-06-01-S537
related:
  - HYP-1981
  - HYP-1982
  - HYP-1998
  - HYP-1999
  - HYP-2001
  - HYP-2009
  - HYP-2012
  - HYP-2016
  - HYP-2017
  - HYP-2018
  - HYP-2019
  - HYP-2020
---

# HYP-2021: Decorated binary danger-lex is the restricted tournament stack for LRC

**Claim.** LRC should be reframed through a small pure tournament quotient,
not through the raw observer-marked class set.  The best current candidate is
a decorated binary danger-lex stack:

```text
Layer 1: binary danger-lex tournament
  blocker > observer > safe

Layer 2: sentinel or left/right defect marks
  remembers which side of the safe arc failed

Layer 3: CRT/coupled-debt labels
  remembers residue-channel arithmetic and higher resonance debt

Layer 4: endpoint-owner protection tournament
  dual static obstruction; a counterexample must create an SCC
```

In this language, LRC becomes:

```text
For every primitive speed set, the decorated BLEX walk exhibits a good-only
source class before the endpoint-owner protection tournament can close into a
realizable coupled-debt SCC.
```

**Evidence.** `lrc_restricted_tournament_mappings_s537.py` compared six
functors on open wall cells for `n=4..7`.  Binary danger-lex had no mixed
good/bad fibers and sharply compressed the standard source quotient:

```text
n=4: standard 8  -> BLEX 4   good-only 1
n=5: standard 22 -> BLEX 7   good-only 2
n=6: standard 50 -> BLEX 11  good-only 4
n=7: standard 96 -> BLEX 17  good-only 6
```

Safe-only had the same sampled counts but forgets blockers entirely; it needs
external arithmetic to explain why the top safe stratum is hit.  Two-sentinel
safe-arc classes were pure and geometrically faithful but much larger, so their
role is boundary memory rather than compression.  CRT-channel classes were
small for composite `n`, but HYP-2017 says raw parity/channel data becomes
vacuous at higher prime-power levels; coupled debt labels are mandatory.

**Interpretation.** Binary danger-lex is the smallest observed fixed-observer
quotient that preserves source purity while retaining runner half-turn
geometry.  The sentinel and endpoint layers repair its main weakness: BLEX
knows how many blockers exist, but not why the safe arc fails.  The CRT-debt
layer repairs the arithmetic weakness: a tiny tournament quotient needs a
decorated memory of unit residues and higher-order resonance.

**External hooks.** This stack matches known outside languages for LRC:
congruence-class proofs for small cases, zonotope covering-radius
discretizations, and tournament Hamiltonian-path/fixed-endpoint structure.
The point is not to replace these languages, but to ask which tournament
isomorphism classes an arithmetic LRC walk can exhibit after each language has
been compressed to its proof-relevant finite shadow.

**Concurrent integration.** Upstream HYP-2018/HYP-2019 landed while this
session was closing.  They identify metric retention as the source of
restricted realizability and name the near-graph/circular-indifference mapping
as a sweet spot.  The BLEX stack should be read as the tournament-directional
partner of that result: near-graph records the symmetric `distance < 1/n`
relation, while BLEX orients the same threshold data around the observer into
source, defect, and protection classes.

A later upstream HYP-2020 claimed the broader restricted quotient-stack
language.  This HYP is therefore the narrow BLEX refinement of that stack:
one specific compressed pure source quotient, with a proposed decoration set
for carrying the missing endpoint/CRT memory.

**Predictions.**

1. The full binary danger-lex class count is a Ferrers/round-tournament
   convolution over blocker and safe blocks.
2. Sentinel danger-lex has pure good fibers with class counts strictly between
   BLEX and two-sentinel safe-arc.
3. A source-avoiding BLEX walk projects to a nontrivial SCC in an
   endpoint-owner protection tournament.
4. For composite `n`, raw CRT majority classes are insufficient, but
   CRT-debt-decorated classes remain small enough to enumerate and may isolate
   the first real obstruction at `n=14` or `n=18`.

**Next tests.**

1. Implement sentinel danger-lex and left/right defect functors.
2. Add score histograms, SCCs, directed cycles, edge flips, and Hamiltonian-path
   counts to the S537 probe.
3. Compute BLEX class-count formulae against round/Ferrers menu counts.
4. Build endpoint-owner protection tournaments for the same speed sets and
   test the predicted BLEX-avoidance/SCC correlation.

**Files.** `04-computation/lrc_restricted_tournament_mappings_s537.py`;
`05-knowledge/results/lrc_restricted_tournament_mappings_s537.out`;
`07-reflections/lrc-restricted-tournament-functor-zoo-s537.md`.
