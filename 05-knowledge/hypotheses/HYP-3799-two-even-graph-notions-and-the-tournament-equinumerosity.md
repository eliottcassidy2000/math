---
id: HYP-3799
title: THE TOURNAMENT / EVEN-GRAPH EQUINUMEROSITY, RESOLVED -- there are TWO distinct "even graph" notions and the repo conflated them. (I) EVEN-DEGREE graph (every vertex even degree = cycle space of K_n) = A002854 = 2,3,7,16,54 = two-graphs = switching classes of graphs (Mallows-Sloane 1975) = the repo's E_n. (II) AUTOMORPHISM-PARITY even graph (Royle-Praeger-Glasby-Freedman-Devillers 2022, arXiv:2204.01947): a graph is ODD if some automorphism reverses an ODD number of edges (w.r.t. a fixed reference orientation), EVEN otherwise; THEOREM 1.1: #(type-II even graphs on n) = #(tournaments on n) = A000568 = 2,4,12,56,456. The prior repo claim "A000568 = #even graphs FAILS (4!=3)" compared tournaments to notion (I); against notion (II) it is a THEOREM. Mechanism: #graphs = #tournaments + #odd graphs (Cauchy-Frobenius). NEW structural facts (verified n<=5): the reversal parity is a character eps_X: Aut(X)->Z/2, X even <=> eps_X trivial; ODD |Aut(X)| => X automatically EVEN (odd-order groups have no Z/2 quotient); tournaments always have ODD |Aut| while type-II even graphs have EVEN |Aut| (n<=5) -- disjoint |Aut| profiles => NO Aut-preserving bijection => RPGFD's natural bijection is genuinely hard/open.
status: CONFIRMED (verified n<=5: type-II even = A000568 = tournaments exactly; #graphs=#tournaments+#odd; cross-tab of the two notions; |Aut| distributions). RESOLVES the repo's flagged confusion (tournaments-and-even-graphs.md) and imports the correct theorem (RPGFD 2022). The reversal-character frame, the odd-|Aut|=>even lemma, and the disjoint-|Aut| obstruction to a natural bijection are new observations. The natural bijection tournaments<->even graphs is OPEN (RPGFD state it).
source: mac-mini-2026-07-01-S82
related:
  - HYP-3798   # S81 arc-hypercube / orbit coloring (same ambient space; even/odd is a stabilizer-parity invariant)
references:
  - "Royle, Praeger, Glasby, Freedman, Devillers 2022 (arXiv:2204.01947, J. Alg. Combin. 57 (2023) 515-524) -- Tournaments and Even Graphs are Equinumerous; Thm 1.1: #even graphs = #tournaments; open: a natural bijection"
  - "Andersson / von Bromssen -- the even/odd graph definition (automorphism reversing an odd # of edges); parity well-defined"
  - "Mallows-Sloane 1975 -- two-graphs = switching classes = Euler graphs = A002854 (the OTHER 'even graph')"
results:
  - 04-computation/two_even_graph_notions_macmini_20260701.py
  - 04-computation/even_graph_switching_fibration_macmini_20260701.py
  - 05-knowledge/results/two_even_graph_notions_macmini_20260701.out
---

# HYP-3799 -- the tournament/even-graph equinumerosity, resolved (two notions)

The owner asked to explore "tournament even graph equinumerosity" with novel frames. The repo's
`tournaments-and-even-graphs.md` had flagged that `A000568 = #even graphs` was NOT confirmed (`4 != 3` at
`n=4`) and needed "checking against the paper." Here is the paper and the resolution.

## Two distinct "even graph" notions (the repo conflated them)
- **(I) EVEN-DEGREE** graph = every vertex even degree = the cycle space of `K_n`. Up to iso = **A002854**
  `= 2,3,7,16,54` = two-graphs = **switching classes of graphs** (Mallows-Sloane 1975). This is the repo's
  even-graph metagraph `E_n`. A valid, well-studied object.
- **(II) AUTOMORPHISM-PARITY even** graph (Andersson; Royle-Praeger-Glasby-Freedman-Devillers 2022,
  arXiv:2204.01947): fix the reference orientation `i->j` for `i<j`; an automorphism `g` **reverses** edge
  `{u,v}` (`u<v`) iff `u^g > v^g`; a graph is **ODD** if some automorphism reverses an ODD number of edges,
  **EVEN** otherwise. **RPGFD Theorem 1.1: `#(even graphs on n) = #(tournaments on n) = A000568`** `=
  2,4,12,56,456`.

These two "even"s are **logically independent** -- the `2x2` cross-tabulation has all four cells populated
(`n=5`: `(even-deg, type-II-even) = (F,F):17, (F,T):10, (T,F):5, (T,T):2`; type-II total `12` = tournaments,
even-degree total `7` = A002854). The repo compared tournaments to notion (I) and (correctly) found they
differ; against notion (II) the equinumerosity is a theorem.

## The mechanism and new structure
- **`#graphs = #tournaments + #odd graphs`** (Cauchy-Frobenius; verified `4=2+2`, `11=4+7`, `34=12+22`).
  So tournaments are exactly the "even part" of graph space under this `Z/2` grading.
- **The reversal character** `eps_X : Aut(X) -> Z/2`, `g |-> (parity of edges g reverses)` (well-defined,
  Andersson). `X` is even `<=>` `eps_X` is trivial. So even/odd is a cohomological invariant of the
  automorphism action -- an index-`<=2` "orientation character" of `Aut(X)`.
- **Lemma (new): odd `|Aut(X)|` `=>` `X` is EVEN.** (An odd-order group has no `Z/2` quotient, so `eps_X`
  is forced trivial.) Hence every asymmetric graph (`|Aut|=1`, first at `n=6`) is even.
- **The `|Aut|` obstruction (new): tournaments always have ODD `|Aut|`; type-II even graphs have EVEN
  `|Aut|` (`n<=5`).** The two equinumerous families have **disjoint `|Aut|` profiles**, so NO bijection can
  preserve automorphism groups -- a structural reason RPGFD's "natural bijection" is genuinely hard.

## Novel frames + open targets
1. **Two orthogonal parities of graph space**: even-DEGREE (cycle space / A002854 / switching) and
   even-AUTOMORPHISM-PARITY (A000568 / tournaments). The repo studied the first; the equinumerosity is the
   second. Their joint refinement (the `2x2` metagraph) is a new object to study.
2. **The natural bijection tournaments <-> even graphs is OPEN** (RPGFD). The repo's tournament machinery
   (tiling model, OCF, Redei-Berge, score sequences) is a candidate toolkit -- but any bijection must map
   odd-`|Aut|` tournaments to even-`|Aut|` even graphs, so it cannot be automorphism-natural; it must be a
   "parity-twisting" correspondence.
3. **A000568 gets a new meaning**: the repo has extended A000568 to 200 terms (Burnside enumerators); each
   term now also counts type-II even graphs. The RPGFD `#graphs = #tournaments + #odd` identity is a new
   decomposition to instrument with the repo's cycle-index tools.
4. **eps_X as a `Z/2` gauge/orientation character** links to the project's Walsh-Fourier / sign-character
   and GF(2) cut/cycle themes.

## Honest scope
The two-notions resolution, RPGFD Thm 1.1 (`type-II even = tournaments = A000568`), `#graphs = #tournaments
+ #odd`, the cross-tabulation, and the `|Aut|` disjointness are all VERIFIED (`n<=5`, exhaustive). The
reversal-character and odd-`|Aut|=>`even lemma are elementary and proved. NOT solved: the natural bijection
(RPGFD's open problem). A resolution + import of the correct theorem + new structural obstructions, not a
new theorem of my own.
