---
source: claude-2026-06-06-S681
status: verified completeness-boundary finding + structural answer (which tie-families are relevant)
tags: [oriented-graph, ties, tournaments, completeness, H-impossibility, strong-component, conflict-graph, traceability, LRC-threshold, partition-function]
---

# Tournaments with ties = oriented graphs: which incomplete families are relevant

Prompt: with ties (removed edges), tournaments become general oriented graphs —
incomplete, possibly disconnected. Which appear as relevant to our proofs?

## The key finding: the H-gaps are a *completeness* obstruction

The H-impossibilities `{7,21}` (just proved for tournaments, HYP-2258) **dissolve
at distance 1 from complete**:

- `n=5`: an orientation of `K₅` minus **one** edge already achieves `H=7`.
- `n=6,7`: **both** `H=7` and `H=21` are achievable with a single tie.
- And **even** `H`-values appear the moment the graph is incomplete (Rédei's
  odd-`H` needs completeness).

So `{7,21}` live *only* on the complete (tournament) stratum. The entire
H-impossibility phenomenon — the converse-of-Rédei gaps, the parity — is a
feature of **completeness**, removed by one tie. The oriented-graph family is
graded by **tie-count** (codistance from complete), and tournaments are the
special dense top where the gaps and Rédei concentrate.

This reframes the gaps cleanly: `{7,21}` are not about the *number* 21, they are
about *completeness forcing the conflict graph to grow* (the pancyclic
odd-cycle blow-up of the H=21 proof). Drop completeness and the forcing is gone.

## Which tie-families are relevant (the answer)

Four families, each tied to a proof mechanism:

**1. Strongly-connected oriented graphs = the atoms.** `H` is multiplicative over
strong components (the partition-function / equidecomposability monoid,
HYP-2183). Ties that **disconnect** factor `H` — or kill it: a Hamiltonian path
needs a **traceable** graph, so a disconnecting tie gives `H=0` (the "vacuum").
The H=21 proof reduced to *strong tournaments* (the indecomposable carriers).
So the relevant incomplete objects on the "structure" side are the
strongly-connected pieces and their condensation.

**2. The conflict graph `Ω` (a general graph; non-edges = ties = vertex-disjoint
odd cycles).** `H = I(Ω,2)`. This is where the count *actually lives*, and the
H=21 proof is entirely about `Ω`'s **tie (independence) structure** — `α₂`
counts ties (disjoint pairs). `Ω` *is* the relevant incomplete object for the
H-impossibility proofs; the tournament is just its source.

**3. One-tie neighbors (tournament minus one edge).** The minimal incomplete
probes: they dissolve the gaps and certify the completeness phenomenon. They are
also the natural **induction neighbors** of tournaments — and the `n→n+2`
source/sink recursion adds *near-tie* vertices (a source/sink is "tied" to the
cyclic core in the sense of contributing no cycles).

**4. LRC threshold-tie graphs (the antipodal matching).** In LRC a runner is
near / far / **exactly at the threshold** `1/n`. The binding antipodal pairs
`{a, n−a}` sit at distance *exactly* `1/n` — these are **ties** (neither near nor
far). At each tight witness `t=j/n` the tie-graph is the antipodal **perfect
matching**, and tight configs are characterized by this tie structure (S679: the
matching = the orbits of the negation involution; an agent has since extended
this to the even-`n` second apex `n/2`, HYP-2259/S700). So the relevant LRC
incomplete object is the **threshold matching**.

## The unifying picture

Ties **stratify** the carrier by codistance from complete. Tournaments are the
dense top, where impossibilities and parity concentrate. Our proofs *descend*
through the strata:

- **strong-component factorization** = disconnecting ties (multiplicativity);
- **the conflict-graph projection** = reading the tournament's odd-cycle ties as
  an independence (general-graph) problem;
- **one-tie induction** = stepping off the complete stratum;
- **LRC threshold matchings** = the ties at the tight wall.

The "useful family" the prompt points at is precisely: **strongly-connected
atoms + the conflict graph `Ω` + the threshold-tie matchings** — the incomplete
objects on which `H` factors, `H` is computed, and LRC tightness is detected.

## Next

1. **Tie-graded H-spectrum:** for each tie-count `t`, what `H`-values appear?
   (`t=0`: the monoid with gaps `{7,21}`; `t≥1`: everything.) Is there a clean
   "gap-closing at `t=1`" theorem?
2. **Traceability threshold:** characterize which disconnecting/tie patterns send
   `H→0` vs factor it — the oriented-graph version of the strong-component law.
3. **LRC via partial-tournament induction:** prove LRC by inducting on ties
   (remove the apex/threshold edges to a simpler oriented graph), the
   partial-tournament analogue of the H=21 strong-component reduction.
