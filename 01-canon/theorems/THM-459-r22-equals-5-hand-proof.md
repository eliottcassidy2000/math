# THM-459: R(2,2) = 5 — structure lemmas + machine-closed case analysis

**Status:** PROVED with a machine-closed final step (Glucose3 UNSAT, 532 CEGAR
clauses, complete verifier; `erdos592_verify_fix_macmini_s1.py`, re-derived by the
S2 persistent SAT verifier). The lemma layer below is PROVED as stated and reduces
the theorem to a bounded case analysis; a fully hand-written closure remains an
open polish task (HYP-2365).
**Source:** mac-mini-2026-06-09-S2 (context: THM-453 D/E)

**Statement.** R(2,2) = 5: every triangle-free graph on the 25 leaves of the
5-ary height-2 tree contains an independent binary subgrid (2 roots, 2 leaves
under each); and 5 is least — a 35-edge triangle-free graph on the 4-ary tree
dominates all its binary subgrids (witness in
`05-knowledge/results/erdos592_verify_fix_macmini_s1.out`).

## Structure lemmas

Rows a ∈ [5], columns [5]; R_a = within-row graph; N_{a'}(a,x) ⊆ [5] = the
cross-neighbourhood of leaf (a,x) in row a′. Let G be triangle-free on [5]² with
no independent binary subgrid (a putative witness).

**L1 (rows).** Each R_a is triangle-free on 5 vertices; hence either R_a ≅ C₅
(the unique triangle-free graph on 5 vertices with independence number 2 —
R(3,3) = 6 criticality) or α(R_a) ≥ 3.

**L2 (cross-neighbourhoods are independent).** For all a ≠ a′ and x:
N_{a'}(a,x) is R_{a'}-independent; dually each fibre {x : y ∈ N_{a'}(a,x)} is
R_a-independent. (Two cross edges sharing an endpoint plus a within-row edge
close a triangle.) In particular |N_{a'}(a,x)| ≤ α(R_{a'}).

**L3 (doubly-dark clique lemma).** For every R_a-independent pair {x₁,x₂} and
every a′ ≠ a, the doubly-dark set D = [5] ∖ (N_{a'}(a,x₁) ∪ N_{a'}(a,x₂))
induces a clique in R_{a'}; by triangle-freeness |D| ≤ 2, and if |D| = 2 then D
is an R_{a'}-edge. Equivalently |N(x₁) ∪ N(x₂)| ≥ 3.
*Proof.* An R_{a'}-independent pair {y₁,y₂} ⊆ D would make the subgrid
(a,{x₁,x₂}; a′,{y₁,y₂}) independent. ∎

**L4 (trace cliques).** Fix a′ and an R_{a'}-independent triple Y. For any other
row a: the columns x with N_{a'}(a,x) ∩ Y = ∅ form an R_a-clique (so ≤ 2 of
them, forming an edge if 2); columns with the same singleton trace into Y
likewise form an R_a-clique.
*Proof.* Two such columns R_a-independent would leave an independent pair of Y
doubly dark (in the empty case any pair of Y; in the equal-singleton-{y} case a
pair of Y ∖ {y}, which exists and is independent since |Y| = 3), contradicting
L3. ∎

**L5 (composition-freeness).** For a < a′ < a″:
B_{aa′} ∘ B_{a′a″} ∩ B_{aa″} = ∅ (three cross edges across three rows close a
triangle). [THM-453 F1's Schur condition, row-resolved.]

## The reduction

By L1 every row is C₅ or carries an independent triple; the regimes clash:
* **C₅-heavy:** cross-neighbourhoods INTO a C₅ row have size ≤ 2 (L2), while L3
  demands |N(x₁) ∪ N(x₂)| ≥ 3 for each of the five diagonal (independent) pairs
  of a C₅ source row — near-disjoint, near-maximal neighbourhood unions must
  coexist across all 10 ordered row pairs.
* **Triple-bearing:** an independent triple Y forces in every other row ≥ 3
  columns with nonempty, pairwise-distinct-or-clique Y-traces (L4), with all
  fibres R_a-independent (L2) — a dense rigid incidence pattern toward Y from
  all four other rows; L5 forbids the three-row compositions such dense
  patterns force.

Both regimes reduce to finite case analyses over the 14 triangle-free graphs on
5 vertices per row. The machine closure certifies no assignment survives. The
lemma layer is rank-free (L3/L5 lift verbatim to the THM-460 miniatures) and is
the skeleton a hand-written closure should follow.

## Lower bound

Q(2,4) is SAT: a 35-edge triangle-free graph on [4]² dominating all its binary
subgrids exists (complete verification). Hence R(2,2) = 5 exactly.
