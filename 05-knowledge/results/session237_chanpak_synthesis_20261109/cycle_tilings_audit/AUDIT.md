# Independent audit: universal cycle carrier and completion tax

Verdict: **PASS WITH TWO MINOR WORDING CORRECTIONS, APPLIED**. No mathematical defect found in Propositions A/B, their all-completion quantifiers, or the stated limitations. Chan–Pak remains OPEN.

## Independent proof review

The construction gives every graph edge its own shared tournament vertex, and every graph vertex a disjoint private string. Distinct prescribed long cycles share at most one vertex, so no unordered pair can receive incompatible prescribed orientations. The private chord closes a triangle and cannot be an edge of a long cycle because the private string has length at least five. All selected long cycles have odd length. This proves the induced corona for **every completion**, at every graph order, without extrapolating from computation.

The corona independent-set identity follows by first choosing its long-cycle part S, then freely choosing any private triangle except the |S| conflicting ones. This is a literal injection of packing sets, hence a coefficientwise lower bound, not merely a numerical inequality at x=2.

For each private fourth vertex a_i, the induced tournament on C_i is strong because it contains a spanning directed cycle. Every vertex v of a strong tournament belongs to a directed triangle: on a shortest directed cycle through v, the chord between v and the second successor either makes a triangle through v or shortens that cycle. Thus there is an extra triangle through a_i inside C_i. It cannot equal the designated triangle of i, any private triangle of j, or a selected long cycle. Extra triangles for different i cannot coincide because their supports would be contained in an intersection of at most one vertex. Therefore there are at least m new singleton packings. All other differences are nonnegative by the original injection. This proves the +m*x tax with no assumption that these new triangles are pairwise disjoint.

The latter is a load-bearing boundary: the audit does **not** strengthen m*x to a binomial packing factor. The count-preserving-reduction obstruction is correctly scoped to selected-carrier inclusions with extras, not arbitrary transformations.

The zero-vertex graph is separately valid using the empty selected carrier. For empty graphs on m>=1 vertices, the same construction and formula work. For complete graphs, incident-cycle intersections are still single edge vertices, so the arc compatibility argument does not deteriorate.

The hereditary-closure corollary is logically valid: induced universality forces every induced-hereditary property of all Omega(T) to hold for all finite graphs. This does not assert that every graph is itself Omega(T); K3 is a genuine hostile example.

## Corrections requested and checked

1. 'For m=1 the construction actually gives H=9' was too broad under the all-completions quantifier. It now says **the lexicographic completion** gives H=9; all sixteen completions have H in {9,11,13,15}.
2. 'Support-forgetting convention' inverted the intended quotient. It now refers to the **cycle-order-forgetting, support-only convention**.

The producer confirmed both changes. They do not change the propositions.

## Independent computation

Reproduce:

    python3 05-knowledge/results/session237_chanpak_synthesis_20261109/cycle_tilings_audit/check.py

The audit imports no producer code. It uses symbolic edge/private vertex labels, memoized graph deletion for independence polynomials instead of packing-subset enumeration, recursive Hamilton continuation counts instead of producer forward DP, and DFS directed-cycle enumeration instead of permutation-by-support enumeration. Checks use explicit exceptions and survive python -O.

Verified:

- All **1,100** labelled graphs of orders zero through five, including the omitted empty boundary; this independently confirms the producer's 1,099 nonempty-order inputs and maximum constructed order **37**.
- All selected cycles and exact corona conflicts, vertex count, compatible prescribed arcs, and distinct extra triangles.
- Corona polynomial via graph deletion versus independent product expansion.
- **112** randomized completions of complete and empty graphs of orders one through seven. These are finite controls, not the basis of the universal completion assertion.
- All **16** K1 completions: independent DFS cycle enumeration, independence-polynomial recursion, Hamilton continuation counts and all-permutation counts agree. H-values are exactly **{9,11,13,15}**. Every completion satisfies the coefficient tax.
- Exact lexicographic controls H(empty)=1, H(K1 construction)=9, H(two isolated input vertices)=81, H(K2 construction)=4617.

Output: `results.json`; optimized-Python output `results.optimized.json` is byte-identical. The producer's tested computational claims match these independently reproduced values. The initial producer used Python assertions. Before final publication these were replaced by explicit RuntimeError checks, and normal/optimized outputs matched. This audit's independent checks are also optimization-safe.

## Provenance and type review

Targeted review confirms the old induced C5/claw predecessors in THM-019 (`01-canon/theorems/THM-019-omega-perfectness.md`). Its broader legacy prose is not needed. THM-029 (`01-canon/theorems/THM-029-alpha1-gap.md`) and THM-343 (`01-canon/theorems/THM-343-H7-impossible.md`) support the extra-cycle warning and the full-carrier K3 obstruction, rather than induced impossibility. THM-4940 (`01-canon/theorems/THM-4940-odd-cycle-packing-counts-in-every-base.md`) credits OCF and the older cycle-count mechanisms. No external novelty claim is certified.

THM-5040 (`01-canon/theorems/THM-5040-sharp-even-graph-chan-pak.md`) uses automorphism edge-reversal evenness, exactly as the producer says, not even degrees. Its theorem realizes 21 whereas tournament H does not. Therefore the stated impossibility of an exact value-preserving map from **all** such even graphs to tournaments is correct. Neither this mismatch nor induced universality rules out a restricted or transformed transfer.

The inclusion-exclusion formula for graph dominating sets is correctly signed and uses unions of closed neighborhoods. The graph cone recurrence and tournament-source/sink invariance are correctly distinguished. The tiling references are contextual coordinate maps, not used to prove the universal carrier proposition. No ternary scalar claim has been silently substituted for binary H.
